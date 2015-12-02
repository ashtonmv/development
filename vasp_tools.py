"""
Suite of commonly used tools for vasp calculations. These tools either cannot be
found in pymatgen, or I thought the implementation in pymatgen could be
improved.
"""

__author__ = 'Michael Ashton'
__maintainer__ = 'Michael Ashton'
__email__ = 'ashtonmv@gmail.com'
__date__ = 'Dec 1, 2015'

import os
from monty.serialization import loadfn


UFILE = loadfn(os.path.join(os.path.expanduser('~'), 'software/U_values.yaml'))
U_VALUES = UFILE['U_values']
D_ELEMENTS = UFILE['d_elements']
F_ELEMENTS = UFILE['f_elements']


def multiply_lattice_parameter(value):
    """
    multiplies your lattice parameter by a float of your choice.

    value = number to multiply your lattice parameter by
    """
    with open('POSCAR', 'r') as readcar:
        readlines = readcar.readlines()
    lattice_parameter_line = readlines[1].split()
    lattice_parameter = float(lattice_parameter_line[0]) * float(value)
    with open('POSCAR', 'w') as writecar:
        writecar.write(readlines[0])
        writecar.write('%s\n' % lattice_parameter)
        for line in readlines[2:]:
            writecar.write(line)


def direct_to_cartesian(inputcar, outputcar):
    """
    This function converts a direct POSCAR (inputcar) to a
    cartesian POSCAR (outputcar).
    """
    with open(inputcar, 'r') as direct_poscar:
        direct_poscar_lines = direct_poscar.readlines()
    name = direct_poscar_lines[0]
    lattice_parameter = direct_poscar_lines[1]
    a_vector = direct_poscar_lines[2].split()
    b_vector = direct_poscar_lines[3].split()
    c_vector = direct_poscar_lines[4].split()
    elements = direct_poscar_lines[5].split()
    stoichiometries = direct_poscar_lines[6].split()
    n_atoms = 0
    for i in range(len(stoichiometries)):
        n_atoms += int(stoichiometries[i])
    direct_atom_lines = []
    cartesian_atom_lines = []
    for i in range(n_atoms):
        direct_atom_lines.append(direct_poscar_lines[i+8].split())
        cartesian_atom_coordinates = []
        for j in range(3):
            cartesian_atom_coordinates.append(float(direct_atom_lines[i][j])
                                              * float(lattice_parameter)*(float(a_vector[j])
                                                                          + float(b_vector[j])
                                                                          + float(c_vector[j])))

        cartesian_atom_lines.append(cartesian_atom_coordinates)

    with open(outputcar, 'w') as cartesian_poscar:
        cartesian_poscar.write('%s' % name)
        cartesian_poscar.write('%s' % lattice_parameter)
        cartesian_poscar.write('%s %s %s\n' % (a_vector[0], a_vector[1],
                                               a_vector[2]))
        cartesian_poscar.write('%s %s %s\n' % (b_vector[0], b_vector[1],
                                               b_vector[2]))
        cartesian_poscar.write('%s %s %s\n' % (c_vector[0], c_vector[1],
                                               c_vector[2]))
        for element in elements:
            cartesian_poscar.write('%s ' % element)
        cartesian_poscar.write('\n')
        for stoichiometry in stoichiometries:
            cartesian_poscar.write('%s ' % stoichiometry)
        cartesian_poscar.write('\n')
        cartesian_poscar.write('cartesian\n')
        for i in range(n_atoms):
            for j in range(3):
                cartesian_poscar.write('%s ' % cartesian_atom_lines[i][j])
            cartesian_poscar.write('\n')


def add_vacuum(delta, cut=0.9):
    """
    Adds vacuum to a POSCAR.

    delta = vacuum thickness in Angstroms
    cut = height above which atoms will need to be fixed. Defaults to 0.9.
    """
    # Fix the POSCAR to put bottom atoms (if they are accidentally above
    # tolerance) at 0.0.

    with open('POSCAR', 'r') as readcar:
        poscar_lines = readcar.readlines()
    n_atoms = int(get_n_atoms())
    atom_lines = []
    for i in range(8, 8+n_atoms):
        atom_lines.append(poscar_lines[i].split())
    atom_line_2s = []
    for atom_line in atom_lines:
        atom_line_2s.append(float(atom_line[2]))
    fixable = False
    addables = []
    for atom_line_2 in atom_line_2s:
        if float(atom_line_2) > cut or\
                (float(atom_line_2) < 0.0 and 1.0 + float(atom_line_2) > cut):
            if float(atom_line_2) < 0.0 and 1.0 + float(atom_line_2) > cut:
                atom_line_2 = float(atom_line_2) + 1.0
            addables.append(atom_line_2)
            fixable = True
    if fixable:
        add_factor = 1.0 - min(addables)
    else:
        add_factor = 0.0
    new_atom_lines = []
    for atom_line in atom_lines:
        new_atom_line_2 = str(float(atom_line[2]) + add_factor)
        if float(new_atom_line_2) >= 1.0:
            new_atom_line_2 = str(float(new_atom_line_2) - 1.0)
        new_atom_lines.append('%s %s %s' % (atom_line[0], atom_line[1],
                                            new_atom_line_2))
    with open('POSCAR', 'w') as writecar:
        for line in poscar_lines[0:8]:
            writecar.write(line)
        for new_atom_line in new_atom_lines:
            writecar.write('%s\n' % new_atom_line)

    # Open files and read in values from POSCAR
    old_poscar = open("POSCAR", "r")
    new_poscar = open("new_POSCAR", "w")
    oldlines = old_poscar.readlines()
    name = oldlines[0].strip()
    lattice_constant = oldlines[1].strip()
    a_lat_par = oldlines[2].split()
    b_lat_par = oldlines[3].split()
    c_lat_par = oldlines[4].split()
    elements = oldlines[5].split()
    stoichiometry = oldlines[6].split()
    coordinate_type = oldlines[7].strip()
    n_atoms = 0
    for item in stoichiometry:
        n_atoms += int(item)
    final_atom_line = n_atoms + 8

    # Elongate c-vector by delta

    save = float(c_lat_par[2])
    c_length = float(c_lat_par[2]) * float(lattice_constant)
    c_length_plus_delta = c_length + float(delta)
    c_lat_par[2] = c_length_plus_delta / float(lattice_constant)
    scalar = c_lat_par[2] / save

    # Create list of atom coordinates and adjust their z-coordinate on the fly

    atoms = []
    for i in range(8, final_atom_line):
        atom = oldlines[i].split()
        atom[2] = float(atom[2]) / scalar
        atoms.append(atom)

    # Write updated values to new_POSCAR, copy it to old_POSCAR, then close
    # files and delete new_POSCAR

    new_poscar.write("%s plus %s\n" % (name, str(delta)))
    new_poscar.write("%s\n" % lattice_constant)
    for item in a_lat_par:
        new_poscar.write("%s " % item)
    new_poscar.write("\n")
    for item in b_lat_par:
        new_poscar.write("%s " % item)
    new_poscar.write("\n")
    for item in c_lat_par:
        new_poscar.write("%s " % item)
    new_poscar.write("\n")
    for item in elements:
        new_poscar.write("%s " % item)
    new_poscar.write("\n")
    for item in stoichiometry:
        new_poscar.write("%s " % item)
    new_poscar.write("\n")
    new_poscar.write("%s\n" % coordinate_type)
    for item in atoms:
        new_poscar.write("%s %s %s\n" % (item[0], item[1], item[2]))

    new_poscar.close()
    old_poscar.close()
    os.remove("POSCAR")

    new_poscar = open("new_POSCAR", "r")
    new_lines = new_poscar.readlines()
    poscar = open("POSCAR", "w")
    for line in new_lines:
        poscar.write(line)
    old_poscar.close()
    poscar.close()
    new_poscar.close()
    os.remove("new_POSCAR")


def write_potcar(types='None'):
    """
    Writes a POTCAR file based on a list of types.

    types = list of same length as number of elements containing specifications
    for the kind of potential desired for each element. If no special potential
    is desired, just enter 'regular', or leave types = 'None'.
    (["pv", "regular", "3"])
    """
    poscar = open("POSCAR", "r")
    lines = poscar.readlines()
    elements = lines[5].split()
    poscar.close()

    if types == 'None':
        types = []
        for i in range(len(elements)):
            types.append('regular')
    potentials = []
    for i in range(len(elements)):
        if types[i] == "regular":
            pass
        else:
            elements[i] += "_%s" % types[i]

        # If specified pseudopotential doesn't exist, try other variations.
        if os.path.exists("/home/mashton/Potentials/POT_GGA_PAW_PBE/%s/POTCAR"
                          % elements[i]):
            pass
        else:
            print "Potential file for %s does not exist. Looking for best"\
                  "variation... " % elements[i]
            if types[i] == 'regular':
                length = 0
            else:
                length = len(types[i]) + 1
                elements[i] = elements[i][:-length]
            elements[i] += "_sv"
            if os.path.exists("/home/mashton/Potentials/POT_GGA_PAW_PBE/%s/"
                              "POTCAR" % elements[i]):
                print "Found one! %s will work." % elements[i]
            else:
                elements[i] = elements[i][:-3]
                elements[i] += "_pv"
                if os.path.exists("/home/mashton/Potentials/POT_GGA_PAW_PBE/%s"
                                  "/POTCAR" % elements[i]):
                    print "Found one! %s will work." % elements[i]
                else:
                    elements[i] = elements[i][:-3]
                    elements[i] += "_3"
                    if os.path.exists(
                        "/home/mashton/Potentials/"
                            "POT_GGA_PAW_PBE/%s/POTCAR" % elements[i]):
                        print "Found one! %s will work." % elements[i]
                    else:
                        elements[i] = elements[i][:-2]
                        if os.path.exists(
                            "/home/mashton/Potentials/"
                                "POT_GGA_PAW_PBE/%s/POTCAR" % elements[i]):

                            print "Found one! %s will work." % elements[i]
                        else:
                            print "I give up looking for a pseudopotential for"\
                                " %s. Do it yourself." % elements[i]

    # Create paths, open files, and write files to POTCAR for each potential.
    for element in elements:
        potentials.append("/home/mashton/Potentials/POT_GGA_PAW_PBE/%s/POTCAR"
                          % element)
    outfile = open("POTCAR", "w")
    for potential in potentials:
        infile = open(potential)
        for line in infile:
            outfile.write(line)
        infile.close()
    outfile.close()


def write_runjob(name, nnodes, nprocessors, pmem, walltime, binary):
    """
    writes a runjob based on a name, nnodes, nprocessors, walltime, and binary.
    Designed for runjobs on the Hennig group_list on HiperGator.
    """
    runjob = open("runjob", "w")
    runjob.write("#!/bin/sh\n")
    runjob.write("#PBS -N %s\n" % name)
    runjob.write("#PBS -o test.out\n")
    runjob.write("#PBS -e test.err\n")
    runjob.write("#PBS -r n\n")
    runjob.write("#PBS -l walltime=%s\n" % walltime)
    runjob.write("#PBS -l nodes=%s:ppn=%s\n" % (nnodes, nprocessors))
    runjob.write("#PBS -l pmem=%s\n" % pmem)
    runjob.write("#PBS -W group_list=hennig\n\n")
    runjob.write("cd $PBS_O_WORKDIR\n\n")
    runjob.write("mpirun ~/bin/%s > job.log\n\n" % binary)
    runjob.write("echo 'Done.'\n")
    runjob.close()


def get_magnetic_moment():
    """
    Returns the final magnetic moment in job.log.
    """
    with open('job.log', 'r') as joblog:
        joblog_lines = joblog.readlines()
    for line in joblog_lines:
        if 'mag=' in line:
            split_line = line.split()
            mag = split_line[9]
            break
    else:
        mag = 'ERROR'
    return mag


def get_polarization():
    """
    Finds and interprets the spontaneous polarization from the OUTCAR.
    """
    with open('OUTCAR', 'r') as outcar:
        outcar_lines = outcar.readlines()
    for line in outcar_lines[::-1]:
        if 'Total electronic dipole moment:' in line:
            electronic_polarization_line = line.split()
            electronic_polarization = electronic_polarization_line[5:8]
            break
    for line in outcar_lines[::-1]:
        if 'Ionic dipole moment' in line:
            ionic_polarization_line = line.split()
            ionic_polarization = ionic_polarization_line[4:7]
            break
    for line in outcar_lines[::-1]:
        if 'volume of cell :' in line:
            cell_volume_line = line.split()
            cell_volume = float(cell_volume_line[4])
            break
    for i in range(len(outcar_lines)):
        if 'length of vectors' in outcar_lines[i]:
            vector_lengths_line = outcar_lines[i+1].split()
            vector_lengths = vector_lengths_line[:3]
            break
    total_polarization = []
    for i in range(3):
        total_polarization.append(
            (float(electronic_polarization[i]) + float(ionic_polarization[i]))
            * 1602.0 * float(vector_lengths[i]) / cell_volume)
        while total_polarization[i] > 1602.0 * float(vector_lengths[i])\
                / cell_volume:
            total_polarization[i] -= 1602.0 * float(vector_lengths[i])\
                / cell_volume
        while total_polarization[i] < -1602.0 * float(vector_lengths[i])\
                / cell_volume:
            total_polarization[i] += 1602.0 * float(vector_lengths[i])\
                / cell_volume
    return total_polarization


def get_toten():
    """
    Returns the final energy in eV from the OUTCAR.
    """
    with open("OUTCAR", "r") as outcar:
        outcar_lines = outcar.readlines()
    for line in outcar_lines[::-1]:
        if "TOTEN" in line:
            toten_line = line.split()
            break
    return float(toten_line[4])


def multiply_poscar(poscar, dimensions):
    """
    Makes a supercell out of a POSCAR. Eliminates need to instantiate a
    Structure object in pymatgen's implementation.
    """
    multiply_factor = float(dimensions[0]) * float(dimensions[1])\
        * float(dimensions[2])
    with open(poscar, 'r') as original_structure:
        original_lines = original_structure.readlines()
    name_line = original_lines[0].split()
    name = name_line[0]
    lattice_parameter_line = original_lines[1].split()
    lattice_parameter = lattice_parameter_line[0]
    a_vector = original_lines[2].split()
    b_vector = original_lines[3].split()
    c_vector = original_lines[4].split()
    elements = original_lines[5].split()
    stoichiometries = original_lines[6].split()
    coordinate_type_line = original_lines[7].split()
    coordinate_type = coordinate_type_line[0]
    start_line = 8
    original_atomic_coordinates = []
    for i in range(len(stoichiometries)):
        original_element_coordinates = []
        for j in range(int(stoichiometries[i])):
            original_element_coordinates.append(
                original_lines[start_line + j].split())
        start_line += int(stoichiometries[i])
        original_atomic_coordinates.append(original_element_coordinates)
    new_stoichiometries = []
    for i in range(len(stoichiometries)):
        new_stoichiometries.append(str(int(stoichiometries[i])
                                   * int(multiply_factor)))

    # Multiply 'a' vector by the first dimension, divide all atomic coordinates
    # by the first dimension, then replicate the structure.
    new_a_vector = []
    for i in range(3):
        new_a_vector.append(str(float(a_vector[i]) * float(dimensions[0])))
    new_atomic_coordinates = []

    for i in range(len(original_atomic_coordinates)):
        new_element_coordinates = []
        for j in range(0, dimensions[0]):
            for k in range(len(original_atomic_coordinates[i])):
                new_coordinates = []
                new_x = str(
                    float(original_atomic_coordinates[i][k][0])
                    / float(dimensions[0]) + float(j) / float(dimensions[0]))

                new_y = str(original_atomic_coordinates[i][k][1])
                new_z = str(original_atomic_coordinates[i][k][2])
                new_coordinates.append(new_x)
                new_coordinates.append(new_y)
                new_coordinates.append(new_z)
                new_element_coordinates.append(new_coordinates)

        new_atomic_coordinates.append(new_element_coordinates)

    # Multiply 'b' vector by the second dimension, divide all atomic coordinates
    # by the first dimension, then replicate the structure.
    new_b_vector = []
    for i in range(3):
        new_b_vector.append(str(float(b_vector[i]) * float(dimensions[1])))
    newer_atomic_coordinates = []

    for i in range(len(new_atomic_coordinates)):
        new_element_coordinates = []
        for j in range(0, dimensions[1]):
            for k in range(len(new_atomic_coordinates[i])):
                new_coordinates = []
                new_x = str(new_atomic_coordinates[i][k][0])
                new_y = str(
                    float(new_atomic_coordinates[i][k][1])
                    / float(dimensions[1]) + float(j) / float(dimensions[1]))

                new_z = str(new_atomic_coordinates[i][k][2])
                new_coordinates.append(new_x)
                new_coordinates.append(new_y)
                new_coordinates.append(new_z)
                new_element_coordinates.append(new_coordinates)

        newer_atomic_coordinates.append(new_element_coordinates)

    # Multiply 'c' vector by the third dimension, divide all atomic coordinates
    # by the third dimension, then replicate the structure.
    new_c_vector = []
    for i in range(3):
        new_c_vector.append(str(float(c_vector[i]) * float(dimensions[2])))
    newest_atomic_coordinates = []

    for i in range(len(newer_atomic_coordinates)):
        new_element_coordinates = []
        for j in range(0, dimensions[2]):
            for k in range(len(newer_atomic_coordinates[i])):
                new_coordinates = []
                new_x = str(newer_atomic_coordinates[i][k][0])
                new_y = str(newer_atomic_coordinates[i][k][1])
                new_z = str(
                    float(newer_atomic_coordinates[i][k][2])
                    / float(dimensions[2]) + float(j) / float(dimensions[2]))

                new_coordinates.append(new_x)
                new_coordinates.append(new_y)
                new_coordinates.append(new_z)
                new_element_coordinates.append(new_coordinates)

        newest_atomic_coordinates.append(new_element_coordinates)

    with open('supercell', 'w') as supercell:
        supercell.write('%s supercell generated by vasp_tools\n' % name)
        supercell.write('%s\n' % lattice_parameter)
        supercell.write('%s\n' % ' '.join(new_a_vector))
        supercell.write('%s\n' % ' '.join(new_b_vector))
        supercell.write('%s\n' % ' '.join(new_c_vector))
        supercell.write('%s\n' % ' '.join(elements))
        supercell.write('%s\n' % ' '.join(new_stoichiometries))
        supercell.write('%s\n' % coordinate_type)
        for i in range(len(newest_atomic_coordinates)):
            for j in range(len(newest_atomic_coordinates[i])):
                supercell.write(
                    '%s\n' % ' '.join(newest_atomic_coordinates[i][j]))


def get_n_atoms():
    """
    Returns the number of atoms in the POSCAR.
    """
    with open('POSCAR', 'r') as poscar:
        poscar_lines = poscar.readlines()
    stoichiometries = poscar_lines[6].split()
    n_atoms = 0
    for stoichiometry in stoichiometries:
        n_atoms += int(stoichiometry)
    return n_atoms


def get_n_electrons():
    """
    Returns the total number of electrons in the POTCAR.
    """
    with open('POTCAR', 'r') as potcar:
        potcar_lines = potcar.readlines()
    n_electrons = 0.0
    zvals = []
    for i in range(len(potcar_lines)):
        if 'ZVAL' in potcar_lines[i]:
            zval_line = potcar_lines[i].split()
            zvals.append(float(zval_line[5]))

    with open('POSCAR', 'r') as poscar:
        poscar_lines = poscar.readlines()
    stoich_line = poscar_lines[6]
    stoichs = stoich_line.split()
    for i in range(len(zvals)):
        n_electrons += zvals[i]*float(stoichs[i])
    return n_electrons


def get_u_values(atoms):
    """
    Obtain Hubbard U values (AFLOWLIB standard values) for each atom in a list.
    """
    return_dict = {}

    u_values = {"Ti": "4.4", "V": "2.7", "Cr": "3.5", "Mn": "4.0", "Fe": "4.6",
                "Co": "5.0", "Ni": "5.1", "Cu": "4.0", "Zn": "7.5", "Ga": "3.9",
                "Nb": "2.1", "Mo": "2.4", "Tc": "2.7", "Ru": "3.0", "Rh": "3.3",
                "Pd": "3.6", "Cd": "2.1", "In": "1.9", "Ta": "2.0", "W": "2.2",
                "Re": "2.4", "Os": "2.6", "Ir": "2.8", "Pt": "3.0", "La": "7.5",
                "Ce": "6.3", "Pr": "5.5", "Gd": "6.0", "Nd": "6.2", "Sm": "6.4",
                "Eu": "5.4", "Tm": "6.0", "Yb": "6.3", "Lu": "3.8", "U": "4.0"}
    d_elements = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                  "Ga", "Nb", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
                  "Ag", "Cd", "In", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
                  "Au"]
    f_elements = ["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
                  "Ho", "Er", "Tm", "Yb", "Lu", "U"]

    for atom in atoms:
        if atom in d_elements:
            l_value = '2'
        elif atom in f_elements:
            l_value = '3'
        else:
            l_value = '-1'
        if atom in u_values:
            return_dict[atom] = [u_values[atom], l_value]
        else:
            return_dict[atom] = ['0.0', l_value]

    return return_dict


def get_volume():
    """
    Returns the volume of the cell in the OUTCAR.
    """
    with open('OUTCAR', 'r') as outcar:
        outcar_lines = outcar.readlines()
    for line in outcar_lines:
        if 'volume' in line:
            volume_line = line
    split_line = volume_line.split()
    return float(split_line[-1])


def get_n_formula_units():
    """
    Returns the number of irreducible formula units in the POSCAR.
    """
    with open('POSCAR', 'r') as poscar:
        p_lines = poscar.readlines()
    stoichiometries = p_lines[6].split()
    float_stoichs = []
    for stoich in stoichiometries:
        float_stoichs.append(float(stoich))
    for i in range(int(min(float_stoichs)), 0, -1):
        for stoich in float_stoichs:
            if not (stoich / i).is_integer():
                break
        else:
            n_formula_units = i
            break

    return n_formula_units


def get_spacing(filename='POSCAR', cutoff=0.95):
    """
    Returns the interlayer spacing for a 2D material.
    """
    lines = open(filename).readlines()
    c_axis = lines[4].split()
    lattice_parameter = lines[1].split()
    split_coords = [line.split() for line in lines[8:8+get_n_atoms()]]
    z_coords = list()
    for coord in split_coords:
        z_coord = float(coord[2])
        if z_coord > cutoff:
            z_coord -= 1
        z_coords.append(z_coord)
    max_height = max([z_coordinate for z_coordinate in z_coords])
    min_height = min([z_coordinate for z_coordinate in z_coords])
    spacing = ((1.0 + min_height) - max_height) * float(c_axis[2])\
        * float(lattice_parameter[0])

    return spacing
