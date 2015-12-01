"""
Suite of commonly used tools for high throughput python scripting.
"""

__author__ = 'Michael Ashton'
__maintainer__ = 'Michael Ashton'
__email__ = 'ashtonmv@gmail.com'
__date__ = 'Dec 1, 2015'

import os


def listdirs(excluded=[]):
    """
    Returns a list of all subdirectories in the current working directory.
    """
    dirs = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d) and 
            d not in excluded]
    return dirs

def write_runjob(name, nnodes, nprocessors, pmem, walltime, command):
    """
    Write a general queue submission script for Hipergator. Very similar to the
    write_runjob in vasp_tools, but does not assume the use of a vasp binary.
    """
    runjob=open("runjob","w")
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
    runjob.write("%s\n\n" % command)
    runjob.write("echo 'Done.'\n")
    runjob.close()


def send_email(subject, content):
    """
    Use Gmail's SMTP server to send yourself an email. Useful for appending at
    the end of time-consuming scripts.

    subject, content (str): subject and content of email to be sent.
    """
    gmail_address = CREDS['gmail_address']
    gmail_password = CREDS['gmail_password']
    message = MIMEText(content)
    message['Subject'] = 'survey.py complete'
    message['From'] = gmail_address
    message['To'] = gmail_address
    try:
        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.ehlo()
        server.starttls()
        server.login(gmail_address, gmail_password)
        server.sendmail(gmail_address, gmail_address, message.as_string())
        server.close()
        print 'email sent'
    except:
        print 'failed to send email'
