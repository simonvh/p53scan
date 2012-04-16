from distutils.sysconfig import *
from distutils.file_util import copy_file
import os
import sys

install = sys.argv[1]
programs = get_special_folder_path("CSIDL_PROGRAMS")
programdir = os.path.join(programs, 'p53scan')
if not(os.path.exists(programdir)):
        os.mkdir(programdir)
	directory_created(programdir)

source = os.path.join(EXEC_PREFIX, 'Scripts',"p53scan_gui.py") 
copy_file(source, os.path.join(programdir, "p53scan.py"))
file_created(os.path.join(programdir, "p53scan.py"))

source = os.path.join(EXEC_PREFIX, 'Scripts',"p63scan_gui.py") 
copy_file(source, os.path.join(programdir, "p63scan.py"))
file_created(os.path.join(programdir, "p63scan.py"))

