#!/usr/bin/env python

from distutils.core import setup
import os, sys
from subprocess import Popen, PIPE

def find_data_files(src_dir, *wildcards):

    src_dir = src_dir.strip()
    while src_dir[-1]=='/':
        src_dir = src_dir[:-1]


    # Find all directory names
    dirList = Popen(['find '+src_dir+' -name "*"'], shell=True, \
            stdout=PIPE).communicate()[0].split()
    
    foundList = []
    for i in range(len(dirList)):
        #if (os.path.isdir(dirList[i])&(dirList[i].strip()!=src_dir)):
        if os.path.isdir(dirList[i]):
            # Find the appropriate files within each directory
            fileList = []
            for wc in wildcards:
                dum = Popen(['ls -1 '+dirList[i]+'/'+wc], shell=True, \
                        stdout=PIPE, stderr=PIPE).communicate()[0].split()
                #dum = Popen(['find '+dirList[i]+' -name "'+wc+'"'], shell=True, \
                        #stdout=PIPE).communicate()[0].split()

                if len(dum)>0:
                    fileList.extend(dum)

            if len(fileList)>0:
                foundList.append((dirList[i], fileList))
                
    return foundList    
   
files = find_data_files('./', '*.*')

python_files = find_data_files('./radmc3dPy', '*.py')[0][1]
moduleNames = []
for i in range(len(python_files)):

    ind1 = python_files[i].strip()[::-1].find('/')
    dum  = python_files[i].strip()[-ind1:-3]
    if dum.strip()!='__init__':
        moduleNames.append('radmc3dPy.'+dum)
    

setup(name='radmc3dPy',
      version='0.25',
      description='Python module for RADMC3D',
      author='Attila Juhasz',
      author_email='juhasz@strw.leidenuniv.nl',
      py_modules=moduleNames)

