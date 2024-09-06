#! /usr/bin/env python3
"""Automatic Quality control and Report Generation for Quantum Chemistry.

QuChemReport is an open-source software for molecular quantum chemistry.
It depends on Openbabel and cclib for the log results manipulation,
orbkit for the molecular orbitals and density discretization,
and NumPy for convenient and fast array manipulation.

QuchemReport can read any logfile results supported by cclib.
After a quality check, it converts the extracted information in JSON.
QuChemreport can merge the results of different logfiles if
they concern the same molecule wuth the same calcualtion parameters,
and there is a Gound state Optimization and/or Frequency calculationi.
Time Dependent and Excited states optimization and also Single Points
calculations can also be treated to generate corresponding results.
The properties are automatically processed and put in a report with a pylatex formatting that contains text data and images.

Required informations :
Ground state should be provided by at least one OPT or Freq
Geometry (atomic positions) has to be present for 3D molecular views (Topo)
Molecular coefficients (MO coeff) and basis sets have to be present in a log for discretization process
TD, Freq, OPT_ES, Freq_ES job types lead a spectrum (UV or IR) process
TD and MO or OPT_ES/Freq_ES and MO launch an electron density differences (EDD) discretization process
Different charges between oxydation states launch a Fukui process and EDD

Conformity check :
All calculations must have the same calculation method (functional and basis set).
All log files must correspond to the same molecular entity / conformer.
That means the same formula and geometry for ground states and oxydation states
Optimized parent excited states like S1, T1 are also accepted.

A restart procedure is possible to pass the discretization process
and use the already available png files.

"""

# Copied from scipy DOCLINES = (__doc__ or '').split("\n")

import os
import sys
import json
import requests
import tempfile
import argparse
import psutil
#from scanlog.scanlog import process_logfile_list

from quchemreport.parser.scanlog import process_logfile_list
from quchemreport.parser import conformity
from quchemreport.visualization import visualization, latex_report

nproc = psutil.cpu_count(logical=False)
availableMem = int(0.9*psutil.virtual_memory().available / 1e+9)

# Arguments parsing
argparser = argparse.ArgumentParser(description='Start QuChemReport')
argparser.add_argument('logfiles', type=str, nargs='+',
                       help='The JSON or log files to include in the report')
argparser.add_argument('--restart', '--r', action='store_const', const=1, default=0,
                       help='Restart mode. Images won\'t be regenerated if they already exist.')
argparser.add_argument('--nproc', '--n', type=int, default=nproc,
                       help='Allocatd cores for discretization process.')
argparser.add_argument('--mem', '--m', type=int, default=availableMem,
                       help='Allocated memory for discretization process (in GB).')
argparser.add_argument('--mode', type=str.lower, default='full', choices=['full', 'si','text'],
                       help='Report verbosity.')
argparser.add_argument('--MEP', action='store_true',
                       help='Generate Molecular Electrostatic Potential images (RAM expensive!)')
argparser.add_argument('--noMO', action='store_true',
                       help='Do not generate Molecular Orbitals images (RAM expensive!)')
argparser.add_argument('--noEDD', action='store_true',
                       help='Do not generate Electronic Density Difference images (RAM expensive!)')
argparser.add_argument('--verbose', '--v', action='store_true',
                       help='Verbosity level')
args = vars(argparser.parse_args())

verbose = args['verbose']
restart = args['restart']
mode = args['mode']

print('Starting QuChemReport...')

tmpdirname = tempfile.mkdtemp()
if verbose:
    print('Creating temporary directory:', tmpdirname)

if verbose:
    if restart:
        print('Requested a restart process, some discretization will not be done if the correct png file exists')
    else:
        print('Requested a full process, all discretizations will done ')

    print('Requested %d processors' %nproc)

    if mode == "full":
        print('Requested a full report')
    else:
        print('Requested a minimalist report for SI')

# Try to read the other arguments as JSON links or log files
# jf is the list of created json
#try:
#    jf = []
#    for api_address in sys.argv[4:]:
#        response = requests.get(api_address)
#        jf.append(response.json()['data'])
#    print("JSON input detected")
#except:
#    print("String could not be converted to JSON. Log files detected")
    # Read a list of log_files, parse them with quality check and return a list of JSON
    # lf is the list of parsed log files in temp directory

# Parsing
print("\nParsing input files...")
lf, jf = process_logfile_list(args['logfiles'], log_storage_path=tmpdirname, verbose=verbose, sparse=False)
# checking that jf is not empty before any processing.
tmpdirname = tempfile.mkdtemp()
if (len(jf)) == 0:
    print("There is no valid JSON file. Program will exit.")
    sys.exit()
print(len(jf), "valid files detected.")

# CONFORMITY: Uniformity tests followed by Parenting and Discretization tests
print('\nStarting conformity tests.')
data, job_types, nres_noES, charges, charge_ref, discret_proc, mo_viz_done, data_for_discretization = conformity.tests(args, jf)

# Generation of a single report for the conformer

# Process discretization and visualization works depending on the job types
print('\nStarting discretization and visualisation process.')
visualization.jobs(args, jf, data)
print('Discretization and visualisation process done.')

# Generate tex and pdf report
print('\nGenerating report.')
latex_report.json2latex(args, jf, data, mode="clean")



