#!/usr/bin/env python3
################################################
#
#subparser 1 MCL_BGC group
#subparser 2 PAN_BGC analysys
################################################
__author__= "Fernando Garcia Guevara"
__credits__= "Fernando Garcia Guevara"
__email__= "garciafertson near gmail.com"
__status__= "Development"

import argparse
import sys
import os
import logging

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')]+sys.path
#print (sys.path)

#import module for running program
import bgc_mcl
from bgc_mcl.run import Run

def print_header():
    bgc_mcl.version()

def phelp():
    print("""
    the program recieve a list with BGC.gbk files
    run PorthoMCL to retrieve clusters_of_orthologue_genes between them (BGCcogs)

    after BGCcogs analysis calulate paired Jaccard index between pairs of BGC.gbk
    and then run MCL to build groups of BGC.gbk files with similar BGCcog content
    
    Ouputs a text file with the filenames of groups constructed

    for more information type:
    
    bgc_mcl.py porthomcl -h

    bgc_mcl.py mcl -h
    """)

##the next section add parsers and subparser to the program

if __name__ == '__main__':
    parser= argparse.ArgumentParser()
    parser.add_argument('--version', action='version', version='bgc_mcl V_%s' % bgc_mcl.__version__)
    ###############################################################################
    ##subcommands
    subparser= parser.add_subparsers(help="sub-comand help:", dest='subparser_name')
    ######################################
    ###### PARSER 1 "porthomcl" subparser#
    ortho_parser=suparser.add_parser('porthomcl',
            description='Recieve BGC.gbk list and build clusters of othologue genes 
            inside the bgc set',
            epilog='''for running call: \n
            bgc_mcl.py --list --folder ''')
    orthoinput= ortho_parser.add_argument_group('input options')
    orthoinput.add_argument('--list','-l',
            dest='list',
            metavar='FILE.TXT',
            help='list of file names to include in porthomcl',
            required=True)
    orthoinput.add_argument('--dir', '-d',
            dest='dir',
            metavar='DIRECTORY',
            help="path to directory with gbk files",
            required=True)
    orthooptions=ortho_parser.add_argument_group('PorthoMCL options')
    orthooptions.add_argument('--threads', "-t",
            dest='threads',
            metavar="INT",
            help="number of threads to employ un PorthoMCL (default 1)",
            default=1)

    ######################################
    ###### PARSER 2 "mcl" subparser######
    mcl_parser = subparser.add_parser( 'mcl',
        description='Recieve BGC.gbk list and build groups MCL based on similar gene content',
        help='Build groups based on similar gene content using MCL clustering',
        epilog='''
    for running call:\n
    bgc_mcl.py find --list --folder BGC_type \n  ''')
    mclinput = mcl_parser.add_argument_group('input options')
    mclinput.add_argument('--bgccogs',
            dest="bgccogs",
            help="file with BGCcogs, to build BGC groups (output 8.all.ort.group from PorthoMCL)"
            metavar="FILE.TXT",
            default="8.all.ort.group")
    mclinput.add_argument('--list',"-l",
            dest='list',
            metavar='FILE.TXT',
            help="list with filenames of Biosynthetic gene clusters in genbank format (BGC.gbk)",
            required=True)
    mclinput.add_argument('--dir','-d',
            dest='dir',
            metavar='DIRECTORY',
            help="path to directory containing the BGC.gbk files "
            required=True)
    mcloptions = mcl_parser.add_argument_group('mcl options')
    mcloptions.add_argument('--inf_lower',
            dest='inf_lower',
            metavar='INT',
            help="set lower limit for Infilation parameter in MCL (default 1)",
            default=2,
            type=float)
    mcloptions.add_argument('--inf_upper',
            dest='inf_upper',
            metavar='INT',
            help='set upper limit for inflation parameter in MCL (default 5)',
            default=6)
    mcloptions.add_argument('--increment',
            dest='increment',
            metavar='INT',
            help='set increment for inflation parameter in MCL scan',
            default=1)
    mclouput = mcl_parser.add_argument_group('output options',
            dest='output',
            metavar="PREFIX",
            help="prefix for MCL groups",
            default="output")

    #check whether --help is needed
    if (len(sys.argv)==1 or sys.argv[1]== '-h' or sys.argv[1]== '--help'):
        phelp()    
    #call Run module passing the arguments in here
    else:
        args=parser.parse_args()
        Run(args).main()

        
