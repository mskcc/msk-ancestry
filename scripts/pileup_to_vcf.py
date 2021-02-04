#!/usr/bin/env python

# Version: 0.1 (2021-02-04)
# Author: Kanika Arora (arorak@mskcc.org)
# Adapted from Conpair scripts (https://github.com/nygenome/Conpair)

import sys
import os
import optparse
import math
import imp
import glob
from collections import defaultdict
from tempfile import NamedTemporaryFile
from shutil import move

desc = """Program to generate genotype matrix for a sample based on pileup file"""
parser = optparse.OptionParser(version='%prog version 0.1 2021-02-04', description=desc)
parser.add_option('-P', '--pileup_file', help='Input pileup file [required]', action='store')
parser.add_option('-C', '--min_cov', help='MIN COVERAGE TO CALL GENOTYPE', default=10, type='int', action='store')
parser.add_option('-O', '--outfile', help='VCF OUTPUT FILE [required]', type='string', action='store')
parser.add_option('-Q', '--min_mapping_quality', help='MIN MAPPING QUALITY', default=10, type='int', action='store')
parser.add_option('-B', '--min_base_quality', help='MIN BASE QUALITY', default=20, type='int', action='store')
parser.add_option('-M', '--marker_file', help='Marker TXT file', action='store')
parser.add_option('-R', '--repository', help='Directory with required modules scripts', default=os.path.dirname(os.path.abspath(sys.argv[0])))
parser.add_option('-F', '--reference_fasta', help='Reference FASTA file [required]', action='store')
parser.add_option('-N', '--sample_name', help='Name of sample in VCF file [required]', type='string', action='store')
(opts, args) = parser.parse_args()

if opts.repository:
    DIR=opts.repository

sys.path.append(DIR)

Marker = imp.load_source('/Marker', DIR + '/Marker.py')

if not opts.pileup_file or not opts.outfile or not opts.marker_file or not opts.reference_fasta or not opts.sample_name:
    parser.print_help()
    sys.exit(1)



MARKER_FILE = opts.marker_file
REFERENCE = opts.reference_fasta
REMOVE_CHR_PREFIX = "no"

if not os.path.exists(MARKER_FILE):
    print('ERROR: Marker file {0} cannot be found.'.format(MARKER_FILE))
    sys.exit(2)

if not os.path.exists(opts.pileup_file):
    print('ERROR: Input pileup file {0} cannot be found.'.format(opts.pileup_file))
    sys.exit(2)

outfile = open(opts.outfile, 'w')

Markers = Marker.get_markers(MARKER_FILE)
COVERAGE_THRESHOLD = opts.min_cov
MMQ = opts.min_mapping_quality
MBQ = opts.min_base_quality
##AA_BB_only = opts.normal_homozygous_markers_only

genotype_likelihoods = Marker.genotype_likelihoods_for_markers(Markers, opts.pileup_file, min_map_quality=MMQ, min_base_quality=MBQ)
D = genotype_likelihoods


G = ['0/0', '0/1', '1/1']

outfile.write("##fileformat=VCFv4.1\n")
outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+opts.sample_name+"\n")
for marker in Markers:
    line_to_output = "\t".join([Markers[marker].chrom, Markers[marker].pos, Markers[marker].id, Markers[marker].ref, Markers[marker].alt])
    line_to_output += "\t100\tPASS\t.\tGT"
    GL = D[marker]
    if GL is None or GL['coverage'] < COVERAGE_THRESHOLD:
        line_to_output += "\t./."
	#continue
    else:
        genotype = G[GL['likelihoods'].index(max(GL['likelihoods']))]
        line_to_output += "\t" + genotype
    outfile.write(line_to_output + "\n")
outfile.close()

#  ===============================
#  Update chromosome prefix if
#  requested
#  ===============================
if REMOVE_CHR_PREFIX == "yes":
    print("Removing 'chr' prefix...")

    with NamedTemporaryFile(delete=False) as tmp_source:
        with open(opts.outfile) as source_file:
            for line in source_file:
                if line.startswith("chr"):
                    tmp_source.write(line[3:])
                else:
                    tmp_source.write(line)

    move(tmp_source.name, source_file.name)

#  ===============================
#  Set output permissions
#  ===============================
try:
    os.chmod(opts.outfile, 0660)
except OSError:
    pass # only the owner can change the file permissions
