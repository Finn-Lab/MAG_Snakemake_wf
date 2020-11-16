#!/usr/bin/env python2

# This file is part of MAG Snakemake workflow.
#
# MAG Snakemake workflow is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# MAG Snakemake workflow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with MAG Snakemake workflow.  If not, see <https://www.gnu.org/licenses/>.


import argparse
import sys

def ren_fasta(args):
	n = 0
	for line in open(args.fasta_file, "rU"):
		if line[0] == ">":
			name = line.strip("\n").replace(">","")
			n += 1
			print ">%s_%i\t%s" % (args.prefix, n, name)
		else:
			print line.strip("\n")


if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Rename multifasta file')
        parser.add_argument('-f', dest='fasta_file', help='Input FASTA file')
	parser.add_argument('-p', dest='prefix', help='Header prefix')
        if len(sys.argv) == 1:
                parser.print_help()
                sys.exit()
        else:
                args = parser.parse_args()
                ren_fasta(args)	

