#!/usr/bin/env python2

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

