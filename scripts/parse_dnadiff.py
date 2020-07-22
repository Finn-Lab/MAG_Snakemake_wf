#!/usr/bin/env python2

import sys
import os

if len(sys.argv) < 5:
    print("usage: python script.py dnadiff.report bestMash.tab label prefix/suffix use_prefix")
    sys.exit()

try:
    if sys.argv[5]:
        if sys.argv[5] == "use_prefix":
            fields = os.path.abspath(sys.argv[1]).split("/")
            for n,i in enumerate(fields):
                if i == "metaspades":
                    store = n
            prefix = os.path.abspath(sys.argv[1]).split("/")[store-1]+"_"
        else:
            prefix = sys.argv[5]+"_"
except:
    prefix = ""

label = sys.argv[3]
ext = sys.argv[4]

if ext == "prefix":
    base = os.path.basename(sys.argv[1]).split(".report")[0].split(label)[-1]
elif ext == "suffix":
    base = os.path.basename(sys.argv[1]).split(".report")[0].split(label)[0]
else:
    print("ERROR! need to indicate if label is either \"prefix\" or \"suffix\"")
    sys.exit()

flag = 0
with open(sys.argv[1], "rU") as f:
    for line in f:
        if "TotalBases" in line:
            cols = line.split()
            lenref = int(cols[1])
            lenquer = int(cols[2])
        if "AlignedBases" in line:
            cols = line.split()
            aliref = cols[1].split("(")[-1].split("%")[0]
            alique = cols[2].split("(")[-1].split("%")[0]
        if "AvgIdentity" in line and flag == 0:
            flag = 1
            cols = line.split()
            ident = cols[1].split("(")[-1].split("%")[0]

linen = 0
with open(sys.argv[2], "rU") as f:
    for line in f:
        cols = line.strip("\n").split("\t")
        query = os.path.basename(cols[0]).split(".fa")[0]
        if query == base:
            dist = cols[2]
            ref = os.path.basename(cols[1])

print("%s%s\t%s\t%i\t%.2f\t%i\t%.2f\t%.2f\t%.4f" % (prefix, base, ref, lenref, float(aliref), lenquer, float(alique), float(ident), float(dist)))
