#!/usr/bin/python3
# ==============================================================================
#                       automod.py
#
# Author: Sandro Valenzuela (sandrolvalenzuelad@gmail.com) 
#
# Please type "python automod.py -h" for usage help
#
# ==============================================================================


import os, re, sys, time
from optparse import OptionParser
from Bio import SeqIO
import pandas as pd
#from math import ceil, floor

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.s1 import *
from modules.s2 import *
from modules.s3 import *
from modules.s4 import *

start_time = time.time()
Ncpu = os.cpu_count()



def fasta2ali(inputseq):
    fasta_sequences = SeqIO.parse(open(inputseq), 'fasta')
    outputAli=open(str("tmp.ali"),'w')
    # "Target pos on sgRNA", "Target context","current AA","potential AA")) #headers

    for fasta in fasta_sequences:
        name, sequence = fasta.description, str(fasta.seq)
        outputAli.write(">P1;%s\nsequence:%s:::::::0.00: 0.00\n%s*\n" % (name,name,sequence))
        outputAli.close()
        break

    return "tmp.ali", len(sequence), name

if __name__ == '__main__':
    parser = OptionParser(usage="Usage: python automod.py -s myseq.fasta -t mytemplate.pdb")
    parser.add_option("-s", "--sequence", dest="inputseq", help="Your sequence to search sgRNA in fasta format")
    parser.add_option("-t", "--template", dest="template", help="Your template in PDB format")
    parser.add_option("-i", "--identity", dest="identity", help="identity(%) of target sequences to consider valid", default=75)
    parser.add_option("-k", "--keepTMP", dest="keepTMP", help="clean tmp files used", default=False, action="store_true")

    (options, args) = parser.parse_args()

    inputseq = options.inputseq
    template = options.template
    identity = int(options.identity)
    keepTMP = options.keepTMP

    if inputseq == "":
        print("Error: No inputseq found")
        sys.exit()

    identity = identity if identity > 0 else 0

    ali_file, slen, input_pdbcode = fasta2ali(inputseq)
    #this function create build_profile.prf
    #step 1
    build_profile(ali_file)


    hits = pd.read_csv("build_profile.prf", comment="#", sep=" ", header=None, skipinitialspace=True)
    hits = hits.rename(columns={
      0:"id", 1:"pdbcode",2:"v3",3:"v4",4:"v5",5:"v6",6:"v7",7:"v8",8:"v9",9:"Laligment",10:"identity",11:"evalue",12:"sequence"})

    hits = hits[hits['identity']>=identity]
    #hits = hits[hits['Laligment']/slen > 0.25]
    knowncodes = list(set(["".join(list(x)[:4]) for x in hits.pdbcode]))

    print("#################################################################################")
    print(hits)
    print("#################################################################################")
    #step2
    align_ali_file = compare(hits, ali_file)

    #step 3
    print("using knowncodes: "+",".join(knowncodes))
    bestpdb_file = buildModel(align_ali_file, knowncodes, input_pdbcode)
    print("#################################################################################")
    
    #step 4
    evalModel(bestpdb_file)

    if not keepTMP:
        os.remove("hits.tree")
        os.remove("build_profile.ali")
        os.remove("build_profile.prf")
        os.remove(ali_file)
        os.remove(align_ali_file)
        os.remove(input_pdbcode+".D00000001")
        os.remove(input_pdbcode+".ini")
        os.remove(input_pdbcode+".rsr")
        os.remove(input_pdbcode+".sch")
        os.remove(input_pdbcode+".V99990001")

    print("Done")
