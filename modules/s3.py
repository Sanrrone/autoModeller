from modeller import *
from modeller.automodel import *


def buildModel(aligned_ali_file, knowns, pdbcode):

  env = Environ()
  a = AutoModel(env, alnfile=aligned_ali_file,
                knowns=knowns, sequence=pdbcode,
                assess_methods=(assess.DOPE,
                                assess.GA341))
  a.starting_model = 1
  a.ending_model = 1
  a.make()

  return pdbcode+".B99990001.pdb"
