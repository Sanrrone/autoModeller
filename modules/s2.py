from modeller import *
import urllib, os.path, sys

def download_pdb(pdbfile):
    url = str('https://files.rcsb.org/download/'+pdbfile)
    print("downloading: "+url)
    urllib.request.urlretrieve(url, 'modules/pdbs/'+pdbfile)


def compare(prfhits, alifile):
    env = Environ()
    aln = Alignment(env)

    match = [ ("".join(list(x)[:4]),"".join(list(x)[4:len(x)])) for x in prfhits.pdbcode]

    with open(alifile) as f:
        lines = f.read() ##Assume the sample file has 3 lines
        first = lines.split('\n', 1)[0]
        alicode = first.split(";")[1]

    for (pdb, chain) in (match):
        filename = pdb+".pdb"
        if not os.path.isfile('modules/pdbs/'+filename):
            download_pdb(filename)

        if len(list(chain))>1:
            m = Model(env, file='modules/pdbs/'+filename, model_segment=('FIRST:'+list(chain)[0], 'LAST:'+list(chain)[1]))
        else:
            m = Model(env, file='modules/pdbs/'+filename, model_segment=('FIRST:'+chain, 'LAST:'+chain))
        
        aln.append_model(m, align_codes=pdb, atom_files='modules/pdbs/'+filename)

 #   for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
 #                                   ((1., 0.5, 1., 1., 1., 0.), False, True),
 #                                   ((1., 1., 1., 1., 1., 0.), True, False)):
#
 #       aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
 #              rr_file='$(LIB)/as1.sim.mat', overhang=30,
 #              gap_penalties_1d=(-450, -50),
 #              gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
 #              dendrogram_file='hits.tree',
 #              alignment_type='tree', # If 'progresive', the tree is not
 #                                     # computed and all structues will be
 #                                     # aligned sequentially to the first
 #              feature_weights=weights, # For a multiple sequence alignment only
 #                                       # the first feature needs to be non-zero
 #              improve_alignment=True, fit=True, #write_fit=write_fit,
 #              write_whole_pdb=False, output='ALIGNMENT QUALITY')
 #   
#
#
#
 #   aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
 #          rr_file='$(LIB)/as1.sim.mat', overhang=30,
 #          gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
 #          gap_gap_score=0, gap_residue_score=0, dendrogram_file='hits.tree',
 #          alignment_type='progressive',
 #          improve_alignment=True, fit=False, write_fit=False,
 #          write_whole_pdb=False, output='QUALITY')


    aln.append(file=alifile, align_codes='all')
    aln_block = len(aln)

    # Structure sensitive variable gap penalty sequence-sequence alignment:
    aln.salign(local_alignment=False, rr_file='${LIB}/blosum62.sim.mat',)

    outputname = str(alifile+"_bestaligns")
    aln.write(file=outputname+".ali", alignment_format='PIR')
    #aln.write(file=outputname+".pap", alignment_format='PAP')

    return outputname+".ali"