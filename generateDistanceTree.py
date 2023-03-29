# imports
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
import os, sys



def generateDistanceTree(alignmentFile):
    if not os.path.isdir("./results/"): os.makedirs("./results/")
    if not os.path.isdir("./results/tree/"): os.makedirs("./results/tree/")
    if os.path.isfile(alignmentFile):
        aln = AlignIO.read(alignmentFile, 'fasta')
        print("Generating distance tree for : ", alignmentFile )
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        constructor = DistanceTreeConstructor(calculator, 'nj')
        tree =constructor.build_tree(aln)
        Phylo.write(tree, "results/tree/"+ alignmentFile.replace(".fasta", "_nj_dist_tree.tree"), "newick")
        print("finish, tree generated")

    else: print("error, file not found, please check your path")

def main():
     if len(sys.argv) == 2:
        generateDistanceTree(sys.argv[1])

     else: print("error, the script take one arugment, the file path ")


if __name__ == '__main__':
    main()