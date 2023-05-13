import Bio as Bio
from Bio import AlignIO
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator,
    DistanceTreeConstructor,
    ParsimonyScorer,
    NNITreeSearcher,
    ParsimonyTreeConstructor,
)
from typing import Union
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.Consensus import (
    bootstrap_trees,
    get_support,
    bootstrap_consensus,
    strict_consensus,
)
from ete3 import Tree as ete3Tree


def read_msa_file(filename: str, format_type: str = "fasta") -> MultipleSeqAlignment:
    with open(filename, "r") as aln:
        alignment = AlignIO.read(aln, format_type)
    return alignment


def build_distance_tree(alignment: MultipleSeqAlignment, method: str = "upgma") -> Tree:
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    print(f"Distance Matrix:\n{distance_matrix}")
    return DistanceTreeConstructor().upgma(distance_matrix)
    constructor = DistanceTreeConstructor(calculator, method)
    return constructor.build_tree(alignment)


def build_parsimony_tree(alignment: MultipleSeqAlignment) -> Tree:
    # Calculates parsimony score of target tree using given alignment
    scorer = ParsimonyScorer()
    # Searches for tree that minimizes parsimony score
    searcher = NNITreeSearcher(scorer)
    return ParsimonyTreeConstructor(searcher).build_tree(alignment)


def get_bootstrapped_trees(
    alignment: MultipleSeqAlignment,
    times: int,
    tree_constructor: str,
) -> list[Tree]:
    if tree_constructor == "distance":
        calculator = DistanceCalculator("identity")
        tree_constructor = DistanceTreeConstructor(calculator, "upgma")
    elif tree_constructor == "parsimony":
        tree_constructor = ParsimonyTreeConstructor(NNITreeSearcher(ParsimonyScorer()))
    return list(bootstrap_trees(alignment, times, tree_constructor))


def draw_phylogenetic_tree(tree: Tree) -> None:
    Phylo.draw(tree)


def analyze_msa(msa_filename: str, num_bootstrap_resamples: int) -> None:
    print(f"\n\n\nAnalyzing {msa_filename}:")

    # Load MSA FASTA file
    msa_fasta = read_msa_file(msa_filename)

    # Build tree objects
    upgma_tree = build_distance_tree(msa_fasta)
    parsimony_tree = build_parsimony_tree(msa_fasta)

    # Draw trees
    draw_phylogenetic_tree(upgma_tree)
    draw_phylogenetic_tree(parsimony_tree)

    # Convert UPGMA tree to Newick format (required for ETE3 Toolkit)
    upgma_tree_newick = ete3Tree(newick=upgma_tree.format(fmt="newick"), format=1)
    upgma_tree_newick.ladderize()

    # Convert Parsimony tree to Newick format (required for ETE3 Toolkit)
    parsimony_tree_newick = ete3Tree(
        newick=parsimony_tree.format(fmt="newick"), format=1
    )
    parsimony_tree_newick.ladderize()

    # Get comparison metrics
    # http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#id1
    # http://etetoolkit.org/docs/latest/reference/reference_tree.html?highlight=compare#ete3.TreeNode.compare
    # http://etetoolkit.org/documentation/ete-compare/
    comparison_metrics = parsimony_tree_newick.compare(upgma_tree_newick)
    print(f"Comparing Parismony Tree to UPGMA Tree:\n{comparison_metrics}")

    # Bootstrap each tree
    upgma_trees = get_bootstrapped_trees(
        alignment=msa_fasta, times=num_bootstrap_resamples, tree_constructor="distance"
    )
    parsimony_trees = get_bootstrapped_trees(
        alignment=msa_fasta, times=num_bootstrap_resamples, tree_constructor="parsimony"
    )

    # Calculate average branch support values for each bootstrapped tree
    # using the `confidence` attribute of the root node
    average_branch_support_upgma = (
        sum([get_support(tree, upgma_trees).root.confidence for tree in upgma_trees])
        / num_bootstrap_resamples
    )
    average_branch_support_parsimony = (
        sum(
            [
                get_support(tree, parsimony_trees).root.confidence
                for tree in parsimony_trees
            ]
        )
        / num_bootstrap_resamples
    )
    print(
        f"Bootstrap Analysis for {num_bootstrap_resamples} resamples:\nAverage Branch Support Value for UPGMA Tree: {average_branch_support_upgma}\nAverage Branch Support for Maximum Parsimony Tree: {average_branch_support_parsimony}"
    )


"""
Unrelated Species
1. Fungus: Saccharomyces cerevisiae (yeast) -- https://www.ncbi.nlm.nih.gov/nuccore/NG_063315.1?report=fasta
2. Plant: Arabidopsis thaliana (thale cress) -- https://www.ncbi.nlm.nih.gov/nuccore/X16077.1?report=fasta
3. Animal:  Drosophila melanogaster (fruit fly)-- https://www.ncbi.nlm.nih.gov/nuccore/NR_133559.1?report=fasta
"""
unrelated_msa_fasta = analyze_msa(
    msa_filename="unrelated/unrelated_msa.fasta", num_bootstrap_resamples=300
)

"""
Closely Related Species
1. Passer domesticus (house sparrow) -- https://www.ncbi.nlm.nih.gov/nuccore/EF462342.2?report=fasta
2. Passer montanus (Eurasian tree sparrow) -- https://www.ncbi.nlm.nih.gov/nuccore/D38344.1?report=fasta
"""
moderately_related_msa_fasta = analyze_msa(
    msa_filename="moderately_related/moderately_related_msa.fasta",
    num_bootstrap_resamples=100,
)


"""
Moderately Related Species
1. Danio rerio (zebrafish) -- https://www.ncbi.nlm.nih.gov/nuccore/KY486501.1?report=fasta
2. Gasterosteus aculeatus (three-spined stickleback) -- https://www.ncbi.nlm.nih.gov/nuccore/EG591101.1?report=fasta
3. Oreochromis niloticus (Nile tilapia) -- https://www.ncbi.nlm.nih.gov/nuccore/XR_003216134.1?report=fasta
4. Salmo salar (salmon) -- https://www.ncbi.nlm.nih.gov/nuccore/XR_006760234.1?report=fasta
5. Gadus morhua (Atlantic cod) -- https://www.ncbi.nlm.nih.gov/nuccore/XR_003975744.1?report=fasta
6. Esox lucius (northern pike) -- https://www.ncbi.nlm.nih.gov/nuccore/XR_003779613.1?report=fasta
"""
related_msa_fasta = analyze_msa(
    msa_filename="related/related_msa.fasta", num_bootstrap_resamples=100
)


"""
All species
"""
all_msa_fasta = analyze_msa(
    msa_filename="all/all_msa.fasta", num_bootstrap_resamples=300
)
