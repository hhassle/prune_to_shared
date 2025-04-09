import os
import re
from glob import glob
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree

def strip_gene_suffix(name):
    """Remove _gene* suffix from taxon name."""
    return re.sub(r'_gene\d+$', '', name)

def get_common_taxa(trees):
    """Find common taxa across all trees (ignoring _gene*)."""
    taxon_sets = []
    for tree in trees:
        taxa = {strip_gene_suffix(term.name) for term in tree.get_terminals()}
        taxon_sets.append(taxa)
    return set.intersection(*taxon_sets)

def prune_tree(tree, common_taxa):
    """Prune tree to only include common taxa."""
    to_remove = []
    for term in tree.get_terminals():
        clean_name = strip_gene_suffix(term.name)
        if clean_name not in common_taxa:
            to_remove.append(term)
    for clade in to_remove:
        tree.prune(clade)
    return tree

def main():
    # Find all .treefile and .tree files
    tree_files = glob("*.treefile") + glob("*.tree")
    
    if not tree_files:
        print("No tree files found.")
        return

    # Read trees
    trees = [Phylo.read(file, "newick") for file in tree_files]
    print(f"Loaded {len(trees)} trees.")

    # Determine shared taxa
    common_taxa = get_common_taxa(trees)
    print(f"Found {len(common_taxa)} shared taxa.")

    # Prune and write pruned trees
    for i, (file, tree) in enumerate(zip(tree_files, trees)):
        pruned = prune_tree(tree, common_taxa)
        output_file = f"{os.path.splitext(file)[0]}_pruned.tree"
        Phylo.write(pruned, output_file, "newick")
        print(f"Pruned tree saved to: {output_file}")

if __name__ == "__main__":
    main()
