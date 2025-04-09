import os
import re
import random
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

def deduplicate_tree(tree):
    """Remove duplicate taxa in a tree, randomly keeping one representative."""
    taxon_to_clades = {}
    for clade in tree.get_terminals():
        base_name = strip_gene_suffix(clade.name)
        taxon_to_clades.setdefault(base_name, []).append(clade)

    for base_name, clades in taxon_to_clades.items():
        if len(clades) > 1:
            keep = random.choice(clades)
            for clade in clades:
                if clade is not keep:
                    tree.prune(clade)

def prune_tree(tree, common_taxa):
    """Prune tree to only include common taxa."""
    to_remove = []
    for term in tree.get_terminals():
        clean_name = strip_gene_suffix(term.name)
        if clean_name not in common_taxa:
            to_remove.append(term)
    for clade in to_remove:
        tree.prune(clade)

def main():
    # Find all .treefile and .tree files
    tree_files = glob("*.treefile") + glob("*.tree")

    if not tree_files:
        print("No tree files found.")
        return

    # Read trees
    trees = [Phylo.read(file, "newick") for file in tree_files]
    print(f"Loaded {len(trees)} trees.")

    # Deduplicate each tree
    for tree in trees:
        deduplicate_tree(tree)

    # Determine shared taxa after deduplication
    common_taxa = get_common_taxa(trees)
    print(f"Found {len(common_taxa)} shared taxa.")

    # Prune and write pruned trees
    for file, tree in zip(tree_files, trees):
        prune_tree(tree, common_taxa)
        output_file = f"{os.path.splitext(file)[0]}_pruned.tree"
        Phylo.write(tree, output_file, "newick")
        print(f"Pruned tree saved to: {output_file}")

if __name__ == "__main__":
    main()
