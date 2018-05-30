#!/usr/bin/env python
"""
Automatic version of grab_gene_descriptions.py

Scrapes NCBI for descriptions matching gene name, then
picks most common. Those with missing or empty string
descriptions are printed at top of output file.

Intended to use on genes in Gene2Products.need-curating.txt
from funannotate annotate formatted as single column,
new-line separated text file.

Outputs 2 column TSV ready for update-gene2products.py

Usage: python auto_gene_descriptions.py <genes.txt> <outfile.txt>

Cam Gilchrist
2018-05-29
"""

import sys
from Bio import Entrez
from collections import Counter

# *Always* tell NCBI who you are
Entrez.email = "cameron.gilchrist@research.uwa.edu.au"

def read_genes(gene_file):
    """Read in list of gene names from \n separated text file and
    return list."""
    genes = []
    with open(gene_file, 'rU') as genefile:
        for gene in genefile:
            gene = gene.strip()
            genes.append(gene)
    return(genes)

def retrieve_descriptions(gene, descriptions, empties):
    """Given single gene name, grab possible descriptions from NCBI
    and prompt user to select one"""

    # Perform ESearch and grab list of IDs
    query = gene + '[Gene Name]'
    handle = Entrez.esearch(db='gene', term=query,
            retmax=30,
            retmode='xml')
    record = Entrez.read(handle)
    idlist = ','.join(record["IdList"])
    handle.close()

    # Ensure you have results, exit if not
    if idlist == '':
        empties.append(gene)
        print('{} has empty description.'.format(gene))
        return

    # Generate summary from UID list
    handle = Entrez.esummary(db='gene', id=idlist)
    record = Entrez.read(handle)
    handle.close()

    # Grab description, counter for unique values
    desc_cnt = Counter()
    doc_sums = record[u'DocumentSummarySet'][u'DocumentSummary']
    for i in range(len(doc_sums)):
        # Use NomenclatureName first if not empty 
        if doc_sums[i][u'NomenclatureName'] != '':
            desc = doc_sums[i][u'NomenclatureName']
            desc_cnt[desc] += 1
        # Otherwise add from OtherDesignations
        else:
            descs = doc_sums[i][u'OtherDesignations'].split('|')
            for desc in descs:
                # Ignore descriptions like 'protein IMPA1'
                if desc == 'protein {}'.format(gene):
                    continue
                desc_cnt[desc] += 1

    desc = desc_cnt.most_common(1)[0][0]
    # Check for empty descriptions
    if desc == '':
        print('{} has empty description.'.format(gene))
        empties.append(gene)
    else:
        descriptions[gene] = desc
        print('{} has {} unique descriptions from {} results. Most common is:\n{}'.format(
            gene, len(desc_cnt), len(doc_sums), desc))

    return(empties)

def print_descriptions(descriptions, empties, outfile):
    """Print descriptions as 2 column TSV for update-gene2products.py"""
    with open(outfile, 'w') as out:
        out.write('Genes with empty descriptions:\n')
        for gene in empties:
            out.write(gene + '\n')
        out.write('\nGenes with non empty descriptions:\n')
        for gene in descriptions:
            out.write('{}\t{}\n'.format(gene, descriptions[gene]))

# Read in genes from file and summarize
genes = read_genes(sys.argv[1])
print('There are {} genes in {}. These are:\n{}\n'.format(
    len(genes), sys.argv[1], ', '.join(genes))
    )

# Fetch descriptions
descriptions = {}
empties = []
for gene in genes:
    retrieve_descriptions(gene, descriptions, empties)

# Write to output file given in second argument
print_descriptions(descriptions, empties, sys.argv[2])
print('\nAll done. Remember to check {} to correct errors or make adjustments!'.format(sys.argv[2]))
