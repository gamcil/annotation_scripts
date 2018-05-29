#!/usr/bin/env python
"""
Fetch descriptions from NCBI given file with gene names.

Intended to use on genes in Gene2Products.need-curating.txt
from funannotate annotate.

Outputs 2 column TSV ready for update-gene2products.py

Usage: python grab_gene_descriptions.py <genes.txt> <outfile.txt>

Cam Gilchrist
2018-05-29
"""

import sys
from Bio import Entrez
from collections import Counter

# *Always* tell NCBI who you are
Entrez.email = "REPLACE WITH YOUR EMAIL ADDRESS"

def read_genes(gene_file):
    """Read in list of gene names from \n separated text file and
    return list."""
    genes = []
    with open(gene_file, 'rU') as genefile:
        for gene in genefile:
            gene = gene.strip()
            genes.append(gene)
    return(genes)

def retrieve_descriptions(gene, descriptions):
    """Given single gene name, grab possible descriptions from NCBI
    and prompt user to select one"""

    # Perform ESearch and grab list of IDs
    query = gene + '[Gene Name]'
    handle = Entrez.esearch(db='gene', term=query,
            retmax=100,
            retmode='xml')
    record = Entrez.read(handle)
    handle.close()
    idlist = ','.join(record["IdList"])

    # Generate summary from UID list
    handle = Entrez.esummary(db='gene', id=idlist)
    record = Entrez.read(handle)
    handle.close()

    # Grab description, counter for unique values
    desc_cnt = Counter()
    doc_sums = record[u'DocumentSummarySet'][u'DocumentSummary']
    for i in range(len(doc_sums)):
        if doc_sums[i][u'NomenclatureName'] != '':
            desc = doc_sums[i][u'NomenclatureName']
        else:
            desc = doc_sums[i][u'OtherDesignations'].split('|')[0]
        desc_cnt[desc] += 1

    # Create list from counter keys for indexing purposes
    desc_list = filter(None, desc_cnt)
    if len(desc_cnt) > 1:
        print('{} has {} unique descriptions from {} results. These are:'.format(
            gene, len(desc_list), len(doc_sums)))
        ans_range = range(len(desc_list))
        for i in ans_range:
            print ('{}: {} [{}/{}]'.format(i+1, desc_list[i], desc_cnt[desc_list[i]], len(doc_sums)))

        # Take user input to accept/reject a description
        while True:
            ans = raw_input('Which do you accept? [{}-{}/N]: '.format(
                min(ans_range)+1, max(ans_range)+1))
            # Check if int or str entered
            try:
                ans = int(ans)-1
                if ans in ans_range:
                    print('Accepting #{}.\n'.format(ans+1))
                    descriptions[gene] = desc_list[ans]
                    break
                else:
                    print('{} is outside acceptable range. Try again.'.format(
                        ans))
            except:
                if ans in ['N', 'n', 'no', 'No']:
                    print('Skipping this gene.\n')
                    break
                else:
                    print('Invalid input, try again.')

    # If there's only one unique description, accept/reject
    elif len(desc_set) == 1:
        desc = next(iter(desc_set))
        print('{} only has one unique description from {} results.'.format(
            gene, len(doc_sums)))
        print('This is:\n{}'.format(desc))

        while True:
            ans = raw_input('Accept? Y/N: ')
            if ans in ['Y', 'y', 'yes', 'Yes']:
                print('Description accepted.\n')
                descriptions[gene] = desc
                break
            elif ans in ['N', 'n', 'no', 'No']:
                print('Skipping this gene.\n')
                break
            else:
                print('Invalid input, try again.')
    return(descriptions)

def print_descriptions(descriptions, outfile):
    """Print descriptions as 2 column TSV for update-gene2products.py"""
    with open(outfile, 'w') as out:
        for gene in descriptions:
            out.write('{}\t{}\n'.format(gene, descriptions[gene]))

# Read in genes from file and summarize
genes = read_genes(sys.argv[1])
print('There are {} genes in {}. These are:\n{}\n'.format(
    len(genes), sys.argv[1], ', '.join(genes))
    )

# Fetch descriptions
descriptions = {}
for gene in genes:
    retrieve_descriptions(gene, descriptions)

# Write to output file given in second argument
print_descriptions(descriptions, sys.argv[2])
print('All done. Remember to check {} to correct errors or make adjustments!'.format(sys.argv[2]))
