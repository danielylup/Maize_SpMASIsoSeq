#!/usr/bin/env python
import re
import sys
from collections import namedtuple

# FindOverlapGene.py: find overlapping genes in a GTF file
# Usage: python FindOverlapGene.py < genes.gtf > overlaps.txt
# The input GTF file should be tab-delimited and contain the following fields:
# seqid, source, feature, start, end, score, strand, frame, and attributes.
# This script is referred and modified from https://gist.github.com/standage/11377530

def is_overlap(range1, range2):
    return (
        range1[0] == range2[0] and 
        range1[2] >= range2[1] and 
          range1[1] <= range2[2]
      )

GeneRange = namedtuple("GeneRange", ["seqid", "start", "end", "strand"])
genes = {}
for line in sys.stdin:
  if line.strip() == "" or line.startswith("#"):
    continue

  line = line.rstrip()
  fields = line.split("\t")
  seqid = fields[0]
  start = int(fields[3])
  end = int(fields[4])
  strand = fields[6]
  attributes = fields[8]

  matches = re.match('gene_id "([^"]+)"', attributes)
  if matches:
    geneid = matches.group(1)
    if geneid not in genes:
      genes[geneid] = GeneRange(seqid, start, end, strand)
      if seqid != genes[geneid].seqid:
          raise ValueError(f"Sequence ID mismatch for gene_id '{geneid}': expected '{genes[geneid].seqid}', found '{seqid}'")
      assert seqid == genes[geneid].seqid
      # Update the end if this feature extends the gene
      if end > genes[geneid].end:
        genes[geneid] = GeneRange(
            genes[geneid].seqid,
            genes[geneid].start,
            end,
            genes[geneid].strand
        )
  else:
    continue

gene_keys = list(genes.keys())
for i in range(len(gene_keys)):
    for j in range(i+1, len(gene_keys)):
      if genes[gene_keys[i]].strand == genes[gene_keys[j]].strand: 
        if gene_keys[i].strand  == gene_keys[j].strand: 
          geneid_i = gene_keys[i]
          geneid_j = gene_keys[j]
          range_i = genes[geneid_i]
          range_j = genes[geneid_j]
          if is_overlap(range_i, range_j):
            print(f"{geneid_i}\t{geneid_j}")
        else:
          continue
