#!/usr/bin/env python

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

import math
from Bio.SeqUtils.CodonUsageIndices import SharpEcoliIndex
from Bio import Data, SeqIO

"""Adapted from https://github.com/Benjamin-Lee/biopython/blob/master/Bio/SeqUtils/CodonUsage.py
Correct CAI Implementation and Multiple Genetic Codes:
Fixes an incorrect implementation of the CAI formula as defined by Sharp and Li. They state "that if a certain codon is never used in the reference set then the CAI for any other gene in which that codon appears becomes zero. To overcome this problem we assign a value of 0.5 to any [codon] that would otherwise be zero." Also adds support for multiple genetic codes.
"""

class CodonAdaptationIndex(object):
    """A codon adaptation index (CAI) implementation.
    Implements the codon adaptation index (CAI) described by Sharp and
    Li (Nucleic Acids Res. 1987 Feb 11;15(3):1281-95).
    """

    def __init__(self):
        self.rscu = {}
        self.genetic_code = 1
        self.index = {}
        self.codon_count = {}
        self.CodonsDict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
        'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
        'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
        'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
        'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
        'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
        'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
        'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
        'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
        'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
        'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

        # this dictionary shows which codons encode the same AA
        self.SynonymousCodons = {}

    # use this method with predefined CAI index
    def set_cai_index(self, index):
        """Sets up an index to be used when calculating CAI for a gene.
        Just pass a dictionary similar to the SharpEcoliIndex in the
        CodonUsageIndices module.
        """
        self.index = index

    def generate_rscu(self, fasta_file):
        """Create an RSCU (Relative Synonymous Codon Usage) table from a FASTA file of CDS sequences.
        Takes a location of a Fasta file containing CDS sequences
        (which must all have a whole number of codons) and generates an RSCU dictionary.
        """
        # first make sure we're not overwriting an existing index:
        if self.index or self.codon_count or self.rscu:
            raise ValueError("an index has already been set or a codon count has been done. cannot overwrite either.")

        # sets up the genetic code dictionary (which is defualted to 1)
        for key, value in Data.CodonTable.unambiguous_dna_by_id[self.genetic_code].forward_table.items():
            try:
                self.SynonymousCodons[value].append(key)
            except KeyError:
                self.SynonymousCodons[value] = [key]

        # count codon occurrences in the file.
        self._count_codons(fasta_file)

        # now to calculate the index we first need to sum the number of times
        # synonymous codons were used all together.
        for aa in self.SynonymousCodons:
            total = 0.0
            rscu_for_amino_acid = []  # RSCU values are CodonCount/((1/num of synonymous codons) * sum of all synonymous codons)
            codons = self.SynonymousCodons[aa]

            for codon in codons:
                total += self.codon_count[codon]

            # calculate the rscu value for each of the codons
            for codon in codons:
                denominator = float(total) / len(codons)
                self.rscu[codon] = self.codon_count[codon] / denominator

    def generate_index(self, *args):
        """Generate a codon usage index from a the instance's RSCU.
        Takes a the instance's RSCU dictionary and generates a codon
        usage index.
        """
        # first make sure an RSCU table is set. Optionally accepts a FASTA file location for backwards compatibility.
        if not self.rscu and not self.index:
            try:
                self.generate_rscu(args[0])
            except IndexError:
                raise IndexError("No RSCU table or index set and no FASTA file location passed")

        if self.genetic_code != 1:
            change_tranlation_table(self.genetic_code)

        for aa in self.SynonymousCodons:
            codons = self.SynonymousCodons[aa]

            # now generate the index W=rscui/rscumax:
            rscu_max = max([self.rscu[codon] for codon in codons])
            for codon in codons:
                self.index[codon] = self.rscu[codon] / float(rscu_max) if self.rscu[codon] != 0 else 0.5 / float(rscu_max)  # uses an RSCU value of 0.5 if the codon is not used in the reference set

    def cai_for_gene(self, dna_sequence):
        """Calculate the CAI (float) for the provided DNA sequence (string).
        This method uses the Index (either the one you set or the one you generated)
        and returns the CAI for the DNA sequence.
        """
        cai_value, cai_length = 0, 0
        # if an RSCU table is given but no index is set or generated, generates an index
        if not self.index and self.rscu:
            self.generate_index()

        # if no index is set or generated, the default SharpEcoliIndex will be used.
        if not self.index and not self.rscu:
            self.set_cai_index(SharpEcoliIndex)
            #print("No index or RSCU set... using default SharpEcoliIndex!")

        if dna_sequence.islower():
            dna_sequence = dna_sequence.upper()

        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i + 3]
            if codon in self.index:
                try:
                    if codon not in ['ATG', 'TGG', 'TGA', 'TAA', 'TAG']:  # these codons are always one; exclude them
                        cai_value += math.log(self.index[codon])
                        cai_length += 1
                except TypeError:
#             elif codon not in ['TGA', 'TAA', 'TAG']:  # some indices may not include stop codons
#                 raise TypeError("illegal codon in sequence: %s.\n%s" % (codon, self.index))

        return math.exp(cai_value / (cai_length - 1.0))

    def _count_codons(self, fasta_file):
        with open(fasta_file, 'r') as handle:

            # make the codon dictionary local
            self.codon_count = self.CodonsDict.copy()

            # iterate over sequence and count all the codons in the FastaFile.
            for cur_record in SeqIO.parse(handle, "fasta"):
                # make sure the sequence is lower case
                if str(cur_record.seq).islower():
                    dna_sequence = str(cur_record.seq).upper()
                else:
                    dna_sequence = str(cur_record.seq)
                for i in range(0, len(dna_sequence), 3):
                    codon = dna_sequence[i:i + 3]
                    if codon in self.codon_count:
                        self.codon_count[codon] += 1
                    else:
                        raise TypeError("illegal codon %s in gene: %s" % (codon, cur_record.id))

    # this just gives the index when the objects is printed.
    def print_index(self):
        """Prints out the index used."""
        for i in sorted(self.index):
            print("%s\t%.3f" % (i, self.index[i]))




