"""
Generate complex fusion-seq amplicons that serve as input for simulate_reads.py

Inputs:
-------
fasta file of the telomere ends (include primer seq)
reference genome

Outputs:
--------
fasta files of telomere fusions with multiple insertions

Notes
-----
fasta sequence must have telomere end on right-hand-side. Make sure to reverse complement reverse strand sequence

"""

import click
import random
from badread import misc, fragment_lengths
import pysam
from scipy.stats import poisson
from scipy.stats import gamma
from sys import stderr
import numpy as np


__version__ = '0.1'


def read_fasta(args):
    fasta = {}
    fin = pysam.FastxFile(args['fasta'])
    for line in fin:
        fasta[line.name] = line.sequence
    if len(fasta) == 0:
        raise ValueError("Empty fasta file")
    return fasta


def generate_insertions(args, ref, fasta, n_seqs, frag_lengths, mu_ins, output_file):
    rev_comps = {k: misc.reverse_complement(a) for k, a in fasta.items()}
    chroms = list(ref.references)
    keys = list(fasta.keys())
    strand = ['reverse', 'forward']

    with open(output_file, 'a') as fasta_file:
        for n in range(n_seqs):
            t = random.choice(keys)
            a = fasta[t].upper()
            b = rev_comps[t]
            if args["mu_ins"] == 0:
                raise ValueError("mu_ins must be > 0")
            blocks = 0
            while not blocks:
                blocks = poisson.rvs(mu=mu_ins, size=1)[0]
                if not blocks:
                    continue

            ins_seqs = []
            tname = f"{t}:0-{len(a)}"
            names = [f">insertion_{blocks}__" + tname]
            blk = 0
            while blk < blocks:
                flen = frag_lengths.get_fragment_length()
                if flen < 15:
                    continue
                c = random.choice(chroms)
                if ref.get_reference_length(c) - flen < 1:
                    continue
                pos = random.randint(1, ref.get_reference_length(c) - flen)

                blk += 1
                seq = ref.fetch(c, pos, pos + flen).upper()
                s = random.choice(strand)
                if s == 'reverse':
                    seq = misc.reverse_complement(seq)
                ins_seqs.append(seq)
                names.append(f"{c}:{pos}-{pos + flen}")
            final_seq = a + "".join(ins_seqs) + b
            names.append(tname)
            final_name = "_".join(names)
            fasta_file.write(f"{final_name}\n")
            fasta_file.write(f"{final_seq}\n")


def generate_insertions_with_false(args, ref, fasta, n_seqs, frag_lengths, mu_ins, output_file):
    chroms = list(ref.references)
    keys = list(fasta.keys())
    strand = ['reverse', 'forward']

    with open(output_file, 'a') as fasta_file:
        for n in range(n_seqs):
            t = random.choice(keys)
            a = fasta[t].upper()
            if args["mu_ins"] == 0:
                raise ValueError("mu_ins must be > 0")
            blocks = 0
            while not blocks:
                blocks = poisson.rvs(mu=mu_ins, size=1)[0]
                if not blocks:
                    continue

            ins_seqs = []
            tname = f"{t}:0-{len(a)}"
            names = [f">insertion_{blocks}__" + tname]

            blk = 0
            while blk < blocks:
                flen = frag_lengths.get_fragment_length()
                if flen < 15:
                    continue
                c = random.choice(chroms)
                pos = random.randint(1, ref.get_reference_length(c) - flen)
                if pos + flen > ref.get_reference_length(c):
                    continue  # happens rarely
                blk += 1
                seq = ref.fetch(c, pos, pos + flen).upper()
                s = random.choice(strand)
                if s == 'reverse':
                    seq = misc.reverse_complement(seq)
                ins_seqs.append(seq)
                ins_seqs.append(ref.fetch(c, pos, pos + flen).upper())
                names.append(f"{c}:{pos}-{pos + flen}")
            final_seq = a + "".join(ins_seqs)
            names.append('False:0-0')
            final_name = "_".join(names)
            fasta_file.write(f"{final_name}\n")
            fasta_file.write(f"{final_seq}\n")

def generate_insertions_with_diff_primers(args, ref, fasta, n_seqs, frag_lengths, mu_ins, output_file):
    rev_comps = {k: misc.reverse_complement(a) for k, a in fasta.items()}
    chroms = list(ref.references)
    keys = list(fasta.keys())
    strand = ['reverse', 'forward']

    with open(output_file, 'a') as fasta_file:
        for n in range(n_seqs):
            t = random.choice(keys)
            a = fasta[t].upper()
            tt = random.choice([primer for primer in keys if primer != t])
            b = rev_comps[tt]
            if args["mu_ins"] == 0:
                raise ValueError("mu_ins must be > 0")
            blocks = 0
            while not blocks:
                blocks = poisson.rvs(mu=mu_ins, size=1)[0]
                if not blocks:
                    continue
            ins_seqs = []
            tname = f"{tt}:0-{len(b)}"
            names = [f">insertion_{blocks}__" + tname]
            blk = 0
            while blk < blocks:
                flen = frag_lengths.get_fragment_length()
                if flen < 15:
                    continue
                c = random.choice(chroms)
                pos = random.randint(1, ref.get_reference_length(c) - flen)
                if pos + flen > ref.get_reference_length(c):
                    continue  # happens rarely
                blk += 1
                seq = ref.fetch(c, pos, pos + flen).upper()
                s = random.choice(strand)
                if s == 'reverse':
                    seq = misc.reverse_complement(seq)
                ins_seqs.append(seq)
                ins_seqs.append(ref.fetch(c, pos, pos + flen).upper())
                names.append(f"{c}:{pos}-{pos + flen}")
            final_seq = a + "".join(ins_seqs) + b
            names.append(tname)
            final_name = "_".join(names)
            fasta_file.write(f"{final_name}\n")
            fasta_file.write(f"{final_seq}\n")


@click.command()
@click.argument('reference')
@click.argument('fasta')
@click.argument('number', type=int)
@click.argument('out')
@click.option("--mu-ins", help="insertion number mean (poisson distribution)", default=3, show_default=True, type=int)
@click.option("--mean-block-len", help="insertion length mean (gamma distribution)", default=150, show_default=True, type=int)
@click.option("--std-block-len", help="insertion length stdev (gamma distribution)", default=150, show_default=True, type=int)
@click.option("--add-diff-primers", help="have different primers in one fusion", is_flag=True)
@click.option("--add-false-tag", help="add fusions with a missing primer at one end", is_flag=True)
@click.option("--diff-primer-fraction", default = 0.1, type=float)
@click.option("--false-tag-fraction", default = 0.1, type=float)
@click.option("--n-reads-mean", help ="mean number of reads in a cluster", default = 12)
@click.option("--n-reads-std", help="standard deviation of the number of reads in a cluster", default = 25)


@click.version_option(__version__)
def generate_fusions(**args):
    ref = pysam.FastaFile(args['reference'])
    fasta = read_fasta(args)

    mean = args['n_reads_mean']
    std_dev = args['n_reads_std']
    alpha = (mean / std_dev) ** 2
    beta = std_dev ** 2 / mean

    unique_n_reads = round(args['number'] / (alpha * beta))

    cluster_sizes = gamma.rvs(a=alpha, scale=beta, size=unique_n_reads)
    cluster_sizes_sorted = sorted(cluster_sizes.tolist(), reverse=True)
    # generate clusters that have the same number of reads
    i = 2
    for c in cluster_sizes_sorted:
        cc = round(c)
        if cc == 0:
            continue

        print(f"Generating {cc} fusions with insertions", file=stderr)
        frag_lengths = fragment_lengths.FragmentLengths(args['mean_block_len'], args['std_block_len'])
        outfile = f"{args['out']}/fusions_{i}_{cc}.fasta"
        i += 1
        n_diff_primers = 0
        n_false_tag = 0

        if args['add_diff_primers']:
            n_diff_primers = round(cc * args['diff_primer_fraction'])
            generate_insertions_with_diff_primers(args, ref, fasta, n_diff_primers, frag_lengths, args['mu_ins'],
                                                  outfile)

        if args['add_false_tag']:
            n_false_tag = round(cc * args['false_tag_fraction'])
            generate_insertions_with_false(args, ref, fasta, n_false_tag, frag_lengths, args['mu_ins'], outfile)

        n_regular = cc - n_diff_primers - n_false_tag
        generate_insertions(args, ref, fasta, n_regular, frag_lengths, args['mu_ins'], outfile)


if __name__ == "__main__":
    generate_fusions()
