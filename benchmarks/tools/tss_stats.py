"""TSS summary statistics comparable to those reported by Ettwiller et al. 2016
(BMC Genomics 17:199, doi:10.1186/s12864-016-2539-z).

usage: tss_stats.py TSS.gff genome.fa annotation.gff

TSS can be GFF/GTF; the annotation needs CDS features. Reported next to each
value is the figure from the paper (E. coli MG1655).
"""
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import cappable_tss as ct

COMP = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    seqs, name, parts = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(parts)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip().upper())
    seqs[name] = "".join(parts)
    return seqs


def base(seqs, seqid, pos, strand, offset):
    """Base at TSS + offset (in transcript direction), 1-based TSS."""
    p = pos + offset if strand == "+" else pos - offset
    b = seqs[seqid][p - 1]
    return b if strand == "+" else b.translate(COMP)


def main(tss_path, fasta, gff):
    seqs = read_fasta(fasta)
    sites = sorted(ct.read_sites(tss_path))
    cds = ct.read_genes(gff, ("CDS",))
    index = ct.GeneIndex(cds)
    n = len(sites)
    plus1, minus1, pair = Counter(), Counter(), Counter()
    intragenic = sense = antisense = leaderless = 0
    codon = Counter()
    for seqid, strand, pos in sites:
        b0, bm1 = base(seqs, seqid, pos, strand, 0), base(seqs, seqid, pos, strand, -1)
        plus1[b0] += 1
        minus1[bm1] += 1
        pair[bm1 + b0] += 1
        genes = index.overlapping(seqid, pos, pos)
        if genes:
            intragenic += 1
            same = [g for g in genes if g.strand == strand]
            if same:
                sense += 1
                g = same[0]
                start = g.start if strand == "+" else g.end
                if pos == start:
                    leaderless += 1
                else:
                    codon[(abs(pos - start)) % 3 + 1] += 1
            else:
                antisense += 1

    def pct(x, total=n):
        return "%5.1f%%" % (100.0 * x / total) if total else "  NA"

    yr = sum(v for k, v in pair.items() if k[0] in "CT" and k[1] in "AG")
    rows = [
        ("TSS", str(n), "16,359 clustered TSS"),
        ("+1 is A or G", pct(plus1["A"] + plus1["G"]), "80%"),
        ("-1 is C or T", pct(minus1["C"] + minus1["T"]), "76%"),
        ("-1/+1 is YR", pct(yr), "57%"),
        ("-1C +1A", pct(pair["CA"]), "18%"),
        ("-1G +1C", pct(pair["GC"]), "0.6%"),
        ("inside a CDS", pct(intragenic), "63% (10,741)"),
        ("  of these, sense", pct(sense, intragenic), "about two-thirds"),
        ("leaderless (at CDS start)", "%s (%d)" % (pct(leaderless), leaderless), "0.4% (82)"),
        ("intragenic sense, codon pos 1", pct(codon[1], sum(codon.values())), "45%"),
        ("intragenic sense, codon pos 2", pct(codon[2], sum(codon.values())), "27%"),
        ("intragenic sense, codon pos 3", pct(codon[3], sum(codon.values())), "27%"),
    ]
    print("%-32s %14s   %s" % ("statistic", "this run", "Ettwiller et al. 2016"))
    for name, value, ref in rows:
        print("%-32s %14s   %s" % (name, value, ref))


if __name__ == "__main__":
    main(*sys.argv[1:4])
