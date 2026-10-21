"""Label-free TSS model prototype: labels come from the Cappable-seq data itself.

1. Candidates: positions with >= MIN_READS read 5' ends in the enriched library.
2. Self-labels from the chemistry:
   positive = strongly enriched (log2 enrichment >= POS_FC), local peak, seen in
              the technical replicate
   negative = depleted processed ends (log2 enrichment <= NEG_FC) plus random
              genomic positions (sequence background)
3. Sequence model: logistic regression on one-hot sequence -45..+5, trained on
   one half of the genome (blocks of 100 kb) and applied to the other half.
4. Evaluation only: Ettwiller et al. 2016 Table S1 (never used for training).
"""
import json
import math
import os
import random
import sys

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import cappable_tss as ct  # noqa: E402

# usage: selftrain.py enriched.bam replicate.bam control.bam genome.fa reference_TSS.gtf out_weights.json
ENR, REP, CTL, FASTA, REF, OUT = sys.argv[1:7]
MIN_READS, POS_FC, NEG_FC, BLOCK = 3, 3.0, -2.0, 100000
UP, DOWN = 45, 5

enr, rep, ctl = (ct.five_prime_from_bam(p) for p in (ENR, REP, CTL))
seq = "".join(l.strip() for l in open(FASTA) if not l.startswith(">")).upper()
COMP = str.maketrans("ACGT", "TGCA")
L = len(seq)


def window(st, p):
    if st == "+":
        w = seq[p - 1 - UP:p - 1 + DOWN]
    else:
        w = seq[p - DOWN:p + UP][::-1].translate(COMP)
    return w if len(w) == UP + DOWN and set(w) <= set("ACGT") else None


def onehot(w):
    x = np.zeros((len(w), 4), dtype=np.float32)
    x[np.arange(len(w)), ["ACGT".index(b) for b in w]] = 1
    return x.ravel()


tot = {name: sum(sum(d.values()) for d in sig.values()) for name, sig in
       (("e", enr), ("r", rep), ("c", ctl))}
cands = []
for (sid, st), d in enr.items():
    r, c = rep.get((sid, st), {}), ctl.get((sid, st), {})
    for p, n in d.items():
        if n < MIN_READS or not (UP + 1 < p < L - UP - 1):
            continue
        fc = math.log2((n / tot["e"]) / ((c.get(p, 0) + 1) / tot["c"]))
        peak = all(d.get(q, 0) <= n for q in range(p - 5, p + 6))
        cands.append((st, p, n, fc, peak, r.get(p, 0)))

pos = [(st, p) for st, p, n, fc, peak, rr in cands if fc >= POS_FC and peak and rr >= 1]
neg_dep = [(st, p) for st, p, n, fc, peak, rr in cands if fc <= NEG_FC]
random.seed(0)
neg_rand = [(random.choice("+-"), random.randint(UP + 2, L - UP - 2)) for _ in range(3 * len(pos))]
print("candidates %d | self-labelled positives %d | depleted negatives %d | random negatives %d"
      % (len(cands), len(pos), len(neg_dep), len(neg_rand)))

ref = {(st, p) for _, st, p in ct.read_sites(REF)}
ref_near = {(st, p + d) for st, p in ref for d in range(-1, 2)}


def half(p):
    return (p // BLOCK) % 2


scores = {}
for test_half in (0, 1):
    X, y = [], []
    for label, group in ((1, pos), (0, neg_dep + neg_rand)):
        for st, p in group:
            if half(p) == test_half:
                continue
            w = window(st, p)
            if w:
                X.append(onehot(w))
                y.append(label)
    model = LogisticRegression(C=0.05, max_iter=2000).fit(np.array(X), np.array(y))
    test = [(st, p) for st, p, *_ in cands if half(p) == test_half]
    Xt = [onehot(w) for w in (window(st, p) for st, p in test) if w]
    keep = [(st, p) for st, p in test if window(st, p)]
    for k, s in zip(keep, model.decision_function(np.array(Xt))):
        scores[k] = s

# Evaluation against Ettwiller (not used above): among candidate positions that
# are local peaks, how well does sequence alone separate published TSS?
peaks = [(st, p, n, fc) for st, p, n, fc, peak, rr in cands if peak and (st, p) in scores]
for lo, hi, label in ((20, 10 ** 9, ">= 20 reads (Ettwiller's range)"), (3, 20, "3-19 reads (below Ettwiller's cut-off)")):
    grp = [(st, p, n, fc) for st, p, n, fc in peaks if lo <= n < hi]
    if lo >= 20:
        y = [int((st, p) in ref_near) for st, p, n, fc in grp]
        s = [scores[(st, p)] for st, p, n, fc in grp]
        auc_seq = roc_auc_score(y, s)
        auc_fc = roc_auc_score(y, [fc for *_, fc in grp])
        auc_both = roc_auc_score(y, [scores[(st, p)] + fc for st, p, n, fc in grp])
        print("%s: %d peaks, %d in Ettwiller | AUROC sequence %.3f, enrichment %.3f, sum %.3f"
              % (label, len(grp), sum(y), auc_seq, auc_fc, auc_both))
    else:
        enriched = [(st, p) for st, p, n, fc in grp if fc >= 1]
        hi_seq = [k for k in enriched if scores[k] > 0]
        print("%s: %d enriched peaks (log2 >= 1); %d (%.0f%%) also have a positive sequence score"
              % (label, len(enriched), len(hi_seq), 100.0 * len(hi_seq) / max(1, len(enriched))))

coef = LogisticRegression(C=0.05, max_iter=2000)
X = [onehot(window(st, p)) for st, p in pos + neg_dep + neg_rand if window(st, p)]
y = [1] * sum(1 for st, p in pos if window(st, p)) + [0] * sum(1 for st, p in neg_dep + neg_rand if window(st, p))
w = coef.fit(np.array(X), np.array(y)).coef_.reshape(UP + DOWN, 4)
json.dump({"positions": list(range(-UP, 0)) + list(range(1, DOWN + 1)), "bases": "ACGT",
           "weights": w.tolist()}, open(OUT, "w"))
print("learned sequence weights, strongest base per position (-45..+5):")
for i in range(UP + DOWN):
    j = int(np.argmax(np.abs(w[i])))
    if abs(w[i, j]) > 0.25:
        print("  %+d  %s %+.2f" % (i - UP if i < UP else i - UP + 1, "ACGT"[j], w[i, j]))
