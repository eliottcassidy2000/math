#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_pattern_linearity_kps.py   (kind-pasteur 2026-06-21, HYP-2724)

Independent test of the KEY structural claim behind the MDS/arc route:
  IF K(n)=D7(n mod 7)/prod n_j depends only on the relation's COEFFICIENT PATTERN
  (not on which elements it links), THEN corr(E)=Sum_{n in Lambda(E)} K(n) is a
  FIXED LINEAR COMBINATION of relation-PATTERN COUNTS -> the wide bound is purely
  COMBINATORIAL (count relations by pattern; AP maximizes the dominant counts).

We test this DATA-DRIVEN (no kernel reconstruction): over a battery of k-sets,
regress the exact corr(E) against the vector of pattern counts. A near-exact fit
(R^2 ~ 1) with stable coefficients == corr is determined by pattern counts.
EXACT measS7; least-squares over rationals->float.
"""
import itertools, random
from fractions import Fraction as Fr
from math import comb, gcd
from functools import reduce

# reuse exact measure machinery
import importlib.util, os
spec = importlib.util.spec_from_file_location(
    "rcm", os.path.join(os.path.dirname(__file__), "lrc_q108_relation_code_mds_kps.py"))
rcm = importlib.util.module_from_spec(spec); spec.loader.exec_module(rcm)
measS7, iid_k = rcm.measS7, rcm.iid_k

def corr(E, p=7):
    return float(measS7(E, p) - iid_k(len(E), p))

def canon_pattern(coefs):
    """canonical coefficient pattern: primitive, sorted by (abs,sign), global sign fixed."""
    g = reduce(gcd, [abs(c) for c in coefs]); coefs = tuple(c//g for c in coefs)
    # fix global sign so the sorted form is canonical under negation
    s = tuple(sorted(coefs)); sn = tuple(sorted(-c for c in coefs))
    return min(s, sn)

def pattern_counts(E, B=2, max_support=4):
    """count primitive relations by (support, canonical-pattern)."""
    E = [int(e) for e in E]; k = len(E)
    counts = {}
    for s in range(2, max_support+1):
        for combo in itertools.combinations(range(k), s):
            for coefs in itertools.product(range(-B, B+1), repeat=s):
                if any(c == 0 for c in coefs):
                    continue
                if sum(c*E[i] for c, i in zip(coefs, combo)) != 0:
                    continue
                pat = canon_pattern(coefs)
                key = (s, pat)
                # dedupe: same (combo, primitive coefs up to sign)
                counts[key] = counts.get(key, 0) + 1
    # each relation counted twice (n and -n) -> halve
    return {k_: v//2 if v % 2 == 0 else v for k_, v in counts.items()}

def main():
    random.seed(7)
    k = 8
    # battery: structured + random, bounded range so measS7 is cheap & exact
    sets = []
    sets.append(tuple(range(8)))                      # AP
    sets.append((0,2,4,6,8,10,12,14))                 # AP step2
    sets.append((0,1,2,4,8,16,32,64))                 # dyadic (cap range)
    sets.append((0,1,3,7,12,20,30,44))                # Sidon
    sets.append((0,1,2,3,40,41,42,43))                # two-block (cap)
    # random k-subsets of [0,24]
    seen=set(sets)
    while len(sets) < 140:
        E = tuple(sorted(random.sample(range(0,25), k)))
        if E[0]!=0:  # pin 0 (observer phase), like the LRC offset convention
            E = tuple(sorted(set([0])|set(E[1:])))
            if len(E)!=k: continue
        if E in seen: continue
        seen.add(E); sets.append(E)

    # collect all patterns that appear, build design matrix
    allpat = {}
    data = []
    for E in sets:
        pc = pattern_counts(E, B=2, max_support=4)
        for key in pc:
            allpat.setdefault(key, len(allpat))
        data.append((E, pc, corr(E)))
    pats = sorted(allpat, key=lambda x:(x[0], x[1]))
    idx = {pat:i for i,pat in enumerate(pats)}
    # design matrix X (rows=sets, cols=patterns), y=corr
    X = [[row[1].get(p,0) for p in pats] for row in data]
    y = [row[2] for row in data]

    # least squares via normal equations (numpy if available, else pure)
    try:
        import numpy as np
        A = np.array(X, float); b = np.array(y, float)
        coef, res, rank, sv = np.linalg.lstsq(A, b, rcond=None)
        pred = A.dot(coef)
        ss_res = float(((b-pred)**2).sum()); ss_tot = float(((b-b.mean())**2).sum())
        r2 = 1 - ss_res/ss_tot
        maxabs = float(np.abs(b-pred).max())
        print(f"#sets={len(sets)}  #patterns={len(pats)}  design rank={rank}")
        print(f"LINEAR FIT corr ~ pattern-counts:  R^2 = {r2:.6f}   max|residual| = {maxabs:.5f}")
        print(f"\nFitted coefficient per pattern (=predicted K(pattern)); |coef|>1e-4:")
        order = sorted(range(len(pats)), key=lambda i:-abs(coef[i]))
        for i in order:
            if abs(coef[i])>1e-4:
                s,pat = pats[i]
                print(f"   support {s}  pattern {str(pat):<18}  K_fit = {coef[i]:+.5f}")
    except ImportError:
        print("numpy not available; reporting pattern inventory only")
        for p in pats: print("  ", p)

    # the leading patterns and AP-maximality of the dominant count
    print("\nDOMINANT support-3 counts across sets (does AP top them?):")
    schur=(3,(-1,1,1)); ap=(3,(-2,1,1)) if (3,(-2,1,1)) in idx else None
    for name,E in [("AP{0..7}",tuple(range(8))),("Sidon",(0,1,3,7,12,20,30,44)),
                   ("random sample",sets[20])]:
        pc=pattern_counts(E,2,4)
        tot3=sum(v for (s,_),v in pc.items() if s==3)
        print(f"   {name:<16} corr={corr(E):+.4f}  total support-3 = {tot3}")
    print("\nINTERPRETATION: R^2~1 ==> corr is (to truncation) a FIXED linear form in")
    print("pattern counts ==> wide bound is COMBINATORIAL; AP maximizes the dominant counts.")
    print("DONE.")

if __name__ == "__main__":
    main()
