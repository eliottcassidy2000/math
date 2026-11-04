#!/usr/bin/env python3
"""fo_rule_modelcheck.py -- HYP-9463 (4), pre-registered test (mac-mini-2026-11-04-oddprimes).

Uniform quantifier-free rules with at most two atoms: phi(x) = B(S(c1 x + d1), S(c2 x + d2)) with fixed integers
c_i in {+-1, ..., +-CMAX}, d_i in {-DMAX, ..., DMAX}, B any Boolean function of two bits (16; the 4 that ignore an
atom cover the one-atom rules).  For each prime p the rule maps a skew set S (circulant tournament on Z_p) to
C_phi(S) = {x in Z_p : phi(x)}.  The rule PASSES at p if (i) C_phi(S) is a symmetric set (C = -C, 0 notin C) for every S,
(ii) the dilation class of C_phi(S) depends only on the dilation class of S, and (iii) the induced map is a bijection from
the t_p dilation classes of skew sets onto the e_p classes of multiplier-even symmetric sets (sign convention of THM-4721).
Reports the number of passing rules per prime and the rules passing at every tested prime.
"""
import sys, itertools, math, time
import numpy as np

def reversal_parity_mult(p, a, C):
    """sign of the dilation x -> a x on the circulant graph (Z_p, C): parity of #{edges {u<v}: au > av mod p}."""
    cnt = 0
    for u in range(p):
        for c in C:
            v = (u + c) % p
            if u < v and (a * u) % p > (a * v) % p: cnt += 1
    return cnt & 1

def setup(p):
    units = [a for a in range(1, p) if math.gcd(a, p) == 1]
    half = list(range(1, (p + 1) // 2))
    # skew sets as bitmasks over 0..p-1
    skew = []
    for bits in itertools.product([0, 1], repeat=len(half)):
        S = 0
        for d, b in zip(half, bits): S |= 1 << (d if b == 0 else p - d)
        skew.append(S)
    skew_arr = np.array([[(S >> x) & 1 for x in range(p)] for S in skew], dtype=np.uint8)
    def dil(mask, a):
        return sum(1 << ((a * x) % p) for x in range(p) if (mask >> x) & 1)
    # skew dilation classes
    cls_skew = {}; nskew = 0
    for S in skew:
        if S in cls_skew: continue
        for a in units: cls_skew[dil(S, a)] = nskew
        nskew += 1
    skew_cls = np.array([cls_skew[S] for S in skew])
    # symmetric sets, their classes, multiplier-evenness
    sym = []
    for bits in itertools.product([0, 1], repeat=len(half)):
        C = 0
        for d, b in zip(half, bits):
            if b: C |= (1 << d) | (1 << (p - d))
        sym.append(C)
    cls_sym = {}; even_cls = set(); nsym = 0
    for C in sym:
        if C in cls_sym: continue
        Cset = [x for x in range(p) if (C >> x) & 1]
        stab = [a for a in units if dil(C, a) == C]
        ev = all(reversal_parity_mult(p, a, Cset) == 0 for a in stab)
        for a in units: cls_sym[dil(C, a)] = nsym
        if ev: even_cls.add(nsym)
        nsym += 1
    pow2 = np.array([1 << x for x in range(p)], dtype=np.int64)
    return dict(p=p, units=units, skew_arr=skew_arr, skew_cls=skew_cls, nskew=nskew, cls_sym=cls_sym,
                even_cls=even_cls, pow2=pow2)

BOOL = [(i, lambda a, b, i=i: ((i >> (2 * a + b)) & 1)) for i in range(16)]   # truth table index i

def test_rule(D, c1, d1, c2, d2, i):
    p = D["p"]; A = D["skew_arr"]
    perm1 = [((c1 * x + d1) % p) for x in range(p)]; perm2 = [((c2 * x + d2) % p) for x in range(p)]
    a1 = A[:, perm1]; a2 = A[:, perm2]
    idx = 2 * a1.astype(np.int64) + a2.astype(np.int64)        # 0..3
    tt = np.array([(i >> j) & 1 for j in range(4)], dtype=np.uint8)
    out = tt[idx]                                               # (nS, p) bits of C_phi(S)
    # symmetric and 0 excluded
    if out[:, 0].any(): return False
    neg = [(-x) % p for x in range(p)]
    if not np.array_equal(out, out[:, neg]): return False
    masks = (out.astype(np.int64) * D["pow2"]).sum(axis=1)
    cls = np.array([D["cls_sym"][int(m)] for m in masks])
    # well-defined on skew classes
    for k in range(D["nskew"]):
        sel = cls[D["skew_cls"] == k]
        if sel.min() != sel.max(): return False
    img = set(cls[np.unique(D["skew_cls"], return_index=True)[1]].tolist())
    return len(img) == D["nskew"] and img == D["even_cls"]

if __name__ == "__main__":
    args = sys.argv[1:]
    CMAX = DMAX = 4
    if args and args[0] == "--box": CMAX = DMAX = int(args[1]); args = args[2:]
    primes = [int(a) for a in args] or [5, 7, 11, 13, 17]
    cs = [c for c in range(-CMAX, CMAX + 1) if c != 0]; ds = list(range(-DMAX, DMAX + 1))
    survivors = None
    for p in primes:
        t0 = time.time(); D = setup(p); passing = []
        atoms = [(c, d) for c in cs for d in ds]
        for (c1, d1), (c2, d2) in itertools.combinations_with_replacement(atoms, 2):
            for i in range(16):
                if test_rule(D, c1, d1, c2, d2, i): passing.append((c1, d1, c2, d2, i))
        print("p=%d: t_p=%d e_p=%d rules tested=%d passing=%d (%.0fs)" % (p, D["nskew"], len(D["even_cls"]),
              16 * len(atoms) * (len(atoms) + 1) // 2, len(passing), time.time() - t0), flush=True)
        if passing[:8]: print("   examples:", passing[:8])
        S = set(passing); survivors = S if survivors is None else survivors & S
        print("   rules passing at every prime so far:", len(survivors), sorted(survivors)[:6], flush=True)
