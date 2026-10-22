# Audit2 (B, C): brute-force connected-cluster sums of log Phi_T(z), Phi_T(z) = sum_L N_L (z/n)^L, by direct
# enumeration of (run-length tuple, set partition of run positions), exact rational arithmetic.
#   (1/p!) sum over ORDERED tuples = sum over MULTISETS of 1/prod(mult_i!) (any fixed ordering);
#   weight mu(pi) = prod_B (-1)^(|B|-1)(|B|-1)!;  Z_pi = sum over labelings f: blocks -> [n] of prod_arcs S_(f f') (integer
#   S; divide by n^L at the end);  excess e = #arcs - #blocks; only connected quotient multigraphs; partitions with an
#   arc inside a block (consecutive positions of one run) have Z = 0 and are skipped (S_vv = 0).
# Options: regular=True prunes blocks of degree 1 (a lone run end; Z = 0 because S 1 = 0) and lets us bound the excess
# (sum_B max(0, deg_B - 2) = 2e once every block has degree >= 2); odd=False keeps only tuples whose runs all have
# even length (a tuple with an odd run sums to zero: reversing that run is a mu-preserving bijection on partitions
# that multiplies Z by (-1)^l; checked explicitly by the validation run).
# Z is computed by generic tensor contraction (numpy einsum, int64), no use of the THM-4757 core/kernel theory.
import sys, itertools
from fractions import Fraction
from math import factorial
from collections import Counter
import numpy as np

LET = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

def multisets(L, odd=True, minlen=1):
    """nondecreasing tuples of positive run lengths summing to L"""
    def rec(rem, lo):
        if rem == 0: yield (); return
        for l in range(lo, rem + 1):
            if not odd and l % 2: continue
            for rest in rec(rem - l, l): yield (l,) + rest
    yield from rec(L, minlen)

def partitions_pruned(pos_deg, prev, regular, emax, L):
    """set partitions of positions 0..P-1 (in order), skipping: block containing a position and its run predecessor
    prev[k]; if regular: final blocks of degree < 2, and (with emax) sum max(0, deg-2) > 2 emax or #blocks < L - emax."""
    P = len(pos_deg)
    blocks = []; bdeg = []; bof = [-1] * P
    budget_cap = None if emax is None else 2 * emax
    Vmin = None if emax is None else L - emax
    def over(d): return d - 2 if d > 2 else 0
    def rec(k, used):
        if k == P:
            if regular and any(d < 2 for d in bdeg): return
            yield [list(b) for b in blocks]; return
        rem_after = P - k - 1
        dq = pos_deg[k]
        # join existing block
        for b in range(len(blocks)):
            if prev[k] >= 0 and bof[prev[k]] == b: continue
            nd = bdeg[b] + dq
            nu = used - over(bdeg[b]) + over(nd)
            if regular and budget_cap is not None and nu > budget_cap: continue
            if regular:
                defic = sum(1 for i, d in enumerate(bdeg) if (d if i != b else nd) < 2)
                if defic > rem_after: continue
            if Vmin is not None and len(blocks) + rem_after < Vmin: continue
            blocks[b].append(k); bdeg[b] = nd; bof[k] = b
            yield from rec(k + 1, nu)
            blocks[b].pop(); bdeg[b] -= dq; bof[k] = -1
        # new block
        if regular:
            defic = sum(1 for d in bdeg if d < 2) + (1 if dq < 2 else 0)
            if defic > rem_after: return
        nu = used + over(dq)
        if Vmin is not None and len(blocks) + 1 + rem_after < Vmin: return
        blocks.append([k]); bdeg.append(dq); bof[k] = len(blocks) - 1
        yield from rec(k + 1, nu)
        blocks.pop(); bdeg.pop(); bof[k] = -1
    yield from rec(0, 0)

def cluster_sums(S, Lmax, regular=True, emax=None, odd=False, Lmin=1, verbose=True):
    """returns dict L -> dict e -> Fraction (coefficient of z^L in the excess-e part of log Phi_T)"""
    n = S.shape[0]
    S64 = np.array(S, dtype=np.int64)
    res = {}
    for L in range(Lmin, Lmax + 1):
        acc = {}; npart = 0
        for tup in multisets(L, odd=odd):
            mult = Counter(tup); wden = 1
            for c in mult.values(): wden *= factorial(c)
            pos = []; prev = []; deg = []
            for j, l in enumerate(tup):
                for i in range(l + 1):
                    prev.append(len(pos) - 1 if i > 0 else -1)
                    deg.append(1 if i in (0, l) else 2)
                    pos.append((j, i))
            arcs_pos = [(k - 1, k) for k in range(len(pos)) if prev[k] >= 0]
            for blocks in partitions_pruned(deg, prev, regular, emax, L):
                bid = [0] * len(pos)
                for t, b in enumerate(blocks):
                    for q in b: bid[q] = t
                arcs = [(bid[a], bid[b]) for a, b in arcs_pos]
                V = len(blocks)
                par = list(range(V))
                def f(a):
                    while par[a] != a: par[a] = par[par[a]]; a = par[a]
                    return a
                for a, b in arcs: par[f(a)] = f(b)
                if len({f(a) for a in range(V)}) != 1: continue
                e = len(arcs) - V
                mu = 1
                for b in blocks: mu *= (-1) ** (len(b) - 1) * factorial(len(b) - 1)
                spec = ','.join(LET[a] + LET[b] for a, b in arcs) + '->'
                Z = int(np.einsum(spec, *([S64] * len(arcs)), optimize='greedy'))
                acc[e] = acc.get(e, Fraction(0)) + Fraction(mu * Z, wden)
                npart += 1
        res[L] = {e: v / n ** L for e, v in acc.items()}
        if verbose:
            print(f"   L={L}: {npart} connected partitions; " +
                  ', '.join(f"e={e}: {float(v):+.10e}" for e, v in sorted(res[L].items())), flush=True)
    return res
