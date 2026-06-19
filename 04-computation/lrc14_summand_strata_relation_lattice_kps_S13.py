#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) ENDGAME — THE SUMMAND-GRAPH STRATA mod C=27=3^3 and the signed cancellation.
kind-pasteur-2026-06-19-S13.

Goal (per task ANGLE):
  For the Freiman dimension-d GAPs and near-AP sets, compute the relation lattice
    Lambda(E) = { m in Z^k : sum_e m_e e = 0 },  E = {e_0=0, e_1, ..., e_{k-1}}.
  Reduce the relations mod C = 2n-1 = 27 = 3^3 (THM-401).  Shells {a, 27-a} = {+a,-a}.
  27 = 3^3 stratifies the residues by gcd(a,27) into: UNIT (gcd 1), gcd-3, gcd-9.
  Classify each relation/generator by:
    (i)   which stratum its support residues live in,
    (ii)  its sign pattern (pos/neg), antipodal +a/-a pairing within a stratum.
  KEY: does the AP's relation lattice live in a DIFFERENT stratum than the GAP's?
       Does the gcd-3 stratum (where V* lives) carry the dangerous near-AP corrections,
       and the unit stratum the harmless ones?

This BUILDS on (does not re-derive): THM-534 L_y duality, caps, HYP-2637 dimension penalty,
THM-401 (modulus 2n-1), HYP-2083 (gcd-stratum, V* in gcd-3 blind spot), HYP-2632/2633 codex
finite character kernel + residue-lift equidistribution.

Engine reused: lrc14_empty_sector_distribution_kps.dist_p / L_y.
"""
import sys, itertools
from fractions import Fraction
from math import gcd
from functools import reduce

sys.path.insert(0, "04-computation")
try:
    from lrc14_empty_sector_distribution_kps import dist_p, g_poly, L_y, consec
except Exception:
    # fallback: inline a minimal copy if import path differs
    def dist_p(E):
        E = sorted(set(E))
        bps = set([Fraction(0), Fraction(1)])
        for e in E:
            if e == 0: continue
            for a in range(0, 7*e+1):
                bps.add(Fraction(a, 7*e))
        bps = sorted(b for b in bps if 0 <= b <= 1)
        p = [Fraction(0)]*7
        for i in range(len(bps)-1):
            lo, hi = bps[i], bps[i+1]
            if hi == lo: continue
            mid = (lo+hi)/2
            hit = set()
            for e in E:
                v = e*mid; v = v - (v.numerator//v.denominator)
                hit.add((v.numerator*7)//v.denominator)
            t = sum(1 for j in range(1,7) if j not in hit)
            p[t] += (hi-lo)
        return p
    def g_poly(k):
        g=[]
        for t in range(7):
            if k==8: val=Fraction((t-1)*(t-2)*(t-4)*(t-5),40)
            elif k in (9,10): val=Fraction(-(t-2)*(t-3)*(t-6),36)
            else: val=Fraction((t-3)*(t-4),12)
            g.append(val)
        return g
    def L_y(E,k):
        p=dist_p(E); g=g_poly(k)
        return sum(p[t]*g[t] for t in range(7)), p
    def consec(k): return list(range(k))

C = 27  # = 2*14 - 1 = 3^3

# ---------------------------------------------------------------------------
# Stratum machinery mod C = 27 = 3^3.
# A nonzero residue a in Z/27 has gcd(a,27) in {1,3,9}.
# UNIT stratum  : gcd=1  (18 residues, (Z/27)^* )
# gcd-3 stratum : gcd=3  (residues 3,6,12,15,21,24 ... i.e. 3*unit-mod-9)  -> 6 residues
# gcd-9 stratum : gcd=9  (residues 9, 18)                                  -> 2 residues
# residue 0 is its own (the trivial shell).
def stratum(a):
    a %= C
    if a == 0:
        return "zero"
    g = gcd(a, C)
    if g == 1: return "unit"
    if g == 3: return "gcd3"
    if g == 9: return "gcd9"
    return "other"

def shell(a):
    """antipodal shell {a, C-a} as a sorted tuple (the +a/-a pair)."""
    a %= C
    b = (-a) % C
    return tuple(sorted({a, b}))

# ---------------------------------------------------------------------------
# Relation lattice Lambda(E) = { m in Z^k : sum_e m_e e = 0 }.
# This is the kernel of the 1xk integer matrix [e_0,...,e_{k-1}] (e_0=0).
# rank of E (the row) is 1 (assuming some e>0), so the lattice has rank k-1.
# We compute an integer basis via the Smith/HNF kernel of the row vector.
def integer_kernel_basis(vec):
    """Integer basis (list of rows) of { m in Z^k : vec . m = 0 }. vec has integer entries."""
    k = len(vec)
    # The kernel of a single linear form sum a_i m_i = 0.
    # Standard basis: for each pair, but cleanest is: kernel of gcd-reduced form.
    # Use the construction via the cofactor / lattice of integer solutions.
    # Reduce: let d = gcd of all a_i. We want all integer m with sum a_i m_i = 0.
    # Build using a unimodular transform: find U (kxk, det +-1) s.t. vec*U = (g,0,...,0).
    a = list(vec)
    # Column HNF of the row: track transform.
    import copy
    U = [[1 if i==j else 0 for j in range(k)] for i in range(k)]  # k x k identity, columns ops
    row = a[:]
    # We perform column operations to zero out all but one entry. U accumulates column ops.
    # row * U_total = reduced. Kernel = columns of U_total that map to zero columns.
    # Do extended-gcd style elimination.
    cols = list(range(k))
    # repeatedly: find pivot = nonzero entry of min abs, reduce others mod it.
    def colop_add(dst, src, mult):
        # column dst += mult * column src   (operates on row entries and U)
        row[dst] += mult * row[src]
        for i in range(k):
            U[i][dst] += mult * U[i][src]
    def colswap(i,j):
        row[i],row[j]=row[j],row[i]
        for r in range(k):
            U[r][i],U[r][j]=U[r][j],U[r][i]
    # Gaussian-like over integers on the single row
    piv = 0
    # bring a nonzero to position 0
    nz = [i for i in range(k) if row[i]!=0]
    if not nz:
        # all zero -> whole space is kernel
        return [[1 if i==j else 0 for j in range(k)] for i in range(k)]
    # reduce to single pivot via gcd
    while True:
        nz = [i for i in range(k) if row[i]!=0]
        if len(nz) <= 1:
            break
        # pick two smallest-abs nonzero
        nz.sort(key=lambda i: abs(row[i]))
        i, j = nz[0], nz[1]
        q = row[j] // row[i]
        colop_add(j, i, -q)
    # now exactly one nonzero column (the pivot). All other columns of U are kernel vectors.
    nz = [i for i in range(k) if row[i]!=0]
    pivcol = nz[0] if nz else 0
    basis = []
    for j in range(k):
        if j == pivcol: continue
        basis.append([U[i][j] for i in range(k)])
    return basis

def relation_reduced(m, E):
    """Given relation m (sum m_e e_e = 0 over Z), its 'shell content':
       the multiset of (residue e_e mod C, coefficient m_e) for e_e != 0, m_e != 0,
       and the stratum classification of the *support residues touched*."""
    touched = {}  # residue -> coefficient
    for me, e in zip(m, E):
        if me == 0: continue
        r = e % C
        touched[r] = touched.get(r, 0) + me
    return touched

def lattice_stratum_profile(E):
    """For the relation lattice of E, classify support residues by stratum, and
       report the sign/antipodal structure. Returns a dict summary."""
    E = sorted(set(E))
    k = len(E)
    vec = E[:]  # the linear form coefficients are the e-values themselves
    basis = integer_kernel_basis(vec)
    # For EACH relation in a reduced LLL-ish basis, record which strata the
    # *nonzero-coefficient support residues* fall in.
    strat_count = {"unit":0,"gcd3":0,"gcd9":0,"zero":0}
    relation_strata = []
    for m in basis:
        residues_with_coef = [(E[i] % C, m[i]) for i in range(k) if m[i] != 0 and E[i] != 0]
        strata = set(stratum(r) for r,_ in residues_with_coef)
        relation_strata.append((tuple(m), residues_with_coef, sorted(strata)))
        for s in strata:
            strat_count[s] = strat_count.get(s,0)+1
    # Also: the FULL set of residues of E by stratum (the "support strata of E")
    E_residues = [e % C for e in E if e != 0]
    E_strat = {"unit":[], "gcd3":[], "gcd9":[], "zero":[]}
    for r in E_residues:
        E_strat[stratum(r)].append(r)
    return {
        "k": k, "E": E, "basis": basis,
        "E_residue_strata": E_strat,
        "relation_strata": relation_strata,
        "strat_count": strat_count,
    }

# ---------------------------------------------------------------------------
# Low-height SHELL relations: enumerate small-coefficient relations sum m_e e = 0 mod C
# (THM-401: relations reduce mod C=27). These are the actual "ghost" terms in the
# correction L_y - L_y^inf. Classify each by stratum + antipodal sign pattern.
def low_height_shell_relations(E, max_coef=2, mod=C):
    """Enumerate relations sum_e m_e * (e mod C) = 0 (mod C), |m_e|<=max_coef, not all zero,
       over the NONZERO elements of E. Classify by stratum and pos/neg shell structure.
       Returns list of (m_vector, touched_residues, strata_set, signpattern)."""
    Ez = [e for e in sorted(set(E)) if e != 0]
    res = [e % mod for e in Ez]
    out = []
    rng = range(-max_coef, max_coef+1)
    for m in itertools.product(rng, repeat=len(Ez)):
        if all(x==0 for x in m): continue
        s = sum(mi*ri for mi,ri in zip(m,res)) % mod
        if s != 0: continue
        touched = [(res[i], m[i]) for i in range(len(Ez)) if m[i]!=0]
        strata = sorted(set(stratum(r) for r,_ in touched))
        # antipodal sign pattern: does it pair +a with -a (same shell, opposite or same sign)?
        shells_used = {}
        for r,coef in touched:
            sh = shell(r)
            shells_used.setdefault(sh, []).append((r,coef))
        # "balanced" = uses an antipodal pair {a, C-a} (true +/- cancellation possible)
        antipodal = any(len(v) >= 2 and len({r for r,_ in v})>=2 for v in shells_used.values())
        out.append((m, touched, strata, antipodal))
    return out

def stratum_histogram(relations):
    """count relations whose support touches purely unit / purely gcd3 / purely gcd9 / mixed."""
    h = {"pure_unit":0,"pure_gcd3":0,"pure_gcd9":0,"mixed":0,"has_gcd3":0,"has_unit":0,"antipodal":0}
    for m, touched, strata, antipodal in relations:
        ss = set(strata)
        if ss == {"unit"}: h["pure_unit"]+=1
        elif ss == {"gcd3"}: h["pure_gcd3"]+=1
        elif ss == {"gcd9"}: h["pure_gcd9"]+=1
        else: h["mixed"]+=1
        if "gcd3" in ss: h["has_gcd3"]+=1
        if "unit" in ss: h["has_unit"]+=1
        if antipodal: h["antipodal"]+=1
    h["total"]=len(relations)
    return h

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("="*78)
    print("LRC(14) SUMMAND-GRAPH STRATA mod C=27=3^3  —  kind-pasteur-S13")
    print("="*78)
    print(f"C = 2n-1 = {C} = 3^3.  Strata by gcd(a,27): unit(gcd1), gcd3, gcd9.")
    print(f"  unit residues  (18): {[a for a in range(1,C) if gcd(a,C)==1]}")
    print(f"  gcd-3 residues (6) : {[a for a in range(1,C) if gcd(a,C)==3]}")
    print(f"  gcd-9 residues (2) : {[a for a in range(1,C) if gcd(a,C)==9]}")
    print()

    # ---- Test families per cluster size k ----
    # AP (consec) vs near-AP (one hole/shift) vs d=2 GAP (Minkowski sum of 2 APs).
    families = {
        8: {
            "AP_consec":        [0,1,2,3,4,5,6,7],
            "nearAP_top_max":   [0,2,3,4,5,6,7,8],   # the verified max non-AP at k=8
            "nearAP_hole":      [0,1,2,3,4,5,6,8],
            "d2_GAP_2x4":       [a+3*b for a in range(4) for b in range(2)],  # {0,1,2,3}+{0,3}
            "d2_GAP_4x2":       [a+4*b for a in range(2) for b in range(4)],  # {0,1}+{0,4,8,12}? trim to 8
        },
        9: {
            "AP_consec":        list(range(9)),
            "nearAP_top_max":   [0,1,2,3,4,5,6,7,9],  # verified k=9 max non-AP (tightest, margin 0.0070)
            "nearAP_hole":      [0,1,2,3,4,5,6,8,9],
            "d2_GAP_3x3":       sorted(set(a+3*b for a in range(3) for b in range(3))),  # {0,1,2}+{0,3,6}=consec! degenerate
            "d2_GAP_3x4":       sorted(set(a+4*b for a in range(3) for b in range(3)))[:9],
        },
        10: {
            "AP_consec":        list(range(10)),
            "nearAP_top_max":   [0,1,2,3,4,5,6,7,8,10],  # verified k=10 max non-AP
            "nearAP_hole":      [0,1,2,3,4,5,6,7,9,10],
        },
    }
    # The tight sporadic V* (HYP-2083) at n=14, |V*|=13 — analyze its summand structure too.
    Vstar = [1,2,3,4,5,6,7,8,9,10,11,13,24]  # NOTE these are speeds, residues mod 27

    for k in [8, 9, 10]:
        g = g_poly(k)
        print("#"*78)
        print(f"## CLUSTER SIZE k={k}")
        print("#"*78)
        for name, E in families[k].items():
            E = sorted(set(E))
            if len(E) != k:
                print(f"  [skip {name}: |E|={len(E)} != {k}  E={E}]")
                continue
            if 0 not in E:
                print(f"  [skip {name}: 0 not in E]")
                continue
            if reduce(gcd, E) != 1:
                # primitive only (scale-invariance THM-531 means non-primitive == its primitive)
                pass
            Lval, p = L_y(E, k)
            excess = len(set(a+b for a in E for b in E)) - (2*k-1)
            prof = lattice_stratum_profile(E)
            relsH2 = low_height_shell_relations(E, max_coef=2)
            hist = stratum_histogram(relsH2)
            print(f"\n  --- {name}: E={E}")
            print(f"      L_y = {float(Lval):.5f}   excess(|E+E|-(2k-1)) = {excess}  (excess=0 <=> full AP)")
            es = prof["E_residue_strata"]
            print(f"      E residues mod 27 by stratum: unit={es['unit']}  gcd3={es['gcd3']}  gcd9={es['gcd9']}")
            print(f"      low-height (|m|<=2) shell relations mod 27: total={hist['total']}")
            print(f"        pure-unit={hist['pure_unit']}  pure-gcd3={hist['pure_gcd3']}  pure-gcd9={hist['pure_gcd9']}  mixed={hist['mixed']}")
            print(f"        touching gcd3={hist['has_gcd3']}  touching unit={hist['has_unit']}  antipodal(+a/-a)={hist['antipodal']}")
        print()

    # ---- V* summand-graph analysis ----
    print("#"*78)
    print("## V* (HYP-2083 tight sporadic, n=14): residues mod 27 and stratum content")
    print("#"*78)
    print(f"  V* = {Vstar}  (speeds = residues mod 27)")
    es = {"unit":[],"gcd3":[],"gcd9":[],"zero":[]}
    for v in Vstar:
        es[stratum(v)].append(v % C)
    print(f"  by stratum: unit={es['unit']}  gcd3={es['gcd3']}  gcd9={es['gcd9']}")
    # HYP-2083: V* misses {12,15} and doubles {3,24}, both gcd-3 stratum. Verify.
    print(f"  gcd-3 residues present in V*: {sorted(set(v%C for v in Vstar if stratum(v)=='gcd3'))}")
    print(f"  gcd-3 residues ABSENT (missed): {sorted(set(a for a in range(1,C) if gcd(a,C)==3) - set(v%C for v in Vstar))}")
    relsV = low_height_shell_relations(Vstar, max_coef=2)
    histV = stratum_histogram(relsV)
    print(f"  low-height shell relations: total={histV['total']} pure-unit={histV['pure_unit']} "
          f"pure-gcd3={histV['pure_gcd3']} mixed={histV['mixed']} antipodal={histV['antipodal']}")
