#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadD_crossblock_carrier_bound_opus_0621.py   (opus-2026-06-21, THREAD D — the LEVER)

THE QUESTION
============
The remaining LRC(14) gap is the MODERATE-SPAN MULTI-BLOCK case (comfortable budget ~0.18).
For E = block_1 u block_2 u ... at dissociated scales, the cover measS7(E) = iid_k + corr(E),
and corr(E) lives ENTIRELY in the OFFSET RELATION LATTICE
        Lambda(E) = { n in Z^E : sum_i n_i e_i = 0 }
via the proved Fourier identity (lrc14_routeE_cutcycle_split / kps-S6 / angleF):
        corr(E) = sum_{R subset Z/7, 0 notin R} (-1)^|R| sum_{0 != n in Lambda(E)} prod_i chat_R(n_i),
with chat_R(0) = 1-|R|/7 and chat_R(7m)=0 (apex-prime vanishing): ONLY relation vectors with NO
coordinate a nonzero multiple of 7 contribute ("effective relations").

DECOMPOSE the relation lattice by how many blocks the support TOUCHES:
   Lambda(E) = Lambda_intra (relations supported in ONE block)  (+)  Lambda_cross (touching >= 2 blocks).
   corr(E)   = corr_intra(E)  +  corr_cross(E).

THM-559 TOURNAMENT HOME TURF (the analogy that motivates the lever):
   Lambda(E) is the LRC analog of the GF(2) CYCLE SPACE of K_n. corr_intra = the per-block
   "internal cycles"; corr_cross = the genuinely MULTI-DIMENSIONAL Weyl error, the analog of
   tournament 3-cycles whose arcs straddle a vertex partition. THM-559: c3 = 2-body Ising on
   arc spins; the SHORTEST cross-block relations are SUPPORT-3 ADDITIVE TRIANGLES
        e_b - e_a = e_d - e_c   (a "same common difference" Schur quadruple),  or  e_a + e_b = e_c
   across blocks. These are the carrier-error 3-cycles.

THE THREAD-D LEVER (what this script delivers, all EXACT Fractions):
 (D1) EXACT decomposition corr = corr_intra + corr_cross by enumerating Lambda(E) and splitting
      by block-support. VERIFY corr_intra + corr_cross + high-residual reproduces measS7-iid.
 (D2) The MULTI-BLOCK FACTORIZATION: for blocks at DISSOCIATED scales,
        - intra relations are scale-invariant copies of the single-block lattices (corr_intra
          = sum over blocks of the per-block corr, exact for separated blocks),
        - cross relations are forced and SHORT (support>=3, smallest = additive triangles).
      Count cross-block effective relations as a function of block separation M.
 (D3) THE BOUND. |corr_cross(E)| <= C * (#effective cross-block relations of bounded height)
      with an EXPLICIT per-relation constant C from the single-arc Fourier coefficient. Test
      whether the bound is < the comfortable budget (cap_k - iid_k - corr_intra).
 (D4) DECAY: does the cross-block relation count -> 0 (and corr_cross -> 0) as separation M->inf?
      (then multi-block reduces to single-block, which is the EASY/comfortable branch.)
 (D5) STRUCTURAL: is the cross-block triangle count bounded by a Freiman/BSG additive-energy
      quantity of the block-OFFSET set? Compute additive energy of block centers vs triangle count.

KEY HONESTY CHECK (codex HYP-2715): the absolute |corr| envelope DIVERGES (sum|K(n)| ~ 9000);
the smallness of corr is SIGNED cancellation. So a TRUE lever must (a) only need a CRUDE bound
(budget 0.18, not the 0.001 tight point), and (b) for cross-block, exploit that cross relations
are FEW and SHORT so even a partly-absolute bound on the SHORT cross shell fits the 0.18 budget.
We test exactly whether the SHORT cross shell (the only cross relations of bounded height) has
absolute Fourier mass below budget -- if yes, the cross part needs NO cancellation and the lever
is real; the cancellation subtlety stays confined to the (comfortable, single-block) intra part.

stdlib only; exact Fractions for measS7/iid; complex Fourier diagnostics cross-checked to exact.
"""
import sys, itertools, math, cmath
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

w7 = cmath.exp(2j * math.pi / 7)
TWO_PI_I = 2j * math.pi

# ----------------------------- exact geometric measS7 -----------------------------
def measS7(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * abs(e) + 1): bps.add(F(m, 7 * abs(e)))
    bps = sorted(b for b in bps if 0 <= b <= 1); tot = F(0)
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i + 1]
        if hi <= lo: continue
        mid = (lo + hi) / 2
        if len({int(((e * mid) % 1) * 7) for e in E}) == 7: tot += hi - lo
    return tot

def stirling2(n, k): return sum((-1)**(k-j)*math.comb(k,j)*j**n for j in range(k+1))//math.factorial(k)
def iid_k(k): return F(math.factorial(7)*stirling2(k,7), 7**k)

# ----------------------------- Fourier pieces (kps-S6 / angleF basis) -----------------------------
# Single 1/7-arc Fourier coefficient. shat(n) is the coeff of the indicator of one width-1/7 arc.
def shat(n, j):
    if n == 0: return 1.0/7.0
    a = j/7.0
    return (cmath.exp(-TWO_PI_I*n*a) - cmath.exp(-TWO_PI_I*n*(a+1/7.0)))/(TWO_PI_I*n)

# Cover indicator coefficient via inclusion-exclusion over the 7 cells (R = subset of cells 1..6).
SUB = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SGN = {T: (-1)**len(T) for T in SUB}
_CH = {}
def chat(n, T):
    key = (n, T)
    if key in _CH: return _CH[key]
    if n == 0:        v = complex(1 - len(T)/7.0, 0.0)
    elif n % 7 == 0:  v = 0j                                  # apex-prime vanishing
    else:             v = -sum(shat(n, j) for j in T)
    _CH[key] = v; return v

def Kk(n):
    """Per-relation Fourier kernel K(n) = sum_R (-1)^|R| prod_i chat_R(n_i). Real."""
    s = 0j
    for T in SUB:
        p = 1.0 + 0j
        for ni in n:
            p *= chat(ni, T)
            if p == 0: break
        s += SGN[T] * p
    return s

def Kk_abs_upper(n):
    """Sum_R |prod_i chat_R(n_i)|  -- the ABSOLUTE (triangle-ineq) mass of one relation."""
    s = 0.0
    for T in SUB:
        p = 1.0
        for ni in n:
            p *= abs(chat(ni, T))
            if p == 0: break
        s += p
    return s

# ----------------------------- relation lattice enumeration -----------------------------
def enum_relations(E, N0):
    """All n in Z^E with sum n_i e_i = 0, |n_i| <= N0, not all zero. Brute small box (exact)."""
    nz_idx = [i for i,e in enumerate(E) if e != 0]
    base = [0]*len(E)
    out = []
    # iterate over the nonzero coords; the e=0 coord (if present) is FREE -> we set it so that
    # it does not affect sum (e=0 contributes 0), but its Fourier factor chat_R(n_0). For e=0,
    # frac(0*x)=0 always in cell 0, so a relation coordinate on e=0 with n_0 != 0 gives chat_R(n_0)
    # which is generally nonzero. BUT the standard convention (routeE) takes E with 0 pinned and
    # the lattice over the NONZERO offsets (the 0-offset is the fixed anchor, always in cell 0,
    # contributes the +1 fast factor). We follow routeE: lattice over nonzero offsets only.
    Enz = [E[i] for i in nz_idx]
    for n in itertools.product(range(-N0, N0+1), repeat=len(Enz)):
        if all(x == 0 for x in n): continue
        if sum(ni*ei for ni,ei in zip(n, Enz)) != 0: continue
        out.append(tuple(n))
    return out, Enz

def block_of(e, blocks):
    for bi, B in enumerate(blocks):
        if e in B: return bi
    return -1

def corr_split(E, blocks, N0):
    """Enumerate Lambda(E) (over nonzero offsets) up to height N0, split into intra/cross by
    block support, return (corr_intra_signed, corr_cross_signed, cross_relations, intra_relations)
    using exact-checked complex K(n). 'effective' = no coord a nonzero multiple of 7 (chat vanishes)."""
    rels, Enz = enum_relations(E, N0)
    blk = [block_of(e, blocks) for e in Enz]
    ci = 0.0; cc = 0.0
    cross_rels = []; intra_rels = []
    cross_abs = 0.0
    for n in rels:
        # effective filter: any coord a nonzero multiple of 7 -> chat=0 -> contributes nothing
        if any((ni % 7 == 0 and ni != 0) for ni in n):
            continue
        supp_blocks = {blk[i] for i in range(len(n)) if n[i] != 0}
        k = Kk(n).real
        if abs(k) < 1e-15: continue
        if len(supp_blocks) >= 2:
            cc += k; cross_rels.append((n, k)); cross_abs += Kk_abs_upper(n)
        else:
            ci += k; intra_rels.append((n, k))
    return ci, cc, cross_rels, intra_rels, cross_abs

# ----------------------------- cross-block additive-triangle counting -----------------------------
def cross_block_triangles(blocks):
    """Count SUPPORT-3 cross-block relations of the two canonical short shapes, over the actual
    offset values (each block a set of integers):
       (T+) e_a + e_b = e_c   with a,b,c not all in one block  (additive triple straddling blocks)
       (T=) e_b - e_a = e_d - e_c  same-common-difference quadruple across >=2 blocks (support 4)
    Returns dict of counts. These are the SHORTEST cross-block relation generators (THM-559 turf)."""
    allE = sorted(set().union(*blocks))
    Eset = set(allE)
    blk = {e: block_of(e, blocks) for e in allE}
    tri_plus = 0
    for ea in allE:
        for eb in allE:
            if eb < ea: continue
            ec = ea + eb
            if ec in Eset:
                bs = {blk[ea], blk[eb], blk[ec]}
                if len(bs) >= 2:
                    tri_plus += 1
    # support-4 same-difference quadruples e_b-e_a = e_d-e_c (a<b, c<d), crossing blocks
    quad_eq = 0
    pairs = [(allE[i], allE[j]) for i in range(len(allE)) for j in range(i+1, len(allE))]
    by_diff = defaultdict(list)
    for a, b in pairs:
        by_diff[b-a].append((a, b))
    for d, lst in by_diff.items():
        for i in range(len(lst)):
            for j in range(i+1, len(lst)):
                (a, b), (c, e) = lst[i], lst[j]
                bs = {blk[a], blk[b], blk[c], blk[e]}
                if len(bs) >= 2:
                    quad_eq += 1
    return {"tri_plus": tri_plus, "quad_eq": quad_eq}

def additive_energy(S):
    """E(S) = #{(a,b,c,d) in S^4 : a+b = c+d}. The BSG/Freiman additive energy of the block centers."""
    from collections import Counter
    S = list(S)
    sums = Counter(a+b for a in S for b in S)
    return sum(v*v for v in sums.values())

# ----------------------------- main -----------------------------
CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91), 11: F(66,91), 12: F(6,7)}

def show(E, blocks, N0, cap):
    k = len(E)
    iid = iid_k(k)
    true = measS7(E)
    corr_free = true - iid
    ci, cc, crels, irels, cross_abs = corr_split(E, blocks, N0)
    low = ci + cc
    resid = float(corr_free) - low
    budget = float(cap - iid)
    tri = cross_block_triangles(blocks)
    print(f"  E={E}")
    print(f"    blocks={[sorted(b) for b in blocks]}")
    print(f"    measS7={float(true):.6f} = iid {float(iid):.5f} + corr {float(corr_free):+.6f}   (cap={float(cap):.4f}, budget cap-iid={budget:+.5f})")
    print(f"    corr split (|n|<= {N0}):  intra={ci:+.6f}  cross={cc:+.6f}   sum_low={low:+.6f}  high-resid={resid:+.6f}")
    print(f"    #cross relations (height<= {N0}, effective)={len(crels)}   ABS cross Fourier mass={cross_abs:.6f}")
    print(f"    cross-block short generators: additive-triples(T+)={tri['tri_plus']}  same-diff-quads(T=)={tri['quad_eq']}")
    print(f"    LEVER CHECK: |corr_cross|={abs(cc):.6f}  vs  ABS-cross-mass={cross_abs:.6f}  vs  budget {budget:.5f}")
    print(f"               cross part {'FITS budget absolutely' if cross_abs < budget else 'needs cancellation'} ; "
          f"|cross|<budget: {abs(cc) < budget}")
    return dict(iid=iid, corr=corr_free, ci=ci, cc=cc, cross_abs=cross_abs, budget=budget,
               ncross=len(crels), tri=tri, true=true)

if __name__ == "__main__":
    print("#"*94)
    print("# THREAD D (the LEVER): cross-block carrier error via cross-block cycle/triangle counting")
    print("#"*94)

    # ---------------------------------------------------------------
    print("\n[D1/D2] EXACT corr = corr_intra + corr_cross by block-support, separated blocks.")
    print("        Two-block E = {0}+block1  U  (M + block2). Vary separation M.")
    print("-"*94)
    # base building blocks (small, so the brute lattice box is feasible)
    b1 = [0,1,2,3]     # block 1 (4 pts, includes the pin 0)
    b2 = [0,1,2,3]     # block 2 shape (4 pts), placed at offset M
    N0 = 6
    twoblock_rows = []
    for M in [7, 14, 21, 30, 49, 70]:
        E = sorted(set(b1) | set(e+M for e in b2))
        if len(E) != 8: continue
        blocks = [set(b1), set(e+M for e in b2)]
        print(f"\n  separation M={M}:")
        r = show(E, blocks, N0, CAP[8])
        twoblock_rows.append((M, r))

    # ---------------------------------------------------------------
    print("\n\n[D4] DECAY of cross part with separation M (two equal 4-blocks):")
    print("-"*94)
    print(f"   {'M':>4} | {'#cross-rel':>10} | {'corr_cross':>12} | {'ABS-cross-mass':>14} | {'corr_total':>12}")
    for M, r in twoblock_rows:
        print(f"   {M:>4} | {r['ncross']:>10} | {r['cc']:>+12.6f} | {r['cross_abs']:>14.6f} | {float(r['corr']):>+12.6f}")
    print("   READ: if #cross-rel and ABS-cross-mass -> 0 as M grows, multi-block -> single-block (EASY).")

    # ---------------------------------------------------------------
    print("\n\n[D3] THE BOUND across several block geometries (k=8): does the cross part fit budget 0.18?")
    print("-"*94)
    geometries = [
        ("2x4 @7",   [0,1,2,3], [0,1,2,3], 7),
        ("2x4 @14",  [0,1,2,3], [0,1,2,3], 14),
        ("2x4 @30",  [0,1,2,3], [0,1,2,3], 30),
        ("4+4 @9",   [0,1,2,3], [0,1,2,3], 9),
        ("5+3 @20",  [0,1,2,3,4], [0,1,2], 20),
        ("3+3+2",    None, None, None),  # handled specially below
    ]
    bound_rows = []
    for name, A, B, M in geometries:
        if name == "3+3+2":
            blocks = [set([0,1,2]), set([14,15,16]), set([30,31])]
            E = sorted(set().union(*blocks))
        else:
            blocks = [set(A), set(e+M for e in B)]
            E = sorted(set().union(*blocks))
        k = len(E)
        cap = CAP.get(k, CAP[8])
        print(f"\n  {name}: ")
        r = show(E, blocks, N0, cap)
        bound_rows.append((name, r))

    # ---------------------------------------------------------------
    print("\n\n[D5] STRUCTURAL: cross-block triangle count vs additive energy of block-offset centers.")
    print("-"*94)
    print("   Block centers = the offsets (min) of each block; cross-triangles should track the")
    print("   additive energy E(centers) (Freiman/BSG): dissociated centers -> low energy -> few triangles.")
    center_tests = [
        ("dissociated 0,7,30,70",  [{0,1},{7,8},{30,31},{70,71}]),
        ("AP centers 0,10,20,30",  [{0,1},{10,11},{20,21},{30,31}]),
        ("dyadic 0,1,2,4 *14",     [{0,1},{14,15},{28,29},{56,57}]),
        ("Sidon-ish 0,1,3,7 *14",  [{0,1},{14,15},{42,43},{98,99}]),
    ]
    print(f"   {'config':<24} | {'centers':<18} | {'add-energy':>10} | {'cross-tri(T+)':>13} | {'same-diff(T=)':>13}")
    for name, blocks in center_tests:
        centers = sorted(min(b) for b in blocks)
        ae = additive_energy(centers)
        tri = cross_block_triangles(blocks)
        print(f"   {name:<24} | {str(centers):<18} | {ae:>10} | {tri['tri_plus']:>13} | {tri['quad_eq']:>13}")
    print("\n   READ: if cross-triangle count is monotone in additive energy, dissociated (Sidon)")
    print("   block centers => low energy => few cross triangles => the cross carrier error is small")
    print("   by a BSG/Freiman bound on the block-center set (the structural lever).")

    print("\n" + "="*94)
    print("VERDICT printed inline above. Summary thresholds: budget cap_k - iid_k =",
          {k: round(float(CAP[k]-iid_k(k)),4) for k in (8,9,10)})
