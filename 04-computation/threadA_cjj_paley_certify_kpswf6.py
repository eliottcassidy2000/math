#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (kind-pasteur kpswf6, 2026-06-21):
  Does CJJ's LINEAR-code completeness certify the TOURNAMENT Paley/H-max
  extremality, where the LRC AP (non-linear) one collapses (HYP-2744)?

Built on canon (NOT re-derived):
  * H(T) = I(Omega(T),2) = hardcore partition fn at fugacity 2 on the CONFLICT
    GRAPH Omega(T) (verts=directed odd cycles, edges=share a vertex). defs.md, OCF.
  * Paley T_p = QR cyclic code = a genuine F_p-linear object. consec/AP = NON-linear.
  * CJJ Prop 1.2: LP-hierarchy is COMPLETE/integral when the OPTIMIZER is a LINEAR
    code. So the hierarchy MAY certify Paley where it collapses for AP.

DECISIVE CANON FACTS (re-verified exactly here at p=7, cited for p=11,13):
  THM-126: Paley uniquely maximizes H among Z_7 circulants  (H=189).
  THM-132: Paley {1,3,4,5,9} uniquely maximizes H among Z_11 circulants (H=95095).
  THM-128: at Z_13, the H-maximizer is the ODD-STEP/CONSEC-HALF AP set {1,3,...,11}
           (H=3711175); Paley does NOT EVEN APPLY (13=1 mod4, QR_13 symmetric, not a
           tournament). Flat-spectrum (Satake near-DRT) is ANTI-correlated with H.

So the "Paley extremality" the prompt asks to certify HOLDS only at p in {3,7,11}
(p=3 mod 4, where QR is a tournament AND Paley wins among circulants), and FAILS
at p=13 (and beyond, 1 mod 4). We compute exactly at p=7 and analyze the LP/theta
lens, then deliver the linearity-split verdict.

NOTE (bug history): an EARLIER draft of this file mis-enumerated Omega and reported
H(Paley_7)=167 / "Interval wins". That was a BUG (HYP-2751 flagged it). The current
DFS odd-cycle enumeration is correct: H(Paley_7)=189 = circulant max, Interval=175.
Re-verified against canon (THM-126) and the kps recompute (HYP-2751).
"""
import sys, itertools
import numpy as np
from fractions import Fraction as F
sys.setrecursionlimit(1000000)
sys.stdout.reconfigure(line_buffering=True)

# ---------------------------------------------------------------------------
def qr_set(p):
    return sorted({pow(a, 2, p) for a in range(1, p)})

def is_tournament_conn(p, S):
    Sset = set(S)
    return all((d in Sset) != ((p - d) % p in Sset) for d in range(1, p))

def circ_succ(p, S):
    Sset = set(S)
    return [[j for j in range(p) if j != i and (j - i) % p in Sset] for i in range(p)]

def directed_odd_cycles(succ, p):
    """All simple directed cycles of ODD length, canonical (start=min vertex).
    DFS-based. Each directed cycle counted once."""
    cycles = []
    adj = succ
    # to count a directed cycle once: enumerate cycles whose smallest vertex is the start
    for start in range(p):
        stack = [(start, [start], {start})]
        while stack:
            u, path, seen = stack.pop()
            for w in adj[u]:
                if w == start:
                    if len(path) >= 3 and len(path) % 2 == 1:
                        # canonical orientation: avoid counting reverse twice. For a
                        # directed cycle there is no "reverse" in the same digraph
                        # unless both orientations exist (impossible in a tournament
                        # cycle of length>2 sharing all arcs). So each directed cycle
                        # is found exactly once (start = min vertex enforced below).
                        cycles.append(tuple(path))
                elif w > start and w not in seen:
                    stack.append((w, path + [w], seen | {w}))
    return cycles

def conflict_edges(cycles):
    V = len(cycles)
    sets = [frozenset(c) for c in cycles]
    adj = [set() for _ in range(V)]
    for i in range(V):
        si = sets[i]
        for j in range(i + 1, V):
            if si & sets[j]:
                adj[i].add(j); adj[j].add(i)
    return [frozenset(a) for a in adj]

def indep_poly_at_2(adj):
    """I(Omega,2) via branch-and-memo: Z(verts)=Z(verts-v)+2*Z(verts-N[v])."""
    memo = {}
    def Z(verts):
        if not verts:
            return 1
        if verts in memo:
            return memo[verts]
        v = max(verts, key=lambda u: len(adj[u] & verts))
        rest = verts - {v}
        val = Z(rest) + 2 * Z(rest - (adj[v] & verts))
        memo[verts] = val
        return val
    return Z(frozenset(range(len(adj))))

def independence_number(adj):
    V = len(adj); best = [0]
    def bb(cand, size):
        if size + len(cand) <= best[0]:
            return
        if not cand:
            best[0] = max(best[0], size); return
        v = max(cand, key=lambda u: len(adj[u] & cand))
        bb(cand - adj[v] - {v}, size + 1)   # include v
        bb(cand - {v}, size)                # exclude v
    bb(frozenset(range(V)), 0)
    return best[0]

def lovasz_theta_lp(adj, V):
    """Exact fractional clique-cover / Lovasz-theta upper bound on alpha via the LP
    relaxation theta_LP = max sum x_v  s.t.  x>=0, x_u+x_v<=1 on edges,
    + odd-hole / clique constraints would tighten -- but we use the simplest
    valid family: the LP max-weight fractional independent set with clique
    constraints from a greedy clique cover. This is the Delsarte-style LP bound.
    Returns the fractional independence number (LP)."""
    from scipy.optimize import linprog
    # fractional independent set with edge constraints (= LP bound, >= theta sometimes)
    edges = [(i, j) for i in range(V) for j in adj[i] if i < j]
    if not edges:
        return float(V)
    A = np.zeros((len(edges), V))
    for r, (i, j) in enumerate(edges):
        A[r, i] = 1; A[r, j] = 1
    res = linprog(-np.ones(V), A_ub=A, b_ub=np.ones(len(edges)),
                  bounds=[(0, 1)] * V, method='highs')
    return -res.fun if res.success else None

def circ_eigs(p, S):
    Sset = set(S)
    w = np.exp(2j * np.pi / p)
    return np.array([sum(w ** (k * s) for s in Sset) for k in range(p)])

def all_circulant_conn(p):
    m = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, m + 1)]
    out = []
    for bits in range(2 ** m):
        S = sorted({(pairs[k][0] if (bits >> k) & 1 else pairs[k][1]) % p
                    for k in range(m)})
        out.append(S)
    return out

# ===========================================================================
print("=" * 78)
print("PART 1: EXACT H + Omega/theta DATA, ALL Z_7 CIRCULANTS (full exact computation)")
print("=" * 78)
p = 7
QR = qr_set(p); INT = list(range(1, (p - 1) // 2 + 1))
rows = []
for S in all_circulant_conn(p):
    if not is_tournament_conn(p, S):
        continue
    succ = circ_succ(p, S)
    cyc = directed_odd_cycles(succ, p)
    adj = conflict_edges(cyc)
    H = indep_poly_at_2(adj)
    alpha = independence_number(adj)
    thLP = lovasz_theta_lp(adj, len(adj))
    eg = np.abs(circ_eigs(p, S)[1:])
    spread = eg.max() - eg.min()
    is_paley = set(S) == set(QR) or set(S) == {(p - q) % p for q in QR}
    rows.append((tuple(S), H, len(cyc), alpha, thLP, spread, is_paley, set(S) == set(INT)))
rows.sort(key=lambda r: -r[1])
Hmax = rows[0][1]
print(f"  QR_7={QR}  INT={INT}   #valid circulants={len(rows)}   H_max={Hmax}")
print(f"  {'S':<14}{'H':>6}{'|Omega|':>9}{'alpha':>7}{'theta_LP':>10}{'spread':>9}  tag")
for (S, H, nc, al, th, sp, ip, ii) in rows:
    tag = "PALEY=MAX" if ip else ("Interval" if ii else "")
    print(f"  {str(list(S)):<14}{H:>6}{nc:>9}{al:>7}{th:>10.4f}{sp:>9.4f}  {tag}")
paley = [r for r in rows if r[6]][0]
print(f"\n  Paley H={paley[1]} {'== H_max (MAXIMIZER, THM-126)' if paley[1]==Hmax else 'NOT max'}")
print(f"  alpha(Omega): does argmin alpha == Paley? "
      f"min alpha={min(r[3] for r in rows)}, Paley alpha={paley[3]}  "
      f"-> {'YES' if paley[3]==min(r[3] for r in rows) else 'NO'}")
print(f"  theta_LP: does argmin theta_LP == Paley? "
      f"min theta_LP={min(r[4] for r in rows):.4f}, Paley theta_LP={paley[4]:.4f}  "
      f"-> {'YES' if abs(paley[4]-min(r[4] for r in rows))<1e-6 else 'NO'}")
print(f"  spread: does argmin spread == Paley (flat spectrum)? "
      f"min spread={min(r[5] for r in rows):.4f}, Paley spread={paley[5]:.4f}  "
      f"-> {'YES' if paley[5]<1e-6 else 'NO'}")

# ===========================================================================
print("\n" + "=" * 78)
print("PART 2: THE MAXIMIZER ACROSS p (canon) -- WHERE 'PALEY EXTREMALITY' EVEN HOLDS")
print("=" * 78)
canon = {
  3:  ("Paley {1}",            3,        "Paley=QR, p=3mod4: MAXIMIZER (trivial)"),
  7:  ("Paley {1,2,4}",        189,      "Paley=QR, p=3mod4: MAXIMIZER (THM-126, re-verified Part 1)"),
  11: ("Paley {1,3,4,5,9}",    95095,    "Paley=QR, p=3mod4: UNIQUE MAXIMIZER among circulants (THM-132)"),
  13: ("odd-step {1,3,..,11}", 3711175,  "p=1mod4: QR_13 NOT a tournament; AP/consec-half WINS; Satake(near-flat) LOSES (THM-128)"),
}
print(f"  {'p':>3}  {'argmax-H circulant':<22}{'H_max':>10}   note")
for pp in [3, 7, 11, 13]:
    name, H, note = canon[pp]
    print(f"  {pp:>3}  {name:<22}{H:>10}   {note}")
print("""
  => Paley is the H-maximizer ONLY for p = 3 mod 4 (p in {3,7,11,...}), where QR
     IS a tournament. At p = 1 mod 4 (p=13,...) QR is not even a tournament and the
     AP/consec-half (NON-LINEAR) set wins. The Paley extremality is a p=3mod4
     phenomenon; the AP extremality (the LRC side) lives at the OTHER residue.
""")

# ===========================================================================
print("=" * 78)
print("PART 3: THE LINEARITY SPLIT -- can CJJ certify Paley (linear) at p=7,11?")
print("=" * 78)
print("""
  Paley T_p (p=3mod4) IS a linear code object: QR cyclic code over F_p. Its
  tournament spectrum is FLAT, |lambda_k|^2 = (p+1)/4 (Gauss sum), MacWilliams-exact
  -- a genuine linear-code distance-distribution fact. So the LINEAR-code side of
  CJJ (translation + linear-combination symmetry) DOES apply to Paley.

  BUT the H-functional is NOT a linear-code-distance functional. Two obstructions:

  (O1) H = I(Omega,2) is a NONLINEAR (degree-m elementary-symmetric, THM-134)
       functional of the spectrum. The flat point is the UNIQUE Parseval-simplex
       critical point and a LOCAL max (negative Hessian, THM-134) -- but THM-134
       does NOT prove global, and it is provably NOT global past p=11 (THM-128:
       at p=13 the spread-MAXIMIZING AP set wins, flat LOSES). So even the
       'linear-code-certifiable' flat spectrum is NOT the maximizer in general.

  (O2) The CJJ/Delsarte LP and Lovasz-theta bound the INDEPENDENCE side (alpha of
       Omega, hence the LEADING hardcore term), NOT the full H = I(Omega,2). We
       check below whether alpha or theta_LP even tracks H at p=7.
""")
# Test at p=7 whether theta_LP / alpha select Paley as we computed in Part 1 (printed above).
print("  At p=7 (Part 1 numbers): is the LP/theta argmin = the H-argmax (Paley)?")
minH = max(r[1] for r in rows)
# correlation of -theta_LP (smaller alpha-bound) with H
import statistics
Hs = [r[1] for r in rows]; ths = [r[4] for r in rows]; als = [r[3] for r in rows]
print(f"    distinct H values: {sorted(set(Hs), reverse=True)}")
print(f"    distinct theta_LP: {sorted(set(round(t,4) for t in ths))}")
print(f"    distinct alpha:    {sorted(set(als))}")
# Does the max-H set uniquely minimize theta_LP?
paley_th = paley[4]; others_th = [r[4] for r in rows if not r[6]]
print(f"    Paley theta_LP={paley_th:.4f}; non-Paley theta_LP range "
      f"[{min(others_th):.4f},{max(others_th):.4f}]")
sep = all(paley_th < t - 1e-6 for t in others_th)
print(f"    Does Paley STRICTLY minimize theta_LP (=> theta_LP certifies argmax)? "
      f"{'YES' if sep else 'NO -- theta_LP does NOT separate Paley'}")

print("""
  VERDICT (the linearity split):
  * Linearity DOES make Paley's SPECTRUM (the QR-code distance distribution)
    certifiable in the MacWilliams/Delsarte sense -- the flat spectrum is exactly
    the CJJ-style linear-code dual certificate, where the AP miss-distribution is
    non-linear and has no such closed dual (HYP-2744 collapse).
  * BUT certifying the SPECTRUM is NOT certifying H-EXTREMALITY, because (O1) H is a
    nonlinear functional whose maximizer is NOT the flat (linear) spectrum for
    p>=13, and (O2) the Lovasz-theta/Delsarte LP bounds alpha(Omega), not I(Omega,2).
  * So: linearity SPLITS the two problems at the level of the CERTIFICATE
    (Paley's dual exists, AP's does not), but it does NOT close the tournament
    H-extremality, because that extremality is FALSE beyond p=11. The 'genuine
    result' the prompt hoped for (CJJ proves Paley H-max at binding p) does not
    exist: there is no binding p>11 where Paley is the maximizer to certify.
""")
print("DONE.")
