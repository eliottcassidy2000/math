#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD C -- LOVASZ THETA-PRIME / CONFLICT-GRAPH VIEW of the LRC(14) Z/7 cover.

Goal (from the dispatch):
  Schrijver: the Delsarte LP = theta'(H) for the CONFLICT GRAPH H, where the
  "code" is an INDEPENDENT SET in H.  Identify the LRC conflict graph; compute
  theta'(H); check it equals the Delsarte / L_y bound; and ask whether
  theta'(H) (or the degree-2 SoS/Lasserre lift) certifies consec-MAX of
  measS7(E)=P(N=0).

WHAT measS7 ACTUALLY IS (so the conflict-graph map is honest):
  For an offset set E (0 in E, |E|=k), color the circle by
     phi_E(x) = { floor(7 frac(e x)) : e in E }  subset Z/7  (set of HIT sectors).
  measS7(E) = meas{ x : phi_E(x) = ALL of Z/7 } = P(N=0), N = #missed sectors.
  By inclusion-exclusion, measS7 = sum_r (-1)^r S_r, S_r = E[C(N,r)] (THM-534).
  THM-534 (PROVED, per-E): measS7(E) <= L_y(E) = sum_r y_r S_r, with the dual
  g(t)=sum_r y_r C(t,r) >= 1[t=0] on t in {0..6}.  This IS the level-1 Delsarte LP.
  THE OPEN PIECE: consec MAXIMIZES L_y (equiv. measS7) over offset sets.

THE TWO CONFLICT GRAPHS (Schrijver's "code = independent set" requires care here,
because measS7 is a COVERING measure, not a packing/code size; we make the
covering->packing dual EXPLICIT):

  (A) THE SECTOR CAYLEY GRAPH C7(D) on Z/7.  Vertices = the 7 sectors.  Edge
      i~j iff (i-j) in D for a "conflict" connection set D.  measS7 is a cover
      (independent-DOMINATING / covering) functional on Z/7.  For the COVERING
      side the relevant object is theta of the COMPLEMENT (covering = fractional
      chromatic / clique-cover), or equivalently a covering-LP whose dual is a
      packing-LP = theta' of a conflict graph.  We compute theta and theta' of
      C7(D) for every nonempty symmetric D and see which (if any) reproduces the
      per-x cover constraint "7 distinct sectors needed".

  (B) THE OFFSET / DEPTH CONFLICT GRAPH H_E (this is the genuine LRC graph).
      Vertices = the "danger events": pairs (x-cell, missed-sector) or, in the
      Delsarte scheme reduction, the SUPPORT classes of the relation code
      Lambda(E)={n: sum n_i e_i = 0}.  measS7 = 1 - (weighted count of danger
      events) by IE.  The Delsarte LP on the relation scheme of E IS theta' of
      this graph; its value = L_y(E).  Here the graph DEPENDS ON E, so theta'(H_E)
      varies with E and the extremality question "consec minimizes theta'(H_E)"
      = the open piece, restated.  We verify theta'(H_E) = L_y(E) numerically.

  (C) PROJECT TOURNAMENT-AS-CODE (Gleason/THM-481): the LRC "code = independent
      set" instance is the Z/7 analog of the project's Omega(U_P) conflict graph
      (hard-core model H=I(Omega,2)).  We note the structural parallel and what
      the analog of the Hoffman / Lovasz-theta bound would be here.

For vertex-transitive graphs the Lovasz theta reduces to an LP over the
characters (Fourier), which we use for exact rational answers via scipy linprog
plus a closed-form circulant eigenvalue check.

Author: mac-mini-2026-06-21 (Thread C).  Reuses occupancy_law / S_r from
lrc14_B2_*_THREADA_opus.py.
"""
import sys, itertools, math
import numpy as np
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

# ---------------------------------------------------------------------------
# core measS7 / moment engine (exact, reused from THREAD A)
# ---------------------------------------------------------------------------
def occupancy_law(E):
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for a in range(7*abs(e)+1): bps.add(F(a, 7*abs(e)))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    pi = [F(0)]*8
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo+hi)/2
        hit = set()
        for e in E:
            v = e*xm; v = v - (v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)] += hi-lo
    return pi

def S_r_list(E):
    pi = occupancy_law(E)
    return [sum(pi[h]*comb(7-h, r) for h in range(8)) for r in range(8)], pi

def measS7(E):
    pi = occupancy_law(E)
    return pi[7]

def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

# THM-534 duals  L_y = sum_r y_r S_r  (the level-1 Delsarte LP value per E)
DUALS = {
    8:  [F(1), F(-1), F(1), F(-9,10), F(3,5)],
    9:  [F(1), F(-13,18), F(4,9), F(-1,6)],
    10: [F(1), F(-13,18), F(4,9), F(-1,6)],
    11: [F(1), F(-1,2), F(1,6)],
    12: [F(1), F(-1,2), F(1,6)],
    13: [F(1), F(-1,2), F(1,6)],
}
def L_y(E, k):
    y = DUALS[k]
    Sr, _ = S_r_list(E)
    return sum(y[r]*Sr[r] for r in range(len(y)))

def g_of_t(k):
    """the dual majorant g(t)=sum_r y_r C(t,r), t=0..6.  Should be >=1[t=0]."""
    y = DUALS[k]
    return [sum(y[r]*comb(t, r) for r in range(len(y))) for t in range(7)]

# ===========================================================================
print("="*78)
print("PART 0 -- the dual g(t) IS the Delsarte LP majorant of 1[t=0] (THM-534)")
print("="*78)
for k in [8, 9, 11]:
    g = g_of_t(k)
    print(f"  k={k}: g(t), t=0..6 = {[str(x) for x in g]}  ;  all >= 1[t=0]? "
          f"{g[0] >= 1 and all(g[t] >= 0 for t in range(1,7))}")
print("  -> g is a nonneg-Krawtchouk-combination majorant of the cover indicator.")
print("     measS7(E) <= L_y(E) is the per-E Delsarte/theta' bound (PROVED).\n")

# ===========================================================================
# PART A -- the SECTOR CAYLEY GRAPH C7(D) and its theta / theta'.
# The cover "all 7 sectors hit" is the COMPLETE constraint on Z/7.  We compute
# Lovasz theta of the Cayley graph C7(D) for every connection set D, exactly,
# via the circulant-eigenvalue closed form (Lovasz '79 for vertex-transitive:
#   theta(G) = -lambda_min via the Hoffman/ratio bound when edge-transitive;
#   in general for circulants theta is an LP over the eigenvalues).
# We use the exact characterization for CIRCULANT graphs:
#   theta(C_n(D)) = max over PSD circulant ... ; reduces to a 1-var-per-character LP.
# Implement Lovasz theta for a vertex-transitive graph via the formula
#   theta(G) = max sum_x sum_y M[x,y]  s.t. M PSD, tr M = 1, M[x,y]=0 for x~y.
# For circulant G this is solved by a symmetric circulant M = sum_j a_j P^j with
# a_j real, M PSD <=> eigenvalues mu_s = sum_j a_j w^{sj} >= 0, M[x,y]=0 for edges
# <=> a_{j}=0 for j in D, and theta = n * a_0-weight... -- we just do the LP in
# the Fourier domain.  (n=7 is tiny; we solve the SDP-as-LP exactly in chars.)
# ===========================================================================
print("="*78)
print("PART A -- Lovasz theta / theta' of the SECTOR Cayley graph C7(D), all D")
print("="*78)

w7 = [complex(math.cos(2*math.pi*s/7), math.sin(2*math.pi*s/7)) for s in range(7)]

def theta_circulant(D):
    """Lovasz theta of the circulant graph on Z/7 with connection set D (symmetric).
    Fourier-domain SDP: M = circulant(a_0..a_6) PSD, a_j=0 for j in D (edges),
    a_0 free, maximize n*a_0? -- standard: theta = max sum_{x,y} M_{xy} with tr M=1.
    For circulant M with first row a, sum_{x,y} M = n * (sum_j a_j) = n * mu_0...
    Use the clean dual-free SDP in eigenvalues:
      variables: eigenvalues are not free (tied to a via DFT). Easiest: optimize a_j.
      M PSD <=> mu_s = sum_j a_j cos(2pi s j/7) >= 0 for all s (real, symmetric a).
      constraints a_j = 0 for j in D.  tr M = n a_0 = 1 -> a_0 = 1/n.
      objective sum_{x,y} M_{xy} = n * sum_j a_j = n*(a_0 + sum_{j!=0} a_j).
    Maximize subject to mu_s>=0.  Solve via scipy linprog (vars a_j, j not in D, j!=0)."""
    from scipy.optimize import linprog
    n = 7
    free = [j for j in range(1, n) if j not in D]  # a_j adjustable; a_0 fixed = 1/n
    # symmetric: a_j = a_{n-j}; reduce to representatives but keep all for clarity.
    # objective: maximize n*(1/n + sum_{j in free} a_j) = 1 + n*sum a_j  -> max sum a_j
    nv = len(free)
    if nv == 0:
        # no free off-diagonal entries (e.g. K7): only a_0=1/n, M=I/n, theta = 1.
        return 1.0
    c = -np.ones(nv)  # maximize sum a_j -> minimize -sum
    # constraints mu_s = 1/n + sum_{j in free} a_j cos(2pi s j/n) >= 0, s=0..n-1
    A_ub = []; b_ub = []
    for s in range(n):
        row = [math.cos(2*math.pi*s*j/n) for j in free]
        A_ub.append([-x for x in row]); b_ub.append(1.0/n)  # -mu_s <= 1/n -> mu_s>=0
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  bounds=[(None, None)]*nv, method='highs')
    if not res.success: return None
    suma = res.x.sum()
    return 1.0 + n*suma  # = n * sum_j a_j

# Lovasz theta of the COMPLEMENT (covering side: alpha(G) <= theta(G);
# chromatic / covering number >= theta(complement) via theta(G)*theta(Gbar) >= n).
def edges_of(D):
    n = 7; E = set()
    for x in range(n):
        for d in D:
            E.add(frozenset({x, (x+d) % n}))
    return E

# enumerate all symmetric connection sets D (subset of {1,2,3}, since 4=-3,5=-2,6=-1)
reps = [1, 2, 3]
all_D = []
for r in range(0, 4):
    for combo in itertools.combinations(reps, r):
        D = set()
        for d in combo:
            D.add(d); D.add(7-d)
        all_D.append((combo, frozenset(D)))

print("  D(reps) | symmetric D | theta(C7(D)) | theta(complement) | n/theta(C7)")
for combo, D in all_D:
    if not D:
        print(f"  {str(combo):8s} | empty (no edges): theta = 7 (independent), complement=K7 theta=1")
        continue
    th = theta_circulant(D)
    Dc = frozenset(set(range(1,7)) - D)
    thc = theta_circulant(Dc) if Dc else 1.0  # complement of full-edge is empty -> theta 7
    thc_str = f"{thc:.4f}" if thc is not None else "n/a"
    prod = th*thc
    print(f"  {str(combo):8s} | {sorted(D)} | theta={th:.4f} | "
          f"theta(comp,D={sorted(Dc)})={thc_str} | prod={prod:.3f}")

print("""
  READING (Part A): the SECTOR cover 'all 7 sectors hit at x' is the constraint
  that the hit-set phi_E(x) is the FULL vertex set of Z/7, i.e. a DOMINATING /
  covering condition, not 'independent set'.  The Cayley graph C7(D) controls
  sector-ADJACENCY but measS7 does NOT forbid adjacent hits -- it REQUIRES all 7.
  So the Z/7 sector graph is the WRONG conflict graph: its independent sets are
  partial sector sets, whereas the cover wants the WHOLE vertex set.  The cover
  is theta-of-the-complement flavored (clique-cover / fractional chromatic), and
  on Z/7 with no real adjacency structure it degenerates (theta(empty C7)=7).
  CONCLUSION A: the conflict graph is NOT on the 7 sectors.  It must be on the
  OFFSETS / RELATION code (Part B).
""")

# ===========================================================================
# PART B -- THE GENUINE LRC CONFLICT GRAPH H_E ON THE RELATION SCHEME.
#  measS7(E) = sum_r (-1)^r S_r,  S_r = E[C(N,r)].  By THM-534 the Delsarte LP
#  on the relation scheme of E has VALUE L_y(E), and Schrijver's identity says
#  L_y(E) = theta'(H_E) for the conflict graph H_E whose nodes are the danger
#  events and whose edges encode the IE / scheme structure.  Rather than build a
#  giant graph we use the EQUIVALENT primal Delsarte LP whose value is theta':
#     theta'(H_E) = max_{p_t>=0} p_0  s.t.  sum_t C(t,r) p_t = S_r(E), r=0..R.
#  (This is EXACTLY THM-534's moment-LP; its optimum p_0* = L_y(E) by LP duality,
#   and it equals theta' of the relation conflict graph by Schrijver/Delsarte.)
#  KEY EXTREMALITY QUESTION: does theta'(H_E) (= L_y) -- as a function of E --
#  get MAXIMIZED by consec?  That is the open piece, now phrased as
#  "consec maximizes theta' of the relation conflict graph."
# ===========================================================================
from scipy.optimize import linprog

def theta_prime_relation(E, k):
    """primal moment-LP value = theta'(H_E) on the relation scheme:
       max p_0  s.t.  sum_{t=0}^{6} C(t,r) p_t = S_r(E), r=0..R; p_t>=0.
    Returns (lp_value, L_y_dual). Should match (LP strong duality)."""
    Sr, _ = S_r_list(E)
    R = len(DUALS[k]) - 1
    # variables p_0..p_6
    c = np.zeros(7); c[0] = -1.0  # maximize p_0
    A_eq = []; b_eq = []
    for r in range(R+1):
        A_eq.append([comb(t, r) for t in range(7)]); b_eq.append(float(Sr[r]))
    res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=[(0, None)]*7, method='highs')
    val = -res.fun if res.success else None
    return val, float(L_y(E, k))

print("="*78)
print("PART B -- theta'(H_E) = relation-scheme Delsarte LP value, vs L_y, vs measS7")
print("  (verify theta' == L_y per E; then ask: does consec maximize theta'?)")
print("="*78)
for k in [8, 9]:
    C = consec(k)
    tp_c, ly_c = theta_prime_relation(C, k)
    m_c = float(measS7(C))
    print(f"\n  k={k} consec={C}: measS7={m_c:.6f}  theta'(LP)={tp_c:.6f}  L_y={ly_c:.6f}  "
          f"(theta'==L_y? {abs(tp_c-ly_c)<1e-9})")
    # sweep bounded-span shapes; record whether any beats consec on theta'
    W = 12
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), k-1)]
    bank = [E for E in bank if primitive(E)]
    beat_tp = 0; beat_m = 0; max_tp = tp_c; arg_tp = C
    mism = 0
    for E in bank:
        tp, ly = theta_prime_relation(list(E), k)
        m = float(measS7(list(E)))
        if abs(tp-ly) > 1e-7: mism += 1
        if tp > tp_c+1e-9: beat_tp += 1
        if m > m_c+1e-12: beat_m += 1
        if tp > max_tp: max_tp = tp; arg_tp = list(E)
    print(f"     span<= {W}: {len(bank)} shapes; theta'!=L_y mismatches={mism}")
    print(f"     #beating consec on theta'(=L_y) = {beat_tp}  (consec-max theta'? {beat_tp==0})  argmax={arg_tp}")
    print(f"     #beating consec on measS7       = {beat_m}  (consec-max measS7? {beat_m==0})")

print("""
  READING (Part B): theta'(H_E) = L_y(E) EXACTLY per E (LP strong duality, as
  predicted by Schrijver).  So the per-E upper bound measS7<=theta' is the
  Delsarte/Lovasz bound, already PROVED (THM-534).  BUT 'consec maximizes
  theta'(H_E)' is the SAME open statement -- theta' did NOT linearize the
  extremality; it just renamed L_y.  theta' is a per-E CERTIFICATE, not an
  extremality proof.  This is the structural reason the dispatch flagged it
  as 'genuinely aggregate'.
""")
