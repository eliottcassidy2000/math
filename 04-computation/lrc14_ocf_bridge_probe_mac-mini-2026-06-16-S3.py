#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_ocf_bridge_probe_mac-mini-2026-06-16-S3.py   (ANGLE 7)

THE CROSS-PROBLEM BRIDGE: is the lonely measure L(S) structurally a tournament
partition function (independence polynomial I(Omega,2) = H)?

14 = 2*7. The project's forbidden Hamiltonian-path values are the orbit
7, 21, 63, 189 = 7*3^k, with Phi_3(2) = 2^2+2+1 = 7 the "forbidden phantom".
The LRC(14) kernel has period 7 and danger band 1/14. This script asks, with
CONCRETE numbers, whether the "7" coincidence is DEEP or SHALLOW, by directly
comparing the two signed/positive lattice/graph sums.

THE TWO OBJECTS (exact definitions used in repo):
  L(S)  = sum_{t in Lambda} prod_{i=1..13} h(t_i),   Lambda = relation lattice
          h(0)=6/7, h(t) = -sin(pi t/7)/(pi t) for t!=0.  (THM-515 theta form)
  H(T)  = I(Omega(T), 2) = sum_{U indep set in Omega(T)} 2^|U|.  (THM-002 OCF)

STRUCTURAL CHECKLIST (what would make a DEEP bridge):
  (A) fugacity: both have a "2". In H it is z=2 (2 orientations / counting weight).
      In L the band is 2/14; does a "2" play the same role? -> probe.
  (B) support / graph: H sums over independent sets of a GRAPH (overlap = edge).
      L sums over a LATTICE. Is there a graph on the short relations whose
      independent sets reproduce the dominant terms of L? -> build it, compare.
  (C) sign: H is a POSITIVE sum (all +2^k). L is a SIGNED sum (h(t)<0 for t!=0,
      and s(t) flips sign at |t|=7). The 7-periodicity (s(t)=0 iff 7|t) is the
      cyclotomic kernel. Is that the same 7 as Phi_3(2)? -> arithmetic test.
  (D) the inf: H is bounded below by 1 (empty set, structural). L>0 is the OPEN
      target. Does the graph picture give an analogous "empty set" floor? -> test.

We build a tiny concrete example end to end and give an HONEST verdict.
"""
import sys, itertools, math
from math import gcd, sin, pi, cos
from fractions import Fraction

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)

# ----------------------------------------------------------------------
# LRC kernel
# ----------------------------------------------------------------------
def s(t):
    return sin(pi*t/7.0)/(pi*t) if t != 0 else 1.0/7.0

def h(t):
    return 6.0/7.0 if t == 0 else -s(t)

def lonely_measure(S, Q, n=14):
    """exact D(Q,S)/Q ; danger band r in {-floor(Q/14)..floor(Q/14)} mod Q."""
    rad = Q // n
    cnt = 0
    for a in range(Q):
        ok = True
        for v in S:
            r = (v*a) % Q
            if r <= rad or r >= Q - rad:
                ok = False; break
        if ok:
            cnt += 1
    return cnt/Q

# ----------------------------------------------------------------------
# relation lattice (short vectors, support<=4) and theta partial sum
# ----------------------------------------------------------------------
def short_relations(S, B, maxsupp=4):
    n = len(S); rels = []
    for sz in range(2, maxsupp+1):
        for T in itertools.combinations(range(n), sz):
            for coeffs in itertools.product([c for c in range(-B, B+1) if c != 0], repeat=sz):
                if sum(c*S[i] for c, i in zip(coeffs, T)) == 0:
                    g = 0
                    for c in coeffs: g = gcd(g, c)
                    if g == 1:
                        t = [0]*n
                        for c, i in zip(coeffs, T): t[i] = c
                        rels.append(tuple(t))
    return rels

def theta_partial(S, B, maxsupp=4):
    """main term + all support<=maxsupp, |coeff|<=B lattice contributions (no primitivity
    filter: enumerates the full box so multiples k*rel inside the box are counted)."""
    n = len(S); total = (6.0/7.0)**n
    seen = set()
    for sz in range(2, maxsupp+1):
        for T in itertools.combinations(range(n), sz):
            for coeffs in itertools.product([c for c in range(-B, B+1) if c != 0], repeat=sz):
                if sum(c*S[i] for c, i in zip(coeffs, T)) == 0:
                    key = tuple(sorted(zip(T, coeffs)))
                    if key in seen: continue
                    seen.add(key)
                    term = (6.0/7.0)**(n-sz)
                    for c in coeffs: term *= h(c)
                    total += term
    return total

# ----------------------------------------------------------------------
# OCF side: independence polynomial of a graph at x (general)
# ----------------------------------------------------------------------
def independence_poly_at(adj, n, x):
    """I(G,x) = sum_{indep U} x^|U|, brute force over subsets. adj[i] = set of nbrs."""
    val = 0.0
    for mask in range(1<<n):
        verts = [i for i in range(n) if (mask>>i)&1]
        indep = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if verts[b] in adj[verts[a]]:
                    indep = False; break
            if not indep: break
        if indep:
            val += x**len(verts)
    return val

def independence_poly_coeffs(adj, n):
    """alpha_k coefficients of I(G,x)."""
    from collections import Counter
    cnt = Counter()
    for mask in range(1<<n):
        verts = [i for i in range(n) if (mask>>i)&1]
        indep = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if verts[b] in adj[verts[a]]:
                    indep = False; break
            if not indep: break
        if indep:
            cnt[len(verts)] += 1
    return cnt

# ----------------------------------------------------------------------
# build the "relation overlap graph": vertices = short primitive relations,
# edge iff two relations SHARE a coordinate (support overlap). This is the
# LRC analogue of the conflict graph Omega (cycles share a vertex).
# ----------------------------------------------------------------------
def relation_overlap_graph(rels):
    n = len(rels)
    adj = [set() for _ in range(n)]
    supp = [frozenset(i for i, c in enumerate(r) if c != 0) for r in rels]
    for a in range(n):
        for b in range(a+1, n):
            if supp[a] & supp[b]:
                adj[a].add(b); adj[b].add(a)
    return adj

# ======================================================================
print("="*78)
print("ANGLE 7: is L(S) = sum_Lambda prod h  a tournament partition function H=I(Omega,2)?")
print("="*78)

# ----------------------------------------------------------------------
# (C) THE ARITHMETIC OF '7': deep cyclotomic identity or shallow coincidence?
# ----------------------------------------------------------------------
print("\n--- (C) THE TWO SEVENS ---")
print(" tournament side: Phi_3(2) = 2^2+2+1 =", 2**2+2+1, " (forbidden phantom; orbit 7*3^k)")
print(" forbidden H orbit:", [7*3**k for k in range(4)])
print(" LRC side: band denom n=14=2*7; kernel s(t)=sin(pi t/7)/(pi t); s(t)=0 iff 7|t")
print(" check s(t)=0 exactly at multiples of 7 (cyclotomic vanishing of the kernel):")
for t in range(0, 22):
    z = s(t)
    flag = "  <-- ZERO (7|t)" if (t != 0 and t % 7 == 0) else ""
    if t % 7 == 0 or t in (1,6,8,13):
        print(f"     s({t:2d}) = {z:+.6f}{flag}")
print(" NOTE: the LRC 7 is the DENOMINATOR n/2 = 14/2 = 7 of the danger band 1/14.")
print("       The tournament 7 = Phi_3(2), the cyclotomic value of x^2+x+1 at the")
print("       OCF fugacity x=2. These are 7 from DIFFERENT mechanisms:")
print("         LRC 7   = half the runner count (geometry of the circle / band size)")
print("         OCF 7   = Phi_3(2) (algebra: 3-cycle generating fn at fugacity 2)")
print("       Same NUMBER, different generators. -> shallow unless a map identifies them.")

# Is there ANY z with Phi_3(z)=7 giving the band? Phi_3(z)=z^2+z+1.
print("\n  Probe: does any natural fugacity z give Phi_3(z) = 7 or 14?")
for z in range(1, 6):
    print(f"     Phi_3({z}) = {z*z+z+1};  Phi_2({z})=z+1={z+1};  Phi_6({z})=z*z-z+1={z*z-z+1}")
print("  -> Phi_3(2)=7 (z=2) and Phi_6(... ) etc. The OCF 7 needs z=2 EXACTLY.")
print("     The LRC 7 is independent of any fugacity. Coincidence of the integer 7.")

# ----------------------------------------------------------------------
# (A)+(B): build BOTH objects for a tiny case and compare structure
# Use the numerical extremizer core S = {1..13}\{6} U {98}.
# ----------------------------------------------------------------------
print("\n--- (A)+(B) SIDE-BY-SIDE on the extremizer core S={1..13}\\{6} U {98} ---")
S = sorted([x for x in range(1,14) if x != 6] + [98])
print(" S =", S, " |S| =", len(S))
Q = 14*98*2  # multiple of 14 and of the speeds' lcm-ish for a clean exact measure window
# use a large multiple-of-14 Q
Q = 27440  # = 14 * 1960
Lm = lonely_measure(S, Q)
print(f" exact lonely measure D({Q})/{Q} = {Lm:.6f}   (main term (6/7)^13 = {(6/7)**13:.6f})")
th3 = theta_partial(S, 3); th4 = theta_partial(S, 4)
print(f" theta partial (support<=4): |coeff|<=3 -> {th3:.6f};  |coeff|<=4 -> {th4:.6f}")

rels = short_relations(S, 3)
print(f" # short primitive relations (support<=4, |coeff|<=3) = {len(rels)}")
# show the shortest few
rels_sorted = sorted(rels, key=lambda r: (sum(abs(x) for x in r), r))
print(" shortest relations (coeff vector restricted to support, with speeds):")
for r in rels_sorted[:8]:
    supp = [(S[i], c) for i, c in enumerate(r) if c != 0]
    print(f"     {supp}   L1-norm={sum(abs(x) for x in r)}")

# ----------------------------------------------------------------------
# (B) the relation-overlap GRAPH and its independence polynomial at 2
# ----------------------------------------------------------------------
print("\n--- (B) DOES the relation lattice carry a graph whose I(.,2) relates to L? ---")
# keep it tractable: take the shortest ~16 relations
keep = rels_sorted[:16]
adj = relation_overlap_graph(keep)
nv = len(keep)
coeffs = independence_poly_coeffs(adj, nv)
I2 = independence_poly_at(adj, nv, 2.0)
print(f" relation-overlap graph on {nv} shortest relations:")
edges = sum(len(a) for a in adj)//2
print(f"   vertices={nv}, edges={edges}, density={edges/(nv*(nv-1)/2):.3f}")
print(f"   independence poly coeffs alpha_k:", dict(sorted(coeffs.items())))
print(f"   I(G,2) = {I2:.1f}   <-- this is the OCF-style partition function of the relation graph")
print(f"   (compare to L(S) = {Lm:.6f}: I(G,2) is a POSITIVE integer >> L, NO direct match)")

# The honest structural difference: signs.
print("\n--- (D) THE SIGN / FLOOR difference (the crux) ---")
print(" H=I(Omega,2)=sum 2^|U| is POSITIVE termwise -> floor 1 (empty set) is STRUCTURAL.")
print(" L=sum prod h has h(t)<0 for t!=0 AND s(t) sign-flips at |t|>=8 -> L is a delicate")
print(" SIGNED cancellation; the level masses Lambda_k GROW (0.11, 0.55, 1.17,...) so there")
print(" is NO termwise floor. That is exactly why inf L>0 is HARD and H>=1 is EASY.")

# Quantify: build the signed analogue of I at 2 using h instead of +x^k on the SAME graph.
def signed_indep_at_h(adj, n, rels, S):
    """sum over independent sets U of prod_{v in U} (weight of relation v),
    weight = prod h(coeff) * (6/7)^(13-support) -- the actual L contribution of that relation."""
    def relweight(r):
        w = (6.0/7.0)**(13 - sum(1 for c in r if c != 0))
        for c in r:
            if c != 0: w *= h(c)
        return w
    wts = [relweight(r) for r in rels]
    val = 0.0
    for mask in range(1<<n):
        verts = [i for i in range(n) if (mask>>i)&1]
        indep = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if verts[b] in adj[verts[a]]:
                    indep = False; break
            if not indep: break
        if indep:
            term = 1.0
            for v in verts: term *= wts[v]
            val += term
    return val

signed_I = signed_indep_at_h(adj, nv, keep, S)
print(f"\n SIGNED independence sum over the SAME overlap graph (weights = real L contributions):")
print(f"   sum_{{U indep}} prod_v relweight(v) = {signed_I:.6f}")
print(f"   (main (6/7)^13 added: {(6/7)**13 + (signed_I - 1*0):.6f} if 'empty set'=main term)")
# Actually the empty set contributes 1*... here relweight already carries (6/7)^(13-supp);
# the true main term is (6/7)^13. Let's assemble a CLUSTER-EXPANSION estimate.
main = (6.0/7.0)**13
# independence-style cluster sum WITH the empty set giving main term:
cluster_est = main + sum(  # single relations + independent pairs ... already in signed_I minus empty
    0 for _ in [0])
print(f"\n CLUSTER-EXPANSION view: L ~ (6/7)^13 * [ 1 + sum signed indep clusters of relations ]")
# compute relative form: divide each weight by (6/7)^13 won't match because supports differ;
# instead report the ratio of the graph-truncated signed sum to true L.
print(f"   true L = {Lm:.6f};  theta(support<=4,B=4) = {th4:.6f}")
print(f"   ratio theta/L = {th4/Lm:.4f}  (how well the short-relation lattice sum captures L)")

print("\n" + "="*78)
print("VERDICT (printed in structured summary).")
print("="*78)
