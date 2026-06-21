#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_challenger56_kpswf2.py  (kind-pasteur 2026-06-21, THREAD 3 of OPEN-Q-108)

THE '56 CHALLENGER SHAPES = 56 TOURNAMENTS ON 6 VERTICES' PROBE.

Frontier (HYP-2723 / HYP-2154): for a k-set E of co-offsets, the carrier error
    corr(E) = measS7(E) - iid_k = sum_{0!=n in Lambda(E)} K(n)
is governed by the relation code Lambda(E)={n in Z^k : sum n_i e_i = 0}.  The
LEADING binding layer is SUPPORT-3: the support-3 primitive relations are exactly
   * 3-APs            e_a - 2 e_b + e_c = 0   (coefs (1,-2,1))
   * Schur triples    e_a + e_b - e_c = 0     (coefs (1,1,-1))
(additive-energy triples; HYP-2154: tightness REQUIRES these support-3 +-1 vectors).

THIS SCRIPT (rigorous, exact):
  (a) Define the SUPPORT-3 RELATION HYPERGRAPH H3(E): vertices = elements of E,
      hyperedges = the (unordered) triples {a,b,c} carrying a primitive support-3
      relation, LABELED by type (AP-middle b, or Schur-sum c).  Enumerate these
      structures up to iso for small k / small range; count distinct 'challenger
      shapes'.  Test which natural parameter gives exactly 56 (or 47).
  (b) Test the DIRECT BIJECTION challenger-shape <-> tournament on 6 vertices.
      A000568(6)=56.  We look for a structural map (relation signs <-> arc orient.)
      and HONESTLY report match / coincidence.
  (c) Cross-check the coarse (A2,A3,dmin) signature counts from the prior script and
      report whether the 47 (k=7) / 56 numbers are robust or signature artifacts.

EXACT arithmetic; integer relation search.
"""
from __future__ import annotations
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import Counter, defaultdict

try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

def banner(t): print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78)

# ----------------------------------------------------------------- A000568
A000568 = {1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880}

# =====================================================================
#  SUPPORT-3 PRIMITIVE RELATIONS of a set E (the leading binding layer)
# =====================================================================
# A support-3 primitive relation on positions (i<j<l) is coefs (p,q,r), all !=0,
# gcd=1, canonical sign (first nonzero >0), with p E_i + q E_j + r E_l = 0.
# With |coef|<=2 (the additive-energy box, matching the prior MDS script B=2) the
# only primitive solutions are, up to the canonical sign and ordering:
#   * (1,1,-1) family   -> Schur triple  E_a + E_b = E_c       (one elt is the SUM)
#   * (1,-2,1) family   -> 3-term AP      E_a + E_c = 2 E_b     (one elt is MIDDLE)
# We classify each support-3 triple by which role each vertex plays.

def support3_relations(E, B=2):
    """Return list of (triple_indices, coefs, kind, special_index).
       kind in {'AP','Schur'}; special_index = AP-middle or Schur-sum vertex.
       |coefs|<=B, primitive, support exactly 3."""
    E = [int(x) for x in E]; k = len(E)
    out = []
    seen = set()
    for (i, j, l) in itertools.combinations(range(k), 3):
        for coefs in itertools.product(range(-B, B+1), repeat=3):
            if any(c == 0 for c in coefs):
                continue
            if coefs[0]*E[i] + coefs[1]*E[j] + coefs[2]*E[l] != 0:
                continue
            g = reduce(gcd, [abs(c) for c in coefs])
            prim = tuple(c//g for c in coefs)
            if prim[0] < 0:
                prim = tuple(-c for c in prim)
            key = ((i, j, l), prim)
            if key in seen:
                continue
            seen.add(key)
            # classify (within |coef|<=2, gcd 1, primitive)
            absset = tuple(sorted(abs(c) for c in prim))
            if absset == (1, 1, 1):
                # Schur: one coef has opposite sign to the other two
                # find the 'sum' vertex = the one whose coef sign differs from majority
                signs = [1 if c > 0 else -1 for c in prim]
                # majority sign vertices are the addends, minority is the sum
                pos = [v for v, s in zip((i, j, l), signs) if s == 1]
                neg = [v for v, s in zip((i, j, l), signs) if s == -1]
                if len(neg) == 1:
                    special = neg[0]
                else:  # canonical sign forces first +, so the lone one is +
                    special = pos[0]
                out.append(((i, j, l), prim, 'Schur', special))
            elif absset == (1, 1, 2):
                # AP: the |coef|==2 vertex is the MIDDLE
                idx2 = [v for v, c in zip((i, j, l), prim) if abs(c) == 2][0]
                out.append(((i, j, l), prim, 'AP', idx2))
            else:
                out.append(((i, j, l), prim, 'other', None))
    return out

# =====================================================================
#  CHALLENGER SHAPE = isomorphism class of the support-3 structure.
#  We encode it as a labeled multi-hypergraph and canonicalize by trying all
#  vertex relabelings (k small).  The label of a hyperedge records, per vertex,
#  its ROLE (AP-middle / AP-end / Schur-sum / Schur-addend).  Two E's are the
#  same challenger shape iff their role-labeled support-3 hypergraphs are iso.
# =====================================================================
def role_hyperedges(E, B=2):
    """Return a frozenset of role-labeled hyperedges on vertex labels 0..k-1.
       Each hyperedge = frozenset of (vertex, role) pairs."""
    rels = support3_relations(E, B)
    edges = set()
    for (trip, prim, kind, special) in rels:
        if kind == 'other':
            continue
        roles = []
        for v in trip:
            if kind == 'AP':
                roles.append((v, 'APmid' if v == special else 'APend'))
            else:  # Schur
                roles.append((v, 'Ssum' if v == special else 'Sadd'))
        edges.add((kind, frozenset(roles)))
    return edges

def _vertex_invariant(edges):
    """Iterated colour-refinement vertex invariant (Weisfeiler-Leman-ish) for the
       role-labeled hypergraph.  Returns dict vertex -> stable colour (hashable)."""
    # initial colour = multiset of (kind, role) memberships
    init = defaultdict(list)
    for (kind, fs) in edges:
        for (v, r) in fs:
            init[v].append((kind, r))
    colour = {v: tuple(sorted(lst)) for v, lst in init.items()}
    for _ in range(len(colour) + 1):
        newc = {}
        for v in colour:
            # gather, per incident edge, the multiset of neighbour colours+roles
            sig = []
            for (kind, fs) in edges:
                members = {vv: rr for (vv, rr) in fs}
                if v in members:
                    others = tuple(sorted((colour[vv], rr) for (vv, rr) in fs if vv != v))
                    sig.append((kind, members[v], others))
            newc[v] = (colour[v], tuple(sorted(sig)))
        # compress to small ints to keep hashable size bounded
        order = {c: i for i, c in enumerate(sorted(set(newc.values())))}
        comp = {v: order[newc[v]] for v in newc}
        if all(comp[v] == colour[v] for v in comp) and set(comp.values()) == set(
                colour[v] for v in colour if isinstance(colour[v], int)):
            colour = comp
            break
        colour = comp
    return colour

def canon_shape(E, B=2):
    """Iso-invariant signature of the role-labeled support-3 hypergraph.  Uses
       colour-refinement (fast, complete for these small sparse hypergraphs) instead
       of brute n! relabeling.  Two E's with the same signature are (WL-)iso; we
       additionally append the sorted edge multiset under the refined colours, which
       is a near-complete invariant for our small cases."""
    E = list(E)
    edges = role_hyperedges(E, B)
    if not edges:
        return ('EMPTY', len(E))
    col = _vertex_invariant(edges)
    # edge signature under refined colours
    esig = tuple(sorted(
        (kind, tuple(sorted((col[v], r) for (v, r) in fs)))
        for (kind, fs) in edges))
    # colour-class multiset (extra invariant)
    csig = tuple(sorted(Counter(col.values()).items()))
    return ('SHAPE', csig, esig)

# =====================================================================
#  measS7 / corr (exact) -- reused from the MDS script for ranking
# =====================================================================
def stirling2(n, kk):
    if kk == 0: return 1 if n == 0 else 0
    S = [[0]*(kk+1) for _ in range(n+1)]; S[0][0] = 1
    for i in range(1, n+1):
        for j in range(1, min(i, kk)+1):
            S[i][j] = j*S[i-1][j] + S[i-1][j-1]
    return S[n][kk]

def iid_k(k, p=7):
    from math import factorial
    return F(factorial(p) * stirling2(k, p), p**k)

def measS7(E, p=7):
    E = sorted(set(int(e) for e in E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, p*e + 1): bps.add(F(m, p*e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0+x1)/2
        secs = {int(((e*xm) % 1)*p) for e in E}
        if len(secs) == p: total += x1-x0
    return total

def corr(E, p=7):
    return measS7(E, p) - iid_k(len(E), p)

# =====================================================================
#  PART A: enumerate challenger shapes over k-subsets of a window
# =====================================================================
def part_A():
    banner("PART A -- challenger shapes = iso classes of the support-3 relation hypergraph")
    print("Window: all k-subsets of {0,...,W}; shape = role-labeled support-3 hypergraph")
    print("up to relabeling.  |coef|<=2 (additive box).  Roles: AP{mid,end}, Schur{sum,add}.")
    print()
    print(f"{'k':>2} {'W':>3} {'#subsets':>9} {'#shapes':>8} {'#nonempty':>10}  A000568 note")
    results = {}
    for (k, W) in [(4,6),(5,7),(5,9),(6,8),(6,10),(6,12),(7,9),(7,11),(8,10),(8,12)]:
        shapes = set()
        nonempty = set()
        cnt = 0
        for E in itertools.combinations(range(0, W+1), k):
            cnt += 1
            sh = canon_shape(list(E))
            shapes.add(sh)
            if sh[0] not in ('EMPTY',):
                nonempty.add(sh)
        note = ''
        for vv, a in A000568.items():
            if len(shapes) == a: note += f" #shapes=A000568({vv})={a}!"
            if len(nonempty) == a: note += f" #nonempty=A000568({vv})={a}!"
        if len(shapes) == 56 or len(nonempty) == 56 or len(shapes) == 47 or len(nonempty) == 47:
            note += '  <-- 56/47 HIT'
        print(f"{k:>2} {W:>3} {cnt:>9} {len(shapes):>8} {len(nonempty):>10}  {note}")
        results[(k,W)] = (cnt, len(shapes), len(nonempty))
    return results

# =====================================================================
#  PART B: the tournament bijection test.
#  A tournament on 6 labeled vertices = an orientation of each of the 15 edges of
#  K_6; 56 = #iso classes (A000568(6)).  We look for a NATURAL 6-vertex object in
#  the support-3 structure whose iso classes number 56, then test a structural map.
#  Candidate object: '6 elements with a SINGLE support-3 triangle/3-cycle of
#  relations' OR 'orientations of a triple-structure'.  We test the most direct:
#  do the support-3 SCHUR triples on a 6-set define a tournament (orient each
#  pair by which is the 'larger' in the additive triple) and does the iso-class
#  count match 56?
# =====================================================================
def tournaments_on_n(n):
    """All tournaments on n labeled vertices as adjacency (i<j -> bit). Return
       canonical-form set (iso classes) by brute relabel.  n<=6 ok-ish (n=6: 32768
       tournaments * 720 perms ~ 23M; we instead count via known A000568)."""
    # We do NOT brute n=6 (too big); use orbit counting only for n<=5 to validate.
    pairs = list(itertools.combinations(range(n), 2))
    idx = {p: i for i, p in enumerate(pairs)}
    def canon(bits):
        best = None
        for perm in itertools.permutations(range(n)):
            b = 0
            for (a, c) in pairs:
                pa, pc = perm[a], perm[c]
                # orientation a->c if bit set; relabel
                if pa < pc:
                    orient = (bits >> idx[(a, c)]) & 1
                    if orient: b |= 1 << idx[(pa, pc)]
                else:
                    orient = (bits >> idx[(a, c)]) & 1
                    if not orient: b |= 1 << idx[(pc, pa)]
            if best is None or b < best: best = b
        return best
    classes = set()
    for bits in range(1 << len(pairs)):
        classes.add(canon(bits))
    return classes

def part_B():
    banner("PART B -- the tournament-on-6 bijection test")
    print("A000568(6) = 56 = # iso classes of tournaments on 6 vertices.")
    # validate the brute tournament counter on small n
    for n in (3, 4, 5):
        c = len(tournaments_on_n(n))
        print(f"  brute tournament iso classes n={n}: {c}  (A000568={A000568[n]}) "
              f"{'OK' if c == A000568[n] else 'MISMATCH'}")
    print()
    print("Candidate bijection 1: '6 elements, support-3 structure = a tournament'.")
    print("  For a 6-set E, orient each pair {a,b} by an additive rule and ask whether")
    print("  the resulting structures, up to iso, number 56.")
    print("  PROBLEM (honest): a tournament needs a DEFINED orientation on ALL C(6,2)=15")
    print("  pairs, but support-3 relations only touch SOME triples -- E's support-3")
    print("  hypergraph is NOT a tournament (it is a 3-uniform partial hypergraph).")
    print("  So a literal 'shape <-> tournament' map cannot be a bijection by type alone.")
    print()
    print("  We therefore test the WEAKER claim: does SOME natural support-3 count = 56?")
    # how many distinct role-labeled support-3 hypergraphs arise on EXACTLY 6 vertices
    # all touched, with at least one edge, over a generous window?
    for W in (10, 12):
        shapes6 = set()
        for E in itertools.combinations(range(0, W+1), 6):
            edges = role_hyperedges(list(E))
            if not edges:
                continue
            used = {v for (_, fs) in edges for (v, _) in fs}
            if len(used) != 6:
                continue  # all 6 vertices participate
            sh = canon_shape(list(E))
            shapes6.add(sh)
        note = '  <-- 56!' if len(shapes6) == 56 else ('  <-- 47!' if len(shapes6) == 47 else '')
        print(f"  6-sets in [0,{W}], all 6 vertices in support-3 edges: "
              f"{len(shapes6)} distinct shapes{note}")
    print()
    print("Candidate bijection 2: the SINGLE-TRIANGLE challenger.  A challenger shape with")
    print("exactly ONE support-3 hyperedge has only finitely many ROLE TYPES:")
    print("  AP triple: 1 way to label (middle distinguished, 2 ends symmetric) ;")
    print("  Schur triple: 1 way (sum distinguished, 2 addends symmetric).")
    print("=> only 2 single-triangle shapes.  Not 56.  The '56' is NOT a single triple.")
    print()
    print("Candidate bijection 3: the support-3 COEFFICIENT SIGNATURES with |coef|<=3.")
    print("Enumerate primitive (p,q,r), gcd=1, canonical sign, that can vanish on 3")
    print("DISTINCT positive reals -- count them by |coef| bound:")
    for B in (2, 3, 4):
        sigs = set()
        for coefs in itertools.product(range(-B, B+1), repeat=3):
            if any(c == 0 for c in coefs):
                continue
            g = reduce(gcd, [abs(c) for c in coefs])
            prim = tuple(c//g for c in coefs)
            if prim[0] < 0:
                prim = tuple(-c for c in prim)
            # realizable on distinct positive E iff signs mixed (not all same sign)
            s = [1 if c > 0 else -1 for c in prim]
            if len(set(s)) == 1:
                continue  # all-same-sign cannot vanish on positives
            sigs.add(prim)
        # up to permutation of the 3 positions (unordered triple):
        upto_perm = set()
        for prim in sigs:
            best = min(itertools.permutations(prim))
            # re-canonicalize sign
            if best[0] < 0 or (best[0] == 0):
                pass
            upto_perm.add(best)
        print(f"  |coef|<= {B}: {len(sigs)} ordered signatures, "
              f"{len(upto_perm)} unordered  {'<-- 56!' if len(sigs)==56 or len(upto_perm)==56 else ''}")
    return

# =====================================================================
#  PART C: reproduce / refine the prior coarse (A2,A3,dmin) signature counts
# =====================================================================
def support_spectrum(E, B=2, max_support=3):
    E = [int(e) for e in E]; k = len(E)
    counts = {s: 0 for s in range(2, max_support+1)}
    seen = set()
    for s in range(2, max_support+1):
        for combo in itertools.combinations(range(k), s):
            for coefs in itertools.product(range(-B, B+1), repeat=s):
                if any(c == 0 for c in coefs): continue
                if sum(c*E[i] for c, i in zip(coefs, combo)) != 0: continue
                g = reduce(gcd, [abs(c) for c in coefs])
                prim = tuple(c//g for c in coefs)
                if prim[0] < 0: prim = tuple(-c for c in prim)
                key = (combo, prim)
                if key in seen: continue
                seen.add(key)
                counts[s] += 1
    nz = [s for s in counts if counts[s] > 0]
    return counts, (min(nz) if nz else None)

def part_C():
    banner("PART C -- the prior coarse (A2,A3,dmin) signature counts: robust or artifact?")
    print("Prior MDS script reported, over k-subsets of [0,k+2]:")
    print("  k=7 -> 47 distinct (A2,A3,dmin) shapes; k=8 -> 55.  Are 47/56 meaningful")
    print("  or signature coincidences?  We recompute AND give the FULL role-hypergraph")
    print("  iso-shape count for the SAME window to compare resolutions.")
    print()
    print(f"{'k':>2} {'window':>10} {'#sub':>6} {'#(A2,A3,dmin)':>14} {'#full-iso-shape':>16}")
    for k in range(4, 9):
        W = k + 2
        coarse = set(); full = set()
        cnt = 0
        for E in itertools.combinations(range(0, W+1), k):
            cnt += 1
            counts, dmin = support_spectrum(list(E), B=2, max_support=3)
            coarse.add((counts.get(2, 0), counts.get(3, 0), dmin))
            full.add(canon_shape(list(E)))
        flag = ''
        if len(coarse) in (47, 56) or len(full) in (47, 56):
            flag = '  <-- 47/56'
        print(f"{k:>2} [0,{W:>2}] {cnt:>9} {len(coarse):>14} {len(full):>16}{flag}")

def main():
    print("LRC(14) OPEN-Q-108  THREAD 3: the '56 challenger shapes' probe (kps-wf2)")
    part_A()
    part_B()
    part_C()
    print("\nDONE.")

if __name__ == "__main__":
    main()
