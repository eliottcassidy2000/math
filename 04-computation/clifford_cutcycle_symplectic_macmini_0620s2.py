#!/usr/bin/env python3
"""
CLIFFORD / STABILIZER  <->  GF(2) CUT/CYCLE   (mac-mini-2026-06-20-S2)

THREAD: Clifford gates {H,S,CNOT} are classically simulable (Gottesman-Knill);
the stabilizer formalism is SYMPLECTIC over GF(2). The repo's tiling space is a
GF(2) hypercube with a CUT (+) CYCLE split. The HYPOTHESIS to test, brutally:

(i)  Is cut/cycle GENUINELY SYMPLECTIC? i.e. are CUT and CYCLE a pair of dual
     LAGRANGIAN subspaces of GF(2)^E under the standard symplectic/bilinear form?
     A Lagrangian L satisfies L = L^perp (self-orthogonal AND maximal, dim=E/2).
     For the stabilizer formalism the stabilizer group <-> a Lagrangian in the
     symplectic phase space GF(2)^{2n}. "Cut cheap / cycle dear" would mirror
     "Clifford cheap / magic dear" ONLY IF the two spaces are symplectic duals.

     CONCRETE FALSIFIABLE CLAIMS:
       (a) dim(CUT) + dim(CYCLE) = |E| = C(n,2)   [orthogonal complements]
       (b) CUT = CYCLE^perp under  B(x,y)=sum x_e y_e mod 2     [duality]
       (c) Is either space ISOTROPIC (self-orthogonal, x.x'=0 for all x,x' in L)?
           A Lagrangian requires BOTH self-orthogonal AND dim=E/2.
       The stabilizer-symplectic analogy is TRUE only if CUT (or CYCLE) is a
       Lagrangian. If neither is, the analogy is DEAD: cut/cycle is the standard
       ORTHOGONAL (not symplectic) homology split, NOT the stabilizer one.

(ii) Is H a stabilizer WEIGHT ENUMERATOR obeying a MacWilliams identity between
     G_n and the even-graph dual E_n?  The cycle space of K_n = the even-degree
     edge sets = exactly the EVEN GRAPHS. So E_n IS the cycle space. Test the
     MacWilliams transform between the cut-space weight enumerator and the
     cycle-space (even-graph) weight enumerator.

(iii) Magic monotone analog: does a "non-Clifford rank" of a tournament (its
     projection onto cycle space = c3 and beyond) equal a cycle-space invariant?
     Test: is the LAST score-determined datum c3 exactly the cut-space content,
     and everything beyond it strictly cycle-space (the THM-555 wall as a
     symplectic-rank statement)?

stdlib only, exact arithmetic.
"""
from itertools import combinations, product
from fractions import Fraction

# ----------------------------------------------------------------------
# GF(2) linear algebra
# ----------------------------------------------------------------------
def rref_gf2(rows, ncols):
    """Return (basis as list of int-bitmasks, rank). rows = list of int bitmasks."""
    rows = [r for r in rows if r]
    basis = []
    for r in rows:
        cur = r
        for b in basis:
            cur = min(cur, cur ^ b)
        if cur:
            basis.append(cur)
            basis.sort(reverse=True)
    # re-reduce to a clean echelon basis
    piv = []
    for r in rows:
        cur = r
        for p in piv:
            if cur ^ p < cur:
                cur ^= p
        if cur:
            piv.append(cur)
            piv.sort(reverse=True)
    return piv, len(piv)

def span(basis):
    """All vectors in the span of a list of int bitmasks (small dims only)."""
    out = {0}
    for b in basis:
        out |= {x ^ b for x in out}
    return out

def dot2(a, b):
    return bin(a & b).count("1") & 1

# ----------------------------------------------------------------------
# K_n edge/cut/cycle spaces over GF(2)
# ----------------------------------------------------------------------
def Kn_edges(n):
    return list(combinations(range(n), 2))  # E = C(n,2) edges

def cut_space_basis(n):
    """Cut space of K_n over GF(2) = row space of the vertex-edge incidence matrix.
       Each vertex v -> the set of edges incident to v (a 'star' = elementary cut).
       Dim = n-1."""
    edges = Kn_edges(n)
    idx = {e: i for i, e in enumerate(edges)}
    rows = []
    for v in range(n):
        mask = 0
        for e in edges:
            if v in e:
                mask |= (1 << idx[e])
        rows.append(mask)
    basis, rk = rref_gf2(rows, len(edges))
    return basis, rk, edges, idx

def cycle_space_basis(n):
    """Cycle space = kernel of incidence matrix = even-degree edge subsets.
       Dim = |E| - (n-1) = C(n,2)-n+1.  Fundamental cycles wrt a spanning tree.
       Use the star tree at vertex 0: tree edges (0,k). Fundamental cycle for a
       non-tree edge (i,j): (0,i)+(0,j)+(i,j)."""
    edges = Kn_edges(n)
    idx = {e: i for i, e in enumerate(edges)}
    tree = {(0, k) for k in range(1, n)}
    rows = []
    for (i, j) in edges:
        if (i, j) in tree:
            continue
        # i,j both >=1 (since (0,k) are tree edges). triangle 0-i-j-0
        mask = (1 << idx[(i, j)])
        mask ^= (1 << idx[tuple(sorted((0, i)))])
        mask ^= (1 << idx[tuple(sorted((0, j)))])
        rows.append(mask)
    basis, rk = rref_gf2(rows, len(edges))
    return basis, rk, edges, idx

# ----------------------------------------------------------------------
# PART (i): IS CUT/CYCLE SYMPLECTIC?
# ----------------------------------------------------------------------
def test_symplectic(n):
    cutB, cutR, edges, idx = cut_space_basis(n)
    cycB, cycR, _, _ = cycle_space_basis(n)
    E = len(edges)
    res = {
        "n": n, "E": E, "dim_cut": cutR, "dim_cycle": cycR,
        "dim_sum": cutR + cycR,
    }
    # (a) complementary dims
    res["complementary_dims"] = (cutR + cycR == E)
    # (b) duality: is CUT = CYCLE^perp ? equivalent to: cut perp cycle AND dims add to E
    # check every basis pair orthogonal
    cross_orth = all(dot2(a, b) == 0 for a in cutB for b in cycB)
    res["cut_perp_cycle"] = cross_orth
    res["cut_eq_cycleperp"] = cross_orth and (cutR + cycR == E)
    # (c) isotropy of each space under standard dot form
    cut_iso = all(dot2(a, b) == 0 for a in cutB for b in cutB)
    cyc_iso = all(dot2(a, b) == 0 for a in cycB for b in cycB)
    res["cut_isotropic"] = cut_iso
    res["cycle_isotropic"] = cyc_iso
    res["cut_is_lagrangian"] = cut_iso and (cutR == E // 2) and (E % 2 == 0)
    res["cycle_is_lagrangian"] = cyc_iso and (cycR == E // 2) and (E % 2 == 0)
    # self-dual?  cut = cut^perp ?  (i.e. cut isotropic & dim = E/2)
    res["cut_selfdual"] = cut_iso and (cutR == E - cutR)
    res["cycle_selfdual"] = cyc_iso and (cycR == E - cycR)
    return res

# ----------------------------------------------------------------------
# PART (ii): MacWilliams between cut-space and cycle-space weight enumerators
# ----------------------------------------------------------------------
def weight_enumerator(basis, E):
    """W_C(x,y) homogeneous: sum_{c in C} x^{E-wt} y^{wt}, returned as dict wt->count."""
    counts = {}
    for c in span(basis):
        w = bin(c).count("1")
        counts[w] = counts.get(w, 0) + 1
    return counts

def macwilliams(counts, E):
    """Apply MacWilliams transform to a weight distribution (dict wt->count).
       W_{C^perp}(x,y) = (1/|C|) W_C(x+y, x-y).
       Returns dict wt->count for the dual (must be integers if C is a code)."""
    size = sum(counts.values())
    # build W_C(x,y) coefficients a_w; compute b_j = (1/size) sum_w a_w * Krawtchouk
    # Krawtchouk K_j(w; E) = sum_l (-1)^l C(w,l) C(E-w, j-l)
    from math import comb
    out = {}
    for j in range(E + 1):
        s = 0
        for w, aw in counts.items():
            Kjw = sum((-1)**l * comb(w, l) * comb(E - w, j - l) for l in range(0, j + 1))
            s += aw * Kjw
        # b_j = s / size
        assert s % size == 0, f"non-integer dual coeff at j={j}: {s}/{size}"
        out[j] = s // size
    return out

def test_macwilliams(n):
    cutB, cutR, edges, idx = cut_space_basis(n)
    cycB, cycR, _, _ = cycle_space_basis(n)
    E = len(edges)
    Wcut = weight_enumerator(cutB, E)
    Wcyc = weight_enumerator(cycB, E)
    dual_of_cut = macwilliams(Wcut, E)
    # compare on nonzero support only (explicit-zero keys are not a mismatch)
    nz_dual = {k: v for k, v in dual_of_cut.items() if v != 0}
    nz_cyc = {k: v for k, v in Wcyc.items() if v != 0}
    match = (nz_dual == nz_cyc)
    return {"n": n, "E": E,
            "W_cut": dict(sorted(Wcut.items())),
            "W_cycle": dict(sorted(Wcyc.items())),
            "MacWilliams(cut)": dict(sorted(dual_of_cut.items())),
            "cut_dual_equals_cycle": match}

# ----------------------------------------------------------------------
# PART (iii): score = cut content ; c3 & beyond = cycle content (THM-555 wall)
# Map a tournament (arc orientation in {0,1}^E on K_n) to:
#   - its score vector  (cut/syndrome content)
#   - its 3-cycle count c3 (claimed last score-determined datum)
# A tournament is a point; the SCORE is a linear (cut-space) functional of arcs.
# Test: does the cut-space projection of the orientation determine the score, and
# does c3 fail to be determined past the cut projection (need cycle data)?
# ----------------------------------------------------------------------
def score_and_c3(orient, n, edges, idx):
    """orient: bitmask, bit set => edge (i,j) oriented i->j (i<j); else j->i."""
    # adjacency: arc[a][b]=1 if a->b
    out = [[0]*n for _ in range(n)]
    for (i, j) in edges:
        if (orient >> idx[(i, j)]) & 1:
            out[i][j] = 1
        else:
            out[j][i] = 1
    score = tuple(sorted(sum(out[v]) for v in range(n)))
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        # count cyclic triangles among a,b,c
        # edges: a-b,a-c,b-c
        e = [(a, b), (a, c), (b, c)]
        # cyclic iff each vertex has out-degree 1 within the triple
        deg = {a: 0, b: 0, c: 0}
        for (u, v) in e:
            if out[u][v]:
                deg[u] += 1
            else:
                deg[v] += 1
        if set(deg.values()) == {1}:
            c3 += 1
    return score, c3

def test_cut_cycle_content(n):
    """Enumerate all 2^E tournaments (small n). Group by score. Within a fixed
       score, is c3 constant (score-determined) ? Beyond c3, is c5 NOT constant?"""
    edges = Kn_edges(n)
    idx = {e: i for i, e in enumerate(edges)}
    E = len(edges)
    by_score = {}
    for orient in range(1 << E):
        sc, c3 = score_and_c3(orient, n, edges, idx)
        by_score.setdefault(sc, set()).add(c3)
    c3_determined = all(len(v) == 1 for v in by_score.values())
    return {"n": n, "num_scores": len(by_score),
            "c3_score_determined": c3_determined,
            "max_c3_spread_within_score": max(len(v) for v in by_score.values())}

# ----------------------------------------------------------------------
# RUN
# ----------------------------------------------------------------------
if __name__ == "__main__":
    print("=" * 72)
    print("  PART (i): IS CUT/CYCLE GENUINELY SYMPLECTIC (Lagrangian pair)?")
    print("=" * 72)
    for n in range(3, 8):
        r = test_symplectic(n)
        print(f"\n n={n}  |E|=C(n,2)={r['E']}")
        print(f"   dim CUT = {r['dim_cut']} (=n-1),  dim CYCLE = {r['dim_cycle']} (=C(n,2)-n+1)")
        print(f"   dim sum = E ?               {r['complementary_dims']}")
        print(f"   CUT _|_ CYCLE (std form)?   {r['cut_perp_cycle']}")
        print(f"   CUT = CYCLE^perp ?          {r['cut_eq_cycleperp']}")
        print(f"   CUT isotropic (self-orth)?  {r['cut_isotropic']}")
        print(f"   CYCLE isotropic?            {r['cycle_isotropic']}")
        print(f"   CUT is LAGRANGIAN?          {r['cut_is_lagrangian']}")
        print(f"   CYCLE is LAGRANGIAN?        {r['cycle_is_lagrangian']}")
        print(f"   CUT self-dual (=cut^perp)?  {r['cut_selfdual']}")
        print(f"   CYCLE self-dual?            {r['cycle_selfdual']}")

    print("\n" + "=" * 72)
    print("  PART (ii): MacWilliams duality  cut-space <-> cycle-space")
    print("=" * 72)
    for n in range(3, 8):
        r = test_macwilliams(n)
        print(f"\n n={n}  E={r['E']}")
        print(f"   W_cut   = {r['W_cut']}")
        print(f"   W_cycle = {r['W_cycle']}")
        print(f"   MacWilliams(W_cut) = {r['MacWilliams(cut)']}")
        print(f"   MacWilliams(cut) == W_cycle ?   {r['cut_dual_equals_cycle']}")

    print("\n" + "=" * 72)
    print("  PART (iii): cut = score (Clifford/simulable), c3 = last cut datum")
    print("=" * 72)
    for n in range(3, 7):
        r = test_cut_cycle_content(n)
        print(f" n={n}: #score classes={r['num_scores']}, "
              f"c3 score-determined? {r['c3_score_determined']} "
              f"(max c3-spread within a score = {r['max_c3_spread_within_score']})")
