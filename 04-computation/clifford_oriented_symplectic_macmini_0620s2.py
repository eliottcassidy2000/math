#!/usr/bin/env python3
"""
FOLLOW-UP: where IS the symplectic structure, if cut/cycle is merely orthogonal?
(mac-mini-2026-06-20-S2, companion to clifford_cutcycle_symplectic_*)

PART (i) showed: CUT = CYCLE^perp under the STANDARD form (orthogonal complements)
but NEITHER is isotropic => NOT a Lagrangian pair => the literal
"stabilizer = Lagrangian" reading of cut/cycle is FALSE.

The stabilizer formalism's symplectic space is GF(2)^{2n} = (X-part, Z-part) with
the form  B((a,b),(a',b')) = a.b' + a'.b  (commutation of Pauli operators).
That requires TWO copies. The tournament has ONE copy (orientation in GF(2)^E).
So the natural question: is there a SECOND copy hiding?

CANDIDATE: the cut-space (scores) and cycle-space (cycles) are dual = the X-part
and Z-part. Then the symplectic form between (cut-coords) and (cycle-coords) is the
NATURAL PAIRING (cut . cycle), which IS the standard form restricted to the
decomposition. Two Lagrangians of GF(2)^E (+) GF(2)^E ... but that's the trivial
doubling. The honest test:

  TEST A (confirm): cut-space and cycle-space are MacWilliams duals (exact, on
          nonzero support). This is the REAL coding-theory content. Cut = a
          [E, n-1] code; Cycle = its dual [E, E-n+1]. Verify W_cyc = MacW(W_cut)
          AND that this is the EVEN-GRAPH enumerator (cycle space = even subgraphs).

  TEST B (symplectic, oriented): On the ORIENTED arc space, define the OCF-relevant
          antisymmetric pairing. A directed 3-cycle is a generator; two directed
          cycles 'link' via shared arcs with sign. Build the directed-cycle
          intersection form  B(C,C') = (# arcs of C that are reversed in C') mod 2
          and ask: is it ANTISYMMETRIC / SYMPLECTIC (nondegenerate alternating)?
          If yes, the cycle space carries a genuine symplectic form -> a real
          stabilizer/Lagrangian story (just on cycle space, not cut+cycle).

  TEST C (magic): the THM-555 wall says c3 is the LAST score(=cut)-determined
          datum; c5, alpha_2 need cycle space. Frame as: define the 'cut-rank' of
          a tournament's OCF data and the 'cycle-rank' (=magic). Check the
          DIMENSION of the cut-determined OCF observables vs total -> a literal
          'stabilizer fraction'.

stdlib, exact.
"""
from itertools import combinations
from math import comb

# ---------- GF(2) ----------
def rref_gf2(rows):
    piv = []
    for r in rows:
        cur = r
        for p in piv:
            if cur ^ p < cur:
                cur ^= p
        if cur:
            piv.append(cur); piv.sort(reverse=True)
    return piv

def span(basis):
    out = {0}
    for b in basis:
        out |= {x ^ b for x in out}
    return out

def dot2(a, b):
    return bin(a & b).count("1") & 1

# ---------- K_n spaces ----------
def Kn(n):
    edges = list(combinations(range(n), 2))
    idx = {e: i for i, e in enumerate(edges)}
    return edges, idx

def cut_basis(n):
    edges, idx = Kn(n)
    rows = []
    for v in range(n):
        m = 0
        for e in edges:
            if v in e:
                m |= 1 << idx[e]
        rows.append(m)
    return rref_gf2(rows), edges, idx

def cycle_basis(n):
    edges, idx = Kn(n)
    rows = []
    for (i, j) in edges:
        if i == 0:
            continue
        m = (1 << idx[(i, j)]) ^ (1 << idx[(0, i)]) ^ (1 << idx[(0, j)])
        rows.append(m)
    return rref_gf2(rows), edges, idx

def wenum(basis):
    c = {}
    for v in span(basis):
        w = bin(v).count("1"); c[w] = c.get(w, 0) + 1
    return c

def macw(counts, E):
    size = sum(counts.values()); out = {}
    for j in range(E + 1):
        s = 0
        for w, aw in counts.items():
            s += aw * sum((-1)**l * comb(w, l) * comb(E - w, j - l) for l in range(j + 1))
        assert s % size == 0
        if s: out[j] = s // size
    return out

# ---------- even-graph enumerator (cycle space = even subgraphs) ----------
def even_subgraph_enum(n):
    """Directly enumerate edge subsets of K_n with all degrees even; weight=#edges.
       This IS the cycle space; confirms cycle_basis spans even subgraphs."""
    edges, idx = Kn(n)
    E = len(edges)
    c = {}
    for s in range(1 << E):
        deg = [0]*n
        for k, (i, j) in enumerate(edges):
            if (s >> k) & 1:
                deg[i] ^= 1; deg[j] ^= 1
        if not any(deg):
            w = bin(s).count("1"); c[w] = c.get(w, 0) + 1
    return c

# ---------- TEST B: oriented directed-cycle intersection form ----------
def directed_cycle_form(n):
    """Generators: all directed 3-cycles of K_n (choose 3 vertices, 2 orientations).
       Represent each as a vector in GF(2)^{2E}? No -- use the arc model:
       arc set A = ordered pairs; a directed 3-cycle (a->b->c->a).
       Pairing B(C,C') = sum over edges {u,v} of [C uses u->v] AND [C' uses v->u]
       i.e. count of edges where the two cycles DISAGREE in direction, mod 2.
       Build the Gram matrix over a basis of directed cycles and test alternating /
       rank (symplectic <=> alternating + full even rank)."""
    edges, idx = Kn(n)
    E = len(edges)
    # represent a directed cycle as a dict edge-index -> direction bit (which endpoint is tail)
    # direction bit: for edge (i,j) with i<j, bit=0 means i->j, bit=1 means j->i.
    def tri_cycles():
        cyc = []
        for a, b, c in combinations(range(n), 3):
            for order in [(a, b, c), (a, c, b)]:
                x, y, z = order
                arcs = [(x, y), (y, z), (z, x)]
                d = {}
                for (u, v) in arcs:
                    e = (u, v) if u < v else (v, u)
                    bit = 0 if u < v else 1
                    d[idx[e]] = bit
                cyc.append(d)
        return cyc
    cyc = tri_cycles()
    def pair(C, Cp):
        s = 0
        for e, b in C.items():
            if e in Cp and Cp[e] != b:
                s ^= 1
        return s
    # Gram on first few independent generators (use a GF2-independent subset by support mask)
    # build support masks (undirected) to pick an independent generating set of cycle space
    masks = []
    chosen = []
    basis = []
    for C in cyc:
        m = 0
        for e in C: m |= 1 << e
        cur = m
        for p in basis:
            if cur ^ p < cur: cur ^= p
        if cur:
            basis.append(cur); basis.sort(reverse=True)
            chosen.append(C); masks.append(m)
    G = [[pair(chosen[i], chosen[j]) for j in range(len(chosen))] for i in range(len(chosen))]
    # tests
    sym = all(G[i][j] == G[j][i] for i in range(len(G)) for j in range(len(G)))
    alt = all(G[i][i] == 0 for i in range(len(G)))
    # rank over GF2
    rows = []
    for i in range(len(G)):
        r = 0
        for j in range(len(G)):
            if G[i][j]: r |= 1 << j
        rows.append(r)
    rk = len(rref_gf2(rows))
    return {"n": n, "ncyc_basis": len(chosen), "form_symmetric": sym,
            "form_alternating": alt, "gram_rank": rk,
            "nondegenerate": rk == len(chosen)}

if __name__ == "__main__":
    print("=" * 70)
    print("  TEST A: cut/cycle MacWilliams duality (exact, nonzero support)")
    print("          and cycle space == even-subgraph (even-graph) enumerator")
    print("=" * 70)
    for n in range(3, 7):
        (cb, edges, idx) = cut_basis(n)
        (cyb, _, _) = cycle_basis(n)
        E = len(edges)
        Wc = wenum(cb); Wcy = wenum(cyb)
        mw = macw(Wc, E)
        nz_cy = {k: v for k, v in Wcy.items() if v}
        even = even_subgraph_enum(n)
        print(f"\n n={n} E={E}")
        print(f"   MacWilliams(W_cut) == W_cycle (nonzero)?   {mw == nz_cy}")
        print(f"   W_cycle == even-subgraph enumerator?       {nz_cy == even}")
        print(f"   (cycle space dim = {len(cyb)} = C(n,2)-n+1 = {E-n+1})")

    print("\n" + "=" * 70)
    print("  TEST B: oriented directed-3-cycle intersection form -- symplectic?")
    print("=" * 70)
    for n in range(3, 7):
        r = directed_cycle_form(n)
        print(f" n={n}: cycle-basis size={r['ncyc_basis']}, "
              f"symmetric={r['form_symmetric']}, alternating(diag=0)={r['form_alternating']}, "
              f"gram_rank={r['gram_rank']}, nondegenerate={r['nondegenerate']}")
