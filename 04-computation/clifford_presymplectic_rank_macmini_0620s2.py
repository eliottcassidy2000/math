#!/usr/bin/env python3
"""
Pin down the directed-cycle alternating form's RANK (mac-mini-2026-06-20-S2).
TEST B found: form is ALTERNATING but DEGENERATE with gram_rank = 0,2,4,6 = 2(n-2)
on a cycle-basis of dimension C(n-1,2). Confirm rank is basis-independent and
identify the radical. If rank = 2(n-2) robustly, the 'symplectic part' of the OCF
cycle space has dimension 2(n-2) -> a candidate 'effective number of magic qubits'.

The form: B(C,C') = sum_e [C,C' use edge e in OPPOSITE directions] mod 2.
We compute it as a bilinear form on the WHOLE cycle space (all 2^{C(n-1,2)} even
subgraphs given a reference orientation), via the matrix in the standard arc basis,
so the rank is intrinsic (rank of the Gram over a full basis).
"""
from itertools import combinations

def rref_gf2(rows):
    piv = []
    for r in rows:
        cur = r
        for p in piv:
            if cur ^ p < cur: cur ^= p
        if cur:
            piv.append(cur); piv.sort(reverse=True)
    return piv

def Kn(n):
    edges = list(combinations(range(n), 2))
    return edges, {e: i for i, e in enumerate(edges)}

def cycle_basis_dir(n):
    """Fundamental directed cycles wrt star tree at 0. Each is a dict edge->dirbit
       (0 = i->j for i<j, 1 = j->i). Triangle 0->i->j->0 for non-tree edge (i,j)."""
    edges, idx = Kn(n)
    B = []
    for (i, j) in edges:
        if i == 0: continue
        # cycle 0 -> i -> j -> 0  (i<j, both>=1)
        d = {}
        # 0->i : edge (0,i), 0<i so dir 0
        d[idx[(0, i)]] = 0
        # i->j : edge (i,j), i<j so dir 0
        d[idx[(i, j)]] = 0
        # j->0 : edge (0,j), 0<j so this is j->0 = reverse => dir 1
        d[idx[(0, j)]] = 1
        B.append(d)
    return B, edges, idx

def form(C, Cp):
    s = 0
    for e, b in C.items():
        if e in Cp and Cp[e] != b:
            s ^= 1
    return s

def gram_rank(n):
    B, edges, idx = cycle_basis_dir(n)
    k = len(B)
    rows = []
    for i in range(k):
        r = 0
        for j in range(k):
            if form(B[i], B[j]): r |= 1 << j
        rows.append(r)
    return k, len(rref_gf2(rows))

if __name__ == "__main__":
    print("directed-cycle alternating form: dim(cycle) vs rank vs radical")
    print(f"{'n':>3} {'dimCyc=C(n-1,2)':>16} {'rank':>6} {'2(n-2)':>8} {'radical':>8} {'rank=2(n-2)?':>13}")
    for n in range(3, 9):
        k, rk = gram_rank(n)
        print(f"{n:>3} {k:>16} {rk:>6} {2*(n-2):>8} {k-rk:>8} {str(rk==2*(n-2)):>13}")
