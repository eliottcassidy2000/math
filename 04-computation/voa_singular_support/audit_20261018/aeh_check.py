#!/usr/bin/env python3
"""Literature cross-check for THM-4671's novelty claim.

Andrews-van Ekeren-Heluani (arXiv:2005.10769, sec. 2.4, second Lemma) state, for Vir_{3,p'} (p' >= 4),
that the lowest piece of ker(pi) is spanned by the degree-(2p'+1) vector
    b(p') = (9-2p')/(3(p'-2)) L_{-5} L_{-2}^{p'-2} + L_{-4} L_{-3} L_{-2}^{p'-3},
and for Ising they exhibit the combination
    w + 256/429 L_{-3}v - 64/429 L_{-1}L_{-2}v - 31/286 L_{-1}^3 v  in lower PBW length,
    w = L_{-5}L_{-2}L_{-2}|0> + 6 L_{-4}L_{-3}L_{-2}|0>.
Here we check (exactly) that
  (1) b(p') == (const) * b_N  modulo J = (x_0^N)_d  with N = p'-1, for p' = 4,5,7,8,10,11;
  (2) b(p') lies in grI for those models (so AEH's element and THM-4671's b_N are the same class);
  (3) AEH's Ising combination equals  w - (64/429) u - (31/286) L_{-1}^3 v  with THM-4671's u
      (v normalized by coefficient 1 on L_{-2}^3), i.e. AEH already used the same lift.
"""
import flint
from alin import rref_vectors, rank, in_span, leading_forms
from virmod import VirVac, partitions_min2, central_charge, poly_mul, poly_d, symbol
from alin import nullspace

Q = flint.fmpq


def sing_vec(vir, N):
    B = partitions_min2(2 * N)
    eq = {}
    for n in (1, 2):
        for b in B:
            for m, cc in vir.act(n, b).items():
                eq.setdefault((n, m), {})[b] = cc
    ker = nullspace(list(eq.values()), B)
    assert len(ker) == 1
    v = ker[0]
    lead = v[(2,) * N]
    return {m: c / lead for m, c in v.items()}


for pp in (4, 5, 7, 8, 10, 11):
    p, q = 3, pp
    N = pp - 1
    c = central_charge(p, q)
    vir = VirVac(c, Q(1))
    v = sing_vec(vir, N)
    W = 2 * N
    # I in weights W..W+3
    I = {W: [v]}
    for D in range(W + 1, W + 4):
        gens = [vir.actv(-1, w) for w in I[D - 1]]
        if D - 2 >= W:
            gens += [vir.actv(-2, w) for w in I[D - 2]]
        I[D] = [row for _, row in rref_vectors(gens, partitions_min2(D))]
    D = W + 3
    B = partitions_min2(D)
    gr = leading_forms(I[D], B, degree=len)
    m1 = (5,) + (2,) * (N - 1)
    m2 = (4, 3) + (2,) * (N - 2)
    m3 = (3, 3, 3) + (2,) * (N - 3)
    baeh = {m1: Q(9 - 2 * pp, 3 * (pp - 2)), m2: Q(1)}
    bN = {m2: Q(1), m3: Q(2 * N - 7, 12)}
    f = {(2,) * N: Q(1)}
    d3f = poly_d(poly_d(poly_d(f)))
    Jg = []
    for s in range(0, 4):
        ds = dict(f)
        for _ in range(s):
            ds = poly_d(ds)
        for mm in partitions_min2(3 - s):
            Jg.append(poly_mul({mm: Q(1)}, ds))
    same_class = (rank(Jg + [bN, baeh], B) == rank(Jg + [bN], B))
    print(f"Vir(3,{pp}) N={N}: AEH b(p') in grI: {in_span(baeh, gr, B)};  b(p') in J: {in_span(baeh, Jg, B)};"
          f"  b(p') and b_N same class mod J: {same_class}")

# (3) Ising combination
p, q, N = 3, 4, 3
vir = VirVac(central_charge(p, q), Q(1))
v = sing_vec(vir, N)
w = {(5, 2, 2): Q(1), (4, 3, 2): Q(6)}
comb = dict(w)


def add(acc, vec, s):
    for m, cc in vec.items():
        nv = acc.get(m, 0) + s * cc
        if nv == 0:
            acc.pop(m, None)
        else:
            acc[m] = nv


add(comb, vir.actv(-3, v), Q(256, 429))
add(comb, vir.actv(-1, vir.actv(-2, v)), Q(-64, 429))
add(comb, vir.actv(-1, vir.actv(-1, vir.actv(-1, v))), Q(-31, 286))
print("Ising: max PBW length of AEH combination (claimed in F^5, i.e. length <= 2):",
      max(len(m) for m in comb))
u = vir.actv(-2, vir.actv(-1, v))
add(u, vir.actv(-3, v), Q(-N))
comb2 = dict(w)
add(comb2, u, Q(-64, 429))
add(comb2, vir.actv(-1, vir.actv(-1, vir.actv(-1, v))), Q(-31, 286))
print("Ising: AEH combination == w - (64/429) u - (31/286) L_{-1}^3 v :", comb == comb2)
print("Ising: AEH right-hand side:", {m: str(c) for m, c in comb.items()})
