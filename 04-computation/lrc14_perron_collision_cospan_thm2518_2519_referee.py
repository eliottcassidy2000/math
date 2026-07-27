#!/usr/bin/env python3
"""Exact finite-tower referee for THM-2518/2519.

The circle grid C_(K*N) is written as q=y+rK, with y in C_K and
r in C_N.  The quotient q -> y is the finite model of x -> N x.
All arithmetic and all reported identities are exact.
"""

from fractions import Fraction


P = 13
L = 2
N = P**L
M = N // P
K = 5
D = K * N


def mean(values):
    values = list(values)
    return sum(values, Fraction(0)) / len(values)


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def check(name, left, right):
    require(left == right, (name, left, right))
    print(f"PASS {name}: {left}")


# A deliberately irregular nonnegative response and a positive future mask.
F = tuple(
    int(((q * q + 7 * q + 3) % 29) < 11)
    + int(((5 * q + 2) % 31) < 7)
    for q in range(D)
)
G = (1, 0, 1, 0, 1)

A = mean(F)
rho = mean(G)


def branch(r):
    return mean(G[y] * F[y + r * K] for y in range(K))


def branch_pair(r, s):
    return mean(
        G[y] * F[y + r * K] * F[y + s * K]
        for y in range(K)
    )


def correlation(d):
    return mean(
        G[q % K] * F[q] * F[(q + d * K) % D]
        for q in range(D)
    )


perron = tuple(
    mean(F[y + r * K] for r in range(N))
    for y in range(K)
)

S = mean(G[q % K] * F[q] for q in range(D))
C = tuple(correlation(d) for d in range(N))

print("=== THM-2518 finite Perron/groupoid identities ===")
check("mean branch = owner-gated first Perron moment", mean(branch(r) for r in range(N)), S)
check("Perron first moment", S, mean(G[y] * perron[y] for y in range(K)))

for d in (0, 1, 17, 83, 168):
    check(
        f"needle d={d} = mean ancestry-sheet pair",
        C[d],
        mean(branch_pair(r, (r + d) % N) for r in range(N)),
    )

check(
    "mean needles = owner-gated squared Perron moment",
    mean(C),
    mean(G[y] * perron[y] * perron[y] for y in range(K)),
)

for d in range(N):
    require(C[d] == C[(-d) % N], ("antipodal", d, C[d], C[(-d) % N]))
print("PASS antipodal needle symmetry: C[d] = C[-d] for all d")


def residue_average(u):
    return mean(C[(u + P * e) % N] for e in range(M))


def shifted_coset_average(q, u):
    return mean(F[(q + (u + P * e) * K) % D] for e in range(M))


B = tuple(residue_average(u) for u in range(P))
for u in range(P):
    check(
        f"fixed last-digit average u={u}",
        B[u],
        mean(
            G[q % K] * F[q] * shifted_coset_average(q, u)
            for q in range(D)
        ),
    )

print("PASS unit residues are exact depth-L first collisions")
for u in range(1, P):
    for e in range(M):
        d = u + P * e
        require(d % P == u, ("last digit", d, u))
        require(
            (P ** (L - 1) * Fraction(d, N)) % 1 == Fraction(u, P),
            ("penultimate separation", d, u),
        )
        require((P**L * Fraction(d, N)) % 1 == 0, ("coalescence", d))


print("=== THM-2519 last-digit conditional variance ===")
# h=P_M F on C_(13K); z -> z mod K is the final digit quotient.
h = tuple(
    mean(F[z + e * (P * K)] for e in range(M))
    for z in range(P * K)
)
qbar = tuple(mean(h[y + u * K] for u in range(P)) for y in range(K))

B0_fibre = mean(G[z % K] * h[z] * h[z] for z in range(P * K))
Ball_fibre = mean(G[y] * qbar[y] * qbar[y] for y in range(K))
variance = mean(
    G[z % K] * (h[z] - qbar[z % K]) ** 2
    for z in range(P * K)
)
dirichlet = mean(
    G[y]
    * mean(
        (h[y + r * K] - h[y + s * K]) ** 2
        for r in range(P)
        for s in range(P)
    )
    / 2
    for y in range(K)
)

check("u=0 correlation is predecessor square mass", B[0], B0_fibre)
check("mean-u correlation is future conditional square", mean(B), Ball_fibre)
check("collision drift = conditional variance", B[0] - mean(B), variance)
check("collision drift = K13 Dirichlet energy", variance, dirichlet)
require(variance > 0, ("positive control", variance))
print(f"PASS hostile control is positive: drift={variance}")


print("=== anchored old-owner / future-owner recurrence averaging ===")
H = tuple(int(((3 * q + 1) % 17) < 8) for q in range(D))
W = tuple(H[q] * G[q % K] for q in range(D))


def owner_edge(d):
    return mean(W[q] * F[(q + d * K) % D] for q in range(D))


def owner_triangle(d, e):
    return mean(
        W[q] * F[(q + d * K) % D] * F[(q + e * K) % D]
        for q in range(D)
    )


u, v = 1, 2
edge_average = mean(owner_edge(u + P * e) for e in range(M))
triangle_average = mean(
    owner_triangle(u + P * e, v + P * f)
    for e in range(M)
    for f in range(M)
)
edge_grid = mean(W[q] * shifted_coset_average(q, u) for q in range(D))
triangle_grid = mean(
    W[q] * shifted_coset_average(q, u) * shifted_coset_average(q, v)
    for q in range(D)
)
check("same-owner edge coset average", edge_average, edge_grid)
check("same-owner triangle coset average", triangle_average, triangle_grid)


print("=== sharp zero-drift future-measurable control ===")
# If F0 factors through the depth-L future, every sibling value is identical.
F0 = tuple(int((q % K) in (0, 1)) for q in range(D))


def correlation0(d):
    return mean(
        G[q % K] * F0[q] * F0[(q + d * K) % D]
        for q in range(D)
    )


C0 = tuple(correlation0(d) for d in range(N))
require(len(set(C0)) == 1, ("zero-drift hostile", len(set(C0))))
print(f"PASS future-measurable hostile: all {N} needles equal {C0[0]}")

print("=== summary ===")
print(f"p={P}; L={L}; N={N}; grid={D}; mean(F)={A}; mean(G)={rho}")
print("all exact finite-tower identities passed")
