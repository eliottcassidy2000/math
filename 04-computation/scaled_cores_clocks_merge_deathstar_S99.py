"""Historical S99 scale-and-clock lens; corrected by MISTAKE-232.

Only the exact Paley spectral identities at the prime moduli 7 and 13 survive.
There is no Paley object at composite modulus 14 and no proved map between
GMC(2) Frobenius and the additive clocks used by THM-2057.
"""


def matmul(left, right):
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def add(*matrices):
    return [
        [sum(matrix[i][j] for matrix in matrices) for j in range(len(matrices[0][0]))]
        for i in range(len(matrices[0]))
    ]


def scale(c, matrix):
    return [[c * entry for entry in row] for row in matrix]


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def ones(n):
    return [[1 for _ in range(n)] for _ in range(n)]


def paley_matrix(p):
    qr = {(x * x) % p for x in range(1, p)}
    return [[int(i != j and (j - i) % p in qr) for j in range(p)] for i in range(p)]


print("STATUS: HISTORICAL / PARTLY REFUTED (MISTAKE-232).")
print("SURVIVES: exact Paley spectra at primes 7 and 13; THM-2057 is independent.")
print("REFUTED: a Paley object at 14 and an exact Frobenius-to-LRC-clock transfer.")
print()
print("=== Exact Paley identities at the two prime moduli ===")

# At p=7, A is a tournament and S=A-A^T is its skew Jacobsthal matrix.
A7 = paley_matrix(7)
S7 = [[A7[i][j] - A7[j][i] for j in range(7)] for i in range(7)]
assert matmul(S7, S7) == add(ones(7), scale(-7, identity(7)))
print("7 (Paley tournament): S^2=J-7I exactly; spectrum 0,(+i√7)^3,(-i√7)^3.")
print("  adjacency nonprincipal spectrum ((-1+i√7)/2)^3, ((-1-i√7)/2)^3.")

# At p=13, A is the strongly regular Paley graph with parameters (13,6,2,3).
A13 = paley_matrix(13)
assert matmul(A13, A13) == add(scale(-1, A13), scale(3, identity(13)), scale(3, ones(13)))
print("13 (Paley graph): A^2+A=3I+3J exactly; spectrum 6,((-1+√13)/2)^6,((-1-√13)/2)^6.")
print("14=2·7 is the runner count only; no Paley graph/tournament at composite modulus 14 is claimed.")
print()
print("HISTORICAL HEURISTIC ONLY:")
print("  Both proofs use words such as 'scale' and 'clock', but no typed map or preserved predicate connects them.")
print("  GMC2: dilation by p and Frobenius in characteristic p produce the exact surviving power Q̄^p.")
print("  LRC : THM-2057 uses additive unit orbits modulo 12a and 14a.")
print("  Exact LRC arithmetic: lcm(12a,14a)=84a; this is not a Frobenius identity.")
print("  The tournament-zeta and transitive-core discussion is archived as motivation, not an LRC certificate.")
