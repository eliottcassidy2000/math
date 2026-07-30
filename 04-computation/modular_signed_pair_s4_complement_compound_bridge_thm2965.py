#!/usr/bin/env python3
"""Exact THM-2965 companion for the two S4 signed-pair complements.

This is deliberately dependency-free.  It checks the representation-theoretic
claim behind the accompanying scratch note:

* W=(C2^3) semidirect S3 has central global flip z and quotient S4;
* exactly two complements of <z> occur, distinguished by flip parity versus
  pair-permutation sign;
* Hom_{H0}(Lambda^5 V,B)=0 but Hom_{H1}(Lambda^5 V,B)=2;
* for H1 the two face/cube incidence orbitals span the Hom space and have the
  stated rank profile.
"""

from fractions import Fraction
from itertools import permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PERMS = tuple(permutations(range(3)))
BITS = tuple(product((0, 1), repeat=3))
FACES = tuple((i, b) for i in range(3) for b in range(2))


def perm_parity(p):
    return sum(p[i] > p[j] for i in range(3) for j in range(i + 1, 3)) % 2


def act_bits(w, x):
    p, d = w
    y = [0, 0, 0]
    for i in range(3):
        y[p[i]] = x[i] ^ d[i]
    return tuple(y)


def act_face(w, face):
    p, d = w
    i, b = face
    return (p[i], b ^ d[i])


def complement_class(x):
    xb = tuple(1 ^ a for a in x)
    return min(x, xb)


OMEGA = tuple(sorted({complement_class(x) for x in BITS}))


def quotient_perm(w):
    return tuple(OMEGA.index(complement_class(act_bits(w, x))) for x in OMEGA)


def cycle_type(q):
    seen = set()
    cycles = []
    for i in range(len(q)):
        if i in seen:
            continue
        j = i
        length = 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = q[j]
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))


def compose_perm(a, b):
    """a after b."""
    return tuple(a[b[i]] for i in range(len(a)))


def generated_perms(generators):
    ident = tuple(range(len(generators[0])))
    seen = {ident}
    frontier = [ident]
    while frontier:
        x = frontier.pop()
        for g in generators:
            y = compose_perm(g, x)
            if y not in seen:
                seen.add(y)
                frontier.append(y)
    return seen


def inverse_perm(q):
    out = [0] * len(q)
    for i, j in enumerate(q):
        out[j] = i
    return tuple(out)


def matrix_rank(rows):
    a = [[Fraction(x) for x in row] for row in rows]
    if not a:
        return 0
    m, n = len(a), len(a[0])
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, m) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][c]
        a[r] = [x / scale for x in a[r]]
        for i in range(m):
            if i == r or not a[i][c]:
                continue
            scale = a[i][c]
            a[i] = [a[i][j] - scale * a[r][j] for j in range(n)]
        r += 1
        if r == m:
            break
    return r


def mat_add(a, b, sign=1):
    return [[a[i][j] + sign * b[i][j] for j in range(len(a[0]))]
            for i in range(len(a))]


def matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))]
            for i in range(len(a))]


def apply_rep_source(w, column):
    # Lambda^5(V) = det(V) tensor V^*.  Pair-block permutations have
    # determinant +1 in dimension six; each within-pair flip contributes -1.
    det6 = -1 if sum(w[1]) % 2 else 1
    return det6, FACES.index(act_face(w, FACES[column]))


def apply_rep_target(w, row):
    # Reordering the three one-dimensional wedge factors contributes sign(sigma).
    wedge_sign = -1 if perm_parity(w[0]) else 1
    return wedge_sign, BITS.index(act_bits(w, BITS[row]))


def equivariant(matrix, group):
    # Check T rho_source(w) = rho_target(w) T entrywise.
    for w in group:
        for col in range(6):
            ss, col2 = apply_rep_source(w, col)
            for row in range(8):
                ts, row2 = apply_rep_target(w, row)
                if ss * matrix[row2][col2] != ts * matrix[row][col]:
                    return False
    return True


W = tuple((p, d) for p in PERMS for d in BITS)
z = ((0, 1, 2), (1, 1, 1))
H0 = tuple(w for w in W if sum(w[1]) % 2 == 0)
H1 = tuple(w for w in W if sum(w[1]) % 2 == perm_parity(w[0]))
require(len(W) == 48 and len(H0) == len(H1) == 24, "group sizes")
require(z not in H0 and z not in H1, "complements avoid central flip")
require(len({quotient_perm(w) for w in H0}) == 24, "H0 quotient is S4")
require(len({quotient_perm(w) for w in H1}) == 24, "H1 quotient is S4")

# H0_ab has order two.  Since W=<z> x H0, complements of <z> are exactly
# graphs of homomorphisms H0 -> <z>, hence there are exactly two.
H0Q = tuple(quotient_perm(w) for w in H0)
commutators = []
for a in H0Q:
    for b in H0Q:
        commutators.append(
            compose_perm(
                compose_perm(compose_perm(a, b), inverse_perm(a)),
                inverse_perm(b),
            )
        )
derived = generated_perms(tuple(commutators))
require(len(derived) == 12 and len(H0Q) // len(derived) == 2, "H0 abelianization C2")

CLASS_SIZES = {
    (1, 1, 1, 1): 1,
    (2, 1, 1): 6,
    (2, 2): 3,
    (3, 1): 8,
    (4,): 6,
}
IRREPS = {
    "1": {(1, 1, 1, 1): 1, (2, 1, 1): 1, (2, 2): 1, (3, 1): 1, (4,): 1},
    "sgn": {(1, 1, 1, 1): 1, (2, 1, 1): -1, (2, 2): 1, (3, 1): 1, (4,): -1},
    "2": {(1, 1, 1, 1): 2, (2, 1, 1): 0, (2, 2): 2, (3, 1): -1, (4,): 0},
    "3": {(1, 1, 1, 1): 3, (2, 1, 1): 1, (2, 2): -1, (3, 1): 0, (4,): -1},
    "3sgn": {(1, 1, 1, 1): 3, (2, 1, 1): -1, (2, 2): -1, (3, 1): 0, (4,): 1},
}


def characters(group):
    source = {}
    target = {}
    census = {}
    for w in group:
        ct = cycle_type(quotient_perm(w))
        census[ct] = census.get(ct, 0) + 1
        fixed_faces = sum(act_face(w, x) == x for x in FACES)
        det6 = -1 if sum(w[1]) % 2 else 1
        fixed_bits = sum(act_bits(w, x) == x for x in BITS)
        wedge_sign = -1 if perm_parity(w[0]) else 1
        cs = det6 * fixed_faces
        ctgt = wedge_sign * fixed_bits
        if ct in source:
            require(source[ct] == cs and target[ct] == ctgt, "class character constant")
        source[ct] = cs
        target[ct] = ctgt
    require(census == CLASS_SIZES, "S4 class census")
    return source, target


def decompose(char):
    out = {}
    for name, irr in IRREPS.items():
        numerator = sum(CLASS_SIZES[c] * char[c] * irr[c] for c in CLASS_SIZES)
        require(numerator % 24 == 0, "integral character multiplicity")
        out[name] = numerator // 24
    return {name: m for name, m in out.items() if m}


def hom_dimension(cs, ct):
    numerator = sum(CLASS_SIZES[c] * cs[c] * ct[c] for c in CLASS_SIZES)
    require(numerator % 24 == 0, "integral Hom dimension")
    return numerator // 24


cs0, ct0 = characters(H0)
cs1, ct1 = characters(H1)
require(decompose(cs0) == {"1": 1, "2": 1, "3": 1}, "H0 source decomposition")
require(decompose(ct0) == {"sgn": 2, "3sgn": 2}, "H0 target decomposition")
require(decompose(cs1) == {"sgn": 1, "2": 1, "3": 1}, "H1 source decomposition")
require(decompose(ct1) == {"1": 1, "sgn": 1, "3": 1, "3sgn": 1}, "H1 target decomposition")
require(hom_dimension(cs0, ct0) == 0, "H0 Hom vanishes")
require(hom_dimension(cs1, ct1) == 2, "H1 Hom has dimension two")

# Face/cube incidence and its opposite.
A = [[int(eps[i] == b) for (i, b) in FACES] for eps in BITS]
Ac = [[1 - A[r][c] for c in range(6)] for r in range(8)]
require(not equivariant(A, H0) and not equivariant(Ac, H0), "H0 incidence hostile")
require(equivariant(A, H1) and equivariant(Ac, H1), "H1 incidence intertwiners")
require(matrix_rank(A) == 4 and matrix_rank(Ac) == 4, "incidence ranks")
require(matrix_rank(mat_add(A, Ac)) == 1, "constant orbital rank")
require(matrix_rank(mat_add(A, Ac, sign=-1)) == 3, "first Walsh orbital rank")

# Global conjugation z pairs complementary cube vertices.  The real trace
# quotient is z-even (Walsh degrees 0+2), whereas the incidence image is
# degrees 0+1.  Hence only its constant line survives trace symmetrization.
Sz = [[0] * 6 for _ in range(6)]
Tz = [[0] * 8 for _ in range(8)]
for col in range(6):
    ss, col2 = apply_rep_source(z, col)
    Sz[col2][col] = ss
for row in range(8):
    ts, row2 = apply_rep_target(z, row)
    Tz[row2][row] = ts
require(matmul(A, Sz) == [[-x for x in row] for row in matmul(Tz, A)],
        "incidence anti-commutes with global conjugation")
omega_reps = OMEGA
Trace = [[int(eps == rep or eps == tuple(1 ^ a for a in rep)) for eps in BITS]
         for rep in omega_reps]
AntiTrace = [[int(eps == rep) - int(eps == tuple(1 ^ a for a in rep)) for eps in BITS]
             for rep in omega_reps]
first_walsh = mat_add(A, Ac, sign=-1)
require(matrix_rank(matmul(Trace, A)) == 1, "trace sees only constant incidence")
require(matrix_rank(matmul(Trace, first_walsh)) == 0, "trace kills first Walsh")
require(matrix_rank(matmul(AntiTrace, A)) == 3, "anti-trace retains first Walsh")

# The group intertwiner is not an intertwiner of exterior-power multiplication.
# For M=diag(2,2,3,3,5,5), a fifth weight omitting a 2 is 900/2=450,
# whereas every balanced third weight is 2*3*5=30.
operator_source_weight = (2 * 2 * 3 * 3 * 5 * 5) // 2
operator_target_weight = 2 * 3 * 5
require(operator_source_weight == 450 and operator_target_weight == 30,
        "operator-intertwining hostile")

# Find a marked (2,3)-generating pair in H0.  The H1 involution is z times
# the H0 involution, while the order-three lift is unchanged.
def qorder(q):
    x = tuple(range(len(q)))
    for n in range(1, 13):
        x = compose_perm(q, x)
        if x == tuple(range(len(q))):
            return n
    raise RuntimeError("order bound")


pair = None
for s0 in H0:
    qs = quotient_perm(s0)
    if qorder(qs) != 2:
        continue
    for t in sorted(set(H0).intersection(H1)):
        qt = quotient_perm(t)
        if qorder(qt) == 3 and len(generated_perms((qs, qt))) == 24:
            pair = (s0, t)
            break
    if pair:
        break
require(pair is not None, "marked (2,3) generators")
s0, t = pair
s1 = (s0[0], tuple(1 ^ b for b in s0[1]))
zt = (t[0], tuple(1 ^ b for b in t[1]))
require(s1 in H1 and quotient_perm(s1) == quotient_perm(s0), "central involution twist")
require(qorder(quotient_perm(s0)) == 2 and qorder(quotient_perm(t)) == 3, "local orders")
face_perm = lambda w: tuple(FACES.index(act_face(w, x)) for x in FACES)
require(qorder(face_perm(s0)) == qorder(face_perm(s1)) == 2, "C2 central twist keeps order")
require(qorder(face_perm(t)) == 3 and qorder(face_perm(zt)) == 6, "C3 central twist doubles order")

print("THM-2965 SIGNED-PAIR S4 COMPLEMENT / COMPOUND BRIDGE")
print(f"W={len(W)} quotient={len({quotient_perm(w) for w in W})} H0={len(H0)} H1={len(H1)} complements=2")
print(f"H0 source={decompose(cs0)} target={decompose(ct0)} Hom={hom_dimension(cs0, ct0)}")
print(f"H1 source={decompose(cs1)} target={decompose(ct1)} Hom={hom_dimension(cs1, ct1)}")
print(f"orbital_ranks A={matrix_rank(A)} Ac={matrix_rank(Ac)} sum={matrix_rank(mat_add(A, Ac))} diff={matrix_rank(mat_add(A, Ac, sign=-1))}")
print("trace_projection rank=1 first_walsh_rank=0 antitrace_rank=3 z_anticommutes=1")
print(f"operator_hostile={operator_source_weight}_vs_{operator_target_weight}")
semantic_s0 = (s0[1], s0[0])
semantic_s1 = (s1[1], s1[0])
semantic_t = (t[1], t[0])
print(f"marked_s0={semantic_s0} marked_s1=z*s0={semantic_s1} common_t={semantic_t}")
print("local_orders=2,3 quotient_generated=24")
print("central_twist_orders C2:2_to_2 C3:3_to_6")
print("PASS")
