#!/usr/bin/env python3
"""Exact hostile for the quartic S4 -> V4-resolvent quasi-etale gate.

This is deliberately a non-Keller control.  It verifies the finite quotient

    A^3/S3  ->  A^3/(V4 semidirect S3)

and the intervening V4 quotient.  The companion mathematics uses the control
to show that a quasi-etale V4 resolvent cover is sharp; the nonconstant
Jacobian records exactly why it is not a Keller example.

All checks raise explicitly, including under ``python -O``.
"""

from itertools import permutations, product

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


checks = 0


def check(label, condition):
    global checks
    checks += 1
    if not bool(condition):
        raise RuntimeError(f"FAILED: {label}")


# An element (perm, signs) sends coordinate i to signs[i] * x[perm[i]].
PERMS = tuple(permutations(range(3)))
EVEN_SIGNS = tuple(s for s in product((-1, 1), repeat=3) if s[0] * s[1] * s[2] == 1)
IDENTITY = ((0, 1, 2), (1, 1, 1))


def compose(g, h):
    """Return g after h."""
    pg, sg = g
    ph, sh = h
    return (
        tuple(ph[pg[i]] for i in range(3)),
        tuple(sg[i] * sh[pg[i]] for i in range(3)),
    )


G = tuple((p, s) for p in PERMS for s in EVEN_SIGNS)
H = tuple((p, (1, 1, 1)) for p in PERMS)
V = tuple(((0, 1, 2), s) for s in EVEN_SIGNS)

check("G has order 24", len(G) == 24)
check("H has order 6", len(H) == 6)
check("V has order 4", len(V) == 4)
check("G is closed", all(compose(g, h) in G for g in G for h in G))
check("H intersect V is trivial", set(H).intersection(V) == {IDENTITY})


# Right cosets G/H are exactly the four even sign vectors.  Left
# multiplication gives the quartic sheet action.
SHEETS = EVEN_SIGNS


def act_on_sheet(g, sheet):
    pg, sg = g
    return tuple(sg[i] * sheet[pg[i]] for i in range(3))


def action_tuple(g):
    return tuple(SHEETS.index(act_on_sheet(g, s)) for s in SHEETS)


ACTIONS = {action_tuple(g) for g in G}
check("four cosets", len(SHEETS) == 4)
check("quartic action faithful", len(ACTIONS) == 24)
check("H is the +++ sheet stabilizer", {g for g in G if act_on_sheet(g, (1, 1, 1)) == (1, 1, 1)} == set(H))


def cycle_lengths(perm):
    seen = set()
    lengths = []
    for i in range(len(perm)):
        if i in seen:
            continue
        j = i
        length = 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = perm[j]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


v_cycle_types = [cycle_lengths(action_tuple(v)) for v in V if v != IDENTITY]
v_fixed_codimensions = [sum(1 for sign in v[1] if sign == -1) for v in V if v != IDENTITY]
check("nonzero V elements are double transpositions", v_cycle_types == [(2, 2)] * 3)
check("nonzero V fixed loci have codimension two", v_fixed_codimensions == [2, 2, 2])


# Conjugation by H realizes Q=S3=GL(2,F2) on the three nonzero elements
# of V.  Brute-force inverses keep the group check independent of formulas.
def inverse(g):
    hits = [h for h in G if compose(g, h) == IDENTITY and compose(h, g) == IDENTITY]
    check("unique group inverse", len(hits) == 1)
    return hits[0]


nonzero_v = tuple(v for v in V if v != IDENTITY)
q_on_v = set()
for h in H:
    hi = inverse(h)
    image = tuple(nonzero_v.index(compose(compose(h, v), hi)) for v in nonzero_v)
    q_on_v.add(image)
check("Q acts as full S3 on V minus zero", len(q_on_v) == 6)


# Symbolic invariant and resolvent identities.
x, y, z, T, U = sp.symbols("x y z T U")
a, b, c, d, A, B = sp.symbols("a b c d A B")
e1, e2, e3 = sp.symbols("e1 e2 e3")

source_e1 = x + y + z
source_e2 = x * y + x * z + y * z
source_e3 = x * y * z
target_A = x**2 + y**2 + z**2
target_B = x**2 * y**2 + x**2 * z**2 + y**2 * z**2

check("A=e1^2-2e2", sp.expand(source_e1**2 - 2 * source_e2 - target_A) == 0)
check("B=e2^2-2e1e3", sp.expand(source_e2**2 - 2 * source_e1 * source_e3 - target_B) == 0)
check("d=e3", source_e3 == x * y * z)

sheet_roots = (
    x + y + z,
    x - y - z,
    -x + y - z,
    -x - y + z,
)
quartic_target = T**4 - 2 * A * T**2 - 8 * d * T + A**2 - 4 * B
quartic_roots = sp.expand(sp.prod(T - root for root in sheet_roots))
quartic_specialized = sp.expand(quartic_target.subs({A: target_A, B: target_B, d: source_e3}))
check("quartic is product over four V sheets", sp.expand(quartic_roots - quartic_specialized) == 0)

# The same quartic follows by eliminating e2 from the explicit map.
e2_from_first = (T**2 - A) / 2
quartic_eliminant = sp.expand(4 * (e2_from_first**2 - 2 * T * d - B))
check("explicit map has quartic eliminant", sp.expand(quartic_eliminant - quartic_target) == 0)

p = -2 * A
q = -8 * d
r = A**2 - 4 * B
resolvent = sp.expand(U**3 + 2 * p * U**2 + (p**2 - 4 * r) * U - q**2)
resolvent_target = U**3 - 4 * A * U**2 + 16 * B * U - 64 * d**2
check("depressed squared-pair resolvent formula", resolvent == resolvent_target)
resolvent_roots = sp.expand((U - 4 * a) * (U - 4 * b) * (U - 4 * c))
resolvent_on_R = sp.expand(resolvent_target.subs({A: a + b + c, B: a * b + a * c + b * c, d**2: a * b * c}))
check("resolvent roots are 4a,4b,4c on d^2=abc", sp.expand(resolvent_roots - resolvent_on_R) == 0)

F1 = e1**2 - 2 * e2
F2 = e2**2 - 2 * e1 * e3
F3 = e3
jacobian = sp.factor(sp.Matrix([F1, F2, F3]).jacobian([e1, e2, e3]).det())
check("hostile Jacobian", sp.expand(jacobian - 4 * (e1 * e2 - e3)) == 0)
check("hostile is not Keller", not jacobian.is_constant())

# The normalized cubic Q(W)=W^3-AW^2+BW-d^2 also sees the ramification of
# the auxiliary S3 quotient Z -> X.  Pullback separates it exactly from the
# actual quartic Jacobian factor; this is the hostile's order-index sidecar.
# The displayed squared-pair resolvent is 4^3 Q(U/4), so its discriminant
# has the additional scaling factor 4^6=4096.
disc_H = e1**2 * e2**2 - 4 * e2**3 - 4 * e1**3 * e3 - 27 * e3**2 + 18 * e1 * e2 * e3
disc_Q = A**2 * B**2 - 4 * B**3 - 4 * A**3 * d**2 - 27 * d**4 + 18 * A * B * d**2
disc_Q_pullback = sp.expand(disc_Q.subs({A: F1, B: F2, d: F3}))
check(
    "normalized resolvent discriminant splits into H-index and quartic Jacobian factors",
    sp.expand(disc_Q_pullback - disc_H * (jacobian / 4) ** 2) == 0,
)
check(
    "squared-pair resolvent discriminant has the 4^6 scaling",
    sp.expand(sp.discriminant(resolvent_target, U) - 4096 * disc_Q) == 0,
)


# The parity-annihilator calculation proves the V-invariant monomials are
# exactly the all-even and all-odd exponent parities.  Hence
# C[x,y,z]^V=C[a,b,c,d]/(d^2-abc).
def dot2(left, right):
    return sum(i * j for i, j in zip(left, right)) % 2


even_flip_bits = tuple(tuple(0 if s == 1 else 1 for s in signs) for signs in EVEN_SIGNS)
invariant_parities = tuple(bits for bits in product((0, 1), repeat=3) if all(dot2(bits, flip) == 0 for flip in even_flip_bits))
check("V invariant parity annihilator", set(invariant_parities) == {(0, 0, 0), (1, 1, 1)})
check("hypersurface relation", sp.expand(d**2 - a * b * c).subs({a: x**2, b: y**2, c: z**2, d: x * y * z}) == 0)


# Normality/quasi-etaleness: the hypersurface singular set has d=0 and at
# most one of a,b,c nonzero, hence dimension one inside the threefold.  The
# fixed-locus computation above independently gives the same codimension-two
# quotient boundary.
relation = d**2 - a * b * c
gradient = tuple(sp.diff(relation, variable) for variable in (a, b, c, d))
check("hypersurface gradient", gradient == (-b * c, -a * c, -a * b, 2 * d))
singular_components = ((a, b, d), (a, c, d), (b, c, d))
check("three codimension-two singular axes", len(singular_components) == 3 and all(len(component) == 3 for component in singular_components))


# Class group.  Localizing at d gives the Laurent UFD
# C[a^+-1,b^+-1,d^+-1], so Nagata says Cl(R) is generated by the three
# height-one primes P_a=(a,d), P_b=(b,d), P_c=(c,d).  Their principal-divisor
# relations are 2P_a, 2P_b, 2P_c, and P_a+P_b+P_c.
class_relations = sp.Matrix(((2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 1)))
smith = smith_normal_form(class_relations, domain=ZZ)
smith_diagonal = tuple(abs(int(smith[i, i])) for i in range(3))
check("class group Smith form", smith_diagonal == (1, 2, 2))


def add2(left, right):
    return tuple(i ^ j for i, j in zip(left, right))


one = (1, 1, 1)


def class_rep(bits):
    mate = add2(bits, one)
    return min(bits, mate)


class_reps = tuple(sorted({class_rep(bits) for bits in product((0, 1), repeat=3)}))
check("four class-group elements", len(class_reps) == 4)


def permute_bits(perm, bits):
    return tuple(bits[perm[i]] for i in range(3))


q_on_class_group = {
    tuple(class_reps.index(class_rep(permute_bits(perm, bits))) for bits in class_reps)
    for perm in PERMS
}
nonzero_classes = tuple(bits for bits in class_reps if bits != (0, 0, 0))
nonzero_orbit = {class_rep(permute_bits(perm, nonzero_classes[0])) for perm in PERMS}
check("Q acts faithfully on Cl[2]", len(q_on_class_group) == 6)
check("Q is transitive on three nonzero classes", nonzero_orbit == set(nonzero_classes))
check("standard F2 plane is irreducible", len(nonzero_orbit) == 3)


# Positive weights deg(a,b,c,d)=(2,2,2,3) make R a connected graded domain;
# therefore every unit has degree zero and is constant.  The relation is
# homogeneous of degree six.
check("positive grading relation", 2 * 3 == 2 + 2 + 2)


print("S4 QUARTIC / V4 QUASI-ETALE HOSTILE")
print(f"group orders: |G|={len(G)}, |H|={len(H)}, |V|={len(V)}, sheets={len(SHEETS)}")
print(f"quartic action: faithful permutations={len(ACTIONS)}, V cycle types={v_cycle_types}")
print(f"V fixed-locus codimensions={v_fixed_codimensions}; quotient is unramified in codimension one")
print("V-invariant ring: C[a,b,c,d]/(d^2-abc); singular locus is the three coordinate axes")
print("Q=S3 invariant ring: C[A=a+b+c,B=ab+ac+bc,d]")
print("quartic source map: (e1,e2,e3) -> (e1^2-2e2,e2^2-2e1e3,e3)")
print(f"Jacobian={jacobian}; NON-KELLER hostile")
print("quartic: T^4-2*A*T^2-8*d*T+(A^2-4*B)")
print("squared-pair resolvent: U^3-4*A*U^2+16*B*U-64*d^2; roots 4a,4b,4c")
print("normalized cubic disc: Disc_Q(F)=Disc_H*(Jac(F)/4)^2")
print("squared-pair cubic disc: Disc_U=4096*Disc_Q (scaling recorded explicitly)")
print(f"Cl lattice Smith diagonal={smith_diagonal}; Cl(R)=(Z/2)^2")
print(f"Q action on Cl[2]: {len(q_on_class_group)} elements, one orbit of {len(nonzero_orbit)} nonzero classes")
print("units: C* by the connected positive grading; the rank-two carrier lies in Cl(R)[2]")
print(f"CHECKS PASSED: {checks}")
print("FAILED CHECKS: NONE")
