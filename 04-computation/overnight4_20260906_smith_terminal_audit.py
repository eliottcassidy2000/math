"""Independent bounded audit of THM-4439 terminal two-jet precision.

No repository imports. Literal integer Smith forms, rational matrix inverses,
reciprocal denominators, and ball enumeration are separate calculation paths.
All verification gates remain active under Python -O.
"""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

GATES = 0


def need(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def valuation(a, p):
    a = int(a)
    if a == 0:
        return None
    result = 0
    while a % p == 0:
        a //= p
        result += 1
    return result


def rational_valuation(q, p):
    if not q:
        return None
    return valuation(q.numerator, p)-valuation(q.denominator, p)


def matrix(nodes):
    degree = 2*len(nodes)
    rows = []
    for x in nodes:
        rows.append([x**j for j in range(degree)])
        rows.append([0]+[j*x**(j-1) for j in range(1, degree)])
    return Matrix(rows)


def terminals(nodes, p):
    if len(nodes) == 1:
        return []
    all_depths = sorted({valuation(x-y, p) for x in nodes for y in nodes if x != y})
    found = set()
    for f in all_depths:
        for x in nodes:
            C = tuple(sorted(y for y in nodes if y == x or valuation(y-x, p) >= f))
            if len(C) < 2:
                continue
            if all(valuation(y-z, p) == f for y in C for z in C if y != z):
                found.add((f, C))
    result = []
    for f, C in sorted(found):
        sums = {sum(valuation(x-y, p) for y in nodes if x != y) for x in C}
        need(len(sums) == 1 and 2 <= len(C) <= p, "terminal constant S and child count")
        S = sums.pop()
        reciprocal_vals = []
        for x in C:
            q = sum((Fraction(1, x-y) for y in nodes if x != y), Fraction(0))
            reciprocal_vals.append(rational_valuation(p**f*q, p))
        finite = [v for v in reciprocal_vals if v is not None]
        need(min(finite) == (1 if p > 2 and len(C) == p else 0),
             "exact complete-cluster reciprocal cancellation budget")
        L = 2*S+max(0, f-int(len(C) == p))
        result.append((f, C, S, L, tuple(reciprocal_vals)))
    need(bool(result), "nonempty terminal list")
    return result


def check(name, p, nodes):
    nodes = tuple(nodes)
    need(len(nodes) == len(set(nodes)), "distinct nodes")
    E = matrix(nodes)
    D = smith_normal_form(E, domain=ZZ)
    diag = [abs(int(D[i, i])) for i in range(D.rows)]
    need(all(diag) and all(diag[i+1] % diag[i] == 0 for i in range(len(diag)-1)),
         "literal integer Smith divisibility")
    exponents = tuple(valuation(x, p) for x in diag)
    tree = terminals(nodes, p)
    predicted = max((entry[3] for entry in tree), default=0)
    need(max(exponents) == predicted, "terminal formula vs full integer Smith")
    inverse = E.inv()
    inverse_loss = max(valuation(int(x.q), p) for x in inverse)
    need(inverse_loss == predicted, "literal rational inverse denominator")
    contributions = []
    for x in nodes:
        S = sum(valuation(x-y, p) for y in nodes if x != y)
        q = sum((Fraction(1, x-y) for y in nodes if x != y), Fraction(0))
        v2q = rational_valuation(2*q, p)
        contributions.append(max(2*S, 2*S-v2q) if v2q is not None else 2*S)
    need(max(contributions) == predicted, "direct reciprocal loss formula")
    determinant_v = 4*sum(valuation(nodes[j]-nodes[i], p)
                          for i in range(len(nodes)) for j in range(i+1, len(nodes)))
    need(sum(exponents) == determinant_v, "confluent determinant normalization")
    print(name, "p", p, "nodes", nodes, "Smith", exponents, "L", predicted,
          "terminals", tree)
    return exponents


def main():
    print("source_sha256", sha256(Path(__file__).read_bytes()).hexdigest())
    cases = [
        ("singleton2", 2, (7,)), ("singleton3", 3, (-3,)),
        ("dyadic_gap0", 2, (0, 1)), ("dyadic_gap1", 2, (0, 2)),
        ("dyadic_gap3", 2, (0, 8)), ("ternary_partial_pair", 3, (0, 3)),
        ("ternary_full_depth0", 3, (0, 1, 2)),
        ("ternary_full_depth1", 3, (0, 3, 6)),
        ("ternary_full_depth2", 3, (0, 9, 18)),
        ("ternary_full_lifted", 3, (0, 9, 45)),
        ("ternary_full_near_outsider", 3, (0, 1, 3, 6)),
        ("ternary_depth2_near_outsider", 3, (0, 3, 9, 18)),
        ("ternary_lifted_outsider", 3, (0, 1, 3, 15)),
        ("quinary_partial_three", 5, (0, 5, 10)),
        ("quinary_full_depth1", 5, (0, 5, 10, 15, 20)),
        ("quinary_full_outsider", 5, (0, 1, 5, 10, 15, 20)),
        ("quinary_full_depth2", 5, (0, 25, 50, 75, 100)),
        ("septenary_full", 7, (0, 7, 14, 21, 28, 35, 42)),
        ("dyadic_twin_A", 2, (0, 8, 16, 40)),
        ("dyadic_twin_B", 2, (0, 8, 24, 32)),
        ("ternary_twin_A", 3, (0, 9, 27, 81)),
        ("ternary_twin_B", 3, (0, 9, 54, 81)),
        ("unequal_terminal_branches", 3, (0, 9, 18, 1, 4)),
        ("signed_translate", 3, (-7, 2, 11)),
    ]
    answers = {name: check(name, p, nodes) for name, p, nodes in cases}
    need(answers["ternary_full_depth1"][-1] == 4, "omitting full-residue subtraction fails")
    need(answers["ternary_partial_pair"][-1] == 3, "subtracting at partial residue set fails")
    A, B = answers["ternary_twin_A"], answers["ternary_twin_B"]
    need(A == (0, 0, 2, 6, 7, 12, 15, 22), "ternary first intermediate profile")
    need(B == (0, 0, 2, 6, 7, 13, 14, 22), "ternary second intermediate profile")
    need(sum(min(13, a) for a in A) == 53 and sum(min(13, b) for b in B) == 54,
         "metric-only largest exponent does not determine intermediate kernels")
    print("named_cases", len(cases))
    print("ternary_twin_kernel_log3_at_precision13", 53, 54)
    print("active_gates", GATES)
    print("PASS THM4439 formula, exact complete-residue exception, and intermediate-factor boundary")


if __name__ == "__main__":
    main()
