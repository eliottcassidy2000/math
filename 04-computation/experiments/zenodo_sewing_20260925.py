"""Exact integral lattice underlying Zenodo 22800071's D3 E block.

Source matrices: deposited README.md, section 'Block conventions'.
Source note: 05-knowledge/results/zenodo_sewing_20260925.md.
This script neither executes Wolfram code nor revalidates amplitude symbols.
"""

from fractions import Fraction as Q
from itertools import product
import json


def check(value, message):
    if not value:
        raise RuntimeError(message)


def mat(a, b, c, d):
    return ((Q(a), Q(b)), (Q(c), Q(d)))


def mul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def transpose(a):
    return tuple(tuple(a[j][i] for j in range(2)) for i in range(2))


def neg(a):
    return tuple(tuple(-x for x in row) for row in a)


def apply(a, x):
    return tuple(sum(row[i] * x[i] for i in range(2)) for row in a)


def mod2(a):
    check(all(x.denominator == 1 for row in a for x in row),
          "cannot reduce a denominator divisible by two")
    return tuple(tuple(int(x) % 2 for x in row) for row in a)


def closure(generators):
    found = {I}
    pending = [I]
    while pending:
        x = pending.pop()
        for g in generators:
            y = mul(g, x)
            if y not in found:
                found.add(y)
                pending.append(y)
    return found


def norm(x):
    a, b = x
    return a * a - a * b + b * b


I = mat(1, 0, 0, 1)
C = mat(Q(-1, 2), Q(-1, 2), Q(3, 2), Q(-1, 2))
F = mat(1, 0, 0, -1)
B = mat(1, -1, 1, 1)
Bi = mat(Q(1, 2), Q(1, 2), Q(-1, 2), Q(1, 2))


def main():
    check(mul(Bi, B) == I, "basis inverse")
    c = mul(Bi, mul(C, B))
    f = mul(Bi, mul(F, B))
    check(c == mat(0, -1, 1, -1), "integral rotation")
    check(f == mat(0, -1, -1, 0), "integral reflection")
    check(mul(mul(c, c), c) == I, "order three")
    check(mul(f, f) == I and mul(mul(f, c), f) == mul(c, c), "dihedral relation")
    group = closure((c, f))
    reduced = {mod2(g) for g in group}
    check(len(group) == 6 and len(reduced) == 6, "faithful D3 reduction")
    full = group | {neg(g) for g in group}
    kernel = {g for g in full if mod2(g) == mod2(I)}
    check(len(full) == 12 and kernel == {I, neg(I)}, "parity kernel")
    check(mod2(c) == ((0, 1), (1, 1)), "Fibonacci color matrix")
    check(mod2(f) == ((0, 1), (1, 0)), "color swap")
    colors = [(1, 0), (0, 1), (1, 1)]
    check([tuple(int(v) % 2 for v in apply(c, color)) for color in colors]
          == [colors[1], colors[2], colors[0]], "three-color rotation")
    H = mat(3, 0, 0, 1)
    check(mul(transpose(C), mul(H, C)) == H, "invariant metric C")
    check(mul(transpose(F), mul(H, F)) == H, "invariant metric F")
    check(mul(transpose(C), C) != I, "identity metric hostile")
    cases = 0
    for x in product(range(-16, 17), repeat=2):
        bx = apply(B, x)
        check(3 * bx[0] ** 2 + bx[1] ** 2 == 4 * norm(x), "hexagonal norm")
        for g in group:
            gx = apply(g, x)
            check(norm(gx) == norm(x), "integral norm preservation")
            gx_mod = tuple(int(v) % 2 for v in gx)
            reduced_result = tuple(int(v) % 2 for v in apply(mod2(g),
                                                            tuple(z % 2 for z in x)))
            check(gx_mod == reduced_result, "reduction intertwines action")
            cases += 1
    minus_c = neg(c)
    check(mul(mul(minus_c, minus_c), minus_c) == neg(I), "parity-C has order six")
    check(mod2(minus_c) == mod2(c), "order six collapses to order three")
    fib = mat(0, 1, 1, 1)
    check(mod2(fib) == mod2(c), "same mod-two Fibonacci action")
    check(mul(mul(fib, fib), fib) == mat(1, 2, 2, 3) != I,
          "integer Fibonacci lift is not order three")
    check(norm((1, 0)) == 1 and norm((3, 0)) == 9, "height-loss hostile")
    check(tuple(x % 2 for x in (1, 0)) == tuple(x % 2 for x in (3, 0)),
          "same color, distinct norm")
    print(json.dumps({"status": "EXACT representation transfer; no amplitude recomputation",
                      "source": "https://doi.org/10.5281/zenodo.22800071",
                      "basis_columns": [[1, 1], [-1, 1]],
                      "integral_C": [[int(x) for x in row] for row in c],
                      "integral_F": [[int(x) for x in row] for row in f],
                      "D3_elements": len(group), "distinct_mod2_actions": len(reduced),
                      "parity_times_D3_elements": len(full), "mod2_kernel_size": len(kernel),
                      "norm_and_intertwining_cases": cases,
                      "source_metric": [3, 1], "norm_lost_example": [1, 9],
                      "negative_C_order": 6, "negative_C_reduction_order": 3,
                      "integer_Fibonacci_cube": [[1, 2], [2, 3]]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
