#!/usr/bin/env python3
"""Dependency-free exact scout for the level-3/level-4 modular quotients.

This is deliberately a scratch companion, not canonical evidence.
"""

from collections import Counter
from itertools import combinations, permutations, product


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def canon_mat(a, n):
    a = tuple(x % n for x in a)
    b = tuple((-x) % n for x in a)
    return min(a, b)


def det(a, n):
    return (a[0] * a[3] - a[1] * a[2]) % n


def mul(a, b, n):
    return canon_mat(
        (
            a[0] * b[0] + a[1] * b[2],
            a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2],
            a[2] * b[1] + a[3] * b[3],
        ),
        n,
    )


def inv(a, n):
    return canon_mat((a[3], -a[1], -a[2], a[0]), n)


def psl(n):
    return tuple(
        sorted(
            {
                canon_mat(a, n)
                for a in product(range(n), repeat=4)
                if det(a, n) == 1 % n
            }
        )
    )


def sl_size(n):
    return sum(
        det(a, n) == 1 % n
        for a in product(range(n), repeat=4)
    )


def closure(gens, n):
    identity = canon_mat((1, 0, 0, 1), n)
    out = {identity}
    todo = [identity]
    while todo:
        x = todo.pop()
        for g in gens:
            y = mul(x, g, n)
            if y not in out:
                out.add(y)
                todo.append(y)
    return frozenset(out)


def order(g, n):
    identity = canon_mat((1, 0, 0, 1), n)
    x = identity
    for k in range(1, 100):
        x = mul(x, g, n)
        if x == identity:
            return k
    raise RuntimeError("order cap")


def canon_point(v, n):
    units = tuple(x for x in range(n) if __import__("math").gcd(x, n) == 1)
    return min(tuple((u * x) % n for x in v) for u in units)


def projective_line(n):
    from math import gcd

    return tuple(
        sorted(
            {
                canon_point((x, y), n)
                for x, y in product(range(n), repeat=2)
                if gcd(gcd(x, y), n) == 1
            }
        )
    )


def act_point(g, v, n):
    return canon_point(
        ((g[0] * v[0] + g[1] * v[1]) % n, (g[2] * v[0] + g[3] * v[1]) % n),
        n,
    )


def action_on_points(group, points, n):
    index = {x: i for i, x in enumerate(points)}
    return {
        g: tuple(index[act_point(g, x, n)] for x in points)
        for g in group
    }


def compose_perm(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def cycle_type(p):
    seen = set()
    sizes = []
    for i in range(len(p)):
        if i in seen:
            continue
        j = i
        size = 0
        while j not in seen:
            seen.add(j)
            size += 1
            j = p[j]
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def census(action):
    return Counter(cycle_type(p) for p in action.values())


def check_action(group, action, n, label):
    identity_perm = tuple(range(len(next(iter(action.values())))))
    identity = canon_mat((1, 0, 0, 1), n)
    require(action[identity] == identity_perm, f"{label}: identity")
    for g in group:
        for h in group:
            require(
                action[mul(g, h, n)] == compose_perm(action[g], action[h]),
                f"{label}: homomorphism",
            )


def orbit(action, base=0):
    return {p[base] for p in action.values()}


def all_two_generated_subgroups(group, n):
    return {
        closure((a, b), n)
        for a in group
        for b in group
    }


def conjugate_subgroup(g, h, n):
    gi = inv(g, n)
    return frozenset(mul(mul(g, x, n), gi, n) for x in h)


def action_on_subgroups(group, subgroups, n):
    index = {h: i for i, h in enumerate(subgroups)}
    return {
        g: tuple(index[conjugate_subgroup(g, h, n)] for h in subgroups)
        for g in group
    }


def action_on_k_subsets(group, base_action, k):
    n = len(next(iter(base_action.values())))
    subsets = tuple(combinations(range(n), k))
    index = {s: i for i, s in enumerate(subsets)}
    return subsets, {
        g: tuple(index[tuple(sorted(p[i] for i in s))] for s in subsets)
        for g, p in base_action.items()
    }


def induced_block_action(group, action, blocks):
    blocks = tuple(frozenset(x) for x in blocks)
    index = {b: i for i, b in enumerate(blocks)}
    out = {}
    for g in group:
        p = action[g]
        out[g] = tuple(index[frozenset(p[i] for i in b)] for b in blocks)
    return out


def find_equivariant_bijection(group, action_a, action_b):
    na = len(next(iter(action_a.values())))
    nb = len(next(iter(action_b.values())))
    require(na == nb, "action sizes")
    for f in permutations(range(nb)):
        if all(
            f[action_a[g][i]] == action_b[g][f[i]]
            for g in group
            for i in range(na)
        ):
            return f
    return None


def fmt_census(c):
    return ", ".join(f"{k}:{c[k]}" for k in sorted(c, reverse=True))


def main():
    S_raw = (0, -1, 1, 0)
    C_raw = (0, -1, 1, 1)

    # Level 3: the four projective points are the four C3 complements to V4.
    G3 = psl(3)
    S3 = canon_mat(S_raw, 3)
    C3 = canon_mat(C_raw, 3)
    T3 = mul(inv(S3, 3), C3, 3)
    require((sl_size(3), len(G3), len(closure((S3, C3), 3))) == (24, 12, 12), "level 3 sizes")
    require((order(S3, 3), order(C3, 3), order(T3, 3)) == (2, 3, 3), "level 3 orders")
    P3 = projective_line(3)
    A3 = action_on_points(G3, P3, 3)
    check_action(G3, A3, 3, "level 3 projective")
    require(len(P3) == 4, "P1(F3)")
    require(len(set(A3.values())) == 12, "level 3 faithful")
    require(len(orbit(A3)) == 4, "level 3 transitive")
    require(census(A3) == Counter({(3, 1): 8, (2, 2): 3, (1, 1, 1, 1): 1}), "A4 census")
    V3 = frozenset(g for g in G3 if order(g, 3) <= 2)
    require(len(V3) == 4, "level 3 V4")
    H3_all = all_two_generated_subgroups(G3, 3)
    C3_complements = tuple(
        sorted(
            (
                h
                for h in H3_all
                if len(h) == 3 and h.intersection(V3) == {canon_mat((1, 0, 0, 1), 3)}
            ),
            key=lambda h: tuple(sorted(h)),
        )
    )
    require(len(C3_complements) == 4, "four C3 complements")
    AC3 = action_on_subgroups(G3, C3_complements, 3)
    check_action(G3, AC3, 3, "level 3 complement")
    require(find_equivariant_bijection(G3, A3, AC3) is not None, "level 3 complement torsor")

    # Level 4: V4 kernel, four S3 complements, and the two six-point actions.
    G4 = psl(4)
    S4 = canon_mat(S_raw, 4)
    C4 = canon_mat(C_raw, 4)
    T4 = mul(inv(S4, 4), C4, 4)
    require((sl_size(4), len(G4), len(closure((S4, C4), 4))) == (48, 24, 24), "level 4 sizes")
    require((order(S4, 4), order(C4, 4), order(T4, 4)) == (2, 3, 4), "level 4 orders")
    red = {g: canon_mat(tuple(x % 2 for x in g), 2) for g in G4}
    G2 = set(red.values())
    require(len(G2) == 6, "level 4 quotient S3")
    I4 = canon_mat((1, 0, 0, 1), 4)
    K4 = frozenset(g for g in G4 if red[g] == canon_mat((1, 0, 0, 1), 2))
    require(len(K4) == 4 and Counter(order(x, 4) for x in K4) == Counter({2: 3, 1: 1}), "kernel V4")

    H4_all = all_two_generated_subgroups(G4, 4)
    complements = tuple(
        sorted(
            (
                h
                for h in H4_all
                if len(h) == 6 and h.intersection(K4) == {I4} and len({red[x] for x in h}) == 6
            ),
            key=lambda h: tuple(sorted(h)),
        )
    )
    require(len(complements) == 4, "four S3 complements")
    Acomp = action_on_subgroups(G4, complements, 4)
    check_action(G4, Acomp, 4, "level 4 complement")
    require(len(set(Acomp.values())) == 24, "complement action faithful S4")
    require(
        census(Acomp)
        == Counter({(3, 1): 8, (2, 1, 1): 6, (4,): 6, (2, 2): 3, (1, 1, 1, 1): 1}),
        "S4 four-complement census",
    )

    edges, Aedge = action_on_k_subsets(G4, Acomp, 2)
    check_action(G4, Aedge, 4, "level 4 edge")
    require(
        census(Aedge)
        == Counter({(2, 2, 1, 1): 9, (3, 3): 8, (4, 2): 6, (1, 1, 1, 1, 1, 1): 1}),
        "edge census",
    )

    P4 = projective_line(4)
    AP4 = action_on_points(G4, P4, 4)
    check_action(G4, AP4, 4, "level 4 projective")
    require(len(P4) == 6 and len(set(AP4.values())) == 24, "P1(Z/4) faithful")
    require(len(orbit(AP4)) == 6 and len(orbit(Aedge)) == 6, "six-actions transitive")
    require(
        census(AP4)
        == Counter(
            {
                (3, 3): 8,
                (2, 2, 2): 6,
                (4, 1, 1): 6,
                (2, 2, 1, 1): 3,
                (1, 1, 1, 1, 1, 1): 1,
            }
        ),
        "projective six census",
    )

    # Stabilizer isomorphism types are the sharp nonisomorphism witness.
    p_stab = tuple(g for g in G4 if AP4[g][0] == 0)
    e_stab = tuple(g for g in G4 if Aedge[g][0] == 0)
    require(Counter(order(g, 4) for g in p_stab) == Counter({4: 2, 2: 1, 1: 1}), "P1 stabilizer C4")
    require(Counter(order(g, 4) for g in e_stab) == Counter({2: 3, 1: 1}), "edge stabilizer V4")
    require(find_equivariant_bijection(G4, AP4, Aedge) is None, "six-actions nonisomorphic")

    # Both six-actions have the same 3-block S3 quotient.
    P2 = projective_line(2)
    P4_to_P2 = {p: canon_point(tuple(x % 2 for x in p), 2) for p in P4}
    p_blocks = tuple(frozenset(i for i, p in enumerate(P4) if P4_to_P2[p] == q) for q in P2)
    AP4_blocks = induced_block_action(G4, AP4, p_blocks)
    matchings = (
        frozenset((edges.index((0, 1)), edges.index((2, 3)))),
        frozenset((edges.index((0, 2)), edges.index((1, 3)))),
        frozenset((edges.index((0, 3)), edges.index((1, 2)))),
    )
    Aedge_blocks = induced_block_action(G4, Aedge, matchings)
    require(find_equivariant_bijection(G4, AP4_blocks, Aedge_blocks) is not None, "same S3 block quotient")
    require(
        {g for g in G4 if AP4_blocks[g] == (0, 1, 2)}
        == {g for g in G4 if Aedge_blocks[g] == (0, 1, 2)}
        == set(K4),
        "same V4 block kernel",
    )

    print("MODULAR LEVEL 3/4 EXACT SCOUT")
    print(f"|PSL2(F3)|={len(G3)}, P1={len(P3)}, generator orders={order(S3,3)},{order(C3,3)},{order(T3,3)}")
    print(f"level3 P1/complement census: {fmt_census(census(A3))}")
    print(f"level3 C3 complements={len(C3_complements)}, V4 kernel={len(V3)}")
    print(f"|PSL2(Z/4)|={len(G4)}, quotient={len(G2)}, kernel={len(K4)}")
    print(f"level4 generator orders={order(S4,4)},{order(C4,4)},{order(T4,4)}")
    print(f"level4 S3 complements={len(complements)}, complement census: {fmt_census(census(Acomp))}")
    print(f"edge census: {fmt_census(census(Aedge))}")
    print(f"P1(Z/4) census: {fmt_census(census(AP4))}")
    print(f"P1 point stabilizer orders={sorted(order(g,4) for g in p_stab)}")
    print(f"edge point stabilizer orders={sorted(order(g,4) for g in e_stab)}")
    print("six-actions equivariant=no; three-block quotients equivariant=yes; block kernel=V4")
    print("PASS")


if __name__ == "__main__":
    main()
