"""Exact one-chord tribonacci grammar and boundary controls.

Run from the repository root:
  python 05-knowledge/results/defect_boundary_20261020/grammar_check.py
  python -O 05-knowledge/results/defect_boundary_20261020/grammar_check.py

Standard library only; all checks remain active under -O.  The inherited
bitset routines provide one independent subset-DP path count; the new
run-composition enumerator is checked as a set of actual vertex paths.
"""
from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations, product
from pathlib import Path

BASE = Path(__file__).resolve().parents[1] / "metagraph_growth_20261020/growth_check.py"
SPEC = spec_from_file_location("inherited_growth", BASE)
OLD = module_from_spec(SPEC)
SPEC.loader.exec_module(OLD)
arc, need, tau = OLD.arc, OLD.need, OLD.tau


def r_tournament(n):
    return OLD.tournament(n, 0, fixed_path=True)


def flip(a, u, v):
    b = list(a)
    b[u] ^= 1 << v
    b[v] ^= 1 << u
    return tuple(b)


def s_tournament(m, r):
    return flip(r_tournament(m + r + 3), m, m + 2)


def j_tournament(m, r):
    n = m + r + 3
    a = tuple(sum(1 << v for v in range(u + 1, n)) for u in range(n))
    return flip(a, m, m + 2)


def sequences(n):
    h, a, f = [1, 1, 1], [1, 0, 1], [0, 1, 2]
    for _ in range(3, n + 1):
        for seq in (h, a, f):
            seq.append(sum(seq[-3:]))
    return h, a, f


H, A, F = sequences(100)


def formula(m, r):
    return (H[m + r + 3] + F[m] * F[r + 1] + F[m + 1] * F[r]
            - H[r] * A[m] - A[r] * H[m])


def hamilton_paths(a):
    n = len(a)
    if not n:
        return {()}
    found = set()

    def visit(path, unused):
        if not unused:
            found.add(path)
            return
        for v in range(n):
            if unused >> v & 1 and (not path or arc(a, path[-1], v)):
                visit(path + (v,), unused ^ (1 << v))

    visit((), (1 << n) - 1)
    return found


def ordinary_runs(n, offset=0):
    """Old March grammar: descending consecutive blocks, no adjacent 1s."""
    if n == 0:
        return {()}
    found = set()

    def visit(left, lengths):
        if not left:
            cuts, blocks = offset, []
            for length in lengths:
                blocks.append(tuple(range(cuts, cuts + length)))
                cuts += length
            found.add(tuple(v for block in reversed(blocks) for v in block))
            return
        for length in range(1, left + 1):
            if length == 1 and lengths and lengths[-1] == 1:
                continue
            visit(left - length, lengths + (length,))

    visit(n, ())
    return found


def special_runs(m, r):
    """New paths using the forward chord; two orders of the hole/run."""
    n, center = m + r + 3, m + 1
    found = set()
    for left in range(1, center + 1):
        for right in range(1, n - center):
            run = tuple(range(center - left, center)) + tuple(
                range(center + 1, center + right + 1))
            middles = []
            if left >= 2:
                middles.append((center,) + run)
            if right >= 2:
                middles.append(run + (center,))
            for top in ordinary_runs(n - center - right - 1, center + right + 1):
                for bottom in ordinary_runs(center - left):
                    for middle in middles:
                        found.add(top + middle + bottom)
    return found


def symmetries(a, reverse=False):
    n = len(a)
    degrees = [row.bit_count() for row in a]
    groups = {}
    for v, d in enumerate(degrees):
        groups.setdefault(d, []).append(v)
    domains = sorted(groups)
    target_groups = [groups.get(n - 1 - d if reverse else d, [])
                     for d in domains]
    if any(len(groups[d]) != len(g) for d, g in zip(domains, target_groups)):
        return []
    answer = []
    for blocks in product(*(permutations(g) for g in target_groups)):
        mapping = [None] * n
        for d, block in zip(domains, blocks):
            for u, v in zip(groups[d], block):
                mapping[u] = v
        if all(arc(a, mapping[u], mapping[v]) == (1 - arc(a, u, v)
                                                 if reverse else arc(a, u, v))
               for u, v in OLD.pairs(n)):
            answer.append(tuple(mapping))
    return answer


def compatible_paths(a, involution):
    n = len(a)
    fixed = [u for u in range(n) if involution[u] == u]
    pairs = [(u, involution[u]) for u in range(n) if u < involution[u]]
    need(len(fixed) == n % 2, "anti-involution fixed vertices")
    count = 0
    for order in permutations(pairs):
        for choices in product((0, 1), repeat=len(pairs)):
            first = tuple(pair[c] for pair, c in zip(order, choices))
            path = first + tuple(fixed) + tuple(involution[u] for u in reversed(first))
            if all(arc(a, u, v) for u, v in zip(path, path[1:])):
                count += 1
    return count


def check_path_cover(edges, n):
    outgoing = Counter(u for u, v in edges)
    incoming = Counter(v for u, v in edges)
    need(max(outgoing.values(), default=0) <= 1, "outdegree exceeds one")
    need(max(incoming.values(), default=0) <= 1, "indegree exceeds one")
    return tuple((incoming[v], outgoing[v]) for v in range(n))


def main():
    print("PROVED formulas / FINITE-EXACT independent controls: one-chord grammar.")
    print("Universe: S_mr = R_(m+r+3) with (m,m+2) flipped; 0-based script labels.")
    print("Inherited R_n tribonacci is March 5 S11, recovered in THM-4731.")
    for n in range(40):
        need(F[n] == sum(H[:n]), "prefix-sum F")
    for n in range(1, 10):
        paths = hamilton_paths(r_tournament(n))
        need(paths == ordinary_runs(n), "March composition bijection")
        need(len(paths) == H[n], "R tribonacci")
        need(sum(p[0] != n - 1 for p in paths) == A[n], "a boundary condition")
        need(sum(p[-1] != 0 for p in paths) == A[n], "dual a condition")
    print("h_0..9 =", H[:10])
    print("a_0..9 =", A[:10], "(a_0=1 is the empty-boundary convention)")
    print("F_0..9 =", F[:10])

    subset_cases, set_cases = 0, 0
    for n in range(3, 13):
        row = []
        for m in range(n - 2):
            r = n - m - 3
            a = s_tournament(m, r)
            count = OLD.hp_subset_dp(a)
            need(count == formula(m, r), "one-chord formula / subset DP")
            row.append(count)
            subset_cases += 1
            if n <= 8:
                actual = hamilton_paths(a)
                old = ordinary_runs(n)
                lost = {p for p in old if (m + 2, m) in tuple(zip(p, p[1:]))}
                gained = special_runs(m, r)
                need(len(lost) == H[r] * A[m] + A[r] * H[m], "lost formula")
                need(len(gained) == F[m] * F[r + 1] + F[m + 1] * F[r], "gain formula")
                need(actual == (old - lost) | gained, "path-set decomposition")
                need(not gained & old, "new and old paths disjoint")
                set_cases += 1
        print("n=%d H(S_mr)=%s" % (n, row))
    print("Independent subset-DP cases:", subset_cases, "; actual path-set cases:", set_cases)

    coefficients = {(0, 0): 1, (1, 0): 4, (0, 1): 4, (2, 0): 3,
                    (1, 1): 4, (0, 2): 3, (2, 1): 2, (1, 2): 2, (2, 2): 1}
    d = (1, -1, -1, -1)
    for m in range(31):
        for r in range(31):
            need(formula(m, r) == formula(r, m), "reflection")
            residue = sum(d[i] * d[j] * formula(m - i, r - j)
                          for i in range(min(3, m) + 1)
                          for j in range(min(3, r) + 1))
            need(residue == coefficients.get((m, r), 0), "bivariate numerator")
            if m >= 3:
                need(formula(m, r) == sum(formula(m - i, r) for i in (1, 2, 3)),
                     "left tribonacci")
            if r >= 3:
                need(formula(m, r) == sum(formula(m, r - i) for i in (1, 2, 3)),
                     "right tribonacci")
    b = [[formula(m, r) for r in range(3)] for m in range(3)]
    det = (b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1])
           - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
           + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]))
    need(det == -4, "bilinear rank")
    print("B =", b, "; det(B) =", det)
    print("Rational bivariate GF and both tribonacci directions: 31 x 31 coefficient controls.")

    symmetry_cases, leaf_cases = 0, 0
    for n in range(3, 21):
        s_classes, merged_s, merged_j = set(), set(), set()
        for m in range(n - 2):
            r = n - m - 3
            a, partner = s_tournament(m, r), j_tournament(m, r)
            need(tau(a) == partner, "literal simple-core partner")
            need(OLD.transpose(a) == s_tournament(r, m), "literal reflection")
            need(len(symmetries(a)) == 1, "S rigidity")
            anti = symmetries(a, True)
            need(bool(anti) == (n <= 5 or m == r), "S self-converse classification")
            need(len(symmetries(partner)) == 3, "leaf automorphisms")
            need(bool(symmetries(partner, True)) == (m == r), "leaf self-converse")
            if n <= 10:
                need(OLD.hp_subset_dp(partner) == 3, "leaf H=3")
                leaf_cases += 1
            symmetry_cases += 1
            c, cr = OLD.canonical(a), OLD.canonical(OLD.transpose(a))
            s_classes.add(c)
            merged_s.add(min(c, cr))
            merged_j.add(min(OLD.canonical(partner),
                             OLD.canonical(OLD.transpose(partner))))
        if n >= 6:
            need(len(s_classes) == n - 2, "distinct ordinary partner classes")
            need(len(merged_s) == (n - 1) // 2, "merged partner count")
            need(len(merged_j) == (n - 1) // 2, "merged leaf count")
    print("Symmetry / literal tau / class controls:", symmetry_cases, "(all n=3..20).")
    print("Independent H=3 leaf subset-DP controls:", leaf_cases)
    for n in range(3, 6):
        blue = []
        for m in range(n - 2):
            a = s_tournament(m, n - m - 3)
            anti = symmetries(a, True)
            blue.append(compatible_paths(a, anti[0]))
        print("Small n=%d blue multiplicities=%s" % (n, blue))
    for m in range(7):
        n = 2 * m + 3
        a, old = s_tournament(m, m), r_tournament(n)
        reflection = tuple(reversed(range(n)))
        new_blue = compatible_paths(a, reflection)
        old_fixed = compatible_paths(old, reflection)
        need(new_blue == old_fixed == H[m + 2], "central blue unchanged")
        print("Central n=%d: H=%d; blue=%d" % (n, formula(m, m), new_blue))

    a = r_tournament(6)
    left = flip(flip(a, 0, 2), 2, 4)
    right = flip(flip(a, 1, 3), 2, 4)
    scores = lambda t: sorted(row.bit_count() for row in t)
    need(scores(left) == scores(right) == [1, 2, 2, 3, 3, 4], "score hostile")
    need(OLD.hp_subset_dp(left) == 25 and OLD.hp_subset_dp(right) == 29,
         "two-defect connectivity hostile")
    print("Range-two hostile: two defects and scores [1,2,2,3,3,4], but H=25 versus 29.")

    # Named vertices c=0,a=1,u=2,v=3,b=4,d=5.
    cover_a = {(1, 2), (2, 4), (0, 3), (3, 5)}
    cover_b = {(1, 2), (2, 5), (0, 3), (3, 4)}
    signature = check_path_cover(cover_a, 6)
    need(signature == check_path_cover(cover_b, 6), "degree-only boundary collision")
    ambient = r_tournament(6)
    for u, v in [(0, 3), (2, 4), (2, 5), (3, 5)]:
        ambient = flip(ambient, u, v)
    need(all(arc(ambient, u, v) for u, v in cover_a | cover_b | {(4, 0)}),
         "common bounded-range ambient")
    one_path = (1, 2, 4, 0, 3, 5)
    need(set(zip(one_path, one_path[1:])) == cover_a | {(4, 0)}, "connected completion")
    need({(0, 3), (3, 4), (4, 0)} <= cover_b | {(4, 0)}, "cycle completion")
    print("Boundary hostile: identical per-vertex in/out degrees; same next arc gives HP versus 3-cycle.")
    print("All controls passed. No general bounded-range automaton or A000568 recurrence is claimed.")


if __name__ == "__main__":
    main()
