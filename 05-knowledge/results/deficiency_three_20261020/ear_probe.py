"""Exact three-black detours for the failed F5 terminal-ear guard."""
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path

spec = spec_from_file_location("ear_black", Path(__file__).with_name("black_check.py"))
K = module_from_spec(spec)
spec.loader.exec_module(K)
B, G, need = K.B, K.G, K.need
F5 = (2, 4, 9, 19, 7)


def path_mask(p, index):
    return sum(1 << index[tuple(sorted(e))] for e in zip(p, p[1:]))


def directed_path_from_mask(rows, mask, pairs):
    n = len(rows)
    if mask.bit_count() != n - 1:
        return None
    adjacency = [[] for _ in range(n)]
    for bit, (u, v) in enumerate(pairs):
        if mask >> bit & 1:
            adjacency[u].append(v)
            adjacency[v].append(u)
    ends = [v for v, nb in enumerate(adjacency) if len(nb) == 1]
    if len(ends) != 2 or any(len(nb) not in (1, 2) for nb in adjacency):
        return None
    p, previous = [ends[0]], -1
    while len(p) < n:
        choices = [v for v in adjacency[p[-1]] if v != previous]
        if len(choices) != 1 or choices[0] in p:
            return None
        previous = p[-1]
        p.append(choices[0])
    p = tuple(p)
    if B.is_path(rows, p):
        return p
    if B.is_path(rows, p[::-1]):
        return p[::-1]
    return None


def search():
    source = B.ordered_join([G.transpose(F5), F5])
    target_path = tuple(range(10))
    target = K.reverse_path(source, target_path)
    need(B.blue(source, target_path), "source standard path must be blue")
    pairs = G.pairs(10)
    index = {e: j for j, e in enumerate(pairs)}
    target_mask = path_mask(target_path, index)
    first = [(p, K.reverse_path(source, p), path_mask(p, index))
             for p in B.hamilton_paths(source) if not B.blue(source, p)]
    last = [(p, K.reverse_path(target, p), path_mask(p, index))
            for p in B.hamilton_paths(target) if not B.blue(target, p)]
    detours = []
    for p1, after1, mask1 in first:
        for reverse3, before3, mask3 in last:
            p2 = directed_path_from_mask(after1, target_mask ^ mask1 ^ mask3, pairs)
            if p2 is None or B.blue(after1, p2):
                continue
            need(K.reverse_path(after1, p2) == before3, "middle labelled match")
            moves = (p1, p2, reverse3[::-1])
            K.verify(source, target_path, moves)
            detours.append(moves)
    print("F5 Hamiltonian paths=", len(B.hamilton_paths(F5)),
          "source_black_paths=", len(first), "target_black_paths=", len(last),
          "pair_tests=", len(first) * len(last), "three_black_detours=", len(detours))
    need(detours, "bounded search found no detour")
    print("first_exact_detour=", detours[0])
    return detours


def insert_middle(p, cut, m, forward):
    rename = lambda v: v if v < 5 else v + m
    middle = tuple(range(5, 5 + m))
    if not forward:
        middle = middle[::-1]
    return tuple(map(rename, p[:cut + 1])) + middle + tuple(map(rename, p[cut + 1:]))


def find_uniform(detours):
    tested = 0
    for moves in detours:
        choices = []
        for j, p in enumerate(moves):
            choices.append([i for i, (u, v) in enumerate(zip(p, p[1:]))
                            if (u < 5 <= v if j != 1 else v < 5 <= u)])
        for cuts in product(*choices):
            tested += 1
            good = True
            for m in (1, 2, 3):
                for code in range(1 << len(G.free_pairs(m))):
                    z = G.tournament(m, code, True)
                    source = B.ordered_join([G.transpose(F5), z, F5])
                    expanded = tuple(insert_middle(p, cut, m, j != 1)
                                     for j, (p, cut) in enumerate(zip(moves, cuts)))
                    try:
                        K.verify(source, tuple(range(m + 10)), expanded)
                    except RuntimeError:
                        good = False
                        break
                if not good:
                    break
            if good:
                print("opaque_middle_candidate=", moves, "insertion_cuts=", cuts,
                      "candidate_tests=", tested)
                return moves, cuts
    print("No opaque-middle template among tested candidates:", tested)
    return None


def audit_uniform(result):
    if result is None:
        return
    moves, cuts = result
    score_sums = set()
    counts = []
    for m in range(7):
        count = 0
        for code in range(1 << len(G.free_pairs(m))):
            z = G.tournament(m, code, True)
            source = B.ordered_join([G.transpose(F5), z, F5])
            expanded = tuple(insert_middle(p, cut, m, j != 1)
                             for j, (p, cut) in enumerate(zip(moves, cuts)))
            score_sums.add(tuple(K.verify(source, tuple(range(m + 10)), expanded)))
            count += 1
        counts.append((m, count))
    print("all_fixed_path_middle_controls=", counts,
          "endpoint_score_sums=", sorted(score_sums))


def family_template(r, m):
    need(r >= 4, "proved family requires r>=4")
    A = tuple(range(r + m))
    a, b, z = r + m, r + m + 1, 2 * r + m - 1
    C = tuple(range(b + 1, z))
    return (A + (b,) + C + (z, a),
            (a, z, b) + A[::-1] + C[::-1],
            C + A + (a, b, z))


def family_audit():
    counts = []
    for r in range(4, 10):
        right = G.tournament(r, 0, True)
        count = 0
        for m in range(7):
            for code in range(1 << len(G.free_pairs(m))):
                middle = G.tournament(m, code, True)
                source = B.ordered_join([G.transpose(right), middle, right])
                moves = family_template(r, m)
                sums = K.verify(source, tuple(range(2 * r + m)), moves)
                need(tuple(sums) == (m + r + 2, 4, r + 1),
                     "uniform family exact endpoint charges")
                need(B.blue(source, tuple(range(2 * r + m))) ==
                     B.blue(middle, tuple(range(m))), "source blue iff middle blue")
                count += 1
        counts.append((r, count))
    tri = G.tournament(3, 0, True)
    small = B.ordered_join([G.transpose(tri), tri])
    rotated = (0, 1, 2, 4, 5, 3)
    need(B.is_path(small, rotated) and B.blue(small, rotated),
         "r=3 hostile must have a blue first move")
    print("Tstar_r_family_controls (r, all middle presentations m<=6)=", counts)
    print("Proved endpoint sums=(m+r+2,4,r+1), n-1=m+2r-1;",
          "r=3 first-move hostile is blue.")


def general_audit():
    counts = []
    for r in range(4, 7):
        eligible = rotation_auto = checked = 0
        for code in range(1 << len(G.free_pairs(r))):
            right = G.tournament(r, code, True)
            if not G.arc(right, r - 1, 0) or not G.arc(right, r - 1, 1):
                continue
            eligible += 1
            rotation = tuple(range(1, r)) + (0,)
            auto = G.induced(right, rotation) == right
            rotation_auto += auto
            for m in range(3):
                middle = G.tournament(m, 0, True)
                source = B.ordered_join([G.transpose(right), middle, right])
                target_path = tuple(range(2 * r + m))
                need(B.blue(source, target_path), "general source must be blue")
                moves = family_template(r, m)
                need(K.edges(target_path) == K.edges(moves[0]) ^
                     K.edges(moves[1]) ^ K.edges(moves[2]), "general edge ledger")
                current, colours = source, []
                for j, path in enumerate(moves):
                    need(B.is_path(current, path), "general directed move")
                    colours.append(B.blue(current, path))
                    if j:
                        score = current[path[0]].bit_count() + current[path[-1]].bit_count()
                        need(score <= 2 * r - 2 <= len(current) - 2,
                             "general endpoint charge bound")
                    current = K.reverse_path(current, path)
                need(colours == [auto, False, False], "rotation is exact first-blue criterion")
                need(current == K.reverse_path(source, target_path), "general labelled target")
                checked += 1
        counts.append((r, eligible, rotation_auto, checked))
    regular = B.cyclic(5)
    need(G.arc(regular, 4, 0) and G.arc(regular, 4, 1), "regular C5 has both guards")
    source = B.ordered_join([G.transpose(regular), regular])
    need(B.blue(source, family_template(5, 0)[0]), "regular C5 first move is blue")
    print("General guarded-endcap controls (r, eligible, rotation-auto, middle cases)=", counts)
    print("Regular cyclic C5 hostile: closing arc and chord hold, but first move is blue.")


def main():
    detours = search()
    audit_uniform(find_uniform(detours))
    family_audit()
    general_audit()
    print("PASS: full labelled matrices, exact edge parity, all three black predicates.")


if __name__ == "__main__":
    main()
