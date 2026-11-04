#!/usr/bin/env python3
"""Named interior repairs and complete bounded arbitrary-birth certificates."""
from pathlib import Path
import hashlib
import itertools
import json
import math
import kernel as K

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
c = K.c
w = c.word


def need(test, message):
    if not test:
        raise RuntimeError(message)


def iso(a, b):
    n = len(a)
    blocks = {d: [i for i, row in enumerate(a) if row.bit_count() == d]
              for d in set(row.bit_count() for row in a)}
    dest = {d: [i for i, row in enumerate(b) if row.bit_count() == d]
            for d in blocks}
    if any(len(blocks[d]) != len(dest[d]) for d in blocks):
        return None
    ds = sorted(blocks)
    for choices in itertools.product(*(itertools.permutations(dest[d]) for d in ds)):
        mapping = [None] * n
        for d, perm in zip(ds, choices):
            for i, j in zip(blocks[d], perm):
                mapping[i] = j
        if all(bool(a[i] >> j & 1) == bool(b[mapping[i]] >> mapping[j] & 1)
               for i in range(n) for j in range(n)):
            return mapping
    return None


def reverse(a):
    return tuple(sum(1 << j for j in range(len(a)) if a[j] >> i & 1)
                 for i in range(len(a)))


def extensions(a, offsets=(2, 4)):
    h = c.tables(a)[0][-1]
    targets = {b: 3*h+b for b in offsets}
    cap = max(targets.values())
    base = K.cut_kernel(a)
    first = []
    hits = {b: [] for b in offsets}
    hist = {}
    for mask in range(1 << len(a)):
        value = K.cut_h(base, mask)
        hist[value] = hist.get(value, 0) + 1
        for b in offsets:
            if value == targets[b]:
                hits[b].append([mask])
        if value <= cap:
            first.append([mask, value])
    for mask, _ in first:
        child = K.add(a, mask)
        data = K.cut_kernel(child)
        for mask2 in range(1 << len(child)):
            value = K.cut_h(data, mask2)
            for b in offsets:
                if value == targets[b]:
                    hits[b].append([mask, mask2])
    return {"H": h, "order": len(a), "targets": targets,
            "total_ordered_two_birth_universe": 1 << (2*len(a)+1),
            "first_masks_tested": 1 << len(a),
            "first_H_histogram": sorted(hist.items()),
            "first_masks_at_most_maximum_target": first,
            "second_masks_tested": len(first) * (1 << (len(a)+1)),
            "pruning": "H of an induced tournament never exceeds H of its extension",
            "hits": hits}


def named_extension(a, names, masks, new_names, root, child):
    source = tuple(a)
    original_names = list(names)
    steps = []
    for mask, name in zip(masks, new_names):
        need(name not in names, "fresh physical name")
        steps.append({"name": name, "old_in_neighbours": [names[i] for i in range(len(a))
                                                          if mask >> i & 1],
                      "in_mask": mask})
        a = K.add(a, mask)
        names = list(names) + [name]
    need(c.induced(a, list(range(len(source)))) == source, "all old arcs survive")
    need(names[:len(source)] == original_names, "all old names survive")
    p = c.profile(a, root)
    return {"rows": a, "names": names, "root": root, "child": child,
            "profile": p, "primitive": math.gcd(p[0], p[1]) == 1,
            "strong": c.strong(a), "steps": steps,
            "old_to_new": list(range(len(source))),
            "two_apex_markings": [{"pair": pair, "word": value[0], "chain": value[1],
                                    "seventeen_two_apex_suffix_queries": w.full_profile(value[0])}
                                   for pair, value in sorted(c.marked_words(a).items())]}


def insert_at(a, names, pair, chain, positions, letters, new_names):
    oldword = ''.join(str((a[pair[0]] >> v & 1) + 2*(a[pair[1]] >> v & 1)) for v in chain)
    holes = dict(zip(positions, zip(letters, new_names)))
    old_vertices = iter(chain)
    old_letters = iter(oldword)
    expanded = []
    expanded_names = []
    name_to_old = {name: i for i, name in enumerate(names)}
    for j in range(len(chain) + len(positions)):
        if j in holes:
            letter, name = holes[j]
        else:
            v = next(old_vertices)
            letter, name = next(old_letters), names[v]
        expanded.append(letter)
        expanded_names.append(name)
    word = ''.join(expanded)
    new_names_all = [names[pair[0]], names[pair[1]]] + expanded_names
    need(len(set(new_names_all)) == len(new_names_all), "insertion names distinct")
    graph = c.graph(word)
    old_to_new = [new_names_all.index(name) for name in names]
    need(c.induced(graph, old_to_new) == tuple(a), "named insertion induced source")
    return graph, new_names_all, old_to_new, word


def recurrence_certificate():
    # The exact inherited transfer has dimension 36. Cayley-Hamilton turns
    # the checked 36 scalar zeros into all future powers of the same matrix.
    rows = []
    terminal = w.suffix_row('2')
    for prefix in ('2', '202'):
        v0 = w.state(prefix+'0')
        v1 = w.step(v0, 0)
        v2 = w.step(v1, 0)
        need(all(z-3*y+2*x == 0 for x, y, z in zip(v0, v1, v2)), "zero-run full-state recurrence")
        controls = []
        for zero_count in (1, 2):
            v0 = w.state(prefix+'0'*zero_count+'3')
            v1 = w.step(v0, 3)
            v2 = w.step(v1, 3)
            residual = tuple(z-3*y+2*x for x, y, z in zip(v0, v1, v2))
            observed = []
            for _ in range(w.D):
                observed.append(w.dot(terminal, residual))
                residual = w.step(residual, 3)
            need(not any(observed), "three-run observed recurrence")
            controls.append(observed)
        initial = [[w.hamilton(prefix+'0'*(m-1)+'3'*r+'2') for r in (1, 2)]
                   for m in (2, 3)]
        rows.append({"prefix": prefix, "full_state_zero_run_recurrence": [2, -3, 1],
                     "terminal_three_run_recurrence": [2, -3, 1],
                     "Cayley_Hamilton_dimension": w.D,
                     "zero_residual_outputs": controls, "initial_values_m2_m3_r1_r2": initial})
    for m in range(2, 12):
        for r in range(1, 9):
            t = (2**r+1)*(2**m+2)-1
            b = 5*(2**r+1)*2**m-1
            need(w.hamilton('2'+'0'*(m-1)+'3'*r+'2') == t, "base two-variable formula")
            need(w.hamilton('202'+'0'*(m-1)+'3'*r+'2') == b, "balanced two-variable formula")
    for m in range(2, 102):
        need(math.gcd(45*2**m-1, 9*2**m+17) == (43 if m % 14 == 13 else 1),
             "plus-two primitive selector")
    need(pow(2, 7, 43) == 42 and pow(2, 14, 43) == 1,
         "order fourteen modulo forty-three")
    return rows


def three_interior_census():
    a = c.graph('20200032')
    h = 239
    hits = {2: [], 4: []}
    tested = 0
    pairs = c.marked_words(a)
    for pair, (oldword, chain) in sorted(pairs.items()):
        for k in (1, 2, 3):
            for positions in itertools.combinations(range(len(oldword)+k), k):
                holes = set(positions)
                for letters in itertools.product('0123', repeat=k):
                    oi, ni = iter(oldword), iter(letters)
                    newword = ''.join(next(ni) if j in holes else next(oi)
                                      for j in range(len(oldword)+k))
                    value = w.hamilton(newword)
                    tested += 1
                    for b in (2, 4):
                        if value == 3*h+b:
                            hits[b].append({"pair": pair, "oldword": oldword,
                                            "old_chain": chain, "newword": newword,
                                            "positions": positions, "letters": letters})
    need(tested == 33948 and len(hits[2]) == 30 and not hits[4], "complete C4 interior menu")
    triangles = ((1, 2, 3), (0, 8, 9))
    for x, y, z in triangles:
        need(a[x] >> y & 1 and a[y] >> z & 1 and a[z] >> x & 1,
             "vertex-disjoint directed triangles")
    need(set(triangles[0]).isdisjoint(triangles[1]), "triangle disjointness")
    return {"H": h, "source_word": '20200032', "source_rows": a,
            "source_names": ['c_4', 'x', 'reset_root', 'reset_child', 'y', 'c_1', 'c_2', 'c_3', 'birth_a', 'birth_b'],
            "root": 2, "child": 3, "disjoint_triangle_witness": triangles,
            "initial_markings": len(pairs), "birth_budget": 3,
            "tested_addresses": tested, "hits": hits,
            "plus_four_seven_prime_factorization": [7, 103]}


def main():
    controls = 0
    for oldword in ('2', '21', '112'):
        a = c.graph(oldword)
        data = K.cut_kernel(a)
        for mask in range(1 << len(a)):
            need(K.cut_h(data, mask) == c.tables(K.add(a, mask))[0][-1], "independent literal subset DP birth")
            controls += 1

    old_receipt_path = ROOT/'05-knowledge/results/orbit_spectrum_session_20261103/coverage/check.json'
    inherited = json.loads(old_receipt_path.read_text())
    raw = next(row['moving_mark_child']['balanced_graph'] for row in inherited['clone_rows'] if row['m'] == 2)
    raw_rows = tuple(raw['rows'])
    pair = (0, 1)
    oldword, chain = c.marked_words(raw_rows)[pair]
    order = list(pair) + list(chain)
    a = c.induced(raw_rows, order)
    names = [raw['names'][i] for i in order]
    need(oldword == '202032' and a == c.graph(oldword), "actual source transport")
    need(c.profile(a, 2) == (59, 17, 17, 59), "source balanced profile")
    source = {"rows": a, "names": names, "word": oldword, "raw_to_word": [order.index(i) for i in range(len(a))],
              "root": 2, "child": 3, "profile": [59, 17, 17, 59],
              "seventeen_two_apex_suffix_queries": w.full_profile(oldword),
              "source_receipt": str(old_receipt_path.relative_to(ROOT))}
    census = extensions(a)
    expected = [[63, 63], [63, 319], [63, 383], [127, 63], [127, 127], [127, 383]]
    need(census['hits'][2] == expected and census['hits'][4] == [], "complete H59 birth menu")
    need(len(census['first_masks_at_most_maximum_target']) == 24, "first-stage survivors")
    positives = [named_extension(a, names, masks, ['insert_a', 'insert_b'], 2, 3) for masks in expected]
    need(all(p['profile'] == (179, 53, 53, 179) and p['strong'] for p in positives), "all plus-two interfaces")

    first, names1, map1, word1 = insert_at(a, names, (0, 1), tuple(range(2, len(a))), (4, 5), '33', ['interior_a', 'interior_b'])
    second, names2, map2, word2 = insert_at(first, names1, (0, 1), tuple(range(2, len(first))), (5,), '0', ['interior_c'])
    need(word1 == '20203332' and word2 == '202030332', "literal two-step program")
    need(c.profile(first, 2) == (179, 53, 53, 179), "first retained interface")
    need(c.profile(second, 2) == (541, 147, 147, 541), "second retained interface")
    need(c.strong(first) and c.strong(second), "strong completed insertion stages")
    path = [dict(rows=first, names=names1, old_to_new=map1, word=word1, positions=[4, 5], letters='33',
                 profile=[179, 53, 53, 179], seventeen_two_apex_suffix_queries=w.full_profile(word1)),
            dict(rows=second, names=names2, old_to_new=map2, word=word2, positions=[5], letters='0',
                 profile=[541, 147, 147, 541], seventeen_two_apex_suffix_queries=w.full_profile(word2))]
    kernel179 = K.cut_kernel(first)
    one179 = {b: [mask for mask in range(1 << len(first)) if K.cut_h(kernel179, mask) == 3*179+b] for b in (2, 4)}
    need(one179 == {2: [638], 4: [124]}, "complete H179 one-birth menu")
    escape = named_extension(first, names1, [638], ['non_chain_birth'], 2, 3)
    need(escape['profile'] == (539, 157, 157, 539) and not escape['two_apex_markings'], "actual exit from two-apex language")
    second_census = extensions(second)
    need(second_census['hits'] == {2: [], 4: []}, "bounded H541 birth obstruction")

    edges = list(itertools.combinations(range(len(a)), 2))
    one_flip = []
    two_flip = []
    for e in edges:
        b = list(a)
        i, j = e
        b[i] ^= 1 << j
        b[j] ^= 1 << i
        if c.tables(b)[0][-1] == 59:
            one_flip.append(e)
    need(one_flip == [], "no neutral single arc flip")
    for e, f in itertools.combinations(edges, 2):
        b = list(a)
        for i, j in (e, f):
            b[i] ^= 1 << j
            b[j] ^= 1 << i
        if c.tables(b)[0][-1] != 59:
            continue
        menu = extensions(b, (4,))
        need(menu['hits'][4] == [], "two-flip retained obstruction")
        two_flip.append({"reversed_arcs": [e, f], "rows": b, "names": names,
                         "strong": c.strong(b), "protected_root_child_arc_retained": bool(b[2] >> 3 & 1),
                         "isomorphism_from_source": iso(a, b),
                         "isomorphism_from_converse": iso(reverse(a), b), "extension_census": menu})
    need(len(two_flip) == 5, "complete neutral two-flip fibre")
    sources = [c.SOURCE, ROOT/'05-knowledge/results/orbit_spectrum_session_20261103/coverage/check.py', old_receipt_path,
               ROOT/'05-knowledge/results/session164_emergent_20261103/nonuniform_operator/NOTE.md',
               ROOT/'05-knowledge/results/two_entry_funnels_20261103/READY.md']
    print(json.dumps({"status": "PASS", "scope": "Complete at stated birth and arc-flip budgets; all-size word formula via exact transfer recurrence.",
        "source": source, "H59_all_one_two_births": census, "H179_positive_extensions": positives,
        "named_interior_path_59_179_541": path, "H179_all_one_births": one179,
        "H539_balanced_exit_from_two_apex": escape, "H541_all_one_two_births": second_census,
        "two_variable_recurrence": recurrence_certificate(),
        "H239_complete_three_interior_births": three_interior_census(),
        "seven_clocks": {"m_residues": [0, 1, 2], "plus_two_mod7": [(45*2**m-1)%7 for m in range(3)],
                          "plus_four_mod7": [(45*2**m+1)%7 for m in range(3)]},
        "neutral_arc_reversal_census": {"one_arc_sets": len(edges), "one_arc_hits": one_flip,
            "two_arc_sets": math.comb(len(edges), 2), "two_arc_hits": two_flip},
        "independent_subset_DP_birth_controls": controls,
        "primality_controls": {str(x): c.prime(x) for x in (59, 179, 181, 541, 1627)},
        "joint_query_scope": "These 17 queries are two-apex word suffix queries; full named graphs are retained also when no two-apex presentation exists.",
        "two_apex_query_suffixes": w.FULL_SUFFIXES,
        "normal_vs_optimized": "No removable assertions carry any validation.",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "kernel_sha256": hashlib.sha256((HERE/'kernel.py').read_bytes()).hexdigest(),
        "sources": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
