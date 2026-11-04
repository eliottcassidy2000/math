#!/usr/bin/env python3
"""Positive-defect charts and their exact bridge from balanced sources."""
from pathlib import Path
import hashlib
import json
import kernel as K

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
c = K.c
w = c.word


def need(test, message):
    if not test:
        raise RuntimeError(message)


def exterior(p):
    n = len(p)
    return tuple(p) + ((1 << n)-1,)


def exterior_data(p, b):
    a = exterior(p)
    j = (len(p), b)
    h = c.tables(a)[0][-1]
    E = sum(c.profile(a, v)[3] for v in j)
    W = K.cut_kernel(a)[2]
    C = W[j[0]][j[1]] + W[j[1]][j[0]]
    return h, E, C


def lift(p, b, k):
    need(k >= 0, "nonnegative literal prefix length")
    n = len(p)
    a = list(p) + [0]*(k+2)
    src, root = n, n+k+1
    a[src] = ((1 << n)-1) | (1 << root)
    a[b] |= 1 << root
    for i in range(n+1, root):
        a[i] = (1 << (n+1))-1
        for j in range(i+1, root):
            a[i] |= 1 << j
    a[root] = (((1 << n)-1) ^ (1 << b)) | sum(1 << i for i in range(n+1, root))
    return tuple(a), root, (src, b)


def main():
    chart = []
    for r in range(1, 5):
        P = c.graph('2'*r)
        p = c.profile(P, 1)
        u = 2**r+1
        need(p == (u, 1, 1, u), "rigid balanced singleton-deletion core")
        data = exterior_data(P, 1)
        need(data == (u, 3*u, u+1), "full exterior joint data")
        for k in range(1, 4):
            word = '3'+'0'*k+'3'*r+'1'
            a = c.graph(word)
            h = 3*u*2**k-1
            expected = (h, u, u, h+u+1)
            need(c.profile(a, 2) == expected, "literal graded chart profile")
            need(c.strong(a), "strong graded word")
            chart.append({"k": k, "r": r, "word": word, "root": 2,
                          "profile": expected, "exterior_data": data})

    core_controls = []
    for word, b in [('2', 1), ('22', 1), ('21', 0), ('112', 0), ('202032', 2)]:
        p = c.graph(word)
        h, s, q, F = c.profile(p, b)
        data = exterior_data(p, b)
        need(data == (h, h+2*F, s+F), "general exterior identities")
        core_controls.append({"word": word, "root": b, "profile": [h, s, q, F],
                              "joint_exterior": data})
    singleton_chain = []
    for word in ('300', '30031', '3003331'):
        profile = c.profile(c.graph(word), 2)
        singleton_chain.append({"word": word, "profile": profile})
    need([x['profile'] for x in singleton_chain] == [(11, 1, 1, 13), (35, 3, 3, 39), (107, 9, 9, 117)],
         "singleton boundary gives two successive plus-two steps")
    need(exterior_data((0,), 0) == (1, 3, 2), "singleton exterior boundary")

    old_source = (102, 124, 120, 113, 33, 64, 16)
    mapping_order = (3, 4, 0, 1, 2, 5, 6)
    need(c.induced(old_source, mapping_order) == c.graph('30031'), "conditional lane physical map")
    prior_base = c.graph('2032')
    before_sealing = c.birth(prior_base, 1)
    before_order = (3, 0, 1, 6, 2, 4, 5)
    need(c.induced(before_sealing, before_order) == c.graph('30031'), "actual first-reset intermediate is graded source35")
    need(c.profile(before_sealing, 1) == (35, 3, 3, 39), "retained old-x profile before sealing")
    sealed = c.birth(before_sealing, 6)
    need(c.profile(sealed, 6) == (59, 17, 17, 59), "second reset birth seals to inherited balanced59")
    hostile_base = c.graph('20032')
    hostile_intermediate = c.birth(hostile_base, 1)
    need(c.profile(hostile_intermediate, 1) == (71, 3, 3, 87), "m3 intermediate does not have graded-chart gap4")
    need(c.profile(c.graph('300031'), 2) == (71, 3, 3, 75), "same first-three coordinates do not identify actual source")
    new_word = '3003331'
    old_to_new = [0, 1, 2, 3, 4, 5, 8]
    need(c.induced(c.graph(new_word), old_to_new) == c.graph('30031'), "graded named insertion")
    named_step = {"source_rows": c.graph('30031'), "source_names": ['a', 'b', 'root', 'u_1', 'u_2', 'w_1', 'v'],
                  "target_rows": c.graph(new_word), "target_names": ['a', 'b', 'root', 'u_1', 'u_2', 'w_1', 'w_2', 'w_3', 'v'],
                  "old_to_new": old_to_new, "root": 2, "in_ports": [0, 1],
                  "source_profile": [35, 3, 3, 39], "target_profile": [107, 9, 9, 117],
                  "new_word_positions": [4, 5], "new_letters": '33',
                  "conditional_source_to_standard": [mapping_order.index(i) for i in range(7)]}

    lifts = []
    for m, r, k in [(2, 1, 1), (3, 3, 1), (6, 3, 1)]:
        word = '202'+'0'*(m-1)+'3'*r+'2'
        p = c.graph(word)
        b = 5*(2**r+1)*2**m-1
        t = (2**r+1)*(2**m+2)-1
        need(c.profile(p, 2) == (b, t, t, b), "actual balanced input")
        a, root, ports = lift(p, 2, k)
        need(c.induced(a, list(range(len(p)))) == p, "lift retains all old arcs")
        h = 3*b*2**k-t
        profile = (h, b, b, h+b+t)
        need(c.profile(a, root) == profile, "literal arbitrary-balanced lift profile")
        need(c.strong(a), "strong arbitrary-balanced lift")
        names = ['apex0', 'apex1', 'balanced_root', 'balanced_child', 'y']
        names += ['clone_'+str(i) for i in range(1, m)]
        names += ['three_'+str(i) for i in range(1, r+1)] + ['last2']
        need(len(names) == len(p), "old names complete")
        names += ['source_port']+['prefix_'+str(i) for i in range(1, k+1)]+['graded_root']
        lifts.append({"m": m, "r": r, "k": k, "old_word": word,
                      "old_balanced_pair": [b, t], "rows": a, "names": names,
                      "old_to_new": list(range(len(p))), "root": root, "in_ports": ports,
                      "profile": profile})

    clocks = []
    for r in (1, 3):
        table = []
        for m in range(3):
            b = 5*(2**r+1)*2**m-1
            t = (2**r+1)*(2**m+2)-1
            table.append([(3*b*2**k-t)%7 for k in range(3)])
        clocks.append({"r": r, "m_rows_mod3": [0, 1, 2], "k_columns_mod3": [0, 1, 2], "H_mod7": table})
    need(all(x for row in clocks[0]['H_mod7'] for x in row), "original balanced family misses seven-lift lane")
    need([(i,j) for i,row in enumerate(clocks[1]['H_mod7']) for j,x in enumerate(row) if x == 0] == [(0,1)],
         "plus-two balanced family opens exact seven clock")
    qrows = []
    for j in range(1, 15):
        h = 261*8**j-23
        need(h%7 == 0 and (h%49 == 0) == (j%7 == 4), "seven and forty-nine clocks")
        q = h//7
        need((261*8**(j+1)-23)//7 == 8*q+23, "seven quotient recurrence")
        qrows.append({"j": j, "H": h, "quotient": q, "divisible_by49": h%49 == 0})
    need(qrows[1]['H'] == 16681 and qrows[1]['quotient'] == 2383 and c.prime(2383), "actual seven-prime control")
    need(qrows[0]['quotient'] == 5*59, "hostile composite quotient")
    boundary_input = c.graph('202032')
    boundary_graph, boundary_root, boundary_ports = lift(boundary_input, 2, 0)
    boundary_profile = c.profile(boundary_graph, boundary_root)
    need(boundary_profile == (101, 0, 59, 236), "zero-prefix boundary is ordinary birth behind a global source")
    gap_controls = []
    for target in (41, 167):
        witness = w.find_count(target)
        need(witness is not None and w.hamilton(witness) == target, "actual target beyond free-u chart")
        gap_controls.append({"H": target, "word": witness, "prime": c.prime(target), "required_u": 7})
    sources = [K.c.SOURCE,
        ROOT/'05-knowledge/results/ready_return_pentagon_20261103/ATLAS.md',
        ROOT/'05-knowledge/results/two_entry_funnels_20261103/READY.md',
        ROOT/'01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md',
        HERE/'NOTE.md', HERE/'kernel.py']
    print(json.dumps({"status": "PASS", "scope": "Scoped symbolic formulas with literal named graph controls; no universal return or coverage claim.",
        "graded_chart": chart, "general_exterior_controls": core_controls,
        "singleton_two_step_path": singleton_chain,
        "previous_construction_before_sealing": {"rows": before_sealing,
            "names": ['c_2', 'x', 'y', 'c_1', 'birth_a', 'birth_b', 'reset_root'],
            "old_root": 1, "standard_vertex_order": before_order,
            "standard_names": ['c_1', 'c_2', 'x', 'reset_root', 'y', 'birth_a', 'birth_b'],
            "profile": [35, 3, 3, 39], "sealed_rows": sealed,
            "sealed_profile": [59, 17, 17, 59],
            "m3_hostile_profiles": [[71, 3, 3, 87], [71, 3, 3, 75]]},
        "named_35_to_107": named_step, "balanced_to_graded_lifts": lifts,
        "seven_clocks": clocks, "seven_quotient_controls": qrows,
        "zero_prefix_hostile": {"input_pair": [59, 17], "rows": boundary_graph,
            "root": boundary_root, "ports": boundary_ports, "profile": boundary_profile,
            "correct_H": 2*59-17, "invalid_positive_prefix_formula": 3*59-17,
            "scope": "k=0 is excluded from G1/G5 and cannot be used to obtain formal carry4"},
        "conditional_free_u_gaps": gap_controls,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "sources": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
