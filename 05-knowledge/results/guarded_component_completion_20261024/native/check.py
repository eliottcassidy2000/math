#!/usr/bin/env python3
"""Ten- and twelve-control native owner stars with a protected charge reservoir.

Reconstructs actual cells, shared arcs, all full native record lifts, complete
pin domains, original normalized capacity, and source/partner host gluing.
No graph enumeration or resource API is changed.
"""
from array import array
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import argparse, csv, gzip, json, sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
INHERITED = RESULTS/'compatible_charged_helpers_20261023/supply'
sys.path.insert(0, str(INHERITED))
import shared_owner as native

CHECKS = 0
U, V = 1965, 2234
SAME = [368, 556, 753, 1369, 2066]
OPPOSITE = [2135, 2397, 2439, 2483, 3187]
TEN = SAME+OPPOSITE


def need(ok, reason):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(reason)


def digest(path):
    b = path.read_bytes()
    if path.suffix in ('.py', '.json', '.tsv', '.md'):
        b = b.replace(b'\r\n', b'\n')
    return sha256(b).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    supply_path = INHERITED/'selected_supply.json'
    cert = json.loads(supply_path.read_text())
    old_path = RESULTS/cert['old_domain']
    old = json.loads(old_path.read_text())['sources']
    nodes_path = RESULTS/'metagraph_order8_20261023/merged_nodes.tsv'
    nodes = list(csv.DictReader(nodes_path.open(), delimiter='\t'))
    merged = {int(z[k]): int(z['merged_id']) for z in nodes for k in ('class_a', 'class_b')}
    owner_path = RESULTS/'metagraph_order8_20261023/owner_map.u16.gz'
    owner = array('H')
    owner.frombytes(gzip.decompress(owner_path.read_bytes()))
    if sys.byteorder != 'little':
        owner.byteswap()
    epath = RESULTS/'metagraph_order8_electrical_20261023/receipt.json'
    electrical = json.loads(epath.read_text())['order8']
    X = electrical['higher_layers'][-1]['lifts'][0]
    Y = electrical['radical_lifts_mod4']
    potentials, moduli = [X]+Y[:5], [256]+[4]*5

    def label(s):
        return ((X[s] % 128) << 4) | sum((Y[j][s] % 2) << (j-1) for j in range(1, 5))

    def charge(ss):
        result = 0
        for s in ss:
            result = native.plus(result, label(s))
        return result

    data = {}
    for z in old:
        choices = {r['target']: r for r in z['exact_all_shared_clocks_target_domain']}
        data[z['source']] = dict(records=[choices[U], choices[V]], fixed=z['fixed_clocks'])
    cores = sorted(data)
    helper_map = {z['source']: z for z in cert['helpers']}
    helpers = sorted(helper_map)
    for s, z in helper_map.items():
        m = z['selected']
        data[s] = dict(records=[m['before'], m['after']], fixed=m['common_clocks'])
    matrices = {s: native.matrix(int(nodes[s]['canon_a'])) for s in data}
    paths = {s: [r['path'] for r in d['records']] for s, d in data.items()}
    costs = {s: sum(a != b for a, b in zip(*paths[s])) for s in data}
    need(len(cores) == 28 and len(helpers) == 53 and len(data) == 81, 'actual protected28 plus new53 universe')
    bank = cert['banks'][cert['main_bank_index']]
    reservoir = bank['source_ids']
    outside = sorted(set(helpers)-set(reservoir))
    need(len(reservoir) == 40 and outside == sorted(TEN+[1952, 2638, 2948]), 'exact protected-reservoir complement')

    # Independently rebuild all40-choose3 residual charge witnesses.
    sections, witnesses = {}, [[] for _ in range(2048)]
    for triple in combinations(reservoir, 3):
        q = charge(triple)
        sections.setdefault(q, triple)
        witnesses[q].append(triple)
    need(set(sections) == set(range(2048)) and sum(map(len, witnesses)) == 9880, 'all residual charges have actual distinct triples')
    for q in range(2048):
        bm = bank['charge_section'][str(q)]
        supplied = tuple(reservoir[i] for i in range(40) if bm >> i & 1)
        need(supplied == sections[q] and charge(supplied) == q, 'protected canonical section independently agrees')

    # Authenticate entire native option reservoir once, with no frame quotient.
    all_words = set()
    for s, d in sorted(data.items()):
        A = matrices[s]
        for bit, r in enumerate(d['records']):
            p, m = r['path'], r['mask']
            need(sorted(p) == list(range(8)) and all(A[a][b] for a, b in zip(p, p[1:])), 'literal native Hamiltonian path')
            need(native.word(A, p) == m and merged[owner[m]] == s, 'literal source owner and normalized word')
            K = native.partner(A, p)
            need(native.word(K, p) == native.FULL ^ m and merged[owner[native.FULL ^ m]] == (U, V)[bit],
                 'literal partner and target class')
            op = [list(row) for row in zip(*A)]
            need(native.word(op, p[::-1]) == native.reflect(m), 'literal converse reversal')
            need(all(p[t] == value for t, value in d['fixed']), 'every inherited source common pin')
            need((X[s]-X[(U, V)[bit]]) % 8 == r['phase'], 'selected original mod8 phase')
            words = native.orbit(m)
            need(len(words) == 4 and not all_words & words, 'capacity-one full F/R option resources')
            all_words |= words
    need(len(all_words) == 648, 'complete81-slot option reservation162units648words')

    # The source176 star needs nonconsecutive shared pairs in BOTH options.
    eligibility = {}
    for s in outside:
        candidates = []
        for clock, (a, b) in enumerate(zip(*paths[s])):
            if a != b and all(abs(p.index(a)-p.index(b)) > 1 for p in paths[s]):
                candidates.append([clock, a, b, matrices[s][a][b]])
        eligibility[s] = candidates
    need([s for s in outside if not eligibility[s]] == [1952, 2638], 'exact two exceptions to the nonconsecutive star seam')
    need(sum(costs[s] for s in TEN) == 61 and costs[1952]+costs[2638] == 11,
         'literal primary control costs61 plus11')

    # Every selected special cell keeps ALL lifts above its chosen normalized
    # units, and the original pins+target margin still select exactly two.
    special = sorted({176, 1532} | set(TEN) | {1952, 2638})
    kernels = {}
    for s in special:
        A, d = matrices[s], data[s]
        lifts, pinned, hp_count, autos = [[], []], [], 0, 0
        for p in permutations(range(8)):
            if all(A[i][j] == A[p[i]][p[j]] for i, j in native.PAIRS):
                autos += 1
            if not all(A[a][b] for a, b in zip(p, p[1:])):
                continue
            hp_count += 1
            m = native.word(A, p)
            t = merged[owner[native.FULL ^ m]]
            for bit, rec in enumerate(d['records']):
                if m in native.orbit(rec['mask']):
                    lifts[bit].append(dict(path=list(p), word=m, target=t))
            if t in (U, V) and all(p[j] == v for j, v in d['fixed']):
                pinned.append((tuple(p), m, t))
        need(set(pinned) == {(tuple(r['path']), r['mask'], r['target']) for r in d['records']},
             'complete two-target native menu under all inherited indexed pins')
        kernels[s] = dict(H=hp_count, automorphism_order=autos, full_canonical_resource_lifts=lifts)
    need(kernels[2439]['automorphism_order'] == 1 and any(len(v) > 1 for v in kernels[2439]['full_canonical_resource_lifts']),
         'self-converse full resource lifts retained without capacity multiplication')

    def build_family(twelve):
        name = 'twelve_two_components' if twelve else 'ten_star'
        controls = sorted(TEN+([1952, 2638] if twelve else []))
        fixed_u = sorted(set(outside)-set(controls))
        # Root176 stays literally0..7. Other root/cell frames use private owners.
        embeddings = {176: list(range(8))}
        next_owner, guards = 8, []

        def attach(s, root, root_clock, helper_clock):
            nonlocal next_owner
            cp, hp = paths[root], paths[s]
            a, b = cp[0][root_clock], cp[1][root_clock]
            c, d = hp[0][helper_clock], hp[1][helper_clock]
            need(a != b and c != d, 'each new shared pin genuinely changes its owner')
            epsilon = int(matrices[root][a][b] != matrices[s][c][d])
            mapped = {c: embeddings[root][(a, b)[epsilon]], d: embeddings[root][(b, a)[epsilon]]}
            embed = []
            for owner_id in range(8):
                if owner_id not in mapped:
                    mapped[owner_id] = next_owner
                    next_owner += 1
                embed.append(mapped[owner_id])
            embeddings[s] = embed
            alpha = [int(abs(p.index(a)-p.index(b)) == 1) for p in cp]
            beta = [int(abs(p.index(c)-p.index(d)) == 1) for p in hp]
            need(all(alpha[v] == beta[v ^ epsilon] for v in range(2)), 'source edge gauge and partner adjacency profiles both agree')
            guards.append(dict(helper=s, root=root, root_clock=root_clock, helper_clock=helper_clock,
                bit_xor=epsilon, root_pair=[a, b], helper_pair=[c, d], root_adjacency=alpha,
                helper_adjacency=beta, global_shared_pair=sorted([embeddings[root][a], embeddings[root][b]])))

        for s in TEN:
            clock, a, b, orientation = eligibility[s][0]
            epsilon = int(matrices[176][2][6] != orientation)
            desired = int(s in OPPOSITE)
            attach(s, 176, 2 if epsilon == desired else 6, clock)
            need(guards[-1]['bit_xor'] == desired, 'five same and five opposite176 controls')
        if twelve:
            embeddings[1532] = list(range(next_owner, next_owner+8))
            next_owner += 8
            attach(1952, 1532, 2, 3)
            attach(2638, 1532, 2, 4)
            need([g['bit_xor'] for g in guards[-2:]] == [0, 1], 'one same and one opposite1532 control')
        roots = sorted({g['root'] for g in guards})
        for s in sorted(set(data)-set(embeddings)):
            embeddings[s] = list(range(next_owner, next_owner+8))
            next_owner += 8
        need(next_owner == 81*8-2*len(controls), 'exact whole source host owner count')
        parent = {g['helper']: g['root'] for g in guards}
        for a, b in combinations(sorted(data), 2):
            ra, rb = parent.get(a, a), parent.get(b, b)
            expected = 2 if ra == rb and ra in roots else 0
            need(len(set(embeddings[a]) & set(embeddings[b])) == expected,
                 'complete host intersection census has only declared component seams')

        def glue(cell_matrices):
            prescribed = {}
            for s, A in cell_matrices.items():
                emb = embeddings[s]
                for i, j in native.PAIRS:
                    a, b = emb[i], emb[j]
                    edge = min(a, b), max(a, b)
                    value = A[i][j] if a < b else A[j][i]
                    need(edge not in prescribed or prescribed[edge] == value, 'ALL common host arcs coherent')
                    prescribed[edge] = value
            return prescribed

        source_host = glue(matrices)

        def arc(u, v):
            edge = min(u, v), max(u, v)
            a = source_host.get(edge, 1)
            return a if u < v else 1-a

        global_records = {}
        for s, d in data.items():
            global_records[s] = []
            for bit, r in enumerate(d['records']):
                p = [embeddings[s][a] for a in r['path']]
                need(all(arc(a, b) for a, b in zip(p, p[1:])), 'actual induced-cell host path')
                need(sum(arc(p[i], p[j]) << k for k, (i, j) in enumerate(native.CHORDS)) == r['mask'],
                     'whole owner embedding preserves normalized native word')
                need(all(p[t] == embeddings[s][v] for t, v in d['fixed']), 'all original pins transported by fixed injective frames')
                global_records[s].append(p)
        for g in guards:
            for rb, hb in product(range(2), repeat=2):
                p, q = global_records[g['root']][rb], global_records[g['helper']][hb]
                actual = p[g['root_clock']] == q[g['helper_clock']]
                need(actual == (hb == rb ^ g['bit_xor']), 'complete actual shared-owner guard table')
                need((p[::-1][7-g['root_clock']] == q[::-1][7-g['helper_clock']]) == actual,
                     'literal converse reflected clocks preserve guard')

        def controlled(core):
            return {g['helper'] for g in guards if int(g['root'] in core) ^ g['bit_xor']}

        patterns = []
        for bits in product(range(2), repeat=len(roots)):
            active_roots = {r for r, b in zip(roots, bits) if b}
            active_control = controlled(active_roots)
            need(len(active_control) == len(controls)//2, 'balanced controls preserve exact helper margin for EVERY root pattern')
            assignment = active_roots | active_control
            # Reservoir and unrelated private cells can use EITHER option: they
            # have no overlaps. Checking each choice verifies both host extensions.
            for private_bit in range(2):
                target_cells = {}
                for s in data:
                    bit = int(s in assignment) if s in set(roots) | set(controls) else private_bit
                    target_cells[s] = native.partner(matrices[s], paths[s][bit])
                glue(target_cells)
            witnesses_sha = sha256()
            for q in range(2048):
                residual = native.minus(q, charge(active_roots | active_control))
                chosen = set(sections[residual])
                need(not chosen & set(controls) and len(chosen) == 3, 'protected reservoir has no forced-control overlap')
                need(charge(active_roots | active_control | chosen) == q, 'every root pattern and totalQ has an exact compatible section')
                witnesses_sha.update(json.dumps(sorted(chosen), separators=(',', ':')).encode())
            patterns.append(dict(root_bits=list(bits), control_V=sorted(active_control),
                                 complete_section_values=2048, section_sha256=witnesses_sha.hexdigest()))

        # Sharp boundary cost of the forced controls, independent of charge.
        boundary = []
        for left, right in product(product(range(2), repeat=len(roots)), repeat=2):
            ll = {r for r, b in zip(roots, left) if b}
            rr = {r for r, b in zip(roots, right) if b}
            moved = sorted(controlled(ll) ^ controlled(rr))
            boundary.append(dict(from_bits=list(left), to_bits=list(right), changed_controls=moved,
                                 forced_slot_cost=len(moved), forced_clock_cost=sum(costs[s] for s in moved)))
        # Exact worst support/clock bound for ANY two-core-slot exchange.
        maximum_support, maximum_clock, bound_witness = 0, 0, None
        for ca, cb in combinations(cores, 2):
            changed_controls = [g['helper'] for g in guards if g['root'] in (ca, cb)]
            support = 2+len(changed_controls)+6
            clocks = costs[ca]+costs[cb]+sum(costs[s] for s in changed_controls)+sum(sorted((costs[s] for s in reservoir), reverse=True)[:6])
            maximum_support = max(maximum_support, support)
            if clocks > maximum_clock:
                maximum_clock, bound_witness = clocks, [ca, cb]
        need(maximum_support == (20 if twelve else 18) and maximum_clock == (120 if twelve else 109),
             'exact core-menu uniform support and source-specific clock upper bound')

        # Explicit returning route, including the coupled two-root exchange.
        initial, final = {1394, 2185}, {1477, 1841}
        need(charge(initial) == charge(final), 'authenticated equal-charge core endpoints')
        core_route = ([initial, {176, 2185}, {1532, 2185}, {176, 1532}, {176, 1841}, final]
                      if twelve else [initial, {176, 2185}, {176, 1841}, final])
        first_reservoir = {650, 676, 760}
        total_q = charge(initial | controlled(initial) | first_reservoir)
        route = []
        for i, core in enumerate(core_route):
            cc = controlled(core)
            free = first_reservoir if i in (0, len(core_route)-1) else set(sections[native.minus(total_q, charge(core | cc))])
            route.append((core, cc | free))
        need(route[0][1] == route[-1][1], 'all literal helper records return exactly')

        def rows(core, active_helpers):
            active = core | active_helpers
            result = {s: [0]*6 for s in list(data)+[U, V]}
            for s in data:
                t = V if s in active else U
                for j, p in enumerate(potentials):
                    value = 2*(p[s]-p[t])
                    result[s][j] += value
                    result[t][j] -= value
            return result

        baseline = rows(*route[0])
        states, events = [], []
        for core, active_helpers in route:
            need(len(core) == 2 and len(active_helpers) == 3+len(controls)//2, 'every actual arrival exact integer margins')
            need(charge(core | active_helpers) == total_q and not active_helpers & set(fixed_u), 'fullQ and all fixed-U controls')
            words, actual_owners = set(), []
            for g in guards:
                rp = global_records[g['root']][int(g['root'] in core)]
                hp = global_records[g['helper']][int(g['helper'] in active_helpers)]
                need(rp[g['root_clock']] == hp[g['helper_clock']], 'every literal dynamic owner equation at arrival')
                actual_owners.append(rp[g['root_clock']])
            for s in data:
                bit = int(s in core or s in active_helpers)
                ww = native.orbit(data[s]['records'][bit]['mask'])
                need(not words & ww, 'no actual arrival resource reuse')
                words |= ww
            need(len(words) == 324, '81 occupied native units324closedwords')
            current = rows(core, active_helpers)
            need(all((current[s][j]-baseline[s][j]) % moduli[j] == 0 for s in current for j in range(6)),
                 'all original physical paired dyadic generator rows, including target diagonals')
            states.append(dict(core_V=sorted(core), helper_V=sorted(active_helpers),
                reservoir_V=sorted(active_helpers & set(reservoir)), shared_owner_values=actual_owners,
                physical_delta_rows={str(s): [current[s][j]-baseline[s][j] for j in range(6)]
                                     for s in current if current[s] != baseline[s]}))
        for (ca, ha), (cb, hb) in zip(route, route[1:]):
            changed = sorted((ca ^ cb) | (ha ^ hb))
            clocks = sum(costs[s] for s in changed)
            need(len(ca ^ cb) == 2 and len(changed) <= maximum_support and clocks <= maximum_clock,
                 'actual event satisfies full support and clock bounds')
            events.append(dict(changed_sources=changed, support=len(changed), owner_clock_occurrences=clocks))
        return dict(name=name, control_count=len(controls), guards=guards, helper_occupancy=3+len(controls)//2,
            fixed_U_helpers=fixed_u, embeddings={str(s): embeddings[s] for s in sorted(data)},
            source_host_order=next_owner, source_host_rule='prescribed induced cell arcs, then lower owner to higher',
            source_host_prescribed_edge_sha256=sha256(json.dumps(sorted((a, b, v) for (a, b), v in source_host.items())).encode()).hexdigest(),
            root_sources=roots, root_patterns=patterns, forced_boundary_costs=boundary,
            uniform_changed_slot_bound=maximum_support, uniform_primary_clock_bound=maximum_clock,
            clock_bound_witness_core_pair=bound_witness,
            returning_route=dict(total_charge=total_q, states=states, events=events),
            all_fixed_global_clocks={str(s): [[j, embeddings[s][v]] for j, v in data[s]['fixed']] for s in sorted(data)})

    families = [build_family(False), build_family(True)]
    # Minimal geometric hostile: a source176/star attachment to either excluded
    # helper cannot preserve all target-host arcs at a changing owner clock.
    hostile = []
    for s in (1952, 2638):
        for clock, (a, b) in enumerate(zip(*paths[s])):
            if a == b:
                continue
            pattern = [int(abs(p.index(a)-p.index(b)) == 1) for p in paths[s]]
            need(any(pattern), 'at least one K partner disagrees with always-reversed176seam')
            hostile.append(dict(helper=s, helper_clock=clock, helper_pair=[a, b], adjacency_profile=pattern))
    output = dict(status='FINITE-EXACT native ten/twelve moving-owner controls with protected completeQ reservoir',
        original_tournament_order=8, slot_count=81, source_helpers=helpers, protected_core_sources=cores,
        reservoir_sources=reservoir, residual_triples=9880, charge_values=2048,
        source176_eligibility={str(s): eligibility[s] for s in outside},
        original_special_kernels={str(s): kernels[s] for s in special}, families=families,
        source176_seam_hostiles=hostile, reserved_normalized_units=162, reserved_closed_words=648,
        input_sha256={p.relative_to(RESULTS).as_posix(): digest(p) for p in
            (supply_path, old_path, nodes_path, owner_path, epath, INHERITED/'shared_owner.py')},
        scope='Local order8 native electrical resources, explicit source and partner host extensions; no global host Hamiltonian observer or all-order supply theorem.', checks=CHECKS)
    target = HERE/'receipt.json'
    if args.write:
        target.write_text(json.dumps(output, indent=2, sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    else:
        need(json.loads(target.read_text()) == output, 'frozen entire native control certificate matches')
    print(f'PASS {CHECKS:,}; '+('wrote receipt' if args.write else 'receipt matches'))
    for family in families:
        print(family['name'], family['uniform_changed_slot_bound'], family['uniform_primary_clock_bound'],
              [(e['support'], e['owner_clock_occurrences']) for e in family['returning_route']['events']])


if __name__ == '__main__':
    main()
