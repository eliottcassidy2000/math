#!/usr/bin/env python3
"""A genuine shared-owner pin, with an independently enumerated charge section.

Two induced eight-owner cells overlap in exactly two vertices of an explicit
fourteen-owner tournament.  The declared indexed-owner equality forces their
bits to be opposite.  The host is not a global Hamiltonian observer.
"""
from array import array
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import argparse, csv, gzip, json, sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
COUNT = 0
U, V = 1965, 2234
PAIRS = list(combinations(range(8), 2))
CHORDS = [(i, j) for i in range(8) for j in range(i+2, 8)]
FULL = (1 << len(CHORDS))-1
REV = [CHORDS.index((7-j, 7-i)) for i, j in CHORDS]


def need(ok, why):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(why)


def matrix(code):
    a = [[0]*8 for _ in range(8)]
    for k, (i, j) in enumerate(PAIRS):
        a[i][j] = (code >> (27-k)) & 1
        a[j][i] = 1-a[i][j]
    return a


def word(a, p):
    return sum(a[p[i]][p[j]] << k for k, (i, j) in enumerate(CHORDS))


def reflect(m):
    return sum(((m >> k) & 1) << REV[k] for k in range(21))


def orbit(m):
    return {m, reflect(m), FULL ^ m, FULL ^ reflect(m)}


def partner(a, p):
    at = {v: i for i, v in enumerate(p)}
    return [[0 if i == j else a[i][j] if abs(at[i]-at[j]) == 1
             else a[j][i] for j in range(8)] for i in range(8)]


def host(cells, size):
    """Consistent cell restrictions, with all unprescribed arcs lower->higher."""
    a = [[0 if i == j else int(i < j) for j in range(size)] for i in range(size)]
    prescribed = {}
    for local, embedding in cells:
        for i, j in PAIRS:
            x, y = embedding[i], embedding[j]
            edge = (min(x, y), max(x, y))
            value = local[i][j] if x < y else local[j][i]
            need(edge not in prescribed or prescribed[edge] == value,
                 'all overlapping cell arcs agree')
            prescribed[edge] = value
            a[edge[0]][edge[1]] = value
            a[edge[1]][edge[0]] = 1-value
    return a


def plus(a, b):
    return ((((a >> 4)+(b >> 4)) % 128) << 4) | ((a ^ b) & 15)


def minus(a, b):
    return ((((a >> 4)-(b >> 4)) % 128) << 4) | ((a ^ b) & 15)


def digest(path):
    b = path.read_bytes()
    if path.suffix in ('.json', '.py', '.md', '.tsv'):
        b = b.replace(b'\r\n', b'\n')
    return sha256(b).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    supply_path = HERE/'selected_supply.json'
    cert = json.loads(supply_path.read_text())
    old_path = RESULTS/cert['old_domain']
    old = json.loads(old_path.read_text())['sources']
    res_path = HERE.parent/'compatibility/resilience.json'
    resilience = json.loads(res_path.read_text())
    source_bank = next(b for b in resilience['banks']
                       if b['helper_count'] == 53 and b['V_occupancy'] == 3)
    helpers = source_bank['single_U_robust_subbank']['source_ids']
    need(len(helpers) == 48 and len(set(helpers)) == 48 and 368 in helpers,
         'exact inherited48 subbank contains the shared helper')
    nodes_path = RESULTS/'metagraph_order8_20261023/merged_nodes.tsv'
    nodes = list(csv.DictReader(nodes_path.open(), delimiter='\t'))
    merged = {int(z[k]): int(z['merged_id']) for z in nodes for k in ('class_a', 'class_b')}
    owner_path = RESULTS/'metagraph_order8_20261023/owner_map.u16.gz'
    owner = array('H')
    owner.frombytes(gzip.decompress(owner_path.read_bytes()))
    if sys.byteorder != 'little':
        owner.byteswap()
    electrical_path = RESULTS/'metagraph_order8_electrical_20261023/receipt.json'
    electric = json.loads(electrical_path.read_text())['order8']
    x = electric['higher_layers'][-1]['lifts'][0]
    ys = electric['radical_lifts_mod4']
    potentials, moduli = [x]+ys[:5], [256]+[4]*5

    def label(s):
        return ((x[s] % 128) << 4) | sum((ys[j][s] % 2) << (j-1) for j in range(1, 5))

    def charge(ss):
        q = 0
        for s in ss:
            q = plus(q, label(s))
        return q

    data = {}
    for z in old:
        choices = {r['target']: r for r in z['exact_all_shared_clocks_target_domain']}
        data[z['source']] = dict(records=[choices[U], choices[V]], fixed=z['fixed_clocks'])
    cores = sorted(data)
    for z in cert['helpers']:
        if z['source'] in helpers:
            sel = z['selected']
            data[z['source']] = dict(records=[sel['before'], sel['after']], fixed=sel['common_clocks'])
    need(len(cores) == 28 and len(data) == 76, 'all old28 and new48 actual slots')
    matrices = {s: matrix(int(nodes[s]['canon_a'])) for s in data}
    embeddings = {176: list(range(8)), 368: [8, 9, 10, 2, 6, 11, 12, 13]}
    a, b = matrices[176], matrices[368]
    need(a[2][6] == 0 and b[4][3] == 1, 'source shared-edge orientation selects opposite framing')
    local_host = host([(a, embeddings[176]), (b, embeddings[368])], 14)
    need(set(embeddings[176]) & set(embeddings[368]) == {2, 6}, 'exactly two shared owners')
    need(all(local_host[i][j]+local_host[j][i] == 1 for i, j in combinations(range(14), 2)),
         'literal14-owner host is a tournament')
    next_owner = 14
    for s in sorted(set(data)-{176, 368}):
        embeddings[s] = list(range(next_owner, next_owner+8))
        next_owner += 8
    need(next_owner == 606, 'private completion of all76 source cells has606 distinct owners')

    # A compact global source host: arcs inside a prescribed cell, else lower->higher.
    prescribed = {}
    for s, embed in embeddings.items():
        aa = matrices[s]
        for i, j in PAIRS:
            u, v = embed[i], embed[j]
            edge = min(u, v), max(u, v)
            value = aa[i][j] if u < v else aa[j][i]
            need(edge not in prescribed or prescribed[edge] == value, 'full606-owner source restrictions agree')
            prescribed[edge] = value

    def arc(u, v):
        edge = min(u, v), max(u, v)
        value = prescribed.get(edge, 1)
        return value if u < v else 1-value

    # Authenticate ALL option records in their fixed global owner frames.
    reservoir, records, costs = set(), {}, {}
    for s, d in sorted(data.items()):
        aa, emb = matrices[s], embeddings[s]
        records[s] = []
        for bit, r in enumerate(d['records']):
            p, m = r['path'], r['mask']
            global_p = [emb[i] for i in p]
            need(sorted(p) == list(range(8)) and all(arc(u, v) for u, v in zip(global_p, global_p[1:])),
                 'literal host-supported Hamiltonian path on its induced eight-owner cell')
            need(word(aa, p) == m == sum(arc(global_p[i], global_p[j]) << k for k, (i, j) in enumerate(CHORDS)),
                 'normalization unchanged by the full owner embedding')
            need(merged[owner[m]] == s and merged[owner[FULL ^ m]] == (U, V)[bit],
                 'actual source/partner class owners unchanged')
            need(all(p[t] == v and global_p[t] == emb[v] for t, v in d['fixed']),
                 'all inherited common-clock owner pins retained')
            rev = list(reversed(global_p))
            need(all(arc(v, u) for u, v in zip(rev, rev[1:])), 'actual converse path in opposite host')
            rm = sum(arc(rev[j], rev[i]) << k for k, (i, j) in enumerate(CHORDS))
            need(rm == reflect(m), 'literal converse word under reversal')
            kk = partner(aa, p)
            need(word(kk, p) == FULL ^ m, 'literal local K-partner word')
            words = orbit(m)
            need(len(words) == 4 and not reservoir & words, 'all full inverse/converse resources are disjoint')
            reservoir |= words
            records[s].append(dict(bit=bit, target=(U, V)[bit], local_path=p, global_path=global_p,
                source_word=m, partner_word=FULL ^ m, normalized_resource=min(words),
                closed_words=sorted(words), global_fixed_clocks=[[t, emb[v]] for t, v in d['fixed']]))
        costs[s] = sum(u != v for u, v in zip(d['records'][0]['path'], d['records'][1]['path']))
    need(len(reservoir) == 608, '76 slots reserve152 normalized alternatives or608 closed words')

    # Full canonical-source lift and complete pin-domain checks for the two shared cells.
    lifts = {}
    for s in (176, 368):
        aa, emb, d = matrices[s], embeddings[s], data[s]
        all_paths, pinned, tile_lifts = [], [], [[], []]
        for p in permutations(range(8)):
            if not all(aa[u][v] for u, v in zip(p, p[1:])):
                continue
            m = word(aa, p)
            pp = [emb[v] for v in p]
            need(m == sum(arc(pp[i], pp[j]) << k for k, (i, j) in enumerate(CHORDS)),
                 'all full source lifts transport without quotient loss')
            all_paths.append(p)
            if all(p[t] == v for t, v in d['fixed']) and merged[owner[FULL ^ m]] in (U, V):
                pinned.append((list(p), m, merged[owner[FULL ^ m]]))
            for bit, r in enumerate(d['records']):
                if m in orbit(r['mask']):
                    tile_lifts[bit].append(dict(local_path=list(p), global_path=pp, source_word=m))
        need({(tuple(p), m, t) for p, m, t in pinned} ==
             {(tuple(r['path']), r['mask'], r['target']) for r in d['records']},
             'complete inherited pin-plus-U/V domain remains exactly two')
        lifts[s] = dict(H=len(all_paths), full_canonical_source_tile_lifts=tile_lifts)
    need(lifts[176]['H'] == 315 and lifts[368]['H'] == 311, 'independent two-cell full path census')

    pair_rows = []
    for cb, hb in product(range(2), repeat=2):
        cp, hp = records[176][cb]['global_path'], records[368][hb]['global_path']
        allowed = cp[2] == hp[1]
        need(allowed == (cb+hb == 1), 'shared indexed-owner equality is EXACT opposite-bit predicate')
        need((cp[::-1][5] == hp[::-1][6]) == allowed,
             'literal converse keeps the equality at the reflected clocks5 and6')
        # Shared edge is nonconsecutive in every option, so both local K's reverse it.
        ak = partner(a, data[176]['records'][cb]['path'])
        bk = partner(b, data[368]['records'][hb]['path'])
        hk = host([(ak, embeddings[176]), (bk, embeddings[368])], 14)
        need(hk[2][6] == 1 and all(hk[embeddings[s][i]][embeddings[s][j]] == kk[i][j]
             for s, kk in ((176, ak), (368, bk)) for i, j in PAIRS),
             'local partner matrices also admit a coherent host restriction')
        pair_rows.append(dict(core_bit=cb, helper_bit=hb, core_clock2=cp[2], helper_clock1=hp[1],
                              allowed=allowed, partner_shared_arc_2_to_6=hk[2][6]))

    # Independent complete quadruple enumeration, specifically for the genuine owner guard.
    section, multiplicity = [{}, {}], [Counter(), Counter()]
    for chosen in combinations(helpers, 4):
        hb = int(368 in chosen)
        q = charge(chosen)
        bm = sum(1 << helpers.index(s) for s in chosen)
        section[hb].setdefault(q, bm)
        multiplicity[hb][q] += 1
    need(sum(sum(h.values()) for h in multiplicity) == 194580, 'all48 choose4 subsets independently enumerated')
    for hb in range(2):
        need(set(section[hb]) == set(range(2048)), 'each actual shared-owner pin value has EVERY residual charge')
        for q, bm in section[hb].items():
            chosen = {helpers[i] for i in range(48) if bm >> i & 1}
            need(len(chosen) == 4 and int(368 in chosen) == hb and charge(chosen) == q,
                 'each concrete signed-section witness retains occupancy, owner pin and fullQ')

    def helper_state(core, total_q):
        hb = 1-int(176 in core)
        bm = section[hb][minus(total_q, charge(core))]
        return {helpers[i] for i in range(48) if bm >> i & 1}

    # A genuine old four-source charge trade, detoured through core176 so the
    # shared owner changes and the SAME helper records must return at the end.
    initial, final = {1394, 2185}, {1477, 1841}
    need(charge(initial) == charge(final) and 176 not in initial | final,
         'actual four-source core trade has equalQ and the shared core bit returns')
    q = 0
    core_route = [initial, {176, 2185}, {176, 1841}, final]
    route = [(set(core), helper_state(core, q)) for core in core_route]
    need(route[0][1] == route[-1][1], 'full helper endpoint records return exactly')

    def rows(core, active_helpers):
        active = core | active_helpers
        out = {s: [0]*6 for s in list(data)+[U, V]}
        for s in data:
            t = V if s in active else U
            for j, p in enumerate(potentials):
                v = 2*(p[s]-p[t])
                out[s][j] += v
                out[t][j] -= v
        return out

    states, events = [], []
    baseline = rows(*route[0])
    for c, h in route:
        need(len(c) == 2 and len(h) == 4 and plus(charge(c), charge(h)) == q, 'route margins and total charge')
        need(int(368 in h) == 1-int(176 in c), 'route literal shared-owner equality')
        currents = rows(c, h)
        need(all((currents[s][j]-baseline[s][j]) % moduli[j] == 0 for s in currents for j in range(6)),
             'all physical paired dyadic generator rows preserved')
        words = set()
        for s in data:
            record = records[s][int(s in c or s in h)]
            need(not words & set(record['closed_words']), 'actual arrival resource injectivity')
            words.update(record['closed_words'])
        need(len(words) == 304, 'every arrival occupies76 four-word resources')
        states.append(dict(core_V=sorted(c), helper_V=sorted(h), total_charge=q,
            shared_owner=records[176][int(176 in c)]['global_path'][2],
            changed_physical_rows={str(s): [currents[s][j]-baseline[s][j] for j in range(6)]
                                  for s in currents if currents[s] != baseline[s]}))
    for (c, h), (cc, hh) in zip(route, route[1:]):
        changed = sorted((c ^ cc) | (h ^ hh))
        cost = sum(costs[s] for s in changed)
        need(len(changed) <= 10 and cost <= 64, 'native coupled event support and actual clock bounds')
        events.append(dict(changed_sources=changed, support=len(changed), owner_clock_occurrences=cost))
    need(route[-1][0] == final and [s['shared_owner'] for s in states] == [2, 6, 6, 2],
         'moving shared owner returns and final protected core records realized')
    bound = sum(sorted((costs[s] for s in cores), reverse=True)[:2]) + \
            sum(sorted((costs[s] for s in helpers), reverse=True)[:8])
    need(bound == 64, 'uniform64 clock-occurrence bound from actual menu costs')

    out = dict(status='FINITE-EXACT native realization of state-dependent opposite helper bit',
        core_source=176, helper_source=368, core_clock=2, helper_clock=1,
        shared_owners=[2, 6], converse_equality_clocks=[5, 6], two_cell_host_adjacency=local_host,
        embeddings={str(s): embeddings[s] for s in sorted(data)},
        two_cell_records={str(s): records[s] for s in (176, 368)},
        all_inherited_pins={str(s): data[s]['fixed'] for s in sorted(data)},
        full_two_cell_kernel_lifts={str(s): lifts[s] for s in (176, 368)},
        all_bit_pairs=pair_rows, helper_source_ids=helpers, helper_V_occupancy=4,
        signed_sections=[dict(forced_helper_bit=b, minimum_witnesses=min(multiplicity[b].values()),
            total_witnesses=sum(multiplicity[b].values()), section_masks=[section[b][q] for q in range(2048)]) for b in range(2)],
        native_route=dict(states=states, events=events),
        slot_count=76, complete_source_host_order=606, reserved_closed_words=608, arrival_closed_words=304,
        source_hashes={p.relative_to(RESULTS).as_posix(): digest(p) for p in
            (supply_path, old_path, res_path, nodes_path, owner_path, electrical_path)},
        support_bound=10, owner_clock_occurrence_bound=64, checks=COUNT,
        scope='Order-eight native resources with explicit shared ownership; induced-cell paths, not a global host observer. All original dyadic aggregate currents, not every ordinary per-port phase.')
    target = HERE/'shared_owner.json'
    if args.write:
        target.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n', encoding='utf-8', newline='\n')
    else:
        need(json.loads(target.read_text()) == out, 'frozen complete shared-owner certificate matches')
    print(f'PASS {COUNT:,}; '+('wrote shared_owner.json' if args.write else 'shared_owner.json matches'))


if __name__ == '__main__':
    main()
