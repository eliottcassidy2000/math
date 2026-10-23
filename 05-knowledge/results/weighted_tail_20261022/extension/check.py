"""Strict consumer for bounded native extensions of frozen source fibres.

Checks every supplied edge, first-arrival boundary, all intermediate cuts,
source ownership, full inherited clocks and exact final weighted gates.
Does not query an unresolved suffix or infer convergence from a bound.
"""
from collections import defaultdict
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import argparse
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WEIGHTED_SHA = '82f051a4112097ceef146e8a59e750b4b787819dd95aca505750c2724189d242'
PREFIX_SHA = '7731a9a9ddf8e5f0478b8d6455e53c6130870f2004353db2e0fefa345bb509c0'
FRONTIER_SHA = '391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f'
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def pin(path, expected):
    raw = path.read_bytes().replace(b'\r\n', b'\n')
    need(sha256(raw).hexdigest() == expected, 'frozen dependency: '+str(path))
    return raw


def integer(n, lo=0):
    need(type(n) is int and n >= lo, 'literal non-Boolean integer')
    return n


def source(n):
    return type(n) is int and n > 0 and n % 2 == 1 and n % 8 != 5


def edge(s):
    """Independent valuation implementation, not the producer's branch rule."""
    need(source(s), 'primitive positive source')
    u = 3*s+1
    a = 0
    while u % 2 == 0:
        u //= 2
        a += 1
    while u % 8 == 5:
        need((u-1) % 4 == 0, 'exact tower strip')
        u = (u-1)//4
        a += 2
    need(source(u), 'primitive positive target')
    return u, a


pin(HERE.parent/'check.py', WEIGHTED_SHA)
spec = importlib.util.spec_from_file_location('native_extension_weighted', HERE.parent/'check.py')
W = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = W
spec.loader.exec_module(W)
G = W.G


def frontier():
    obj = json.loads(pin(ROOT/'session18_native_annulus/frontier.json', FRONTIER_SHA))
    V = obj['frontier']; ranks = obj['frontier_ranks']
    need(type(V) is list and len(V) == len(set(V)) == 627, 'exact 627-point frontier')
    need(type(ranks) is list and len(ranks) == 627, 'rank bank cardinality')
    rr = {}
    for pair in ranks:
        need(type(pair) is list and len(pair) == 2, 'rank pair')
        s, r = pair
        integer(s, 1); integer(r)
        need(s not in rr, 'unique rank address')
        rr[s] = r
    need(set(rr) == set(V) and rr[1] == 0, 'grounding root and rank domain')
    for s in V:
        need(source(s), 'frontier primitive source')
        if s != 1:
            target, _ = edge(s)
            need(target in rr and rr[target] < rr[s], 'strict grounded frontier rank')
    need(edge(1) == (1, 2), 'native root loop clock retained')
    return set(V), rr


def charge(N, s):
    """Independent two-coordinate charge formula inherited from native packet."""
    x = y = 0
    while s < 4*N:
        row = (0, -3, 3) if s < 2*N else (-1, -1, 2)
        x += row[s % 3]; y += row[(s+1) % 3]
        s = 4*s+1
    need((x-y) % 3 == 0 and (x+2*y) % 3 == 0, 'charge integrality')
    return (x-y)//3, -(x+2*y)//3


def check_witness(N, D, row, witness, V):
    node, *payload = row
    b = tuple(payload[:4])
    count, lo, hi, s, clock = payload[4:]
    for x in (node, count, s):
        integer(x, 1)
    for x in (lo, hi, clock):
        integer(x)
    need(len(b) == 4 and all(type(x) is int for x in b) and any(b), 'nonzero linked original block')
    need(source(s) and s < 256*N and lo <= clock <= hi, 'source carrier and clock range')
    need(type(witness) is dict and witness['source'] == s and witness['node'] == node, 'source witness binding')
    u = s; full_word = []
    for _ in range(D):
        need(u not in V, 'pending witness never previously arrived')
        u, a = edge(u); full_word.append(a)
    need(u == node and u not in V and sum(full_word) == clock, 'independent full prefix clock')
    expected = {'source': s, 'mod15': s % 15, 'mod35': s % 35, 'node': node,
                'P_depth': D, 'integer_clock': clock, 'edge_clocks': full_word,
                'arrived': False, 'old': s < 4*N, 'charge': list(charge(64*N, s))}
    need(G.enc(witness) == G.enc(expected), 'all witness coordinates retained')
    need(s % 15 % 5 == s % 35 % 5, 'compatible mod15/mod35 guards')


def validate_paths(old, proof, depth, V):
    need(type(proof) is dict and set(proof) == {'N','D','owners'}, 'literal case schema')
    integer(proof['N'], 1); integer(proof['D'])
    need(proof['N'] == old['N'] and proof['D'] == old['D'], 'case source binding')
    owners = proof['owners']
    need(type(owners) is list and len(owners) == len(old['pending']), 'complete original owner bank')
    bank = {r[0]: r for r in old['pending']}
    seen = set(); answers = {}
    for item in owners:
        need(type(item) is dict and set(item) == {'owner','path'}, 'owner/path schema')
        u = item['owner']; integer(u, 1)
        need(u in bank and u not in seen, 'unique original owner; no credit recycling')
        seen.add(u)
        path = item['path']
        need(type(path) is list and 0 < len(path) <= depth, 'bounded nonempty original pending path')
        points = [(u, 0)]
        x = u; clock = 0
        for triple in path:
            need(type(triple) is list and len(triple) == 3, 'literal edge triple')
            a, b, c = triple
            for z in triple: integer(z, 1)
            need(a == x and x not in V, 'path continuity and first-arrival stopping')
            need((b,c) == edge(a), 'literal native edge with full integer clock')
            x = b; clock += c; points.append((x, clock))
        need(len(path) == depth or x in V, 'unresolved path cannot stop early')
        answers[u] = points
    need(seen == set(bank), 'every source owner appears exactly once')
    return answers


def owner_state(old, points, cut, V):
    ports = {t:(0,0,0,0) for t in V}
    ports.update({r[0]:tuple(r[1:5]) for r in old['arrived_ports']})
    pending = defaultdict(lambda: (0,0,0,0))
    ownership = defaultdict(list)
    grounded = []
    final_rows = []
    for row in old['pending']:
        owner, *p = row
        path = points[owner]
        idx = min(cut, len(path)-1)
        node, clock = path[idx]
        b = G.scale(tuple(p[:4]), 1-2*(clock % 2))
        count, lo, hi, witness, witness_clock = p[4:]
        final_rows.append([owner,node,clock,count,lo+clock,hi+clock,witness,witness_clock+clock])
        if node in V:
            ports[node] = G.add(ports[node], b)
            grounded.append(owner)
        else:
            pending[node] = G.add(pending[node], b)
            ownership[node].append(owner)
    # Zero fibres remain in the lossless owner bank; only their energy shadow
    # is erased after all four signed coordinates have vanished at one node.
    state = G.State(tuple(sorted(ports.items())), tuple((u,b) for u,b in sorted(pending.items()) if any(b)))
    G.audit_state(state)
    return state, grounded, dict(ownership), final_rows


def advance(state, V):
    """Second path: synchronous aggregate transport, independent of owner replay."""
    ports = dict(state.ports); pending = defaultdict(lambda:(0,0,0,0))
    for u, b in state.pending:
        v, clock = edge(u)
        value = G.scale(b, 1-2*(clock % 2))
        if v in V:
            ports[v] = G.add(ports[v], value)
        else:
            pending[v] = G.add(pending[v], value)
    return G.State(tuple(sorted(ports.items())), tuple((u,b) for u,b in sorted(pending.items()) if any(b)))


def run_case(old, proof, depth, V):
    points = validate_paths(old, proof, depth, V)
    witnesses = {w['node']:w for w in old['witnesses']}
    need(set(witnesses) == {r[0] for r in old['pending']}, 'complete inherited witness bank')
    for row in old['pending']:
        check_witness(old['N'], old['D'], row, witnesses[row[0]], V)
    history = []; previous = None
    for cut in range(depth+1):
        state, grounded, ownership, owner_rows = owner_state(old, points, cut, V)
        if previous is not None:
            need(G.state_hash(advance(previous, V)) == G.state_hash(state), 'all-cut independent synchronous replay')
        energy = G.physical(v for _,v in state.ports)
        source_counts = {r[0]: r[5] for r in old['pending']}
        owners_pending = set(source_counts)-set(grounded)
        need(sum(source_counts[u] for u in grounded)+sum(source_counts[u] for u in owners_pending) == sum(source_counts.values()),
             'literal source-count partition at every cut')
        need(set(u for owners in ownership.values() for u in owners) == owners_pending, 'no lost pending owner')
        history.append({'D':old['D']+cut, 'increment':energy, 'grounded_owners':len(grounded),
                        'grounded_original_sources':sum(source_counts[u] for u in grounded),
                        'remaining_owners':len(owners_pending), 'remaining_owner_nodes':len(ownership),
                        'nonzero_pending_nodes':len(state.pending), 'state_sha256':G.state_hash(state)})
        previous = state
    need(history[0]['increment'] == old['history'][-1][1], 'original cut exact binding')
    cert = W.optimize(state)
    W.audit_certificate(state, cert)
    angular = G.projected_certificate(state)['integer_bound']
    U = G.absolute_envelope(state)
    earlier = old['early_cuts_pass'] and all(r['increment'] <= 153*old['N'] for r in history)
    future = cert['integer_bound'] <= 153*old['N']
    need((earlier and future) == (old['N'] >= 4096), 'positive declared scales and retained N7 hostile')
    return {'N':old['N'], 'D_before':old['D'], 'D_after':old['D']+depth, 'budget':153*old['N'],
            'source_rows_sha256':old['source_rows_sha256'], 'original_prefix_edges_sha256':old['prefix_edges_sha256'],
            'original_owner_rows_sha256':sha256(G.enc(old['pending'])).hexdigest(),
            'extension_paths_sha256':sha256(G.enc(proof)).hexdigest(), 'history':history,
            'grounded_owners':grounded, 'remaining_owner_fibres':[[u,ids] for u,ids in sorted(ownership.items())],
            'owner_clock_rows_schema':['original_node','current_node','extra_clock','source_count','min_full_clock','max_full_clock','witness','witness_full_clock'],
            'owner_clock_rows':owner_rows, 'arrived_ports':state.ports, 'pending':state.pending,
            'angular':angular, 'absolute_envelope':U, 'weighted':cert,
            'all_earlier_cuts_pass':earlier, 'all_future_cuts_pass':future}


def rejected(fn):
    try:
        fn()
    except (ValueError,TypeError,KeyError):
        need(True, 'hostile rejected')
    else:
        need(False, 'hostile accepted')


def controls(old, proof, V):
    variants = []
    x = deepcopy(proof); x['owners'][0]['path'][0][2] += 1; variants.append(x)
    x = deepcopy(proof); x['owners'][0]['owner'] = x['owners'][1]['owner']; variants.append(x)
    x = deepcopy(proof); x['owners'].pop(); variants.append(x)
    x = deepcopy(proof); x['owners'][0]['path'] = []; variants.append(x)
    x = deepcopy(proof); x['owners'][0]['path'][0][0] = float(x['owners'][0]['path'][0][0]); variants.append(x)
    x = deepcopy(proof); x['owners'][0]['path'][0][2] = True; variants.append(x)
    x = deepcopy(proof); x['N'] = float(x['N']); variants.append(x)
    x = deepcopy(proof); x['extra'] = 1; variants.append(x)
    for x in variants:
        rejected(lambda x=x:validate_paths(old,x,16,V))
    # A literal first-arrival path may not be extended beyond its boundary.
    x = deepcopy(proof)
    grounded = next(v for v in x['owners'] if v['path'][-1][1] in V)
    last = grounded['path'][-1][1]; target, a = edge(last)
    grounded['path'].append([last,target,a])
    rejected(lambda:validate_paths(old,x,16,V))
    need(source(1) and not source(True) and not source(5), 'primitive source type controls')


def run():
    original = json.loads(pin(ROOT/'future_tail_blocks_20261021b/receipt.json', PREFIX_SHA))
    V, ranks = frontier()
    raw = (HERE/'paths.json').read_bytes().replace(b'\r\n',b'\n')
    data = json.loads(raw)
    need(type(data) is dict and set(data) == {'kind','depth','cases'}, 'packet schema')
    need(data['kind'] == 'ordinary native extension with immutable original owners', 'synchronous state type')
    integer(data['depth'], 1)
    need(data['depth'] == 16 and type(data['cases']) is list and len(data['cases']) == len(original['cases']), 'declared complete bounded universe')
    controls(original['cases'][1],data['cases'][1],V)
    cases = [run_case(old,proof,16,V) for old,proof in zip(original['cases'],data['cases'])]
    return {'status':'FINITE-EXACT sixteen-edge native extension and grounded inherited children; uniform 153N OPEN',
            'weighted_code_sha256_lf':WEIGHTED_SHA, 'original_prefix_receipt_sha256_lf':PREFIX_SHA,
            'frontier_sha256_lf':FRONTIER_SHA, 'frontier_size':len(V), 'maximum_frontier_rank':max(ranks.values()),
            'paths_sha256_lf':sha256(raw).hexdigest(), 'cases':cases, 'checks':CHECKS}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true'); args = parser.parse_args()
    result = run(); raw = G.enc(result)
    if args.write:
        (HERE/'receipt.json').write_bytes(raw)
    else:
        need((HERE/'receipt.json').read_bytes().replace(b'\r\n',b'\n') == raw, 'exact frozen receipt comparison')
    print(json.dumps({'checks':CHECKS,'sha256_lf':sha256(raw).hexdigest(),
                      'cases':[[r['N'],r['D_after'],r['history'][-1]['grounded_owners'],r['history'][-1]['grounded_original_sources'],
                                r['history'][-1]['remaining_owners'],r['history'][-1]['nonzero_pending_nodes'],r['weighted']['integer_bound']] for r in result['cases']]}))
