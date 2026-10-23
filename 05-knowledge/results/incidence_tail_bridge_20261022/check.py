"""Search-free capacity audit of linked optional-arrival cycle certificates.

Discovery greedily chooses triangles, but no optimality is claimed. The saved
certificate names literal pending nodes, and verification checks every edge's
total capacity independently. Frozen native prefixes retain sources and clocks.
"""
from collections import defaultdict
from hashlib import sha256
from importlib.util import spec_from_file_location, module_from_spec
from itertools import combinations, product
from pathlib import Path
import argparse
import gzip
import json
import sys

HERE = Path(__file__).resolve().parent
SPEC = spec_from_file_location('linked_geometry', HERE.parent/'tail_geometry_20261021c/check.py')
G = module_from_spec(SPEC)
sys.modules[SPEC.name] = G
SPEC.loader.exec_module(G)
CHECKS = 0


def need(ok, msg):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(msg)


def encode(x):
    return (json.dumps(x, sort_keys=True, separators=(',', ':'))+'\n').encode()


def certify(state, excluded=()):
    """Integral triangle packing; excluded edges already belong to clusters."""
    G.audit_state(state)
    blocked = {tuple(sorted(e)) for e in excluded}
    bank = dict(state.pending)
    vertices = sorted(bank)
    raw = {(u, v): G.cross(bank[u], bank[v])
           for u, v in combinations(vertices, 2) if (u, v) not in blocked}
    capacity = {e: abs(k) for e, k in raw.items()}
    cycles = []
    for i, u in enumerate(vertices):
        neighbors = [v for v in vertices[i+1:] if capacity.get((u, v), 0)]
        for v, w in combinations(neighbors, 2):
            edges = ((u, v), (u, w), (v, w))
            if all(capacity.get(e, 0) for e in edges) and raw[edges[0]]*raw[edges[1]]*raw[edges[2]] < 0:
                weight = min(capacity[e] for e in edges)
                for e in edges:
                    capacity[e] -= weight
                cycles.append({'nodes': [u, v, w], 'weight': weight})
    # Triangle exhaustion need not exhaust negative cycles. Keep the same
    # residual edge ledger and look for longer fundamental negative cycles.
    adjacency = {u: set() for u in vertices}
    for (u, v), amount in capacity.items():
        if amount:
            adjacency[u].add(v); adjacency[v].add(u)
    for _ in range(1000):
        signs = {}; parent = {}; found = None
        for root in vertices:
            if root in signs:
                continue
            signs[root] = 1; parent[root] = None; queue = [root]
            for u in queue:
                for v in sorted(adjacency[u]):
                    e = tuple(sorted((u, v)))
                    wanted = signs[u]*(1 if raw[e] > 0 else -1)
                    if v not in signs:
                        signs[v] = wanted; parent[v] = u; queue.append(v)
                    elif signs[v] != wanted:
                        left = []; x = u
                        while x is not None:
                            left.append(x); x = parent[x]
                        right = []; x = v
                        while x not in left:
                            right.append(x); x = parent[x]
                        found = left[:left.index(x)+1]+list(reversed(right))
                        break
                if found:
                    break
            if found:
                break
        if not found:
            break
        edges = [tuple(sorted((u, v))) for u, v in zip(found, found[1:]+found[:1])]
        weight = min(capacity[e] for e in edges)
        for e in edges:
            capacity[e] -= weight
            if not capacity[e]:
                u, v = e; adjacency[u].remove(v); adjacency[v].remove(u)
        cycles.append({'nodes': found, 'weight': weight})
    return {'state_sha256': G.state_hash(state), 'excluded_edges': sorted(map(list, blocked)),
            'cycles': cycles, 'saving': 2*sum(c['weight'] for c in cycles)}


def audit(state, certificate):
    G.audit_state(state)
    need(certificate['state_sha256'] == G.state_hash(state), 'literal source-state commitment')
    bank = dict(state.pending)
    need(type(certificate['saving']) is int and certificate['saving'] >= 0, 'literal total saving')
    blocked = {tuple(e) for e in certificate['excluded_edges']}
    need(all(len(e) == 2 and e[0] < e[1] and all(type(u) is int and u in bank for u in e) for e in blocked), 'excluded owner edges')
    need(len(blocked) == len(certificate['excluded_edges']), 'duplicate excluded edge')
    spent = defaultdict(int)
    for row in certificate['cycles']:
        nodes, weight = row['nodes'], row['weight']
        need(type(weight) is int and weight > 0, 'positive literal capacity')
        need(len(nodes) >= 3 and len(set(nodes)) == len(nodes) and all(type(u) is int and u in bank for u in nodes), 'simple owned cycle')
        sign = 1
        for u, v in zip(nodes, nodes[1:]+nodes[:1]):
            e = tuple(sorted((u, v)))
            need(e not in blocked, 'cluster edge cannot be spent again')
            value = G.cross(bank[u], bank[v])
            need(value != 0, 'zero tie is not an interaction edge')
            sign *= 1 if value > 0 else -1
            spent[e] += weight
        need(sign == -1, 'negative cycle sign')
    for (u, v), total in spent.items():
        need(total <= abs(G.cross(bank[u], bank[v])), 'no edge-capacity recycling')
    saving = 2*sum(c['weight'] for c in certificate['cycles'])
    need(saving == certificate['saving'], 'exact capacity saving')
    return saving


def controls():
    # All signed weights of magnitudes1,2; all optional assignments to3ports.
    # Pair terms occur only at a shared nonzero destination, with actual sign.
    for ks in product((-2, -1, 1, 2), repeat=3):
        packing = 2*min(map(abs, ks)) if ks[0]*ks[1]*ks[2] < 0 else 0
        upper = sum(map(abs, ks))-packing
        best = -10**9
        for labels in product(range(-3, 4), repeat=3):
            value = 0
            for (i, j), k in zip(((0, 1), (1, 2), (0, 2)), ks):
                if labels[i] and abs(labels[i]) == abs(labels[j]):
                    value += k*(1 if labels[i]*labels[j] > 0 else -1)
            need(value <= upper, 'optional destination/sign triangle inequality')
            best = max(best, value)
        need(best == upper, 'triangle bound sharp on each scalar sign/magnitude control')
    vectors = ((-3, 2, 2, -3), (2, 0, -2, -2), (3, 0, -1, 3))
    state = G.State(tuple((i, (0, 0, 0, 0)) for i in range(3)), tuple((i+3, v) for i, v in enumerate(vectors)))
    cert = certify(state)
    need(audit(state, cert) == 36, 'mixed-channel cycle saving')
    upper = G.absolute_envelope(state)-cert['saving']
    need(upper == 21, 'mixed gate')
    choices = (None,)+tuple((t, s) for t in range(3) for s in (-1, 1))
    best = max(G.assignment(state, labels) for labels in product(choices, repeat=3))
    need(best == 18, 'independent full physical mixed optimum')
    # For one signed sum, 4Q(g_lambda)=37lambda²-2+73/lambda².
    # AMGM >84 gives every scalar Lorentz gate >21, including unsampled lambda.
    need(4*37*73 > 86**2, 'uniform all-positive-lambda separation by exact square comparison')
    bad = {**cert, 'cycles': cert['cycles']*2, 'saving': cert['saving']*2}
    rejected = False
    try:
        audit(state, bad)
    except ValueError:
        rejected = True
    need(rejected, 'double credit hostile rejected')
    for bad in ({**cert, 'saving': float(cert['saving'])},
                {**cert, 'cycles': [{**cert['cycles'][0], 'nodes': [3.0, 4, 5]}]},
                {**cert, 'excluded_edges': [[3.0, 4]]}):
        rejected = False
        try:
            audit(state, bad)
        except ValueError:
            rejected = True
        need(rejected, 'literal owner and saving aliases rejected')
    # A positive rank-two Gram example where triangle packing alone misses19.
    # Explicit primal/dual witnesses, no optimizer is used in this audit.
    planar = ((-4, -5), (-3, -8), (3, -7), (-7, -1), (-5, -5))
    kernel = {(i, j): G.beta(planar[i], planar[j]) for i, j in combinations(range(5), 2)}
    triangle_cover = {(0, 3), (1, 3), (2, 4)}
    for tri in combinations(range(5), 3):
        es = tuple(combinations(tri, 2))
        if kernel[es[0]]*kernel[es[1]]*kernel[es[2]] < 0:
            need(bool(set(es) & triangle_cover), 'triangle packing dual cover')
    need(sum(abs(kernel[e]) for e in triangle_cover) == 48, 'triangle upper dual48')
    packing = (((0, 2, 3), 27), ((1, 3, 4), 1), ((2, 3, 4), 20),
               ((0, 2, 3, 4), 6), ((1, 2, 3, 4), 13))
    spent = defaultdict(int)
    for cycle, weight in packing:
        sign = 1
        for u, v in zip(cycle, cycle[1:]+cycle[:1]):
            e = tuple(sorted((u, v))); spent[e] += weight
            sign *= 1 if kernel[e] > 0 else -1
        need(sign == -1, 'larger-cycle negative sign')
    need(all(amount <= abs(kernel[e]) for e, amount in spent.items()), 'larger-cycle capacity ledger')
    spin = (-1, -1, -1, 1, -1)
    loss = sum(abs(k) for (i, j), k in kernel.items() if spin[i]*spin[j]*k < 0)
    need(loss == sum(w for _, w in packing) == 67, 'all-cycle primal/dual equality67')
    return {'mixed_vectors': vectors, 'mixed_envelope': 57, 'cycle_bound': upper,
            'exact_optional_assignment_maximum': best, 'all_lambda_lower_expression': '(sqrt(2701)-1)/2',
            'rank_two_vectors': planar, 'triangle_packing_optimum': 48, 'larger_cycle_packing_optimum': 67}


def run(write=False, refresh=False):
    control = controls()
    raw = (HERE.parent/'future_tail_blocks_20261021b/receipt.json').read_bytes().replace(b'\r\n', b'\n')
    need(sha256(raw).hexdigest() == G.PREFIX_SHA, 'frozen source/guard/clock receipt')
    frontier = (HERE.parent/'session18_native_annulus/frontier.json').read_bytes().replace(b'\r\n', b'\n')
    need(sha256(frontier).hexdigest() == G.FRONTIER_SHA, 'fixed terminal universe')
    V = sorted(json.loads(frontier)['frontier'])
    cases = []
    saved_path = HERE/'packings.json.gz'
    saved = None if write else json.loads(gzip.decompress(saved_path.read_bytes()))
    packings = []
    for row in json.loads(raw)['cases']:
        arrived = {r[0]: tuple(r[1:]) for r in row['arrived_ports']}
        state = G.State(tuple((t, arrived.get(t, (0, 0, 0, 0))) for t in V),
                        tuple((r[0], tuple(r[1:5])) for r in row['pending']))
        U = G.absolute_envelope(state)
        need(U == row['envelope'][-1], 'inherited envelope reproduced')
        # Same frozen source registers and actual clocks; no new future suffix.
        witnesses = {w['node']: w for w in row['witnesses']}
        for p in row['pending']:
            node, _, _, _, _, count, lo, hi, source, clock = p
            w = witnesses[node]
            need(w['source'] == source and w['integer_clock'] == sum(w['edge_clocks']) == clock
                 and w['P_depth'] == row['D'] and not w['arrived']
                 and w['mod15'] == source % 15 and w['mod35'] == source % 35
                 and count > 0 and lo <= clock <= hi, 'retained literal source, guards and integer clocks')
        cert = certify(state) if write else saved[len(cases)]
        gain = audit(state, cert)
        angular = G.projected_certificate(state, 1, 1)['integer_bound']
        need(angular <= U-gain, 'native pure-new exact angular comparison')
        cases.append({'N': row['N'], 'D': row['D'], 'state_sha256': G.state_hash(state),
                      'U': U, 'cycle_saving': gain, 'cycle_bound': U-gain,
                      'angular_bound': angular, 'number_of_cycles': len(cert['cycles']),
                      'longer_cycles': sum(len(c['nodes']) > 3 for c in cert['cycles']),
                      'source_rows_sha256': row['source_rows_sha256'],
                      'prefix_edges_sha256': row['prefix_edges_sha256']})
        packings.append(cert)
    if write:
        saved_path.write_bytes(gzip.compress(encode(packings), mtime=0))
    result = {'status': 'PROVED scoped linked cycle-capacity bound; FINITE-EXACT frozen prefixes; uniform153N OPEN',
              'prefix_receipt_sha256_lf': G.PREFIX_SHA, 'controls': control, 'cases': cases,
              'packings_sha256_lf': sha256(encode(packings)).hexdigest(), 'checks': CHECKS}
    output = encode(result)
    if write or refresh:
        (HERE/'receipt.json').write_bytes(output)
    else:
        need((HERE/'receipt.json').read_bytes().replace(b'\r\n', b'\n') == output, 'exact frozen result replay')
    print(json.dumps({'checks': CHECKS, 'sha256': sha256(output).hexdigest(), 'cases': cases}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--refresh', action='store_true', help='refresh result metadata using existing checked packing')
    args = parser.parse_args()
    run(args.write, args.refresh)
