"""Authenticated lagged cancellation at every finite native-depth cut.

The output is a relaxation by optional pulses, not a new synchronous state.
No unresolved suffix is computed or assumed to terminate. Discovery explores
at most16 native edges per supplied pending node; the consumer only checks
the selected paths and source-disjoint ownership.
"""
from hashlib import sha256
from importlib.util import spec_from_file_location, module_from_spec
from itertools import product
from pathlib import Path
import argparse
import json
import sys

HERE = Path(__file__).resolve().parent
SPEC = spec_from_file_location('pulse_geometry', HERE.parent/'tail_geometry_20261021c/check.py')
G = module_from_spec(SPEC); sys.modules[SPEC.name] = G; SPEC.loader.exec_module(G)
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def encode(x):
    return (json.dumps(x, sort_keys=True, separators=(',', ':'))+'\n').encode()


def edge(s):
    need(type(s) is int and s > 0 and s % 2 == 1 and s % 8 != 5, 'native primitive source')
    a = 2 if s % 8 == 1 else 1
    need((3*s+1) % (2**a) == 0, 'native first halving')
    u = (3*s+1)//(2**a)
    while u % 8 == 5:
        u = (u-1)//4; a += 2
    need(u > 0 and u % 2 == 1 and u % 8 != 5, 'native primitive target')
    return u, a


def independent_edge(s):
    # Inverse-form check of the removed 4u+1 tower. This path does not use edge.
    first = 4 if s % 8 == 1 else 2
    u = (3*s+1)//first; removed = 0
    while (u-5) % 8 == 0:
        old = u; u = (u-1)//4; removed += 1
        need(old == 4*u+1, 'literal tower inversion')
    return u, (2 if first == 4 else 1)+2*removed


def produce(state, frontier):
    bank = dict(state.pending); used = set(); selected = []
    for u, b in state.pending:
        x = u; clock = 0; path = []
        for length in range(1, 17):
            y, a = edge(x); path.append([x, y, a]); x = y; clock += a
            if x in frontier:
                break
            if x in bank and x != u and all(bank[x][i]+(-1)**clock*b[i] == 0 for i in range(4)):
                if x not in used and u not in used:
                    selected.append({'late_node': u, 'early_node': x, 'native_delay': length,
                                     'integer_clock': clock, 'path': path})
                    used.update((u, x))
                break
    return {'input_state_sha256': G.state_hash(state), 'pulses': selected}


def consume(state, certificate, frontier):
    G.audit_state(state)
    need(type(frontier) in (set, frozenset) and all(type(t) is int for t in frontier), 'literal stopped frontier set')
    need(frontier == {t for t, _ in state.ports}, 'frontier equals complete marked terminal universe')
    need(certificate['input_state_sha256'] == G.state_hash(state), 'original source-state hash')
    bank = dict(state.pending); used = set(); remove = set()
    for row in certificate['pulses']:
        u, v = row['late_node'], row['early_node']
        need(type(u) is int and type(v) is int and u != v and u in bank and v in bank, 'literal distinct pending owners')
        need(u not in used and v not in used, 'source-disjoint pulse ownership')
        used.update((u, v)); remove.add(u)
        need(type(row['native_delay']) is int and row['native_delay'] == len(row['path']) > 0, 'native delay retained')
        need(type(row['integer_clock']) is int and row['integer_clock'] > 0, 'literal positive integer clock')
        x = u; clock = 0
        for source, target, a in row['path']:
            need(all(type(value) is int and value > 0 for value in (source, target, a)), 'literal native path integers')
            need(source == x and source not in frontier and target not in frontier, 'first-hit safe prefix')
            need((target, a) == edge(source) == independent_edge(source), 'independent literal native transition')
            clock += a; x = target
        need(x == v and clock == row['integer_clock'], 'actual summed integer clock and endpoint')
        need(all(bank[v][i]+(-1)**clock*bank[u][i] == 0 for i in range(4)), 'separate old/new opposite charges')
    relaxed = G.State(state.ports, tuple((u, b) for u, b in state.pending if u not in remove))
    G.audit_state(relaxed)
    return relaxed


def controls():
    # Independent arbitrary-tail model. Infinite means neither arrives.
    for delay, clock, tau, sign in product(range(1, 7), range(1, 9), (None, 1, 2, 7), (-1, 1)):
        early = (2, -1, -3, 4); late = tuple(-(-1)**clock*x for x in early)
        for cut in range(16):
            original = tuple(sign*((early[i] if tau is not None and cut >= tau else 0)
                                  + ((-1)**clock*late[i] if tau is not None and cut >= tau+delay else 0))
                             for i in range(4))
            pulse = tuple(sign*x if tau is not None and tau <= cut < tau+delay else 0 for x in early)
            need(original == pulse, 'all four coordinates at every native cut, including nonarrival')
    # Arrival starts differ even when the eventual sum is zero.
    need(G.delta((0, 0, 1, 0)) == 1, 'a live pulse cannot be erased early')
    # u=1865 ->1399 has native delay1 but integer clock2.
    need(edge(1865) == (1399, 2), 'clock/delay hostile')
    state = G.State(((1, (0, 0, 0, 0)),), ((1399, (0, 0, -1, -2)), (1865, (0, 0, 1, 2))))
    good = {'input_state_sha256': G.state_hash(state), 'pulses': [
        {'late_node': 1865, 'early_node': 1399, 'native_delay': 1, 'integer_clock': 2, 'path': [[1865, 1399, 2]]}]}
    need(len(consume(state, good, {1}).pending) == 1, 'synthetic legal pulse')
    for bad in ({**good, 'pulses': good['pulses']*2},
                {**good, 'pulses': [{**good['pulses'][0], 'integer_clock': 1}]},
                {**good, 'pulses': [{**good['pulses'][0], 'integer_clock': 2.0}]},
                {**good, 'pulses': [{**good['pulses'][0], 'path': [[1865, 1399.0, 2.0]]}]},
                {**good, 'pulses': [{**good['pulses'][0], 'native_delay': True}]}):
        rejected = False
        try:
            consume(state, bad, {1})
        except ValueError:
            rejected = True
        need(rejected, 'ownership and clock corruption rejected')
    boundary_state = G.State(((11, (0, 0, 0, 0)),), ((7, (0, 0, 1, 0)), (17, (0, 0, -1, 0))))
    crossing = {'input_state_sha256': G.state_hash(boundary_state), 'pulses': [
        {'late_node': 7, 'early_node': 17, 'native_delay': 2, 'integer_clock': 2,
         'path': [[7, 11, 1], [11, 17, 1]]}]}
    for false_boundary in (set(), {11.0}, {11}):
        rejected = False
        try:
            consume(boundary_state, crossing, false_boundary)
        except ValueError:
            rejected = True
        need(rejected, 'boundary omission, alias or path crossing rejected')


def run(write):
    controls()
    raw = (HERE.parent/'future_tail_blocks_20261021b/receipt.json').read_bytes().replace(b'\r\n', b'\n')
    need(sha256(raw).hexdigest() == G.PREFIX_SHA, 'complete frozen source/guard/clock receipt')
    front = (HERE.parent/'session18_native_annulus/frontier.json').read_bytes().replace(b'\r\n', b'\n')
    need(sha256(front).hexdigest() == G.FRONTIER_SHA, 'fixed frontier')
    V = sorted(json.loads(front)['frontier']); frontier = set(V)
    certs = [] if write else json.loads((HERE/'paths.json').read_bytes())
    cases = []; exports = []
    for row in json.loads(raw)['cases']:
        ports = {a[0]: tuple(a[1:]) for a in row['arrived_ports']}
        state = G.State(tuple((t, ports.get(t, (0, 0, 0, 0))) for t in V),
                        tuple((a[0], tuple(a[1:5])) for a in row['pending']))
        # Recheck source guards/clocks inherited by every selected owner fibre.
        witnesses = {w['node']: w for w in row['witnesses']}
        for p in row['pending']:
            node, _, _, _, _, count, lo, hi, s, clock = p; w = witnesses[node]
            need(w['source'] == s and w['integer_clock'] == sum(w['edge_clocks']) == clock
                 and w['mod15'] == s % 15 and w['mod35'] == s % 35 and w['P_depth'] == row['D']
                 and not w['arrived'] and count > 0 and lo <= clock <= hi, 'immutable source registers and prior clocks')
        cert = produce(state, frontier) if write else certs[len(cases)]
        if write:
            certs.append(cert)
        relaxed = consume(state, cert, frontier)
        original = G.projected_certificate(state, 1, 1)['integer_bound']
        angular = G.projected_certificate(relaxed, 1, 1)['integer_bound']
        absolute = G.absolute_envelope(relaxed)
        need(angular <= original and angular <= absolute, 'positive pure-new controls improve')
        cases.append({'N': row['N'], 'D': row['D'], 'budget': 153*row['N'], 'pulses': len(cert['pulses']),
                      'original_blocks': len(state.pending), 'relaxed_blocks': len(relaxed.pending),
                      'input_state_sha256': G.state_hash(state), 'relaxation_sha256': G.state_hash(relaxed),
                      'original_angular': original, 'pulse_angular': angular, 'pulse_absolute': absolute,
                      'selected_native_edges': sum(p['native_delay'] for p in cert['pulses']),
                      'all_earlier_cuts_pass': row['early_cuts_pass'],
                      'source_rows_sha256': row['source_rows_sha256'], 'prefix_edges_sha256': row['prefix_edges_sha256']})
        exports.append({'N': row['N'], 'D': row['D'], 'kind': 'optional pulse relaxation, NOT synchronous native state',
                        'ports': relaxed.ports, 'pending': relaxed.pending,
                        'pulse_early_nodes': [p['early_node'] for p in cert['pulses']],
                        'relaxation_sha256': G.state_hash(relaxed)})
    result = {'status': 'PROVED all-future delayed-pulse replacement; FINITE-EXACT selected native paths; uniform153N and entry OPEN',
              'prefix_receipt_sha256_lf': G.PREFIX_SHA, 'frontier_sha256_lf': G.FRONTIER_SHA,
              'paths_sha256_lf': sha256(encode(certs)).hexdigest(), 'cases': cases}
    if write:
        (HERE/'paths.json').write_bytes(encode(certs)); (HERE/'relaxations.json').write_bytes(encode(exports))
        (HERE/'receipt.json').write_bytes(encode(result))
    else:
        need((HERE/'receipt.json').read_bytes().replace(b'\r\n', b'\n') == encode(result), 'frozen result replay')
        need((HERE/'relaxations.json').read_bytes().replace(b'\r\n', b'\n') == encode(exports), 'typed relaxation replay')
    print(json.dumps({'checks': CHECKS, 'sha256': sha256(encode(result)).hexdigest(), 'cases': cases}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true')
    run(parser.parse_args().write)
