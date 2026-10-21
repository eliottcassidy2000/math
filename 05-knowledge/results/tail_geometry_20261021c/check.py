"""Linked pending-block geometry. Consumes one frozen finite-prefix receipt.

No native edge, future suffix, ROOT completion, or new source census is run.
All certificates bind the complete ordered current-port/pending-node state.
"""
from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key, lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import isqrt
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
PREFIX_SHA = '7731a9a9ddf8e5f0478b8d6455e53c6130870f2004353db2e0fefa345bb509c0'
FRONTIER_SHA = '391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f'
CHECKS = 0


def need(ok, msg):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(msg)


def integer(x, lower=0):
    if type(x) is not int or x < lower:
        raise ValueError('literal integer domain')
    return x


def vec(x, n=4):
    if type(x) is not tuple or len(x) != n or any(type(a) is not int for a in x):
        raise ValueError('literal integer vector')
    return x


def Q(v):
    return v[0]**2 - v[0]*v[1] + v[1]**2


def beta(v, w):
    return 2*v[0]*w[0] - v[0]*w[1] - v[1]*w[0] + 2*v[1]*w[1]


def delta(v):
    return Q(v[2:]) + beta(v[:2], v[2:])


def cross(v, w):
    return beta(v[2:], w[2:]) + beta(v[:2], w[2:]) + beta(w[:2], v[2:])


def add(v, w):
    return tuple(a+b for a, b in zip(v, w))


def scale(v, s):
    return tuple(s*a for a in v)


def enc(x):
    return (json.dumps(x, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()


@dataclass(frozen=True)
class State:
    # Rows are (literal terminal/current-node identifier, four-vector).
    ports: tuple
    pending: tuple


def audit_state(state):
    if type(state) is not State:
        raise ValueError('State required')
    ids = []
    for bank in (state.ports, state.pending):
        if type(bank) is not tuple:
            raise ValueError('tuple bank required')
        current = []
        for row in bank:
            if type(row) is not tuple or len(row) != 2:
                raise ValueError('marked row required')
            integer(row[0]); vec(row[1]); current.append(row[0])
        if len(set(current)) != len(current):
            raise ValueError('duplicate source-state identifier')
        ids.append(set(current))
    if not state.ports or ids[0] & ids[1]:
        raise ValueError('nonempty terminal universe, disjoint pending nodes')
    return state


def state_hash(state):
    audit_state(state)
    return sha256(enc([state.ports, state.pending])).hexdigest()


def physical(ports):
    return sum(Q(add(v[:2], v[2:]))-Q(v[:2]) for v in ports)


def assignment(state, choices):
    """Full physical increment for absent or (terminal ID, future sign) blocks."""
    audit_state(state)
    if type(choices) is not tuple or len(choices) != len(state.pending):
        raise ValueError('complete choice tuple')
    ports = dict(state.ports)
    for (_, v), choice in zip(state.pending, choices):
        if choice is None:
            continue
        if type(choice) is not tuple or len(choice) != 2:
            raise ValueError('marked future choice')
        t, sign = choice
        integer(t)
        if type(sign) is not int or sign not in (-1, 1) or t not in ports:
            raise ValueError('literal legal terminal/sign')
        ports[t] = add(ports[t], scale(v, sign))
    return physical(ports.values())


def support(state, v):
    return max(abs(cross(c, v)) for _, c in state.ports)


def absolute_envelope(state):
    audit_state(state)
    b = [v for _, v in state.pending]
    return (physical(v for _, v in state.ports)
            + sum(max(0, delta(v)+support(state, v)) for v in b)
            + sum(abs(cross(v, w)) for v, w in combinations(b, 2)))


def half(v):
    return 0 if v[1] > 0 or (v[1] == 0 and v[0] > 0) else 1


def angle_cmp(v, w):
    if half(v) != half(w):
        return half(v)-half(w)
    det = v[0]*w[1]-v[1]*w[0]
    return -1 if det > 0 else 1 if det < 0 else (v > w)-(v < w)


def zonotope(vectors):
    """Exact planar signed-sum norm maximum; arbitrary integer directions."""
    if type(vectors) is not tuple:
        raise ValueError('tuple generator bank')
    generators = []
    for v in vectors:
        vec(v, 2)
        if any(v):
            generators.append(v if half(v) == 0 else scale(v, -1))
    generators.sort(key=cmp_to_key(angle_cmp))
    z = tuple(-sum(v[i] for v in generators) for i in range(2))
    boundary = [z]
    for v in generators:
        z = add(z, scale(v, 2)); boundary.append(z)
    boundary += [scale(v, -1) for v in boundary]
    best = max(boundary, key=lambda v: (Q(v), v))
    return Q(best), best, len(boundary)


def projection(v, p, q):
    """Integer numerators G,F with common denominator 2pq; lambda=p/q."""
    vec(v); integer(p, 1); integer(q, 1)
    g = tuple(2*p*p*v[i]+(p*p+q*q)*v[i+2] for i in range(2))
    f = tuple(2*p*p*v[i]+(p*p-q*q)*v[i+2] for i in range(2))
    return g, f, 2*p*q


def projected_certificate(state, p=1, q=1):
    audit_state(state); integer(p, 1); integer(q, 1)
    pending = [v for _, v in state.pending]
    projected = tuple(projection(v, p, q)[0] for v in pending)
    maximum, vertex, candidates = zonotope(projected)
    base = physical(v for _, v in state.ports)
    linear = sum(support(state, v) for v in pending)
    denominator = (2*p*q)**2
    return {'state_sha256': state_hash(state), 'lambda': [p, q],
            'base': base, 'linear_support': linear,
            'pending_numerator_max': maximum, 'denominator_squared': denominator,
            'maximizing_boundary_point': list(vertex), 'boundary_candidates': candidates,
            'integer_bound': base+linear+maximum//denominator}


def audit_projected(state, certificate):
    if type(certificate) is not dict or type(certificate.get('lambda')) is not list or len(certificate['lambda']) != 2:
        raise ValueError('projected certificate shape')
    p, q = certificate['lambda']; integer(p, 1); integer(q, 1)
    expected = projected_certificate(state, p, q)
    if enc(certificate) != enc(expected):
        raise ValueError('source-bound projected certificate mismatch')
    return expected['integer_bound']


def perimeter_majorant(state, p=1, q=1):
    """Rational radii above exact norm sums; usable without radical arithmetic."""
    audit_state(state); integer(p, 1); integer(q, 1)

    def radius(v):
        n = Q(v); r = isqrt(n)
        return r+(r*r != n)

    block = [projection(v, p, q) for _, v in state.pending]
    ports = [projection(v, p, q) for _, v in state.ports]
    P = sum(radius(g) for g, _, _ in block)
    R = sum(radius(f) for _, f, _ in block)
    G = max(radius(g) for g, _, _ in ports)
    F = max(radius(f) for _, f, _ in ports)
    den = 2*p*q
    numerator = P*P+2*G*P+2*F*R
    return {'lambda': [p, q], 'common_radius_denominator': den,
            'radius_numerators': {'pending_g': P, 'pending_f': R, 'port_g': G, 'port_f': F},
            'tail_numerator': numerator, 'denominator_squared': den*den,
            'integer_bound': physical(v for _, v in state.ports)+numerator//(den*den)}


def pair_local(state, v, w):
    """Exact local maximum for two blocks, with a maximizing marked choice."""
    audit_state(state); vec(v); vec(w)
    dv, dw, k = delta(v), delta(w), cross(v, w)
    fv = [(abs(cross(c, v)), t, cross(c, v)) for t, c in state.ports]
    fw = [(abs(cross(c, w)), t, cross(c, w)) for t, c in state.ports]
    tv = sorted(fv, reverse=True)[:2]; tw = sorted(fw, reverse=True)[:2]
    best, witness = 0, [None, None]

    def offer(value, choices):
        nonlocal best, witness
        if value > best:
            best, witness = value, choices

    for i, (d, top) in enumerate(((dv, tv), (dw, tw))):
        m, t, a = top[0]; choices = [None, None]
        choices[i] = [t, 1 if a >= 0 else -1]
        offer(d+m, choices)
    # Equal destination: every sign choice, including negative k.
    for (t, c) in state.ports:
        a, b = cross(c, v), cross(c, w)
        for sv, sw in product((-1, 1), repeat=2):
            offer(dv+dw+sv*a+sw*b+sv*sw*k, [[t, sv], [t, sw]])
    # Different destinations: the best two supports in each list suffice.
    for (a, t, fa), (b, u, fb) in product(tv, tw):
        if t != u:
            offer(dv+dw+a+b, [[t, 1 if fa >= 0 else -1], [u, 1 if fb >= 0 else -1]])
    baseline = max(0, dv+tv[0][0])+max(0, dw+tw[0][0])+abs(k)
    return {'maximum': best, 'gain': baseline-best, 'witness': witness}


def pairing_certificate(state):
    """Greedy positive-gain disjoint pairs; no matching optimality claim."""
    audit_state(state)
    local = lru_cache(None)(lambda v, w: pair_local(state, v, w))
    candidates = []
    for i, (u, v) in enumerate(state.pending):
        for z, w in state.pending[:i]:
            record = local(v, w)
            if record['gain'] > 0:
                candidates.append((record['gain'], u, z, record))
    candidates.sort(key=lambda row: row[:3], reverse=True)
    used = set(); pairs = []
    for gain, u, z, record in candidates:
        if u not in used and z not in used:
            used.update((u, z)); pairs.append({'nodes': [u, z], **record})
    gain = sum(r['gain'] for r in pairs)
    return {'state_sha256': state_hash(state), 'pairs': pairs, 'gain': gain,
            'integer_bound': absolute_envelope(state)-gain}


def audit_pairing(state, certificate):
    audit_state(state)
    if type(certificate) is not dict or type(certificate.get('pairs')) is not list:
        raise ValueError('pair certificate shape')
    bank = dict(state.pending); used = set(); gain = 0
    for row in certificate['pairs']:
        if type(row) is not dict or type(row.get('nodes')) is not list or len(row['nodes']) != 2:
            raise ValueError('two source nodes required')
        u, z = row['nodes']; integer(u); integer(z)
        if u == z or u in used or z in used or u not in bank or z not in bank:
            raise ValueError('disjoint actual source-node owners')
        expected = {'nodes': [u, z], **pair_local(state, bank[u], bank[z])}
        if enc(row) != enc(expected):
            raise ValueError('local exact maximum mismatch')
        used.update((u, z)); gain += row['gain']
    expected = {'state_sha256': state_hash(state), 'pairs': certificate['pairs'],
                'gain': gain, 'integer_bound': absolute_envelope(state)-gain}
    if enc(certificate) != enc(expected):
        raise ValueError('full partition bound mismatch')
    return expected['integer_bound']


def rejects(fn):
    try:
        fn()
    except (ValueError, TypeError, KeyError):
        return True
    return False


def controls():
    # Independent full-sign enumeration; no angular ordering in this consumer.
    alphabet = tuple((a, b) for a in range(-1, 2) for b in range(-1, 2))
    for n in range(5):
        for vs in combinations_with_replacement(alphabet, n):
            brute = max(Q(tuple(sum(s*v[i] for s, v in zip(signs, vs)) for i in range(2)))
                        for signs in product((-1, 1), repeat=n))
            need(zonotope(vs)[0] == brute, 'planar zonotope versus complete signs')
    ports = ((0, (1, 0, 0, 1)), (1, (0, 1, -1, 0)))
    allv = tuple(product(range(-1, 2), repeat=4))
    choices = (None, (0, 1), (0, -1), (1, 1), (1, -1))
    base = physical(v for _, v in ports)
    for v, w in product(allv, repeat=2):
        state = State(ports, ((2, v), (3, w)))
        values = [assignment(state, cs) for cs in product(choices, repeat=2)]
        local = pair_local(state, v, w)
        need(base+local['maximum'] == max(values), 'exact optional terminal/sign pair maximum')
        need(local['gain'] >= 0, 'partition refinement never worsens U')
        for p, q in ((1, 1), (1, 2), (2, 1)):
            cert = projected_certificate(state, p, q)
            need(max(values) <= cert['integer_bound'], 'mixed linked projection dominates all assignments')
            g, f, d = projection(v, p, q)
            need(Q(g)-Q(f) == d*d*delta(v), 'exact rational hyperbolic identity')
            gg, ff, _ = projection(w, p, q)
            need(beta(g, gg)-beta(f, ff) == d*d*cross(v, w), 'linked polarization identity')
    v = (-5, 0, 8, 0)
    initial = State(((0, (0, 0, 0, 0)),), ((1, v), (2, v)))
    later = State(((0, v),), ((2, v),))
    need(delta(v) == -16 and max(delta(v), delta(scale(v, -1))) < 0,
         'indefinite vertex-only maximum misses optional nonarrival')
    need(projected_certificate(initial, 2, 1)['integer_bound'] == 0,
         'linked negative cone is annihilated by a rational projection')
    need(projected_certificate(initial)['integer_bound'] == 36 and absolute_envelope(initial) == 32,
         'generic projection has no universal dominance over U')
    need(projected_certificate(later, 2, 1)['integer_bound'] == 16,
         'sound projection need not be horizon monotone; retain earlier certificates')
    need(perimeter_majorant(initial, 2, 1)['integer_bound'] == 0,
         'computable coalescence majorant also retains linked negative cone')
    one = State(((0, (0, 0, 0, 0)),), ((1, (0, 0, 1, 0)), (2, (0, 0, -1, 0))))
    need(assignment(one, ((0, 1), None)) == 1, 'different pending nodes cannot be erased by their sum')
    small = State(((0, (0, 0, 1, 0)), (1, (0, 0, 0, 1))),
                  ((2, (-1, 0, -1, -1)), (3, (-1, 1, 1, -1))))
    need(pair_local(small, small.pending[0][1], small.pending[1][1])['gain'] == 3,
         'strict terminal/sign/absence conflict gain')
    cert = projected_certificate(small); audit_projected(small, cert)
    pair = pairing_certificate(small); audit_pairing(small, pair)
    need(rejects(lambda: audit_projected(small, {**cert, 'integer_bound': float(cert['integer_bound'])})), 'float result alias rejected')
    need(rejects(lambda: projected_certificate(small, True, 1)), 'Boolean boost rejected')
    need(rejects(lambda: projection((0, False, 1, 0), 1, 1)), 'Boolean coordinate rejected')
    need(rejects(lambda: audit_state(State(small.ports, ((True, (0, 0, 0, 0)),)))), 'Boolean source node rejected')
    need(rejects(lambda: audit_pairing(small, {**pair, 'pairs': pair['pairs']*2})), 'pair-owner reuse rejected')
    changed = State(small.ports, ((4, small.pending[0][1]), small.pending[1]))
    need(rejects(lambda: audit_projected(changed, cert)), 'unchanged vectors with different source marks rejected')
    need(rejects(lambda: assignment(small, ((0, True), None))), 'Boolean future sign rejected')


def load_native():
    path = HERE.parent/'future_tail_blocks_20261021b/receipt.json'
    raw = path.read_bytes().replace(b'\r\n', b'\n')
    need(sha256(raw).hexdigest() == PREFIX_SHA, 'frozen independently audited finite-prefix receipt')
    front = (HERE.parent/'session18_native_annulus/frontier.json').read_bytes().replace(b'\r\n', b'\n')
    need(sha256(front).hexdigest() == FRONTIER_SHA, 'same stopped frontier')
    V = tuple(json.loads(front)['frontier'])
    for t in V: integer(t, 1)
    need(len(set(V)) == len(V) == 627, 'complete terminal universe')
    cases = []
    for row in json.loads(raw)['cases']:
        N, D = row['N'], row['D']; integer(N, 1); integer(D)
        arrived = {r[0]: tuple(r[1:]) for r in row['arrived_ports']}
        need(set(arrived) <= set(V), 'literal arrived ports lie in fixed frontier')
        state = State(tuple((t, arrived.get(t, (0, 0, 0, 0))) for t in sorted(V)),
                      tuple((r[0], tuple(r[1:5])) for r in row['pending']))
        audit_state(state)
        need(absolute_envelope(state) == row['envelope'][-1], 'original U independently recovered')
        witnesses = {w['node']: w for w in row['witnesses']}
        for r in row['pending']:
            node, _, _, _, _, count, lo, hi, s, clock = r
            for x in (node, count, lo, hi, s, clock): integer(x)
            w = witnesses[node]
            need(count > 0 and lo <= clock <= hi and s < 256*N, 'source count, clock range, original carrier')
            need(w['source'] == s and w['integer_clock'] == clock and w['P_depth'] == D
                 and not w['arrived'] and sum(w['edge_clocks']) == clock,
                 'same literal source and full prefix clock in pinned handoff')
            need(w['mod15'] == s % 15 and w['mod35'] == s % 35 and w['mod15'] % 5 == w['mod35'] % 5,
                 'joint source registers preserved')
        certs = [projected_certificate(state, p, q) for p, q in ((1, 1), (1, 2), (2, 1))]
        for cert in certs: audit_projected(state, cert)
        pair = pairing_certificate(state); audit_pairing(state, pair)
        perimeter = perimeter_majorant(state)
        need(perimeter['integer_bound'] >= certs[0]['integer_bound'], 'rational perimeter relaxation dominates angular gate')
        best = min([row['envelope'][-1], pair['integer_bound']]+[c['integer_bound'] for c in certs])
        need(all(not any(v[:2]) for _, v in state.pending), 'finite controls happen to have zero pending old channel')
        need(certs[0]['integer_bound'] <= row['envelope'][-1], 'pure-new coalesced angular bound improves U')
        need((best <= 153*N) == (N >= 4096), 'positive controls and N7 hostile remain separated')
        cases.append({'N': N, 'D': D, 'budget': 153*N, 'pending_blocks': len(state.pending),
                      'source_rows_sha256': row['source_rows_sha256'], 'prefix_edges_sha256': row['prefix_edges_sha256'],
                      'state_sha256': state_hash(state), 'U': row['envelope'][-1],
                      'projected': certs, 'pairing': pair, 'perimeter_majorant': perimeter, 'selected_bound': best,
                      'selected_gain': row['envelope'][-1]-best,
                      'all_earlier_cuts_pass': row['early_cuts_pass']})
    return cases


def run():
    controls()
    cases = load_native()
    return {'status': 'PROVED linked mixed-channel majorant and cluster bound; FINITE-EXACT frozen-prefix consumers; uniform stopping OPEN',
            'prefix_receipt_sha256_lf': PREFIX_SHA, 'frontier_sha256_lf': FRONTIER_SHA,
            'cases': cases, 'checks': CHECKS}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true'); args = parser.parse_args()
    result = run(); raw = enc(result); path = HERE/'receipt.json'
    if args.write:
        path.write_bytes(raw)
    else:
        need(path.read_bytes().replace(b'\r\n', b'\n') == raw, 'saved exact receipt comparison')
    print(json.dumps({'checks': CHECKS, 'sha256_lf': sha256(raw).hexdigest(),
                      'cases': [[r['N'], r['D'], r['U'], r['selected_bound'], r['selected_gain']] for r in result['cases']]}))
