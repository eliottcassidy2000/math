"""Exact affine certificate of order-ten endpoint-three black detours.

Standard library only. No tournament census and no imported research kernels.
Every check remains active under python -O. See README.md for exact scope.
"""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
N, D = 10, 12
PAIRS = list(combinations(range(N), 2))
CUT_PAIRS = [(1, 8), (2, 8), (2, 7)]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def parity(x):
    return x.bit_count() & 1


class Affine:
    """GF(2) elimination with exact XOR provenance, not a SAT black box."""
    def __init__(self, equations):
        self.equations = equations
        self.basis = {}
        self.contradiction = None
        for index, (mask, value) in enumerate(equations):
            proof = 1 << index
            while mask:
                pivot = mask.bit_length() - 1
                if pivot not in self.basis:
                    self.basis[pivot] = (mask, value, proof)
                    break
                m, v, p = self.basis[pivot]
                mask, value, proof = mask ^ m, value ^ v, proof ^ p
            if not mask and value:
                self.contradiction = proof
                break

    def proof_indices(self, proof):
        return [i for i in range(len(self.equations)) if proof >> i & 1]

    def proof_value(self, proof):
        mask = value = 0
        for i in self.proof_indices(proof):
            mask ^= self.equations[i][0]
            value ^= self.equations[i][1]
        return mask, value

    def certificate(self):
        if self.contradiction is not None:
            need(self.proof_value(self.contradiction) == (0, 1), 'false contradiction')
            return {'contradiction': self.proof_indices(self.contradiction)}
        need(len(self.basis) == D, 'unexpected positive-dimensional failure')
        units, assignment = [], 0
        for k in range(D):
            m, c, p = self.basis[k]
            for j in range(k):
                if m >> j & 1:
                    mm, cc, pp = units[j]
                    m, c, p = m ^ mm, c ^ cc, p ^ pp
            need(m == 1 << k and self.proof_value(p) == (m, c), 'false unit equation')
            units.append((m, c, p))
            assignment |= c << k
        need(all(parity(m & assignment) == c for m, c in self.equations), 'bad solution')
        return {'assignment': assignment, 'units': [self.proof_indices(p) for m, c, p in units]}


def source_forms(incoming):
    """T(u,v) = parity(mask & x) XOR constant, in twelve free middle bits."""
    out, seen, index, variables = {}, set(), 0, []
    for i, j in combinations(range(1, 9), 2):
        if (i, j) in seen:
            continue
        orbit = {(i, j), (9-j, 9-i)}
        seen.update(orbit)
        if j == i+1:
            form = (0, 1)
        else:
            form = (1 << index, 0)
            variables.append(sorted(orbit))
            index += 1
        for u, v in orbit:
            out[u, v], out[v, u] = form, (form[0], form[1] ^ 1)
    for i in range(1, 9):
        f = int(i not in incoming)
        out[0, i], out[i, 0] = (0, f), (0, f ^ 1)
        f = int(9-i in incoming)
        out[9, i], out[i, 9] = (0, f), (0, f ^ 1)
    out[9, 0], out[0, 9] = (0, 1), (0, 0)
    need(index == D and len(out) == N*(N-1), 'source dimension')
    for u, v in PAIRS:
        need(out[u, v] == out[9-v, 9-u], 'symbolic source anti-reversal')
        need(out[u, v][0] == out[v, u][0] and out[u, v][1] != out[v, u][1], 'tournament')
    need(all(out[i, i+1] == (0, 1) for i in range(9)), 'source word directed')
    need(sum(out[0, i][1] for i in range(1, 10)) == 6, 'first degree')
    need(sum(out[9, i][1] for i in range(9)) == 3, 'last degree')
    return out, variables


def equation(form, value):
    return form[0], form[1] ^ value


def path_edges(word):
    need(tuple(sorted(word)) == tuple(range(10)), 'word must be a permutation')
    return {tuple(sorted(p)) for p in zip(word, word[1:])}


def move_forms(forms, word):
    edges = path_edges(word)
    return {p: (m, c ^ int(tuple(sorted(p)) not in edges)) for p, (m, c) in forms.items()}


def anti_equations(forms, word):
    rev = dict(zip(word, word[::-1]))
    equations = []
    for u, v in PAIRS:
        m, c = forms[u, v]
        mm, cc = forms[rev[v], rev[u]]
        equations.append((m ^ mm, c ^ cc))
    return equations


def realize(forms, assignment):
    return {(u, v) for (u, v), (m, c) in forms.items() if parity(m & assignment) ^ c}


def direct_replay(initial, words):
    now, witnesses = set(initial), []
    for q in words:
        need(tuple(sorted(q)) == tuple(range(10)), 'direct permutation')
        retained = set(zip(q, q[1:]))
        need(retained <= now, 'direct directedness')
        r = dict(zip(q, q[::-1]))
        failures = [(u, v, r[v], r[u]) for u, v in sorted(now) if (r[v], r[u]) not in now]
        need(failures, 'direct actual blackness')
        witnesses.append(failures[0])
        now = {(v, u) for u, v in now - retained} | retained
    p = set(zip(range(9), range(1, 10)))
    target = {(v, u) for u, v in initial - p} | p
    need(now == target, 'direct literal target')
    return witnesses


def ledger(n, incoming, arcs):
    a,b,x,y,c,z = 0,1,2,n-3,n-2,n-1
    d,w = tuple(range(3,n-2)),tuple(range(2,n-2))
    if (b,c) in arcs:
        if y not in incoming:
            return 'open_release', ((b,)+w+(c,z,a), (z,)+w+(c,b,a), (a,z)+w+(c,b))
        if c not in incoming:
            return 'open_penultimate', ((a,b,c,z,x)+d, (b,c,a,z,x)+d, (b,x)+d+(c,a,z))
        interior = tuple(range(3,n-3))
        return 'open_endpair', ((a,b,c,z,x)+interior+(y,), (b,c,z,x)+interior+(a,y), (b,x)+interior+(a,y,c,z))
    if c in incoming:
        if y not in incoming:
            h = next(u for u in incoming if u != c)
            k = n-1-h
            left,right = tuple(range(1,k)),tuple(range(k,n-2))
            return 'closed_movable_cut', ((a,)+left+(z,)+right+(c,), left+(z,)+right+(a,c), left+right+(a,c,z))
        if (x,c) in arcs:
            q = (a,x)+d+(c,z,b)
            s = (a,b,x,c)+d[::-1]+(z,)
            v = (b,z)+d+(c,x,a)
            return 'closed_endpair_two_block', (q,s[::-1],v)
        interior = tuple(range(3,n-3))
        q1 = (a,x)+interior+(y,c,z,b)
        if (y,x) in arcs:
            return 'closed_endpair_yx', (q1, interior+(a,x,y,c,z,b), interior+(a,b,x,y,c,z))
        return 'closed_endpair_xy', (q1, (z,)+interior+(y,x,b,a,c), (b,z)+interior+(y,x,a,c))
    if (c,x) in arcs:
        return 'closed_epsilon_zero', ((c,)+w+(b,z,a), (z,c)+w+(b,a), (a,z,b)+w+(c,))
    if x not in incoming:
        return 'closed_epsilon_one', ((c,b,z,a)+w, (z,c,b,a)+w, (a,z,b)+w+(c,))
    if 3 not in incoming:
        interior = tuple(range(3,n-3))
        return 'closed_second', ((a,b,x)+interior+(z,y,c), (a,x)+interior+(z,y,b,c), (a,x)+interior+(y,b,c,z))
    interior,e = tuple(range(3,n-4)),n-4
    return 'closed_second_pair', ((a,b,x)+interior+(z,e,y,c), (a,x)+interior+(z,e,y,b,c), (a,x)+interior+(e,y,b,c,z))



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'receipt.json')
    parser.add_argument('--certificates', type=Path, default=HERE/'certificates.json')
    parser.add_argument('--repo', type=Path, help='optionally verify frozen inherited sources')
    args = parser.parse_args()
    pins = json.loads((HERE/'source_pins.json').read_text())
    if args.repo:
        for p, h in pins['sources'].items():
            need(hashlib.sha256((args.repo/p).read_bytes()).hexdigest() == h, 'source hash changed: '+p)
    exceptions = json.loads((HERE/'exceptions.json').read_text())
    indexed = {(tuple(r['incoming']), r['mask']): r for r in exceptions}
    need(len(indexed) == 3, 'three distinct exception repairs')
    cells, positive, counts, direction_checks = [], set(), Counter(), 0
    for incoming in combinations(range(2, 9), 2):
        original, variables = source_forms(incoming)
        for bits in product(range(2), repeat=3):
            guards = [equation(original[p], bit) for p, bit in zip(CUT_PAIRS, bits)]
            g = Affine(guards)
            need(g.contradiction is None and len(g.basis) == 3, 'nonempty disjoint rank-three cell')
            sample = sum(bit * original[p][0] for p, bit in zip(CUT_PAIRS, bits))
            need(all(parity(m & sample) == c for m, c in guards), 'cell sample')
            branch, words = ledger(10, incoming, realize(original, sample))
            now, proofs = original, []
            for step, q in enumerate(words):
                for p in zip(q, q[1:]):
                    failure = Affine(guards + [equation(now[p], 0)])
                    need(failure.contradiction is not None, 'ledger legality not uniform')
                    need(failure.proof_value(failure.contradiction) == (0, 1), 'legality proof')
                    direction_checks += 1
                space = Affine(guards + anti_equations(now, q))
                cert = space.certificate()
                if 'assignment' in cert:
                    key = (incoming, cert['assignment'])
                    need(key in indexed, 'unrepaired exceptional presentation')
                    positive.add(key)
                    counts['singleton_blue_systems'] += 1
                else:
                    counts['inconsistent_blue_systems'] += 1
                proofs.append(cert)
                now = move_forms(now, q)
            need(now == move_forms(original, tuple(range(10))), 'symbolic literal target')
            cells.append({'incoming': incoming, 'cell': bits, 'branch': branch, 'words': words, 'proofs': proofs})
    need(positive == set(indexed), 'exactly the three claimed exceptions')
    replay = []
    for record in exceptions:
        incoming = tuple(record['incoming'])
        F, _ = source_forms(incoming)
        initial = realize(F, record['mask'])
        rows = [sum(1 << v for u, v in initial if u == i) for i in range(10)]
        need(rows == record['rows'], 'exception source rows')
        need(sorted(9-i for i in incoming) == record['outgoing'], 'reflected ports')
        target = realize(move_forms(F, tuple(range(10))), record['mask'])
        target_rows = [sum(1 << v for u, v in target if u == i) for i in range(10)]
        need(target_rows == record['partner_rows'], 'exception target rows')
        witnesses = direct_replay(initial, record['words'])
        replay.append({'incoming': incoming, 'mask': record['mask'], 'anti_failures': witnesses})
    # Sharp cyclic-cut hostile from the inherited theorem, not a failure of reachability.
    hostile_rows = [190, 44, 760, 176, 994, 704, 907, 770, 559, 11]
    hostile = {(u, v) for u in range(10) for v in range(10) if hostile_rows[u] >> v & 1}
    cut = tuple(range(1, 10)) + (0,)
    need(set(zip(cut, cut[1:])) <= hostile, 'hostile cut directed')
    r = dict(zip(cut, cut[::-1]))
    need(all((r[v], r[u]) in hostile for u, v in hostile), 'hostile cut actually blue')
    # The same three words as the end-pair exception also repair this hostile.
    hostile_witnesses = direct_replay(hostile, exceptions[2]['words'])
    payload = {'order': 10, 'dimension': D, 'cut_pairs': CUT_PAIRS, 'variables': variables, 'cells': cells}
    certificate_text = json.dumps(payload, indent=2) + '\n'
    args.certificates.write_text(certificate_text)
    own_hashes = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ['check.py', 'exceptions.json', 'source_pins.json']}
    receipt = {'status': 'PASS: complete exact affine certificate',
        'claim': 'Every marked closed blue order-ten source with endpoint scores (6,3) has a literal three-black K detour.',
        'filters': 'None: no isomorphism, source/target mixedness, or strong-partner filter.',
        'source_universe_size': 21 * 2**D, 'complete_source_census_performed': False,
        'cell_representatives_constructed': len(cells),
        'port_pairs': 21, 'cells': len(cells), 'sources_per_cell': 2**9,
        'blue_systems': sum(counts.values()), **dict(counts),
        'symbolic_direction_checks': direction_checks,
        'literal_exception_replays': replay,
        'cyclic_cut_hostile_repaired': hostile_witnesses,
        'source_hashes': own_hashes, 'inherited_repository_commit': pins['repository_commit'],
        'certificates_sha256': hashlib.sha256(certificate_text.encode()).hexdigest()}
    need(len(cells) == 168 and counts == {'singleton_blue_systems': 3, 'inconsistent_blue_systems': 501}, 'expected exact reduction')
    args.output.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
