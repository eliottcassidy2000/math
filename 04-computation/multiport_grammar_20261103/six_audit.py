#!/usr/bin/env python3
"""Independent strong-six coverage, signature inverse and named44 audit.

No producer implementation is imported.  The producer JSON supplies candidate
class representatives, whose complete, disjoint labelled orbits are verified
against an independent enumeration of every labelled tournament.  Atom counts
come from full permutations/cuts; query states are independent signature
counters under a strictly larger Hamilton envelope.
"""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path
import hashlib
import json
import math

GATES = 0


def require(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def graphs(n):
    edges = tuple(combinations(range(n), 2))
    for bits in range(1 << len(edges)):
        rows = [0] * n
        for k, (u, v) in enumerate(edges):
            if (bits >> k) & 1:
                rows[u] |= 1 << v
            else:
                rows[v] |= 1 << u
        yield tuple(rows)


def strong(rows):
    n = len(rows)
    for start in range(n):
        reached, queue = {start}, [start]
        for u in queue:
            for v in range(n):
                if (rows[u] >> v) & 1 and v not in reached:
                    reached.add(v)
                    queue.append(v)
        if len(reached) != n:
            return False
    return True


def relabel(rows, labels):
    n = len(rows)
    # Send old vertex i to labels[i], unlike the producer's pullback convention.
    target = [0] * n
    for i in range(n):
        for j in range(n):
            if (rows[i] >> j) & 1:
                target[labels[i]] |= 1 << labels[j]
    return tuple(target)


@lru_cache(None)
def tensor(rows):
    n = len(rows)
    h = p = 0
    F, G, E, S = ([0] * n for _ in range(4))
    J = [[0] * n for _ in range(n)]
    for path in permutations(range(n)):
        good = [bool((rows[path[i]] >> path[i + 1]) & 1)
                for i in range(n - 1)]
        if all(good):
            h += 1
            E[path[-1]] += 1
            S[path[0]] += 1
        for cut in range(n + 1):
            if not all(good[:max(0, cut - 1)]) or not all(good[cut:]):
                continue
            p += 1
            if cut:
                F[path[cut - 1]] += 1
            if cut < n:
                G[path[cut]] += 1
            if 0 < cut < n:
                J[path[cut - 1]][path[cut]] += 1
    return h, p, F, G, J, E, S


def literal_signature(rows, ports):
    h, p, F, G, J, E, S = tensor(tuple(rows))
    W = [v for v in range(len(rows)) if (ports >> v) & 1]
    return (h, p, sum(F[v] for v in W), sum(G[v] for v in W),
            sum(J[u][v] for u in W for v in W),
            sum(E[v] for v in W), sum(S[v] for v in W))


def multiply(a, b):
    h, p, f, g, j, e, s = a
    H, P, F, G, J, E, S = b
    return (h * H, p * P, f * H + p * F, g * P + h * G,
            j * H + h * J + e * S + g * F, h * E, s * H)


def odd(n):
    return n // (n & -n)


def attach(rows, ports, clock):
    n = len(rows)
    root = n + clock
    block = ((1 << n) - 1) << clock
    answer = [sum(1 << j for j in range(i + 1, clock)) | block
              for i in range(clock)]
    answer.extend((row << clock) | ((1 << root) if (ports >> i) & 1 else 0)
                  for i, row in enumerate(rows))
    answer.append(((1 << clock) - 1) | (block ^ (ports << clock)))
    return tuple(answer)


def endpoint_signature(rows, ports):
    """Independent last-vertex recurrence for the larger actual attachments."""
    n = len(rows)
    size = 1 << n
    predecessors = tuple(sum(1 << u for u in range(n) if (rows[u] >> v) & 1)
                         for v in range(n))

    def endings(pred):
        table = [[0] * n for _ in range(size)]
        for v in range(n):
            table[1 << v][v] = 1
        for mask in range(1, size):
            if not mask & (mask - 1):
                continue
            for v in range(n):
                if not (mask >> v) & 1:
                    continue
                old = mask ^ (1 << v)
                available = old & pred[v]
                while available:
                    bit = available & -available
                    available -= bit
                    table[mask][v] += table[old][bit.bit_length() - 1]
        return table

    end, start = endings(predecessors), endings(rows)
    H = [sum(row) for row in end]
    H[0] = 1
    E = [sum(row[v] for v in range(n) if (ports >> v) & 1) for row in end]
    S = [sum(row[v] for v in range(n) if (ports >> v) & 1) for row in start]
    full = size - 1
    return (H[full], sum(H[m] * H[full ^ m] for m in range(size)),
            sum(E[m] * H[full ^ m] for m in range(size)),
            sum(H[m] * S[full ^ m] for m in range(size)),
            sum(E[m] * S[full ^ m] for m in range(size)), E[full], S[full])


def main():
    parser = ArgumentParser()
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_suffix('.json'))
    parser.add_argument('--tree', action='store_true',
                        help='also independently audit the complete scalar tree receipt')
    args = parser.parse_args()
    folder = Path(__file__).parent
    receipt_bytes = (folder / 'atom_six_probe.json').read_bytes()
    receipt = json.loads(receipt_bytes)
    require(receipt['shift_inverse'] is not None, 'full shifted receipt required')
    bank = set()
    census = []
    raw_marked_six = 0
    six_sig = set()
    for item in receipt['class_census']:
        n = item['n']
        labelled = set(graphs(n))
        actual_strong = {rows for rows in labelled if strong(rows)}
        orbit_union, strong_union = set(), set()
        strong_classes = 0
        for rep in item['representatives']:
            rows = tuple(rep['rows'])
            require(rows in labelled, 'actual representative tournament')
            orbit = {relabel(rows, labels) for labels in permutations(range(n))}
            require(min(orbit) == rows, 'independent canonical representative')
            require(not orbit_union.intersection(orbit), 'class orbits are disjoint')
            orbit_union.update(orbit)
            require(len(orbit) * rep['automorphisms'] == math.factorial(n),
                    'independent stabilizer orbit mass')
            connected = strong(rows)
            require(connected == rep['strong'], 'independent representative strongness')
            if connected:
                strong_classes += 1
                strong_union.update(orbit)
                for ports in range(1 << n):
                    sig = literal_signature(rows, ports)
                    bank.add(sig)
                    if n == 6:
                        raw_marked_six += 1
                        six_sig.add(sig)
        require(orbit_union == labelled, 'complete labelled tournament orbit cover')
        require(strong_union == actual_strong, 'complete independent strong orbit cover')
        require(len(actual_strong) == item['strong_labelled_mass'],
                'independent labelled strong count')
        census.append({'n': n, 'classes': len(item['representatives']),
                       'strong_classes': strong_classes,
                       'labelled': len(labelled), 'strong_labelled': len(actual_strong)})
    require(census[-1] == {'n': 6, 'classes': 56, 'strong_classes': 35,
                          'labelled': 32768, 'strong_labelled': 22320},
            'full sixth-order census')
    require(raw_marked_six == 2240 and len(six_sig) == 1950,
            'all sixth-order marked classes and signatures')
    require(min(a[0] for a in six_sig) == 15 and
            min(a[2] for a in six_sig if a[2]) == 15,
            'independent sixth-order lower statistics')
    require(not [a for a in six_sig if a[2] and a[0] - a[4] in (-35, -29)
                 and odd(a[2]) == 63], 'neither target has a single-six supplier')
    supplied = {tuple(a['sigma']) for a in receipt['bank']}
    require(bank == supplied and len(bank) == len(receipt['bank']) == 2120,
            'complete independently reconstructed signature bank')
    for atom in receipt['bank']:
        require(literal_signature(tuple(atom['rows']), atom['ports']) == tuple(atom['sigma']),
                'each selected representative independently recounted')
    atoms = tuple(sorted(bank))
    zero = tuple(a for a in atoms if a[2] == 0)
    marked = tuple(a for a in atoms if a[2] > 0)
    M = max(a[0] for a in atoms)
    require(M == 45, 'Hamilton maximum')
    by_h = defaultdict(lambda: defaultdict(list))
    for a in atoms:
        by_h[a[0]][a[1:3]].append(a)

    @lru_cache(None)
    def inverse(h, f, first_marked):
        if h < 1 or f < h:
            return ()
        answer = Counter()
        for H, groups in by_h.items():
            if h % H:
                continue
            tail_h = h // H
            for (P, F), bucket in groups.items():
                if first_marked and F == 0:
                    continue
                if F and (h, f) == (H, F):
                    answer.update(bucket)
                remainder = f - F * tail_h
                if remainder <= 0 or remainder % P:
                    continue
                tails = inverse(tail_h, remainder // P, False)
                for b, number in tails:
                    for a in bucket:
                        answer[multiply(a, b)] += number
        return tuple(answer.items())

    def h_values(bound):
        seen, queue = {1}, [1]
        for h in queue:
            for H in by_h:
                value = h * H
                if value <= bound and value not in seen:
                    seen.add(value)
                    queue.append(value)
        return sorted(seen)

    def prefixes(c, R):
        # Preserve one literal (h,p) as well as the quotient's (h,odd(p)).
        out = {(1, 1): (1, 1, True)}
        queue = [(1, 1)]
        for h, op in queue:
            ph, pp, _ = out[h, op]
            for a in zero:
                target = h * a[0], op * odd(a[1])
                if abs(c) % target[0] or R % target[1] or target in out:
                    continue
                out[target] = (ph * a[0], pp * a[1], False)
                queue.append(target)
        return out

    decisions = []
    total_queries = 0
    for test in receipt['shift_inverse']:
        c, R = test['c'], test['R']
        answers = Counter()
        states = prefixes(c, R)
        queries = []
        require(set(states) == {(s['prefix_h'], s['prefix_odd_p']) for s in test['searches']},
                'independent unmarked-prefix quotient')
        for (hp, op), (ph, pp, empty) in sorted(states.items()):
            cc, rr = c // hp, R // op
            # Twice the largest bank-product bound, plus |c|: strictly larger.
            B = 2 * M * M + abs(cc)
            cap = max(max(a[2] for a in marked),
                      max(((a[2] + a[1]) * B - a[1] * cc) // a[0] for a in marked))
            hs = h_values(B)
            f = rr
            while f <= cap:
                values = Counter()
                for h in hs:
                    if h > f:
                        break
                    values.update(dict(inverse(h, f, True)))
                for sig, multiplicity in values.items():
                    if sig[0] - sig[4] != cc:
                        continue
                    full = sig if empty else (
                        ph * sig[0], pp * sig[1], pp * sig[2], ph * sig[3],
                        ph * sig[4], ph * sig[5], 0)
                    require(full[0] - full[4] == c and odd(full[2]) == R,
                            'prefix reconstruction of accepted signature')
                    answers[full] += multiplicity
                matching = [s for s in test['searches'] if
                            (s['prefix_h'], s['prefix_odd_p'], s['f']) == (hp, op, f)]
                if matching:
                    s = matching[0]
                    tight_B = max(M * M, abs(cc))
                    tight_cap = max(max(a[2] for a in marked), max(
                        ((a[2] + a[1]) * tight_B - a[1] * cc) // a[0] for a in marked))
                    require(s['proved_h_bound'] == tight_B and s['proved_f_bound'] == tight_cap,
                            'independent tight receipt bounds')
                    require(sum(number for sig, number in values.items() if sig[0] <= tight_B)
                            == s['candidates'], 'independent normalized word multiplicities')
                queries.append({'prefix': [hp, op], 'f': f, 'h_bound': B,
                                'f_bound': cap, 'candidate_words': sum(values.values())})
                total_queries += 1
                f *= 2
        expected = Counter(tuple(w['sigma']) for w in test['witnesses'])
        require(answers == expected, 'independent complete target decision ' + test['address'])
        decisions.append({'address': test['address'], 'accepted': sum(answers.values()),
                          'signature_count': len(answers), 'queries': queries})
    require([(d['address'], d['accepted']) for d in decisions] == [('24', 0), ('44', 1)],
            'all-shifts decisions')

    tree_audit = None
    if args.tree:
        tree_bytes = (folder / 'six_tree_probe.json').read_bytes()
        tree = json.loads(tree_bytes)
        tree_decisions = []
        wider_tree_queries = 0
        for test in tree['tests']:
            c, R = test['c'], test['R']
            answers = Counter()
            states = prefixes(c, R)
            require(set(states) == {tuple(s['prefix']) for s in test['searches']},
                    'tree prefix quotient ' + test['address'])
            for (hp, op), (ph, pp, empty) in sorted(states.items()):
                cc, rr = c // hp, R // op
                B = 2 * M * M + abs(cc)
                cap = max(max(a[2] for a in marked), max(
                    ((a[2] + a[1]) * B - a[1] * cc) // a[0] for a in marked))
                hs = h_values(B)
                f = rr
                while f <= cap:
                    values = Counter()
                    for h in hs:
                        if h > f:
                            break
                        values.update(dict(inverse(h, f, True)))
                    for sig, multiplicity in values.items():
                        if sig[0] - sig[4] != cc:
                            continue
                        full = sig if empty else (
                            ph * sig[0], pp * sig[1], pp * sig[2], ph * sig[3],
                            ph * sig[4], ph * sig[5], 0)
                        require(full[0] - full[4] == c and odd(full[2]) == R,
                                'tree accepted target reconstruction')
                        answers[full] += multiplicity
                    matching = [s for s in test['searches']
                                if (tuple(s['prefix']), s['f']) == ((hp, op), f)]
                    if matching:
                        s = matching[0]
                        tight_B = max(M * M, abs(cc))
                        tight_cap = max(max(a[2] for a in marked), max(
                            ((a[2] + a[1]) * tight_B - a[1] * cc) // a[0] for a in marked))
                        require(s['c'] == cc and s['h_bound'] == tight_B and
                                s['f_bound'] == tight_cap, 'tree tight bounds')
                        require(sum(number for sig, number in values.items()
                                    if sig[0] <= tight_B) == s['candidates'],
                                'tree independent candidate multiplicities')
                    wider_tree_queries += 1
                    f *= 2
            require(answers == Counter(tuple(w['sigma']) for w in test['witnesses']),
                    'independent complete tree successor ' + test['address'])
            for witness in test['witnesses']:
                rows, W = tuple(witness['rows']), witness['ports']
                require(endpoint_signature(rows, W) == tuple(witness['sigma']),
                        'tree witness independently recounted')
            tree_decisions.append({'address': test['address'],
                                   'accepted_words': sum(answers.values()),
                                   'accepted_signatures': len(answers)})
        declared = {s['address']: (s['c'], s['R']) for s in tree['states']}
        accepted = {''} | {d['address'] for d in tree_decisions if d['accepted_words']}
        require(set(declared) == accepted and declared[''] == (-5, 7),
                'complete declared scalar state set')
        require(len(tree_decisions) == 2 * len(declared), 'both children of every state tested')
        for address, (c, R) in declared.items():
            for digit in (2, 4):
                child = [t for t in tree['tests'] if t['address'] == address + str(digit)]
                require(len(child) == 1 and (child[0]['c'], child[0]['R']) == (3 * c + digit, 3 * R),
                        'actual arithmetic successor')
        require(tree['complete_tree'] is True and not tree['frontier'] and
                tree['completed_depth'] < tree['requested_depth_cap'],
                'empty frontier before chosen cap')
        require(accepted == {'', '2', '4', '22', '42', '44', '422', '442', '444',
                             '4222', '4442', '4444', '42222', '44444'},
                'entire fourteen-state scalar tree')
        tree_audit = {
            'source_sha256': hashlib.sha256((folder / 'atom_six_probe_tree.py').read_bytes()).hexdigest(),
            'receipt_sha256': hashlib.sha256(tree_bytes).hexdigest(),
            'states': sorted(accepted), 'decisions': tree_decisions,
            'wider_queries': wider_tree_queries,
            'producer_queries': sum(len(t['searches']) for t in tree['tests']),
            'producer_candidate_visits': sum(s['candidates'] for t in tree['tests'] for s in t['searches']),
            'complete_tree': True,
        }

    # Independently pay the actual old-name-preserving adjacency rewrite.
    chain5 = tuple(sum(1 << j for j in range(i + 1, 5)) for i in range(5))
    initial = tuple(sum(1 << j for j in range(i + 1, 7)) for i in range(7))
    target = (66, 68, 121, 67, 75, 91, 0)
    ports = (1 << 2) | (1 << 4) | (1 << 6)
    flips = ((0, 2), (0, 3), (0, 4), (0, 5), (1, 3),
             (1, 4), (1, 5), (3, 4), (3, 5), (4, 5))
    rows = list(initial)
    states = [(tuple(rows), ports)]
    for u, v in flips:
        require((rows[u] >> v) & 1 and not (rows[v] >> u) & 1, 'actual old arc reversed')
        rows[u] ^= 1 << v
        rows[v] ^= 1 << u
        states.append((tuple(rows), ports))
    ports ^= 1 << 2
    states.append((tuple(rows), ports))
    require(tuple(rows) == target and ports == 80, 'literal eleven-step endpoint')
    require(literal_signature(target, ports) == (17, 448, 252, 105, 46, 17, 4),
            'independent seven-vertex target signature')
    require(literal_signature(initial, 84) == (1, 128, 84, 21, 12, 1, 0),
            'reblocked source signature')
    rewrite = []
    for step, (rows, W) in enumerate(states):
        require(rows[6] == 0 and all((rows[v] >> 6) & 1 for v in range(6)) and (W >> 6) & 1,
                'terminal marked singleton survives every step')
        sig = literal_signature(rows, W)
        observed = []
        for t in (3, 4, 5):
            graph = attach(rows, W, t - 2)
            actual = endpoint_signature(graph, 1 << (len(graph) - 1))
            require(strong(graph), 'actual intermediate full graph strong')
            require(actual[0] == sig[0] - sig[4] + (1 << (t - 2)) * sig[2],
                    'actual intermediate Hamilton count')
            require(actual[2] - actual[0] == sig[4], 'actual root defect')
            if step == 0:
                require(graph == attach(chain5, 21, t), 'reblocking preserves every named arc')
            if step == len(states) - 1:
                require(actual[0] == 3 * (21 * 2 ** t - 11) + 4,
                        'literal uniform44 positive')
            observed.append([t, actual[0], actual[2], actual[2] - actual[0]])
        rewrite.append({'step': step, 'suffix_sigma': sig, 'observed': observed})
    empty_prefix_count = endpoint_signature(attach(target, 80, 0), 0)[0]
    require(empty_prefix_count == 219 and 3 * (21 * 2 ** 2 - 11) + 4 == 223,
            'clock-two bulk boundary hostile')

    output = {
        'status': 'INDEPENDENT ANALYTIC/SOURCE/FINITE-EXACT AUDIT PASS',
        'scope': 'Complete <=6 atom quotient and all-shifts24/44 scalar decisions; named strong11-flip44 rewrite for t>=3. No native or reusable closure claim.',
        'supplier_source_sha256': hashlib.sha256((folder / 'atom_six_probe.py').read_bytes()).hexdigest(),
        'supplier_receipt_sha256': hashlib.sha256(receipt_bytes).hexdigest(),
        'gates': GATES, 'class_census': census, 'raw_six_port_cases': raw_marked_six,
        'six_signatures': len(six_sig), 'bank_signatures': len(bank),
        'bank_sha256': hashlib.sha256(json.dumps(sorted(bank), separators=(',', ':')).encode()).hexdigest(),
        'wider_queries': total_queries, 'decisions': decisions,
        'named_adjacency_rewrite': rewrite,
        'tree_audit': tree_audit,
        'bulk_at_original_clock_two': {'actual': empty_prefix_count, 'required': 223,
                                       'outside_nonempty_prefix_domain': True},
    }
    content = json.dumps(output, sort_keys=True, indent=2) + '\n'
    args.output.write_text(content)
    print(json.dumps({'gates': GATES, 'wider_queries': total_queries,
                      'decisions': [(d['address'], d['accepted']) for d in decisions],
                      'tree_states': len(tree_audit['states']) if tree_audit else None,
                      'receipt_sha256': hashlib.sha256(content.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
