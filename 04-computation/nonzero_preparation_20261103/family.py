#!/usr/bin/env python3
"""Complete one-variable-run inverse for the explicit ready funnel bank.

The prefix-length bound is proved from the requested positive dyadic slope.
No chosen prefix/suffix word-length cap is used. Actual source rearming and
all-family acceptance are different questions from this restricted inverse.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from surgery_ready import BANK, core_data, ready_count

CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def step(state, bit):
    T, B, eta, L, Z = state
    return (2*T-bit, B if bit == 0 else B+T-eta, bit, L+bit*Z, 2*Z)


def state(word):
    s = (1, 1, 1, 1, 1)
    for b in word:
        s = step(s, b)
    return s


def affine(data, prefix, bit, suffix):
    """Exact (A,B) for p bit^k q, k>=1, as A*2^k+B."""
    T, B, eta, L, Z = state(prefix)
    if bit == 0:
        pairs = [(T, 0), (0, B), (0, L), (Z, 0)]
    else:
        pairs = [(T-1, 1), (T-1, B-T+2-eta), (Z, L-Z), (Z, 0)]
    eta = bit
    # Coordinates T,B,L,Z, each coefficient of X followed by constant.
    for b in suffix:
        t, v, l, z = pairs
        pairs = [(2*t[0], 2*t[1]-b),
                 v if b == 0 else (v[0]+t[0], v[1]+t[1]-eta),
                 (l[0]+b*z[0], l[1]+b*z[1]), (2*z[0], 2*z[1])]
        eta = b
    _, B, L, _ = pairs
    D = 1-(prefix[0] if prefix else bit)
    C = data['G']+data['p2']
    return (C*L[0]+data['h']*B[0],
            C*L[1]+data['h']*B[1]-data['e']*D-data['C']-data['F']-data['b'])


def bits(v):
    return tuple((v >> k) & 1 for k in range(v.bit_length()))


def inverse(A, B, bank=BANK):
    need(A > 0, 'positive family slope')
    hits, queries, prefixes = [], 0, 0
    for idx, (gap, rows, ports) in enumerate(bank):
        data = core_data(rows, ports)
        C, h = data['G']+data['p2'], data['h']
        n = 0
        while C * (1 << n) <= A:
            for prefix in product((0, 1), repeat=n):
                T, _, _, _, Z = state(prefix)
                prefixes += 1
                for bit in (0, 1):
                    divisor = C*Z+h*(T-bit)
                    need(divisor >= C*(1 << n), 'proved prefix length bound')
                    if A % divisor:
                        continue
                    value = A//divisor-bit
                    if value < (1-bit):
                        continue
                    suffix = bits(value)
                    found = affine(data, prefix, bit, suffix)
                    need(found[0] == A, 'slope factorization')
                    for k in (1, 2, 3, 5):
                        need(ready_count(data, prefix+(bit,)*k+suffix)
                             == found[0]*(1 << k)+found[1], 'literal word recurrence')
                    queries += 1
                    if found[1] == B:
                        hits.append({'bank_index': idx, 'gap': gap, 'core': rows,
                                     'ports': ports, 'prefix': prefix, 'bit': bit,
                                     'suffix': suffix, 'A': A, 'B': B})
            n += 1
    return {'A': A, 'B': B, 'prefixes': prefixes, 'divisor_queries': queries, 'hits': hits}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--depth', type=int, default=3)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_suffix('.json'))
    args = parser.parse_args()
    need(0 <= args.depth <= 6, 'declared finite requested family set')
    controls = 0
    for length in range(7):
        for word in product((0, 1), repeat=length):
            T, B, eta, L, Z = state(word)
            from surgery_ready import tail_coefficients
            need(L == tail_coefficients(word, 2)[0], 'binary L recurrence')
            need(B == tail_coefficients(word, 1)[2], 'binary B recurrence')
            need(Z == 1 << length and T >= 1, 'state positivity')
            controls += 1
    frontier = [('', 756, -83)]
    requested = []
    for depth in range(args.depth+1):
        nxt = []
        for address, A, B in frontier:
            result = inverse(A, B)
            result['address'] = address
            requested.append(result)
            if depth < args.depth:
                for d in (2, 4):
                    nxt.append((address+str(d), 3*A, 3*B+d))
        frontier = nxt
    # Positive inherited shared family, with repeated one-run clock.
    positive = inverse(90, 5)
    need({r['gap'] for r in positive['hits']} == {1, 2}, 'inherited common-family positive')
    receipt = {'scope': 'complete one-variable-run template inverse in the explicit 34-entry bank; not all ready families',
               'requested_depth': args.depth, 'state_controls': controls,
               'requested': requested, 'positive_control': positive, 'checks': CHECKS,
               'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
               'dependency_sha256': sha256(Path(__file__).with_name('surgery_ready.py').read_bytes()).hexdigest()}
    raw = (json.dumps(receipt, indent=2, sort_keys=True)+'\n').encode()
    args.output.write_bytes(raw)
    print(json.dumps({'checks': CHECKS, 'families': len(requested),
                      'accepted': {r['address']: sorted({h['gap'] for h in r['hits']})
                                   for r in requested if r['hits']},
                      'sha256': sha256(raw).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()
