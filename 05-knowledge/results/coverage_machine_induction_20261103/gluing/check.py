#!/usr/bin/env python3
"""Exact two-port cyclic-gluing transfer and independent small controls.

No repository imports. Runtime checks remain enabled under python -O.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial, gcd
import json


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def skeletons(k):
    parts = [(a, b) for a in range(-1, k) for b in range(-1, k)
             if (a, b) != (-1, -1)]
    ans = []
    for size in range(k+1):
        for subset in combinations(parts, size):
            used = [v for part in subset for v in set(part) if v >= 0]
            if len(used) == len(set(used)):
                ans.append(tuple(sorted(subset)))
    return tuple(sorted(ans))


STATES = skeletons(2)


def split(parts):
    parts = tuple(sorted(parts))
    return tuple(p for p in parts if p != (-1, -1)), parts.count((-1, -1))


def kernel(skeleton, m):
    """Insert a,b into atomic old paths; old->a->b->old uniformly."""
    parts = list(skeleton) + [(-1, -1)]*m
    a, b = len(parts), len(parts)+1
    ans = Counter()
    for incoming, outgoing, ab in product([None]+list(range(a)),
                                          [None]+list(range(a)), range(2)):
        if ab and incoming is not None and incoming == outgoing:
            continue  # old path -> a -> b -> same old path is a cycle.
        nxt, prev = {}, {}
        for u, v in ((incoming, a), (b, outgoing), (a if ab else None, b)):
            if u is not None and v is not None:
                nxt[u] = v
                prev[v] = u
        endpoints = parts + [(-1, -1), (-1, -1)]
        output, reached = [], set()
        for head in range(a+2):
            if head in prev:
                continue
            tail = head
            reached.add(tail)
            while tail in nxt:
                tail = nxt[tail]
                need(tail not in reached, 'no component cycle')
                reached.add(tail)
            output.append((endpoints[head][0], endpoints[tail][1]))
        need(len(reached) == a+2, 'every atomic component occurs')
        ans[split(output)] += 1
    return ans


def shift(poly, power, multiplier=1):
    return Counter({(s, d+power): multiplier*v for (s, d), v in poly.items()})


def subtract(a, b):
    out = Counter(a)
    for key, value in b.items():
        out[key] -= value
    return Counter({key: value for key, value in out.items() if value})


def operator():
    ops = []
    for state in STATES:
        A = kernel(state, 0)
        B = subtract(kernel(state, 1), shift(A, 1))
        twiceC = subtract(subtract(kernel(state, 2), shift(A, 2)), shift(B, 1, 2))
        need(all(v % 2 == 0 for v in twiceC.values()), 'integral second derivative')
        C = Counter({key: value//2 for key, value in twiceC.items()})
        need(all(v > 0 and 0 <= degree <= 2
                 for block in (A, B, C) for (_, degree), v in block.items()),
             'positive degree-at-most-two differential matrices')
        ops.append((A, B, C))
    return ops


def apply(vector, ops):
    ans = Counter()
    for state, blocks in zip(STATES, ops):
        for (old, m), value in vector.items():
            if old != state:
                continue
            for j, block in enumerate(blocks):
                if m < j:
                    continue
                falling = factorial(m)//factorial(m-j)
                for (target, power), weight in block.items():
                    ans[target, m-j+power] += value*falling*weight
    return +ans


def literal(T, ports=(0, 1)):
    """Independent permutation/cut enumeration of all unordered covers."""
    slot = {v: i for i, v in enumerate(ports)}
    n, covers = len(T), set()
    for word in permutations(range(n)):
        for cuts in range(1 << (n-1)):
            parts, start = [], 0
            for i in range(n-1):
                if cuts >> i & 1:
                    parts.append(word[start:i+1])
                    start = i+1
                elif not (T[word[i]] >> word[i+1] & 1):
                    break
            else:
                parts.append(word[start:])
                covers.add(tuple(sorted(parts)))
    ans = Counter()
    for cover in covers:
        ans[split([(slot.get(p[0], -1), slot.get(p[-1], -1)) for p in cover])] += 1
    return ans


@lru_cache(None)
def H(T):
    n = len(T)
    if n == 0:
        return 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, value in enumerate(dp[mask]):
            choices = T[v] & ~mask
            while value and choices:
                bit = choices & -choices
                choices -= bit
                dp[mask | bit][bit.bit_length()-1] += value
    return sum(dp[-1])


def G(T):
    n = len(T)
    return tuple(row | (1 << n) for row in T) + (1 << (n+1), (1 << n)-1)


def scalar(vector):
    return sum(v for (s, m), v in vector.items() if len(s)+m == 1)


def unmarked(vector):
    out = Counter()
    for (s, m), v in vector.items():
        out[len(s)+m] += v
    return out


def slot(poly, r):
    return sum(v*factorial(r)//factorial(r-k) for k, v in poly.items() if k <= r)


def unmarked_G(poly):
    ans = Counter()
    for j, v in poly.items():
        for k, coefficient in ((j-1, j*(j-1)), (j, j*(j+2)),
                               (j+1, 2*j+1), (j+2, 1)):
            if coefficient:
                ans[k] += v*coefficient
    return +ans


def transitive(n):
    return tuple(sum(1 << j for j in range(i+1, n)) for i in range(n))


def triangle(A, B, C):
    blocks = (A, B, C)
    starts = (0, len(A), len(A)+len(B))
    result = []
    for i, block in enumerate(blocks):
        j = (i+1) % 3
        cross = ((1 << len(blocks[j]))-1) << starts[j]
        result.extend((row << starts[i]) | cross for row in block)
    return tuple(result)


def fixed_divisor(poly):
    return gcd(*(factorial(k)*v for k, v in poly.items()))


def ca(values):
    M = len(values)
    return [(i*i*values[i]+i*(values[(i+1) % M]-values[(i-1) % M])) % M
            for i in range(M)]


def history_coefficients(length):
    """h_j=sum_r a[j,r]F_r, from the adjoint slot difference rule."""
    rows, current = [], Counter({1: 1})
    for j in range(length):
        need(current[j+1] == factorial(j), 'triangular wrapper-history diagonal')
        rows.append(current)
        child = Counter()
        for r, value in current.items():
            child[r+1] += r*value
            child[r] += r*r*value
            if r > 1:
                child[r-1] -= r*value
        current = Counter({r: v for r, v in child.items() if v})
    return rows


def power_coefficients(poly):
    out, falling = Counter(), Counter({0: 1})
    for k in range(max(poly)+1):
        for degree, value in falling.items():
            out[degree] += poly.get(k, 0)*value
        nxt = Counter()
        for degree, value in falling.items():
            nxt[degree+1] += value
            nxt[degree] -= k*value
        falling = Counter({degree: value for degree, value in nxt.items() if value})
    return Counter({degree: value for degree, value in out.items() if value})


def jets(poly):
    coefficients = power_coefficients(poly)
    V = [sum(v*j**d for d, v in coefficients.items()) % 49 for j in range(7)]
    D = [sum(d*v*j**(d-1) for d, v in coefficients.items() if d) % 7 for j in range(7)]
    need(V[0] == 0, 'nonempty slot polynomial zero constant')
    return V, D


def jet_value(state, r):
    V, D = state
    a, j = divmod(r, 7)
    return (V[j]+7*a*D[j]) % 49


def jet_G(state):
    V, D = state
    W = [(j*j*V[j]+j*(jet_value(state, j+1)-jet_value(state, j-1))) % 49 for j in range(7)]
    E = [(2*j*V[j]+j*j*D[j]+V[(j+1) % 7]-V[(j-1) % 7]
          +j*(D[(j+1) % 7]-D[(j-1) % 7])) % 7 for j in range(7)]
    return W, E


def jet_product(A, B):
    V, D = A
    W, E = B
    return ([v*w % 49 for v, w in zip(V, W)],
            [(d*w+v*e) % 7 for v, d, w, e in zip(V, D, W, E)])


def parity_restore(state, epsilon):
    V, D = state
    W, E = [0]*7, [0]*7
    E[0] = D[0] if epsilon == -1 else 0
    for j in (1, 2, 3):
        W[j], E[j] = V[j], D[j]
        W[7-j] = epsilon*(V[j]-7*D[j]) % 49
        E[7-j] = -epsilon*D[j] % 7
    return W, E


def labelled(n):
    pairs = list(combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        rows = [0]*n
        for i, (a, b) in enumerate(pairs):
            u, v = (a, b) if bits >> i & 1 else (b, a)
            rows[u] |= 1 << v
        yield tuple(rows)


def birth(T, root):
    n = len(T)
    return tuple(row | ((1 << n) if i == root else 0) for i, row in enumerate(T)) + (((1 << n)-1) ^ (1 << root),)


def main():
    ops = operator()
    need(len(STATES) == 18 and len(skeletons(3)) == 88, 'two/three-port skeleton counts')
    checks = 0
    for s in STATES:
        for m in range(7):
            out = apply(Counter({(s, m): 1}), ops)
            need(out == kernel(s, m), 'universal multiplicity control')
            need(unmarked(out) == unmarked_G(Counter({len(s)+m: 1})), 'unmarked projection')
            checks += 1
    rows, ca_checks, escape_checks, jet_checks, jet_product_checks = [], 0, 0, 0, 0
    history_rows = history_coefficients(12)
    for n in (2, 3, 4):
        for T in labelled(n):
            initial, child = literal(T), G(T)
            out = apply(initial, ops)
            actual = literal(child)
            need(out == actual, 'literal complete named cover transfer')
            need(scalar(out) == H(child), 'independent Hamilton subset DP')
            p = unmarked(initial)
            need(H(child) == H(T)+slot(p, 2), 'cyclic count identity')
            for r in (1, 2, 3):
                need(slot(unmarked(out), r) == r*r*slot(p, r)+r*(slot(p, r+1)-slot(p, r-1)), 'slot-profile difference identity')
            need(initial[((0, 0), (1, 1)), n-2] == 1, 'all-singletons coefficient one')
            d = fixed_divisor(p)
            dp = fixed_divisor(unmarked(out))
            need(dp % d == 0, 'fixed slot divisor survives G')
            need(gcd(*(slot(p, r) for r in range(1, n+1))) == d,
                 'binomial-basis fixed divisor from n consecutive samples')
            weighted = gcd(*(factorial(r-1)*slot(p, r) for r in range(1, n+1)))
            history, current = [], p
            for j in range(12):
                history.append(current[1])
                need(current[1] == sum(v*slot(p, r) for r, v in history_rows[j].items()),
                     'triangular history identity')
                current = unmarked_G(current)
            need(gcd(*history[:n]) == weighted == gcd(*history),
                 'finite exact wrapper-content control')
            for ell in (2, 3, 5, 7, 11):
                hits = [j for j in range(ell-1) if history[j] % ell]
                need(bool(hits) == bool(d % ell), 'prime escape iff not in fixed slot divisor')
                if ell > n:
                    need(any(history[j] % ell for j in range(n)), 'large-prime bounded escape')
                escape_checks += 1
            for M in (2, 3, 7, 49):
                residue = [slot(p, i) % M for i in range(M)]
                need(all(slot(p, i+M) % M == residue[i] for i in range(M)),
                     'integer-coefficient slot polynomial has modular period M')
                need(ca(residue) == [slot(unmarked(out), i) % M for i in range(M)],
                     'radius-one modular cellular automaton')
                ca_checks += M
            current = p
            for _ in range(4):
                state = jets(current)
                need([jet_value(state, r) for r in range(49)] == [slot(current, r) % 49 for r in range(49)],
                     'thirteen-coordinate jet reconstruction')
                need(parity_restore(state, (-1)**n) == state, 'seven-coordinate parity jet reconstruction')
                nxt = unmarked_G(current)
                need(jet_G(state) == jets(nxt), 'exact derivative wrapper transfer')
                need([jet_value(jet_G(state), r) for r in range(49)] == ca([jet_value(state, r) for r in range(49)]),
                     'jet transfer equals all49-cell automaton')
                jet_checks += 1
                current = nxt
            # The separately constructed uniform context T -> singleton multiplies F by r.
            U = tuple(row | (1 << n) for row in T) + (0,)
            need(jet_product(jets(p), jets(Counter({1: 1}))) == jets(unmarked(literal(U))),
                 'ordinal product jet Leibniz rule against literal covers')
            jet_product_checks += 1
            rows.append({'rows': T, 'H': H(T), 'F2': slot(p, 2), 'wrapped_H': H(child), 'fixed_slot_divisor': d, 'wrapped_fixed_slot_divisor': dp, 'wrapper_content': weighted})
    c3 = (2, 4, 1)
    adaptive = birth(birth(c3, 0), 3)
    ap = unmarked(literal(adaptive))
    c3p = unmarked(literal(c3))
    c3_adaptive = birth(c3, 0)
    need(H(c3_adaptive) == 5, 'adaptive birth escapes C3 fixed divisor')
    fixed_controls = []
    for T in (c3, adaptive, (10, 4, 1, 54, 7, 23)):
        p = unmarked(literal(T))
        need(fixed_divisor(p) == 3, 'three inherited actual fixed-divisor hostiles')
        q = unmarked_G(p)
        need(q[1] % 3 == 0 and fixed_divisor(q) % 3 == 0,
             'gluing does not reset fixed slot divisor')
        fixed_controls.append({'rows': T, 'H': H(T), 'F2': slot(p, 2),
                               'fixed_divisor': 3, 'wrapped_H': q[1]})
    context_controls = []
    for n in (2, 3):
        for T in labelled(n):
            U = triangle(T, transitive(2), transitive(1))
            d = fixed_divisor(unmarked(literal(T)))
            dp = fixed_divisor(unmarked(literal(U)))
            need(dp % d == 0, 'uniform context divisibility')
            context_controls.append({'rows': T, 'context_H': H(U), 'input_divisor': d, 'output_divisor': dp})
    # Two seventeen-query twins inherited from the old observer's ordinal-sink law.
    twins = [{'order': n, 'H': H(transitive(n)), 'wrapped_H': H(G(transitive(n)))} for n in (2, 3)]
    need([x['wrapped_H'] for x in twins] == [5, 9], 'G separates old positive-birth twins')
    formal_f = lambda r: r**7 + 4*r**5 + 4*r**3 + 18*r
    formal_d = gcd(*(formal_f(r) for r in range(1, 8)))
    formal_e = gcd(*(factorial(r-1)*formal_f(r) for r in range(1, 8)))
    need((formal_d, formal_e) == (9, 27), 'fixed and wrapper contents can differ on formal polynomials')
    serialized = []
    ix = {s: i for i, s in enumerate(STATES)}
    for i, blocks in enumerate(ops):
        for derivative, block in enumerate(blocks):
            for (s, power), value in sorted(block.items()):
                serialized.append([i, ix[s], derivative, power, value])
    print(json.dumps({'status': 'PASS', 'scope': 'all-order two-port cyclic-gluing operator; finite literal controls, no spectrum coverage claim',
                      'skeletons': STATES, 'operator_entries': serialized,
                      'multiplicity_and_projection_controls': checks, 'literal_sources': len(rows),
                      'modular_ca_cell_controls': ca_checks, 'moduli': [2, 3, 7, 49],
                      'prime_escape_controls': escape_checks,
                      'mod49_jet_wrapper_controls': jet_checks,
                      'mod49_jet_literal_ordinal_controls': jet_product_checks,
                      'wrapper_history_rows': [dict(sorted(row.items())) for row in history_rows],
                      'formal_not_claimed_realizable_valuation_hostile': {'power_coefficients': [0,18,0,4,0,4,0,1], 'fixed_divisor': formal_d, 'wrapper_content': formal_e},
                      'literal_rows': rows, 'old_observer_twins': twins,
                      'fixed_divisor_control': {'C3_F2': slot(c3p, 2), 'C3_wrapped_H': H(G(c3)), 'C3_adaptive_birth_H': H(c3_adaptive)},
                      'fixed_divisor_hostiles': fixed_controls, 'uniform_context_controls': context_controls,
                      'adaptive_probe': {'rows': adaptive, 'H': H(adaptive), 'F2': slot(ap, 2)}}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
