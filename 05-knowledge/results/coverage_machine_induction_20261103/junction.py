#!/usr/bin/env python3
"""Independent small controls for native, cyclic and singular reset junctions."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
from math import gcd, isqrt
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def validate(T):
    n = len(T)
    need(all(0 <= x < 1 << n and not (x >> i & 1) for i, x in enumerate(T)), 'row bounds')
    need(all(((T[i] >> j) ^ (T[j] >> i)) & 1 for i in range(n) for j in range(i)), 'tournament')


@lru_cache(None)
def endpoints(T):
    n = len(T)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v, count in enumerate(dp[S]):
            rest = T[v] & ~S
            while rest and count:
                bit = rest & -rest
                rest -= bit
                dp[S | bit][bit.bit_length() - 1] += count
    return tuple(map(tuple, dp))


def H(T):
    return sum(endpoints(T)[-1]) if T else 1


def sub(T, names):
    return tuple(sum(((T[v] >> w) & 1) << j for j, w in enumerate(names)) for v in names)


def profile(T, root):
    dp = endpoints(T)
    hs = [sum(x) for x in dp]
    hs[0] = 1
    full = (1 << len(T)) - 1
    F = sum(dp[S][root] * hs[full ^ S] for S in range(1, full + 1))
    reverse = tuple(full ^ row ^ (1 << i) for i, row in enumerate(T))
    s = endpoints(reverse)[-1][root]
    q = hs[full ^ (1 << root)]
    return H(T), s, q, F


def F2(T):
    hs = [sum(x) for x in endpoints(T)]
    hs[0] = 1
    full = (1 << len(T)) - 1
    return sum(hs[S] * hs[full ^ S] for S in range(full + 1))


def slot_values(T):
    hs = [sum(x) for x in endpoints(T)]
    hs[0] = 1
    size = 1 << len(T)
    state = [1] + [0] * (size - 1)
    result = []
    for r in range(1, len(T) + 1):
        nxt = [0] * size
        for S in range(size):
            A = S
            while True:
                nxt[S] += hs[A] * state[S ^ A]
                if not A:
                    break
                A = (A - 1) & S
        state = nxt
        result.append(state[-1])
    return result


def ear(T, incoming):
    n = len(T)
    return tuple(row | (((incoming >> i) & 1) << n) for i, row in enumerate(T)) + (((1 << n) - 1) ^ incoming,)


def G(T):
    n = len(T)
    return ear(ear(T, (1 << n) - 1), 1 << n)


def native(T, P):
    need(sorted(P) == list(range(len(T))) and all(T[a] >> b & 1 for a, b in zip(P, P[1:])), 'native legal spine')
    retained = {frozenset((a, b)) for a, b in zip(P, P[1:])}
    U = list(T)
    for a, b in combinations(range(len(T)), 2):
        if frozenset((a, b)) not in retained:
            U[a] ^= 1 << b
            U[b] ^= 1 << a
    validate(tuple(U))
    return tuple(U)


def paths(T):
    return [P for P in permutations(range(len(T))) if all(T[a] >> b & 1 for a, b in zip(P, P[1:]))]


def auts(T):
    n = len(T)
    return {p for p in permutations(range(n)) if all(((T[i] >> j) & 1) == ((T[p[i]] >> p[j]) & 1) for i in range(n) for j in range(i))}


def strong(T):
    full = (1 << len(T)) - 1
    for v in range(len(T)):
        reach = 1 << v
        while True:
            nxt = reach
            for i in range(len(T)):
                if reach >> i & 1:
                    nxt |= T[i]
            if nxt == reach:
                break
            reach = nxt
        if reach != full:
            return False
    return bool(T)


# The inherited static observer basis; values are independently recounted.
QUERIES = ((), (1,), (2,), (3,), (0,1), (1,1), (2,1), (3,1),
           (0,2), (1,2), (2,2), (3,2), (2,0,1), (3,0,1), (0,2,1), (2,2,1), (3,2,1))


def queries(T, ports):
    answer = []
    for word in QUERIES:
        U = T
        for mask in word:
            U = ear(U, sum(1 << ports[i] for i in range(2) if mask >> i & 1))
        answer.append(H(U))
    return answer


def R(z):
    h, s, q, F = z
    return h + F - s, q, q, 2 * F - q


def N(z):
    h, s, q, F = z
    return h + F - s, h - s, h, h + F


def delta(z):
    h, s, q, F = z
    return F - h - q + s


def sealed_formula(z, k):
    h, s, q, F = z
    a, d = h - s, F - q
    return a + 3 * q + (3 * 2**k - 1) * d, a + q + (2**k - 1) * d


def literal_seal(T, root, k, initial_ports=None):
    ports = list(initial_ports or (root, next(v for v in range(len(T)) if v != root)))
    need(root in ports and len(set(ports)) == 2, 'initial ports')
    steps = []
    for letter in 'R' * k + 'NR':
        oldroot = root
        child = len(T)
        root_slot = ports.index(root)
        T = ear(T, 1 << root)
        ports[1 - root_slot] = child
        if letter == 'N':
            root = child
        need(root in ports, 'named root retained')
        steps.append(dict(letter=letter, root=oldroot, newborn=child, ports=list(ports)))
    return T, root, tuple(ports), steps


def all_tournaments(n):
    edges = list(combinations(range(n), 2))
    for bits in range(1 << len(edges)):
        T = [0] * n
        for j, (a, b) in enumerate(edges):
            u, v = (a, b) if bits >> j & 1 else (b, a)
            T[u] |= 1 << v
        yield tuple(T)


def main():
    source = (2, 4, 1)
    wrapper = G(source)
    spine = (4, 0, 1, 2, 3)
    reset = native(wrapper, spine)
    need(wrapper == (10, 12, 9, 16, 7), 'wrapper rows')
    need(reset == (6, 20, 24, 3, 9), 'reset rows')
    need(H(wrapper) == H(reset) == 15, 'neutral native')
    need(strong(wrapper) and strong(reset), 'strong junction')
    need(len(auts(wrapper)) == 3 and len(auts(reset)) == 5, 'different automorphism groups')
    need(auts(wrapper) & auts(reset) == {tuple(range(5))}, 'intersection identity')
    X, Y = queries(wrapper, (3, 4)), queries(reset, (3, 4))
    need(gcd(*X) == 3 and gcd(*Y) == 1, 'two-port content reset')
    need(all(profile(reset, r) == (15, 3, 5, 21) for r in range(5)), 'regular root profile')
    slots_before, slots_after = slot_values(wrapper), slot_values(reset)
    need(slots_before == [r**5+8*r**3+6*r for r in range(1,6)], 'wrapper slot polynomial')
    need(slots_after == [r**5+10*r**3+4*r for r in range(1,6)], 'native slot polynomial')
    need(gcd(*slots_before)==3 and gcd(*slots_after)==15, 'uniform slot divisor worsens while named divisor improves')
    examples = []
    for k in (1, 2, 3):
        T, root, ports, steps = literal_seal(reset, 4, k, (3, 4))
        b, t = 48 * 2**k + 11, 16 * 2**k + 1
        need(profile(T, root) == (b, t, t, b), 'all-order formula control')
        need(gcd(b, t) == 1 and strong(T), 'primitive strong return')
        need(gcd(*queries(T, ports)) == 1, 'final pair primitive')
        if k <= 2:
            need(len(paths(T)) == b, 'independent permutation recount')
        examples.append(dict(k=k, rows=T, root=root, ports=ports, profile=(b,t,t,b), steps=steps))
    need(all(107 % d for d in range(2, isqrt(107) + 1)), '107 prime')
    need(203 == 7 * 29 and all(29 % d for d in range(2, isqrt(29) + 1)), '203 seven-prime')
    # A deliberately small source bank, not a new scalar-support scan.
    small = list(all_tournaments(2)) + list(all_tournaments(3)) + [(10,4,1,6), wrapper, reset]
    seals = wrapper_checks = 0
    for T in small:
        V = G(T)
        h, f = H(T), F2(T)
        need(profile(V, len(T)+1) == (h+f,h,h,h+f), 'balanced cyclic seal')
        wrapper_checks += 1
        for r in range(len(T)):
            z = profile(T, r)
            need(delta(z) >= 0 and delta(R(z)) == delta(z) and delta(N(z)) == 0, 'defect laws')
            for k in (1,2,3):
                U, root, pp, _ = literal_seal(T, r, k)
                b, t = sealed_formula(z, k)
                expected = (b,t,t,b)
                need(profile(U, root) == expected, 'literal general seal')
                need(gcd(b,t) == gcd(delta(z), 2**k*(z[3]-z[2])+z[2]), 'exact primitive clock')
                seals += 1
    # Complete native/common-automorphism control through order four.
    common_controls = 0
    for n in (2,3,4):
        for T in all_tournaments(n):
            A = auts(T)
            for P in paths(T):
                U = native(T,P)
                need(A & auts(U) == {tuple(range(n))}, 'common automorphism theorem')
                common_controls += 1
    B = (42,4,1,22,7,94,31)
    P = (4,1,2,0,5,6,3)
    V = native(B,P)
    need(H(B) == H(V) == 31 and (F2(B),F2(V)) == (644,596), 'H31 neutral cover tuning')
    need((H(G(B)),H(G(V))) == (675,627), 'gluing sees native source change')
    def ratio(b,t):
        return t * pow(b,-1,11) % 11 if b % 11 else 'infinity'
    def C3(b,t):
        return 61*b-2*t,9*b
    def C5(b,t):
        return 421*b-2*t,33*b
    entries=[]
    for k in range(1,11):
        b,t=48*2**k+11,16*2**k+1
        b,t=C3(b,t)
        initial=(b,t)
        word=['C3']
        need(gcd(b,t)==1 and t%3==0 and b%3!=0,'arithmetic entry guard')
        if ratio(b,t)==7:
            b,t=C3(b,t)
            word.append('C3')
        need((3*b-2*t)%11!=0,'C5 exact entry guard')
        b,t=C5(b,t)
        word.append('C5')
        need(gcd(b,t)==1 and t%3==0 and b%3!=0 and ratio(b,t)==0,'actual controller state entry')
        entries.append(dict(k=k,after_first_C3=initial,initial_ratio=ratio(*initial),word=word,entry=(b,t)))
    need(entries[1]['initial_ratio']==10,'arithmetic guard alone does not enter eight-state controller')
    print(json.dumps(dict(status='FINITE-EXACT controls; proofs in RESET_JUNCTION.md',
        universe=dict(wrapper_sources=len(small), rooted_seal_controls=seals,
                      common_automorphism_controls=common_controls, direct_permutation_orders=[8,9]),
        junction=dict(source=source,wrapper=wrapper,spine=spine,reset=reset,ports=[3,4],
                      observer_before=X,observer_after=Y,automorphism_orders=[3,5],
                      slot_values_before=slots_before,slot_values_after=slots_after,slot_divisors=[3,15]),
        examples=examples, controller_entries=entries,
        cover_tuning=dict(source=B,path=P,target=V,h=31,F2=[644,596],glued_counts=[675,627])),
        indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
