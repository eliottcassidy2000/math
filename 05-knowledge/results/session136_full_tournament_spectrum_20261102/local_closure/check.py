#!/usr/bin/env python3
"""Standalone positive-word local constructor and hostile controls (stdlib only)."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path


def need(x, message):
    if not x:
        raise RuntimeError(message)

# Compressed actual words; all repeats have nonnegative integer exponents.
def seq(*nodes):
    return ('seq', tuple(nodes))


def repeat(node, k):
    need(isinstance(k, int) and k >= 0, 'nonpositive-word exponent')
    return ('rep', node, k)


def length(node):
    if node == 'A':
        return 1
    if node == 'B':
        return 2
    if node[0] == 'rep':
        return node[2] * length(node[1])
    return sum(map(length, node[1]))


def expand(node):
    if node == 'A':
        return '0'
    if node == 'B':
        return '10'
    if node[0] == 'rep':
        return expand(node[1]) * node[2]
    return ''.join(map(expand, node[1]))

# Affine (a,b,g,h) means (d,t) -> (d+a*t+b,g*t+h).
IDENTITY = (0, 0, 1, 0)


def compose(f, q, modulus):
    a, b, g, h = f
    A, B, G, H = q
    return tuple(v % modulus for v in
                 (A+a*G, B+a*H+b, g*G, g*H+h))


def power(f, exponent, modulus):
    ans = tuple(v % modulus for v in IDENTITY)
    while exponent:
        if exponent & 1:
            ans = compose(f, ans, modulus)
        f = compose(f, f, modulus)
        exponent >>= 1
    return ans


def affine(node, modulus):
    if node == 'A':
        return tuple(v % modulus for v in (0, 0, 2, 0))
    if node == 'B':
        return tuple(v % modulus for v in (-1, 0, 4, 1))
    if node[0] == 'rep':
        return power(affine(node[1], modulus), node[2], modulus)
    answer = tuple(v % modulus for v in IDENTITY)
    for child in node[1]:
        answer = compose(affine(child, modulus), answer, modulus)
    return answer


def action(node, state, modulus):
    d, t = state
    a, b, g, h = affine(node, modulus)
    return ((d+a*t+b) % modulus, (g*t+h) % modulus)


def order2(n):
    if n == 1:
        return 1
    value = 2 % n
    e = 1
    while value != 1:
        value = 2*value % n
        e += 1
        need(e <= n, 'odd-modulus order failed')
    return e


def odd_macros(n, L, omit_reset=False, reverse_commutator=False):
    e = order2(n)
    K = n*e
    P = repeat('A', L)
    X = 'B' if omit_reset else seq(P, 'B', P)
    Ai = repeat('A', e-1)
    Xi = repeat(X, K-1)
    C = seq(X, repeat('A', (-2) % e))
    D = seq(Ai, X, Ai)
    Ci = repeat(C, n-1)
    Di = repeat(D, n-1)
    E = seq(C, D)
    J = seq(Ci, Di, C, D)
    if reverse_commutator:
        J = seq(Di, Ci, D, C)
    return E, J, e, K


def construct(M, z, omit_reset=False, reverse_commutator=False):
    need(M >= 1, 'bad modulus')
    z %= M
    n, a = M, 0
    while n % 2 == 0:
        n //= 2
        a += 1
    D = 1 << a
    e = order2(n)
    L = e * max(1, (a+e-1)//e)
    r0 = max(1, (a+1)//2)
    rr = ((9*z-1) * pow(3, -1, D)) % D
    r = r0 + (rr-r0) % D
    initial = seq(repeat('B', r), repeat('A', L))
    j = k = 0
    if n > 1:
        E, J, _, K = odd_macros(n, L, omit_reset, reverse_commutator)
        d, t = action(initial, (0, 0), n)
        j = (z-d-t) % 3 if n % 3 == 0 else 0
        prefix = seq(initial, repeat(E, j))
        d, t = action(prefix, (0, 0), n)
        g = math.gcd(3, n)
        delta = 3*pow(8, -1, n) % n
        need((z-d-t) % g == 0, '3-primary targeting failed')
        if n//g > 1:
            k = (((z-d-t) % n)//g * pow(delta//g, -1, n//g)) % (n//g)
        result = seq(prefix, repeat(J, k))
    else:
        K = 1
        result = initial
    w = 2*L+2
    bound = 2*r+L+(2+n*(n-1))*(2*w+3*e)
    need(length(result) <= bound, 'claimed word-length bound failed')
    return result, dict(M=M, target=z, a=a, n=n, r=r, L=L, e=e,
                        j=j, k=k, length=length(result), bound=bound)


def word_state(word):
    c = t = ell = 0
    for bit in word:
        if bit == '0':
            c, t, ell = c+t+ell, 2*t+ell, 0
        else:
            t, ell = 2*t+ell, 1
    return c, t, ell


def adjacency(word):
    # Apex is vertex 0, chain is 1,...,len(word).
    n = len(word)+1
    out = [0]*n
    for i in range(1, n):
        for j in range(i+1, n):
            out[i] |= 1 << j
        if word[i-1] == '1':
            out[0] |= 1 << i
        else:
            out[i] |= 1
    return out


def paths(out):
    n = len(out)
    size = 1 << n
    dp = [[0]*n for _ in range(size)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, size):
        rest = size-1-mask
        for i, count in enumerate(dp[mask]):
            if count:
                targets = out[i] & rest
                while targets:
                    bit = targets & -targets
                    targets -= bit
                    dp[mask | bit][bit.bit_length()-1] += count
    return sum(dp[-1])


def strong(out):
    if len(out) == 1:
        return True
    full = (1 << len(out))-1
    for start in range(len(out)):
        visited = 1 << start
        todo = visited
        while todo:
            bit = todo & -todo
            todo -= bit
            new = out[bit.bit_length()-1] & ~visited
            visited |= new
            todo |= new
        if visited != full:
            return False
    return True


def f(n):
    if n < 2:
        return 1
    if n % 2:
        return f(n//2)
    h = n//2
    return f(h)+h+h % 2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bound', type=int, default=160)
    args = parser.parse_args()
    receipt = {'status': 'PASS', 'c_residue_modulus_bound': args.bound}
    digest = hashlib.sha256()
    count = expanded = dp_constructor = maximum_length = 0
    short_examples = []
    for M in range(1, args.bound+1):
        for z in range(M):
            node, data = construct(M, z)
            d, t = action(node, (0, 0), M)
            need((d+t) % M == z, f'constructor failed M={M},z={z}')
            need(t % (1 << data['a']) == 0, 'dyadic invariant lost')
            need(action(seq(node, repeat('A', data['L'])), (0, 0), M) == (d,t),
                 'unbounded-fiber loop changed residue')
            digest.update(json.dumps(data, sort_keys=True).encode())
            count += 1
            maximum_length = max(maximum_length, data['length'])
            if data['length'] <= 1000:
                word = expand(node)
                need(word.startswith('1') and word.endswith('0'), 'not strong word')
                c0, t0, ell = word_state(word)
                need((c0-t0) % M == d and t0 % M == t and ell == 0,
                     'expanded/compressed mismatch')
                expanded += 1
                if len(word) <= 10:
                    need(paths(adjacency(word)) == 1+2*c0, 'constructor graph mismatch')
                    dp_constructor += 1
                    if len(short_examples) < 12:
                        short_examples.append(dict(M=M,z=z,word=word,H=1+2*c0))
    receipt.update(constructor_residues_checked=count,
                   expanded_words_checked=expanded,
                   direct_constructor_graphs=dp_constructor,
                   largest_constructor_bit_length=maximum_length,
                   constructor_digest=digest.hexdigest(), short_examples=short_examples)
    # Deliberately mixed 2/3/high odd-primary moduli; no separate CRT tests.
    large = [2**a * 3**b * q for a,b,q in [(6,3,5),(8,2,7),(1,4,11),(0,6,1),
                                                (10,0,1),(3,0,101),(4,1,127)]]
    selected = 0
    for M in large:
        for z in [0,1,2,3,6,10,12,M//2,M-2,M-1]:
            node,data = construct(M,z)
            d,t=action(node,(0,0),M)
            need((d+t)%M==z%M, 'selected mixed modulus failed')
            selected += 1
    receipt.update(selected_mixed_moduli=large, selected_mixed_targets=selected)
    # Independent graph-level universe, not using affine A/B evaluation.
    graph_count = strong_count = 0
    for s in range(10):
        for bits in itertools.product('01', repeat=s):
            word=''.join(bits)
            out=adjacency(word)
            c,t,ell=word_state(word)
            need(paths(out)==1+2*c, 'binary graph identity failed')
            want = (s == 0) or (word[0]=='1' and word[-1]=='0')
            need(strong(out)==want, 'strongness criterion failed')
            graph_count += 1
            strong_count += want
    receipt.update(binary_graph_word_length_bound=9,
                   binary_graphs_checked=graph_count,strong_binary_graphs=strong_count)
    affine_cases=0
    for n in range(3,64,2):
        e=order2(n)
        E,J,_,K=odd_macros(n,e)
        need(affine(repeat('B',K),n)==tuple(v%n for v in IDENTITY),
             'odd B order bound failed')
        for d in range(n):
            for t in range(n):
                need(action(J,(d,t),n)==((d+3*pow(8,-1,n))%n,t),
                     'commutator formula failed')
                if n%3==0:
                    D,T=action(E,(d,t),n)
                    need((D+T-d-t)%3==1, 'ternary correction failed')
                affine_cases+=1
    receipt['odd_affine_states_checked']=affine_cases
    candidates=[f(u)+u+u%2 for u in range(6,13)]
    need(13 not in candidates, 'inherited H13 hostile failed')
    out13=[22,12,24,17,2]
    need(paths(out13)==13, 'unrestricted H13 witness failed')
    receipt['one_apex_13_inverse_candidates']=candidates
    receipt['unrestricted_13_adjacency']=out13
    mutations=[]
    for key in ['omit_reset','reverse_commutator']:
        witness=None
        for M in range(2,25):
            for z in range(M):
                node,data=construct(M,z,**{key:True})
                d,t=action(node,(0,0),M)
                if (d+t)%M != z:
                    witness=dict(mutation=key,M=M,target=z,observed=(d+t)%M)
                    break
            if witness:
                break
        need(witness is not None,'mutation escaped hostile universe')
        mutations.append(witness)
    receipt['rejected_mutations']=mutations
    print(json.dumps(receipt,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
