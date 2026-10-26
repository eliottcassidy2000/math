#!/usr/bin/env python3
"""Exact clock relations versus literal integer traces; standard library only."""
from itertools import product
from math import gcd, isqrt, prod
import json

GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def step(n):
    return n//2 if n % 2 == 0 else (3*n+1)//2


def arrows(M, r):
    return tuple((e, s) for e in (0, 1) for s in range(M)
                 if (2*s-3**e*r-e) % M == 0)


def trace(n, M, L):
    out = [n % M]
    for unused in range(L):
        e = n % 2
        n = step(n)
        out.append((e, n % M))
    return tuple(out)


def paths(M, r, L):
    out = {(r,)}
    for unused in range(L):
        out = {p+((e, s),) for p in out for e, s in arrows(M, p[-1][1] if len(p)>1 else r)}
    return out


def lift(path, M):
    a, Q = path[0], M
    for j in range(1, len(path)):
        candidates = [b for b in (a, a+Q) if trace(b, M, j) == path[:j+1]]
        need(len(candidates) == 1, 'one-bit lift')
        a = candidates[0]
        Q *= 2
    return a, Q


def affine(w):
    a, b, den = 1, 0, 1
    for e in w:
        a, b, den = 3**e*a, 3**e*b+e*den, 2*den
    return a, b, den


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def E(n, B):
    return all((n % p or n == p) and ((n+2) % p or n+2 == p)
               for p in range(2, B+1) if prime(p))


def main():
    census = []
    for M in range(1, 19):
        for r in range(M):
            ar = arrows(M, r)
            need(len(ar) == 2, 'two outgoing labelled arrows')
            expected = {(n % 2, step(n) % M) for n in (r, r+M)}
            need(set(ar) == expected, 'direct one-step existential relation')
            if M % 2:
                need(len({s for e, s in ar}) == (1 if r == (M-1)//2 else 2), 'odd parallel-arrow boundary')
            else:
                need({e for e, s in ar} == {r % 2} and len({s for e, s in ar}) == 2, 'even parity and successor split')
        for L in range(8):
            Q = M*2**L
            source = {}
            for n in range(Q):
                p = trace(n, M, L)
                need(p not in source, 'literal trace source injectivity')
                source[p] = n
            graph = set().union(*(paths(M, r, L) for r in range(M)))
            need(set(source) == graph and len(graph) == Q, 'all paths equal direct trace image')
            for p, a in source.items():
                aa, qq = lift(p, M)
                need((aa, qq) == (a, Q), 'constructive lift versus full source enumeration')
                need(trace(a+Q, M, L) == p, 'positive-source representative')
                w = tuple(e for e, s in p[1:])
                A, b, den = affine(w)
                n = a
                for unused in range(L):
                    n = step(n)
                need(A*a+b == den*n, 'affine evaluator on exact cylinder')
            census.append([M, L, Q])
    covering_pairs = 0
    noncovering_pairs = 0
    for N in range(1, 65):
        for M in range(1, N+1):
            if N % M:
                continue
            all_unique = True
            for b in range(N):
                images = [(e, s % M) for e, s in arrows(N, b)]
                target = arrows(M, b % M)
                counts = [images.count(edge) for edge in target]
                need(all(edge in target for edge in images), 'reduction is graph morphism')
                unique = counts == [1, 1]
                need(unique == ((N//M) % 2 == 1), 'cover iff ratio odd at each fine source')
                all_unique &= unique
            covering_pairs += all_unique
            noncovering_pairs += not all_unique
    for M in range(1, 13):
        for L in range(1, 7):
            Q = M*2**L
            for n in range(Q):
                need((step(n+Q)-step(n)) % (Q//2) == 0, 'finite-horizon deterministic map')
    for L in range(1, 8):
        for w in product((0, 1), repeat=L):
            for cut in range(L+1):
                A, b, Q = affine(w)
                Au, bu, Qu = affine(w[:cut])
                Av, bv, Qv = affine(w[cut:])
                need((A, b, Q) == (Av*Au, Av*bu+Qu*bv, Qu*Qv), 'ordered affine composition')
    for n in range(1, 201):
        short, raw, elapsed = n, n, 0
        for L in range(1, 13):
            e = short % 2
            short = step(short)
            for unused in range(1+e):
                raw = raw//2 if raw % 2 == 0 else 3*raw+1
                elapsed += 1
            need(raw == short, 'actual ordinary-time transport')
    need(trace(3, 7, 1) == (3, (1, 5)) and trace(10, 7, 1) == (3, (0, 5)), 'same endpoints distinct native durations')
    need(set((e, s) for e, s in arrows(7, 3)) == {(0, 5), (1, 5)}, 'parallel arrows retain labels')
    # Generic existential quotient composition is not automatically sound.
    f = {'a':'b', 'b':'b', 'bp':'c', 'c':'c'}
    q = {'a':'A', 'b':'B', 'bp':'B', 'c':'C'}
    relation = {(q[x], q[f[x]]) for x in f}
    need(('A','B') in relation and ('B','C') in relation and (q['a'], q[f[f['a']]]) == ('A','B'), 'generic phantom composition hostile')
    for M in range(1, 14):
        for L in range(1, 9):
            Q = M*2**L
            zeros = (0,)+((0, 0),)*L
            ones = (M-1,)+((1, M-1),)*L
            need(lift(zeros, M) == (0, Q) and lift(ones, M) == (Q-1, Q), 'signed boundary cylinders')
            need(trace(Q, M, L) == zeros and trace(Q-1, M, L) == ones, 'positive finite boundary realizations')
            need(step(0) == 0 and step(-1) == -1, 'actual nonpositive fixed sources')
    for M in range(1, 9):
        for L in range(5):
            Q = M*2**L
            for a in range(Q):
                for H in (0, 1, 5, Q-1, Q, Q+3, 2*Q+1):
                    count = (H-a)//Q-(-a)//Q
                    need(count == sum(n % Q == a for n in range(1,H+1)), 'positive-height exact floor count')
                    if Q >= H:
                        need(count <= 1, 'height uniqueness')
    need(E(3, 3) and not E(9, 3) and 3 % 6 == 9 % 6, 'prime equality does not descend')
    sieve_counts = []
    for B, H in ((3, 22), (5, 46), (7, 118)):
        M = prod(p for p in range(2, B+1) if prime(p))
        L, Q = 3, 8*M
        seen_twins = []
        for a in range(Q):
            p = trace(a, M, L)
            guard = int(3 <= a <= H and E(a, B))
            direct = [n for n in range(3, H+1) if trace(n, M, L) == p and prime(n) and prime(n+2)]
            need(guard == len(direct), 'composed complete-height twin-start path certificate')
            if guard:
                seen_twins.append(a)
        need(seen_twins == [n for n in range(3, H+1) if prime(n) and prime(n+2)], 'all twin starts accounted')
        sieve_counts.append({'B':B,'H':H,'M':M,'L':L,'twins':seen_twins})
    for M in range(1, 101):
        m, n, determinant = 2, 1, 1
        for letter in 'A'+'D'*(M-1)+'UD':
            if letter == 'A':m,n,determinant=2*m+n,m,-determinant
            elif letter == 'U':m,n=2*m-n,m
            else:m,n=m+2*n,n
        need((m,n,determinant) == (16*M+2,4*M+1,-1), 'recovered every-modulus branch family')
        need((m-2) % M == (n-1) % M == 0 and gcd(m,n)==1 and (m+n)%2 == 1, 'root residues and Euclid guard')
        if M % 2 and M > 1:need(1 % M != -1 % M, 'full matrix determinant retains sign')
    print(json.dumps({'status':'FINITE-EXACT controls; all-parameter proofs in README.md',
                      'gates':GATES,'clock_moduli':[1,18],'path_lengths':[0,7],
                      'source_trace_classes':sum(row[2] for row in census),
                      'refinement_N_max':64,'covering_divisor_pairs':covering_pairs,
                      'noncovering_divisor_pairs':noncovering_pairs,
                      'height_sieve_compositions':sieve_counts,
                      'scope':'Finite paths and guarded finite-height realization; no infinite positive-source theorem.'},sort_keys=True,indent=2))


if __name__ == '__main__':
    main()
