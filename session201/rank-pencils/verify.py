"""Exact source-pencil verifier. Standard library, no optimization-sensitive asserts."""
import argparse
import hashlib
import heapq
import itertools
import json
import math
from pathlib import Path

T = (1, 0, 2, 1)
B = (3, 1, 3, 1)
Z0 = (1, 0, 2, 0)
Q = ((0, 0, 0, 1), (0, 1, 0, -1), (0, 0, 1, -1), (1, -1, -1, 1))
QPOS = ((1, 1, 2, 2), (1, 1, 1, 1), (2, 1, 4, 3), (2, 1, 3, 3))
PROBES = (T, (1, 0, 4, 1), B, (5, 2, 6, 3))


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def step(s, c):
    D, E, U, d = s
    if c == 'O':
        return D, 0, 2 * U, d
    if c == 'Z':
        return D - E + U, d, 2 * U - d, d
    if c == 'J':
        return D - E + U, D - E, D + U, D
    raise ValueError(c)


def word(s, letters):
    for c in letters:
        s = step(s, c)
    return s


def dagger(s):
    return ''.join({'O': 'Z', 'Z': 'O', 'J': 'J'}[c] for c in s[::-1])


def K(v, w):
    D, E, U, d = v
    A, B_, C, e = w
    return D*e + d*A + E*B_ + U*C - E*e - d*B_ - U*e - d*C + d*e


def positive(v):
    D, E, U, d = v
    return D-E-d, U-D+E-d, E, d-E


def raw(y):
    a, f, e, p = y
    return a+2*e+p, e, a+f+2*e+2*p, e+p


def minors(w, q):
    return tuple(w[i]*q[j]-w[j]*q[i] for i in range(4) for j in range(i+1, 4))


def collapsed(w, q):
    return not any(minors(w, q))


def probe_data(w, q):
    ms = [K(v, w) for v in PROBES]
    hs = [K(v, q) for v in PROBES]
    ds = [ms[0]*hs[i]-ms[i]*hs[0] for i in range(1, 4)]
    return ms, hs, ds


def det(A):
    n = len(A)
    total = 0
    for p in itertools.permutations(range(n)):
        sign = (-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        total += sign*math.prod(A[i][p[i]] for i in range(n))
    return total


def bank(cap, decorated=False, prune=False):
    seed = (B, Z0) if decorated else T
    height = (lambda x: x[0][2]) if decorated else (lambda x: x[2])
    if height(seed) > cap:
        return {}
    out = {seed: (0, '')}
    heap = [(height(seed), seed)]
    while heap:
        _, s = heapq.heappop(heap)
        cost, text = out[s]
        for c in 'OZJ':
            q = (step(s[0], c), step(s[1], c)) if decorated else step(s, c)
            need(height(q) > height(s), 'height must strictly increase')
            if height(q) > cap or (prune and collapsed(*q)):
                continue
            candidate = cost+(c == 'J'), text+c
            if q not in out:
                out[q] = candidate
                heapq.heappush(heap, (height(q), q))
            elif candidate < out[q]:
                out[q] = candidate
    return out


def sieve(n):
    p = [True]*(n+1)
    p[0:2] = [False, False]
    for i in range(2, math.isqrt(n)+1):
        if p[i]:
            p[i*i:n+1:i] = [False]*((n-i*i)//i+1)
    return p


def image(S, R, N, primes):
    full, balanced = {}, {}
    tests = rank_one_pairs = 0
    for v, (cv, _) in S.items():
        for (w, q), (cw, _) in R.items():
            m, h = K(v, w), K(v, q)
            need(3 <= h < m < 2*h and m % 2 and h % 2, 'physical joint cone')
            if m > N:
                continue
            tests += 1
            dep = collapsed(w, q)
            rank_one_pairs += dep
            cost = cv+cw
            if primes[m]:
                need(not dep, 'rank-one balanced prime')
                balanced[m] = min(cost, balanced.get(m, N))
                full[m] = min(cost, full.get(m, N))
            k = 1
            while True:
                n = m+2*(2**k-1)*h
                if n > N:
                    break
                if primes[n]:
                    need(not dep, 'rank-one Mersenne prime')
                    full[n] = min(cost, full.get(n, N))
                k += 1
    return full, balanced, tests, rank_one_pairs


def graph(letters):
    rows, root = [2, 0], 0
    for c in letters:
        n = len(rows)
        rows.append(0)
        for v in range(n):
            if c != 'O' and v == root:
                rows[n] |= 1 << v
            else:
                rows[v] |= 1 << n
        if c == 'J':
            root = n
    return rows


def count_graph(rows):
    n = len(rows)
    full = (1 << n)-1

    def count(starts, cycles):
        dp = [{} for _ in range(full+1)]
        for v in starts:
            dp[1 << v][v] = 1
        for mask in range(1, full+1):
            for v, c in dp[mask].items():
                todo = rows[v] & (full ^ mask)
                while todo:
                    bit = todo & -todo
                    todo -= bit
                    nxt = bit.bit_length()-1
                    row = dp[mask | bit]
                    row[nxt] = row.get(nxt, 0)+c
        return sum(c for v, c in dp[full].items() if not cycles or rows[v] & 1)

    return count(range(n), False), count([0], True)


def main(N):
    need(N >= 171, 'cap must include literal controls')
    need(det(PROBES) == -4 and det(Q) == -1, 'basis determinant')
    PQ = tuple(tuple(sum(v[k]*Q[k][j] for k in range(4)) for j in range(4)) for v in PROBES)
    need(det(PQ) == 4, 'measurement determinant')
    need(PQ == ((1,-1,1,0),(1,-1,3,-2),(1,0,2,0),(3,-1,3,0)), 'response matrix')
    for i in range(4):
        for j in range(4):
            a = tuple(int(i == k) for k in range(4))
            b = tuple(int(j == k) for k in range(4))
            need(K(raw(a), raw(b)) == QPOS[i][j], 'positive bilinear matrix')
            need(QPOS[i][j] >= (3 if i >= 2 and j >= 2 else 0), 'cone lower bound')

    cap = math.isqrt(4*(N-1))
    S, R, Rp = bank(cap), bank(cap, True), bank(cap, True, True)
    need(Rp == {r: data for r, data in R.items() if not collapsed(*r)}, 'exact prefix-pruned inventory and costs')
    rank_one = 0
    factors = []
    for (w, q), (cost, letters) in R.items():
        need(w == word(B, letters) and q == word(Z0, letters), 'common-word certificate')
        ms, hs, ds = probe_data(w, q)
        need(collapsed(w, q) == (not any(ds)), 'rank equivalent to physical probes')
        if not collapsed(w, q):
            need(any(abs(d)*6*w[2]**2 >= ms[0]*ms[i+1] for i, d in enumerate(ds) if d), 'quantitative separation')
            need(all(d % 2 == 0 for d in ds), 'even probe determinants')
            continue
        rank_one += 1
        g = math.gcd(*w)
        p = tuple(x//g for x in w)
        a = math.gcd(*q)
        need(q == tuple(a*x for x in p), 'integral primitive factorization')
        need(g % 2 == 1 and g >= 3 and g < 2*a < 2*g, 'strict factor bounds')
        need(min(positive(p)) >= 0 and p[3] >= 1, 'primitive cone only')
        for c in 'OZJ':
            need(collapsed(step(w, c), step(q, c)), 'absorbing rank one')
        for v in PROBES:
            z = K(v, p)
            need(z >= 3*v[3]*p[3] >= 3, 'positive integer factor')
            for k in range(1, 6):
                L = 2*(2**k-1)
                need(K(v, w)+L*K(v, q) == (g+L*a)*z, 'target factor')
        if len(factors) < 4:
            factors.append(dict(word=letters, g=g, a=a, primitive=p))

    primes = sieve(N)
    direct = bank(N-1)
    dF, dG = {}, {}
    for s, (cost, _) in direct.items():
        F, G = s[0]-s[1]+s[2], s[0]+2*s[2]
        if 5 <= F <= N and primes[F]:
            dF[F] = min(cost, dF.get(F, N))
        if 5 <= G <= N and primes[G]:
            dG[G] = min(cost, dG.get(G, N))
    f, g, tests, discarded_pairs = image(S, R, N, primes)
    fp, gp, pruned_tests, unused = image(S, Rp, N, primes)
    need(f == fp == dF and g == gp == dG, 'complete prime image and minimum J costs')
    need(unused == 0, 'rank-one states survived guard')

    controls = []
    for alpha, beta, k, target, dep in [('', 'OJOOJ', 1, 171, True), ('', 'ZOJ', 2, 139, False), ('O', '', 2, 51, False)]:
        v, w, q = word(T, alpha), word(B, beta), word(Z0, beta)
        need(collapsed(w, q) == dep, 'literal rank control')
        m, h = K(v, w), K(v, q)
        need(m+2*(2**k-1)*h == target, 'literal paired target')
        letters = alpha+dagger(beta)+'O'+'Z'*k+'J'
        rows = graph(letters)
        need(count_graph(rows) == (target, 1), 'literal C1 graph count')
        controls.append(dict(alpha=alpha, beta=beta, k=k, m=m, h=h, target=target, rank=1 if dep else 2, graph_word=letters, vertices=len(rows), H=target, C=1))
    # Equal first boundary need not give the same rank or joint slope.
    x = (word(B, 'OJOJ'), word(Z0, 'OJOJ'))
    y = (word(B, 'ZZOJ'), word(Z0, 'ZZOJ'))
    need(x[0] == y[0] and collapsed(*x) and not collapsed(*y), 'same w different rank')
    return dict(status='PASS', target_cap=N, complete_forward_U_cap=N-1,
                direct_sources=len(direct), paired_height_cap=cap, paired_sources=len(S),
                synchronized_states=len(R), rank_one_states=rank_one, retained_states=len(Rp),
                complete_pair_tests=tests, pruned_pair_tests=pruned_tests,
                removed_rank_one_pairs=discarded_pairs, prime_C1_targets=len(f),
                prime_balanced_targets=len(g), prime_minimum_J_costs_preserved=True,
                rank_one_examples=factors, literal_controls=controls,
                dependencies={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((Path(__file__).parent/'dependencies').glob('*.md'))},
                verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scope='Finite complete height-bounded validation; all-order theorems require the written proof; no universal prime coverage.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cap', type=int, default=5001)
    args = parser.parse_args()
    print(json.dumps(main(args.cap), indent=2, sort_keys=True))
