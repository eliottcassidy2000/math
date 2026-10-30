#!/usr/bin/env python3
"""Exact controls; no orbit-to-ROOT oracle, floating-point rank, or assert."""
import argparse
from fractions import Fraction as F
from itertools import product, permutations
import json
import hashlib
from math import gcd
from pathlib import Path
from random import Random

gates = 0


def need(ok, msg):
    global gates
    gates += 1
    if not ok:
        raise RuntimeError(msg)


def val2(n):
    need(n != 0, "valuation domain")
    return (abs(n) & -abs(n)).bit_length() - 1


def affine(word):
    P, Q, B = 1, 1, 0
    for a in word:
        P, B, Q = 3 * P, 3 * B + Q, Q * 2**a
    return P, Q, B


def replay(n, word):
    for a in word:
        z = 3 * n + 1
        need(val2(z) == a, "exact supplied valuation")
        n = z // 2**a
    return n


def weight(n):
    return F(4, (n + 1) * (n + 3))


def mass_interval(r, M, count):
    """Exact partial sum and decreasing-function integral tail enclosure."""
    lo = sum((weight(r + M*t) for t in range(count)), F(0))
    a = r + M * count
    return lo, lo + weight(a) + F(4, M * (a + 1))


def shadows():
    words = [(1,), (1, 2), (1, 1, 2), (1, 1, 1, 2),
             (1, 1, 1, 2, 1, 1, 4)]
    # Every expanding positive word of length <= 5 on valuations {1,2,3}.
    words += [w for d in range(1, 6) for w in product(range(1, 4), repeat=d)
              if 3**d > 2**sum(w)]
    words = sorted(set(words), key=lambda w: (len(w), w))
    rows, instances = [], 0
    for word in words:
        P, Q, B = affine(word)
        gg = gcd(B, P-Q)
        b, D = B//gg, (P-Q)//gg
        need(D > 0 and b % 2 == D % 2 == 1, "negative odd anchor")
        phases = []
        for j in range(1, 9):
            M = 2 * Q**j
            r = (-b * pow(D, -1, M)) % M
            need(r % 2 == 1, "odd residue")
            h = (D*r + b)//M
            need(h*M == D*r+b, "exact source phase")
            e = F(2*h*P**j-b, D)
            need(e.denominator == 1, "endpoint integrality")
            for t in (0, 1, 2, 7):
                x = r + M*t
                y = replay(x, word*j)
                need(y == e + 2*P**j*t, "endpoint progression")
                need(y > x, "full block expands at positive source")
                instances += 1
            if M > b:
                need(1 <= h <= D, "bounded phase once modulus clears anchor")
                need(r == (h*M-b)//D, "canonical residue")
                need(e >= F(2*P**j-b, D), "independent endpoint floor")
                phases.append(h)
            if M >= 2*b:
                lo, hi = mass_interval(r, M, 24)
                upper = F(16*D*D + 8*D, M*M)
                need(hi <= upper, "certified original-source tail bound")
                need(lo >= F(4, (M+1)*(M+3)), "mass lower bound")
                L = 2*P**j
                alpha = F(h, D)
                for t0 in (0, 1, 3, 17):
                    _, tail_hi = mass_interval(r+M*t0, M, 16)
                    Y = L*(alpha+t0)
                    weak_bound = 16*(D+1)*F(L, M*M)
                    need(Y*tail_hi <= weak_bound, "transported weak-L1 endpoint tail")
            r2 = (-b * pow(D, -1, 2*Q**(j+1))) % (2*Q**(j+1))
            need(r2 % M == r, "source-preserving refinement")
        if word in [(1,), (1, 2), (1, 1, 2), (1, 1, 1, 2),
                    (1, 1, 1, 2, 1, 1, 4)]:
            period = 1
            while pow(Q, period, D) != 1 % D:
                period += 1
            rows.append(dict(word=word, P=P, Q=Q, anchor=f"-{b}/{D}",
                             phase_period=period, first_phases=phases))
    # An all-depth estimate may not keep the fixed-word constant while changing words.
    hostile = []
    for a in (3, 7, 15, 31):
        w = (1,)*a
        P, Q, B = affine(w)
        need(F(B, P-Q) == 1, "same anchor under word powers")
        hostile.append(dict(length=a, unreduced_D=P-Q, reduced_D=1))
    # Contracting positive anchor: root atom does NOT decay with precision.
    for j in range(1, 9):
        need(replay(1, (2,)*j) == 1, "positive-anchor hostile")
        need(weight(1) == F(1, 2), "nondecaying root mass")
    return dict(words=len(words), literal_instances=instances, anchors=rows,
                reduction_controls=hostile)


def paths(m):
    P = [r*m+c for r in range(m)
         for c in (range(m) if r % 2 == 0 else reversed(range(m)))]
    Q = [r*m+c for c in range(m)
         for r in (range(m) if c % 2 == 0 else reversed(range(m)))]
    return P, Q


def edges(P):
    return {tuple(sorted((x, y))) for x, y in zip(P, P[1:])}


def path_matrix(P):
    A = [[0]*len(P) for _ in P]
    for x, y in zip(P, P[1:]):
        A[x][y], A[y][x] = 1, -1
    return A


def native(K, P):
    need(all(K[x][y] == 1 for x, y in zip(P, P[1:])), "legal native path")
    A = path_matrix(P)
    return [[-x+2*a for x, a in zip(row, ar)] for row, ar in zip(K, A)]


def rank_mod(A, p):
    A = [[x % p for x in row] for row in A]
    r = 0
    for c in range(len(A[0]) if A else 0):
        piv = next((i for i in range(r, len(A)) if A[i][c]), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        z = pow(A[r][c], -1, p)
        A[r] = [(x*z) % p for x in A[r]]
        for i in range(r+1, len(A)):
            if A[i][c]:
                z = A[i][c]
                A[i] = [(x-z*y) % p for x, y in zip(A[i], A[r])]
        r += 1
    return r


def prime_tournament(K):
    """Closure of each pair under all distinguishing outside vertices."""
    n = len(K)
    full = (1 << n)-1
    rows = [sum(1 << j for j in range(n) if row[j] == 1) for row in K]
    for a in range(n):
        for b in range(a+1,n):
            S, pending = (1 << a) | (1 << b), [b]
            while pending and S != full:
                v = pending.pop()
                new = (rows[v] ^ rows[a]) & (full ^ S)
                S |= new
                while new:
                    bit = new & -new
                    pending.append(bit.bit_length()-1)
                    new -= bit
            if S != full:
                return False
    return True


def prime_examples():
    samples = []
    for m in (9, 11):
        n = m*m
        P,Q = paths(m)
        EP = edges(P)
        rng = Random(200+n)
        for attempt in range(100):
            K = [[0]*n for _ in range(n)]
            for x in range(n):
                for y in range(x+1,n):
                    K[x][y] = rng.choice((-1,1))
                    K[y][x] = -K[x][y]
            for x,y in zip(P,P[1:]):
                K[x][y],K[y][x] = 1,-1
            for x,y in zip(Q,Q[1:]):
                a = 1 if tuple(sorted((x,y))) in EP else -1
                K[x][y],K[y][x] = a,-a
            K1 = native(K,P)
            K2 = native(K1,Q)
            if all(prime_tournament(G) for G in (K,K1,K2)):
                break
        else:
            raise RuntimeError('sample prime search exhausted')
        # Gershgorin applied to K^T K is an exact independent norm certificate.
        gram_rows = [sum(abs(sum(K[k][i]*K[k][j] for k in range(n))) for j in range(n)) for i in range(n)]
        need(max(gram_rows) <= (10*m+4)**2, 'exact sampled norm certificate')
        digest = hashlib.sha256(bytes(K[i][j]+1 for i in range(n) for j in range(i+1,n))).hexdigest()
        samples.append(dict(N=n, seed=200+n, attempt=attempt, all_three_prime=True,
                            gram_row_bound=max(gram_rows), source_sha256=digest))
    # Hostile controls for the module detector, and an exhaustive separate path.
    for n in (3,4):
        ps = [(i,j) for i in range(n) for j in range(i+1,n)]
        for mask in range(1<<len(ps)):
            K = [[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(ps):
                K[i][j] = 1 if mask >> k & 1 else -1
                K[j][i] = -K[i][j]
            bad = False
            for s in range(1, (1<<n)-1):
                vs = [i for i in range(n) if s >> i & 1]
                if len(vs) < 2:
                    continue
                if all(len({K[x][v] for v in vs}) == 1 for x in range(n) if x not in vs):
                    bad = True
            need(prime_tournament(K) == (not bad), 'module closure versus exhaustive subsets')
    return samples


def grid_cases():
    rows = []
    for m in (3, 5, 7, 9, 11, 17, 25):
        n, s = m*m, m-2
        P, Q = paths(m)
        EP, EQ = edges(P), edges(Q)
        need(P[0] == Q[0] == 0 and P[-1] == Q[-1] == n-1, "named ports")
        need(len(EP & EQ) == 2*(m-1), "shared boundary edges")
        # Deterministic unspecified-edge completion. No primality asserted here.
        K = [[0 if x == y else (1 if x < y else -1) for y in range(n)] for x in range(n)]
        for x, y in zip(P, P[1:]):
            K[x][y], K[y][x] = 1, -1
        for x, y in zip(Q, Q[1:]):
            sign = 1 if tuple(sorted((x, y))) in EP else -1
            if sign == 1:
                need(K[x][y] == 1, "shared orientation agreement")
            K[x][y], K[y][x] = sign, -sign
        K1, K2 = native(K, P), None
        K2 = native(K1, Q)
        AP, AQ = path_matrix(P), path_matrix(Q)
        need(all(K2[x][y]-K[x][y] == 2*(AQ[x][y]-AP[x][y])
                 for x in range(n) for y in range(n)), "two-move identity")
        need(list(map(sum, K2)) == list(map(sum, K)), "exact scores preserved")
        H = EP ^ EQ
        degrees = [sum(x in e for e in H) for x in range(n)]
        need(len(H) == 2*(m-1)**2, "edge count")
        need(max(degrees) <= 4 and all(d % 2 == 0 for d in degrees), "Eulerian support")
        inn = [(r, c) for r in range(1, m-1) for c in range(1, m-1)]
        D = [[(AQ[r*m+c][u*m+v]-AP[r*m+c][u*m+v]) * (-1)**(r*c+u*v)
              for u, v in inn] for r, c in inn]
        J = lambda a, b: 1 if b == a+1 else (-1 if a == b+1 else 0)
        expected = [[J(r,u)*(c == v)-(r == u)*J(c,v) for u,v in inn] for r,c in inn]
        need(D == expected, "diagonal gauge to Kronecker difference")
        ranks = [rank_mod(D, p) for p in (3, 5)]
        need(ranks == [s*s-s]*2, "exact finite-field rank")
        rows.append(dict(m=m, N=n, toggles=len(H), interior_rank=ranks[0]))
    # Elementary union-bound constants for simultaneous prime completions.
    module_bound = F(1,16)**2/(1-F(1,16)) + F(1,2**19)/(1-F(1,2**19))
    need(3*module_bound < F(1,64), "three prime completions: union budget")
    return dict(grids=rows, three_module_failure_bound=str(3*module_bound))


def carrier_controls():
    count = 0
    for n in (3, 4):
        pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
        for mask in range(2**len(pairs)):
            B = [[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(pairs):
                B[i][j] = 1 if mask >> k & 1 else -1
                B[j][i] = -B[i][j]
            for P in permutations(range(n)):
                if not all(B[x][y] == 1 for x,y in zip(P,P[1:])):
                    continue
                K1 = native(B,P)
                for Q in permutations(range(n)):
                    if not all(K1[x][y] == 1 for x,y in zip(Q,Q[1:])):
                        continue
                    K2 = native(K1,Q)
                    H = edges(P)^edges(Q)
                    decoded = [[B[i][j]*(-1 if tuple(sorted((i,j))) in H else 1)
                                for j in range(n)] for i in range(n)]
                    need(decoded == K2, "exact marked-carrier decoder")
                    rev = list(reversed(Q))
                    opposite = [[-x for x in row] for row in K1]
                    need(native(opposite,rev) == [[-x for x in row] for row in K2], "converse covariance")
                    count += 1
    return dict(all_two_move_histories_orders_3_4=count)


def compile_program(program):
    pins, H = {}, set()
    for time, P in enumerate(program):
        for x,y in zip(P,P[1:]):
            pair = tuple(sorted((x,y)))
            sign = (-1 if (time + (pair in H)) % 2 else 1) * (1 if x < y else -1)
            if pair in pins and pins[pair] != sign:
                return None
            pins[pair] = sign
        H ^= edges(P)
    return pins


def program_controls():
    cases, inconsistent = 0, 0
    for n, depths in ((3,(1,2,3)), (4,(1,2))):
        pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
        paths_ = list(permutations(range(n)))
        bases = []
        for mask in range(1<<len(pairs)):
            B = [[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(pairs):
                B[i][j] = 1 if mask >> k & 1 else -1
                B[j][i] = -B[i][j]
            bases.append(B)
        for depth in depths:
            for program in product(paths_,repeat=depth):
                pins = compile_program(program)
                count = 0
                for B in bases:
                    K, valid = B, True
                    for P in program:
                        if not all(K[x][y] == 1 for x,y in zip(P,P[1:])):
                            valid = False
                            break
                        K = native(K,P)
                    count += valid
                predicted = 0 if pins is None else 2**(len(pairs)-len(pins))
                need(count == predicted, 'compiler iff and exact labelled completion count')
                cases += 1
                inconsistent += pins is None
    for n in (64,81,96,121,128,256):
        for r in range(1,n//32+1):
            q1 = F(n)*F(2)**(2*r-n//4)
            q2 = F(n)*F(2)**(2*r+1-n//2)
            need(q1 < 1 and q2 < 1, 'program union series converge')
            need((r+1)*(q1*q1/(1-q1)+q2/(1-q2)) < F(1,64), 'program prime budget')
    P,Q = paths(5)
    need(compile_program((P,tuple(reversed(P)))) is None, 'reverse-path hostile')
    need(compile_program((P,Q,Q,P)) is not None, 'four-move exact return program')
    return dict(programs=cases,inconsistent=inconsistent,
                literal_universes='all labelled sources n3 and n4; depths1..3 at n3,1..2 at n4')


def archaeology_controls():
    # The second occurrence of 200 is an orbit count, not a norm estimate.
    colors = range(1<<13)
    perms = [[(pow(3,k,13)*x+b)%13 for x in range(13)] for k in range(3) for b in range(13)]
    regular = 0
    for mask in colors:
        stab = sum(all(((mask>>x)&1) == ((mask>>perm[x])&1) for x in range(13)) for perm in perms)
        regular += stab == 1
    need(regular == 7800 and regular//39 == 200, "local13 regular binary orbit count")
    poly = [1]
    for multiplicity in (4,3,2,2,2,1,1,1):
        new = [0]*(len(poly)+multiplicity)
        for i, c in enumerate(poly):
            for j in range(multiplicity+1):
                new[i+j] += c
        poly = new
    need(poly[4] == poly[12] == 200 and sum(poly)-1 == 4319, "physical Gamma-bank count")
    return dict(regular_binary_orbits=200, physical_bank_size=4319, depth_census=poly)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('receipt.json'))
    args = parser.parse_args()
    out = dict(scope='FINITE-EXACT; proofs and open boundaries in companion notes',
               collatz=shadows(), tournaments=grid_cases(), carrier=carrier_controls(),
               archaeology=archaeology_controls(), prime_samples=prime_examples(),
               path_programs=program_controls())
    out['gates'] = gates
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
