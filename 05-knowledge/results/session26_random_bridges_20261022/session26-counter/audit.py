from itertools import combinations, permutations, product
from math import comb, factorial
from fractions import Fraction
import json
from counter import count_hamilton


def check(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def core(k, mask):
    Q = [[0]*k for _ in range(k)]
    for bit, (u,v) in enumerate(combinations(range(k), 2)):
        Q[u][v] = (mask >> bit) & 1
        Q[v][u] = 1-Q[u][v]
    return Q


def graph(Q, word):
    k, n = len(Q), len(Q)+len(word)
    G = [[0]*n for _ in range(n)]
    for u in range(k):
        for v in range(k): G[u][v] = Q[u][v]
    for i, mask in enumerate(word):
        for u in range(k):
            G[u][k+i] = (mask >> u) & 1
            G[k+i][u] = 1-G[u][k+i]
        for j in range(i): G[k+j][k+i] = 1
    return G


def subset_dp(G):
    n = len(G)
    if not n: return 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for used in range(1, 1 << n):
        for v in range(n):
            if dp[used][v]:
                for u in range(n):
                    if not used >> u & 1 and G[v][u]:
                        dp[used | (1 << u)][u] += dp[used][v]
    return sum(dp[-1])


def profile(Q):
    k = len(Q)
    p = [0]*(k+1)
    for sigma in permutations(range(k)):
        bad = sum(not Q[u][v] for u,v in zip(sigma, sigma[1:]))
        for r in range(bad+1, k+1):
            p[r] += comb(k-1-bad, r-1-bad)
    return p


def surjections(m,q):
    if q < 0: return 0
    return sum((-1)**(q-j)*comb(q,j)*j**m for j in range(q+1))


def expected(Q,m):
    return sum(Fraction(p,4**r)*(surjections(m,r+1)+4*surjections(m,r)+4*surjections(m,r-1))
               for r,p in enumerate(profile(Q)) if r)


def one_apex(word):
    h,a,p = 1,0,0
    for b in word:
        if not b: h += a+p
        a,p = 2*a+b,b
    return h


def main():
    totals = {}
    means = 0
    for k, max_m in [(0,8),(1,7),(2,5),(3,3),(4,2)]:
        tested = 0
        for mask in range(1 << comb(k,2)):
            Q = core(k,mask)
            for m in range(max_m+1):
                mass = 0
                for word in product(range(1 << k), repeat=m):
                    h = count_hamilton(Q, word)
                    check(h == subset_dp(graph(Q,word)), ('DP mismatch',k,mask,word,h))
                    if k == 1: check(h == one_apex(word), ('one-apex mismatch',word))
                    mass += h
                    tested += 1
                if k:
                    check(Fraction(mass,2**(k*m)) == expected(Q,m), ('mean mismatch',k,mask,m))
                    means += 1
        totals[k] = tested
    # All binary words through length 12 check independently inherited 4-coordinate recurrence.
    recurrence_tests = 0
    for m in range(13):
        for word in product(range(2), repeat=m):
            check(count_hamilton([[0]],word) == one_apex(word), ('recurrence',word))
            recurrence_tests += 1
    # Fixed known five-vertex mask-41 witness; derive a marked transitive deletion.
    witness = core(5,41)
    check(subset_dp(witness) == 13, 'known H13 witness')
    found = None
    for K in combinations(range(5),2):
        B = [v for v in range(5) if v not in K]
        for order in permutations(B):
            if all(witness[u][v] for i,u in enumerate(order) for v in order[i+1:]):
                Q = [[witness[u][v] for v in K] for u in K]
                word = [sum(witness[u][b] << j for j,u in enumerate(K)) for b in order]
                check(count_hamilton(Q,word) == 13, 'H13 counter')
                found = dict(original_core=K, original_base_order=order,Q=Q,word=word,H=13)
                break
        if found: break
    check(found is not None,'transitive triple missing')
    with open('witness13.json','w') as f: json.dump(found,f,indent=2); f.write('\n')
    # Orientation-sensitive and first-versus-last hostiles.
    check(count_hamilton([[0]],[1,0]) == 3, '10 triangle')
    check(count_hamilton([[0]],[0,1]) == 1, '01 transitive')
    check(count_hamilton([[0]],[1,0,1]) == 3, 'bad exit may later become good elsewhere')
    check(count_hamilton([[0,1],[0,0]],[1]) == 1, 'fixed-word core orientation A')
    check(count_hamilton([[0,0],[1,0]],[1]) == 3, 'fixed-word core orientation B')
    # Each invalid input must be rejected even in python -O.
    invalid = [([[1]],[]), ([[0,0],[0,0]],[]), ([[0,1]],[]), ([[0]],[-1]), ([[0]],[2]), ([],[1]), ([[0]],[True]), ([[0.0]],[])]
    for Q,word in invalid:
        try: count_hamilton(Q,word)
        except ValueError: pass
        else: raise RuntimeError(('accepted invalid input',Q,word))
    # Longer-word controls are not an exhaustive spectrum enumeration.
    long_word = [3,1,2,0]*25
    long_h = count_hamilton([[0,1],[0,0]],long_word)
    check(long_h > 0 and long_h <= factorial(2)*3**100,'long bound')
    result = dict(status='PASS',subset_dp_cases_by_k=totals,subset_dp_total=sum(totals.values()),
                  exact_mean_pairs=means,one_apex_recurrence_cases=recurrence_tests,
                  invalid_inputs=len(invalid),witness13=found,long_word_length=100,long_word_H=long_h)
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__ == '__main__': main()
