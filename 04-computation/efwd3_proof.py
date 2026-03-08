#!/usr/bin/env python3
"""
efwd3_proof.py — Prove E[fwd^3] = A(n) + (6/n)*t3.

PROOF STRATEGY:
fwd = sum_{i=0}^{n-2} X_i where X_i = 1_{P[i]->P[i+1] in T}

fwd^3 = (sum X_i)^3
Since X_i^r = X_i for binary variables:
fwd^3 = sum X_i + 3*sum_{i!=j} X_i*X_j + 6*sum_{i<j<k} X_i*X_j*X_k
     = fwd + 3(fwd^2 - fwd) + 6*C(fwd,3)*6/6   (wrong expansion)

Actually: (sum X)^3 = sum X_i^3 + 3*sum_{i<j}(X_i^2*X_j + X_i*X_j^2) + 6*sum_{i<j<k} X_i*X_j*X_k
= sum X_i + 3*sum_{i<j}(X_i*X_j + X_i*X_j) + 6*sum_{i<j<k} X_i*X_j*X_k
= sum X_i + 6*sum_{i<j} X_i*X_j + 6*sum_{i<j<k} X_i*X_j*X_k

But sum_{i<j} X_i*X_j = (fwd^2 - fwd)/2
And sum_{i<j<k} X_i*X_j*X_k = C(fwd,3) (since X_i binary)

Hmm, let me just use: fwd^3 = fwd*(fwd-1)*(fwd-2) + 3*fwd*(fwd-1) + fwd
  = fwd(fwd-1)(fwd-2) + 3*fwd(fwd-1) + fwd

So E[fwd^3] = E[fwd(fwd-1)(fwd-2)] + 3*E[fwd(fwd-1)] + E[fwd]

E[fwd] = (n-1)/2
E[fwd(fwd-1)] = sum_{i<j} E[X_i*X_j] (factorial moment)
E[fwd(fwd-1)(fwd-2)] = 6*sum_{i<j<k} E[X_i*X_j*X_k]

So the t3-dependent part comes from E[X_i*X_j*X_k].

By THM-089: non-adjacent X_i, X_j are uncorrelated.
For triples: if all three are pairwise non-adjacent, then E[X_i*X_j*X_k] = 1/8.
If exactly one pair is adjacent: need Cov(X_i, X_{i+1}) contribution.
If two pairs adjacent (consecutive triple i, i+1, i+2): need 3-point correlation.

Let me compute E[X_i*X_{i+1}*X_{i+2}] directly and see how it depends on t3.

Author: opus-2026-03-07-S46c
"""
from itertools import permutations, combinations
from math import comb, factorial
from fractions import Fraction

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def count_3cycles(adj, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t3 += 1
    return t3

print("=" * 60)
print("DECOMPOSITION OF E[fwd^3]")
print("=" * 60)

for n in [4, 5, 6]:
    m_vals = n*(n-1)//2
    seen = set()
    data = []

    for bits in range(1 << m_vals):
        adj = tournament_from_bits(n, bits)

        # Compute moments of X_i products directly
        total = factorial(n)

        # E[X_i X_{i+1} X_{i+2}] for consecutive triple
        # This is Pr[P[i]->P[i+1], P[i+1]->P[i+2], P[i+2]->P[i+3]]
        # For fixed i, this averages over all perms

        F_key = []
        for P in permutations(range(n)):
            fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
            F_key.append(fwd)
        F_key = tuple(sorted(F_key))

        if F_key in seen:
            continue
        seen.add(F_key)

        t3 = count_3cycles(adj, n)

        # Compute triple correlations
        # E[X_i X_{i+1} X_{i+2}] for consecutive triples
        consec_triple_sum = Fraction(0)
        for P in permutations(range(n)):
            for i in range(n-3):
                if adj[P[i]][P[i+1]] and adj[P[i+1]][P[i+2]] and adj[P[i+2]][P[i+3]]:
                    consec_triple_sum += 1

        # Normalize: (n-3) positions * n! perms
        E_consec_triple = consec_triple_sum / (total * (n-3)) if n > 3 else Fraction(0)

        # Also compute:
        # E[X_i X_{i+1} X_j] for non-adjacent j (|j-i|>=2 and |j-(i+1)|>=2)
        adj_nonadj_sum = Fraction(0)
        adj_nonadj_count = 0
        for P in permutations(range(n)):
            for i in range(n-2):
                for j in range(n-1):
                    if j != i and j != i+1 and abs(j-i) >= 2 and abs(j-(i+1)) >= 2:
                        if adj[P[i]][P[i+1]] and adj[P[i+1]][P[i+2]] if i+1 == j-1 else True:
                            pass
                        # Actually just check X_i, X_{i+1}, X_j
                        if adj[P[i]][P[i+1]] and adj[P[i+1]][P[i+2]] and adj[P[j]][P[j+1]]:
                            if abs(j - i) >= 2 and abs(j - (i+1)) >= 2:
                                adj_nonadj_sum += 1
                                adj_nonadj_count += 1

        # Count directed 4-paths (forward 3-step paths)
        fwd4path = 0
        for combo in combinations(range(n), 4):
            for perm in permutations(combo):
                if adj[perm[0]][perm[1]] and adj[perm[1]][perm[2]] and adj[perm[2]][perm[3]]:
                    fwd4path += 1

        data.append({
            't3': t3,
            'E_consec': E_consec_triple,
            'fwd4path': fwd4path,
        })

    print(f"\nn={n}: {len(data)} F-classes")
    print(f"  {'t3':>3} {'E[X_i X_{i+1} X_{i+2}]':>28} {'#fwd4path':>10}")
    for d in sorted(data, key=lambda x: x['t3']):
        print(f"  {d['t3']:>3} {str(d['E_consec']):>28} {d['fwd4path']:>10}")

    # Check: is E_consec linear in t3?
    items = [(d['t3'], d['E_consec']) for d in data]
    t3_vals = sorted(set(t for t, _ in items))
    t3_to_e = {t: None for t in t3_vals}
    for t, e in items:
        if t3_to_e[t] is None:
            t3_to_e[t] = e
        elif t3_to_e[t] != e:
            t3_to_e[t] = 'AMBIGUOUS'

    if all(v != 'AMBIGUOUS' for v in t3_to_e.values()) and len(t3_vals) >= 2:
        B = (t3_to_e[t3_vals[1]] - t3_to_e[t3_vals[0]]) / (t3_vals[1] - t3_vals[0])
        A = t3_to_e[t3_vals[0]] - B * t3_vals[0]
        ok = all(t3_to_e[t] == A + B*t for t in t3_vals)
        print(f"  E[consec triple] = {A} + {B}*t3  {'[EXACT]' if ok else '[APPROX]'}")

    # fwd4path linear in t3?
    items2 = [(d['t3'], d['fwd4path']) for d in data]
    t3_to_f4 = {t: None for t in t3_vals}
    for t, f in items2:
        if t3_to_f4[t] is None:
            t3_to_f4[t] = f
        elif t3_to_f4[t] != f:
            t3_to_f4[t] = 'AMBIGUOUS'

    if all(v != 'AMBIGUOUS' for v in t3_to_f4.values()) and len(t3_vals) >= 2:
        B2 = Fraction(t3_to_f4[t3_vals[1]] - t3_to_f4[t3_vals[0]], t3_vals[1] - t3_vals[0])
        A2 = t3_to_f4[t3_vals[0]] - B2 * t3_vals[0]
        ok2 = all(t3_to_f4[t] == A2 + B2*t for t in t3_vals)
        print(f"  #fwd4path = {A2} + {B2}*t3  {'[EXACT]' if ok2 else '[APPROX]'}")
    else:
        print(f"  #fwd4path: NOT determined by t3 alone")
