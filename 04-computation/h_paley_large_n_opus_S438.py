#!/usr/bin/env python3
"""
Is Paley the H-maximiser for LARGER n?  (tests THM-1820's "Paley beaten for large n")
opus-2026-07-20-S438.

Exhaustive n<=7 (companion script) shows Paley(7) is the UNIQUE H-maximiser. Here we probe
n=11,13,19 by enumerating CIRCULANT (rotational) tournaments on Z/p: a connection set
S subset (Z/p)* with S disjoint-union -S=(Z/p)*, i->j iff (j-i) mod p in S. Paley = S=QR.
All 2^((p-1)/2) circulants are regular; if a NON-Paley circulant beats Paley(p) we have a
concrete "Paley is not the H-maximiser at n=p" witness (among the natural quasi-random family).
Also random-regular control via random circulant + relabelling has no extra power (circulant
H is relabelling-invariant), so we additionally sample random NON-circulant regular tournaments.
"""
import itertools, random

def ham_paths_circulant(S, p):
    """H of the circulant tournament on Z/p with connection set S (i->j iff (j-i)%p in S)."""
    Sset = set(S)
    adj = [[1 if (j-i) % p in Sset else 0 for j in range(p)] for i in range(p)]
    return ham_paths(adj, p)

def ham_paths(adj, n):
    size = 1 << n
    dp = [[0]*n for _ in range(size)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c:
                av = adj[v]
                for u in range(n):
                    if not (mask >> u) & 1 and av[u]:
                        dp[mask | (1 << u)][u] += c
    return sum(dp[size-1])

def circulants(p):
    """all connection sets S with, for each antipodal pair {a,-a}, exactly one chosen."""
    half = (p-1)//2
    reps = list(range(1, half+1))          # one rep per pair {a, p-a}
    for choice in itertools.product([0,1], repeat=half):
        S = []
        for r, c in zip(reps, choice):
            S.append(r if c else (p - r))
        yield tuple(sorted(S)), choice

def qr_set(p):
    return tuple(sorted(set((x*x) % p for x in range(1, p))))

for p in (7, 11, 13):
    has_paley = (p % 4 == 3)
    qr = qr_set(p) if has_paley else None
    results = []
    for S, choice in circulants(p):
        h = ham_paths_circulant(S, p)
        results.append((h, S, has_paley and S == qr))
    results.sort(reverse=True)
    hmax = results[0][0]
    n_circ = len(results)
    print(f"\n p={p}  ({n_circ} circulant tournaments){'' if has_paley else '  [no Paley: p=1 mod 4]'}")
    if has_paley:
        hpaley = next(h for h, S, isq in results if isq)
        rank_paley = 1 + sum(1 for h, S, isq in results if h > hpaley)
        n_beating = sum(1 for h, S, isq in results if h > hpaley)
        print(f"   H(Paley(Z/{p})) = {hpaley}")
        print(f"   max H over circulants = {hmax}   (Paley rank = {rank_paley} of {n_circ})")
        print(f"   #circulants beating Paley = {n_beating}")
    else:
        print(f"   max H over circulants = {hmax}  (no Paley tournament exists)")
    for h, S, isq in results[:4]:
        print(f"      H={h}  S={S}{'  <-- PALEY (QR)' if isq else ''}")
    if has_paley:
        if n_beating > 0:
            h, S, isq = results[0]
            print(f"   *** PALEY BEATEN at p={p}: circulant S={S} has H={h} > {hpaley} ***")
        else:
            print(f"   Paley is the (or tied-for) TOP circulant at p={p}.")

# random non-circulant regular control at p=11 (does a non-circulant regular beat Paley(11)?)
print("\n" + "="*60)
print("random NON-circulant regular control at n=11 (vs Paley(11))")
print("="*60)
p = 11
hpaley11 = ham_paths_circulant(qr_set(11), 11)
def random_regular_11(rng, tries=2000):
    """random regular tournament on 11 vertices via random round-robin edge orientation with
       score-balancing by local swaps (Eulerian-orientation style)."""
    n = 11
    # start from circulant Paley then randomize by 3-cycle reversals that preserve scores
    S = set(qr_set(11))
    adj = [[1 if (j-i) % n in S else 0 for j in range(n)] for i in range(n)]
    for _ in range(tries):
        a, b, c = rng.sample(range(n), 3)
        # if a->b->c->a is a directed 3-cycle, reverse it (preserves all scores)
        if adj[a][b] and adj[b][c] and adj[c][a]:
            adj[a][b]=adj[b][c]=adj[c][a]=0
            adj[b][a]=adj[c][b]=adj[a][c]=1
    return adj
rng = random.Random(999)
best = hpaley11; beat = 0
for _ in range(3000):
    adj = random_regular_11(rng, tries=300)
    h = ham_paths(adj, 11)
    if h > best: best = h; beat += 1
print(f" H(Paley(11)) = {hpaley11}")
print(f" best over 3000 random regular tournaments on 11 vertices = {best}  (#beating Paley = {beat})")
