#!/usr/bin/env python3
"""rado_tournament_converse_swap_pfaffian_s690.py — "a loop of the input causes a
swap of the two outputs; consider the Rado graph as a tournament."

THE TWO CLUES, MADE PRECISE.

(I) A relation is a function on ORDERED pairs R(a,b). The transposition (a,b)↦(b,a)
is the deck transformation of the 2-fold cover (ordered → unordered pairs) — "a loop
of the input". For a GRAPH (symmetric) the loop FIXES the output R(a,b)=R(b,a); for a
TOURNAMENT (antisymmetric) it SWAPS/flips it: a→b XOR b→a. So a tournament is the
ODD / twisted section of the double cover, a graph the EVEN section — the ±√ branch
structure (loop around the branch point swaps the two roots). The generic EVEN object
is the RADO graph; the generic ODD object is the generic countable tournament
(Fraïssé limit = a.s. the iid-arc random tournament).

(II) The "loop" as a global involution = the CONVERSE T↦T* (reverse every arc). On the
skew-adjacency S (S_ij=+1 if i→j else −1) the converse is S ↦ −S, so the Pfaffian
            Pf(−S) = (−1)^{n/2} Pf(S).
THM-203: H(T*) = H(T) (converse-invariant). THM-174: det(I+2A)=Pf(S)², H²≡Pf² (mod 8),
Q=(H²−Pf²)/8 ∈ ℤ≥0. Therefore, for n ≡ 2 (mod 4) (Pf negated), the converse loop
            SWAPS THE TWO OUTPUTS   H+Pf  ⟷  H−Pf,
the two roots of x² − 2H x + 8Q (product H²−Pf² = 8Q, discriminant 4·Pf²). The branch
point — where the two outputs MERGE — is Pf = 0. So the user's sentence is literally
the degree-2 FTA monodromy of (H ± Pf) under the converse, with H the swap-EVEN part
and Pf the swap-ODD part (the √ that the loop acts on). For n ≡ 0 (mod 4) the loop
FIXES Pf (no swap).

Self-converse tournaments (T*≅T, the FIXED POINTS of the loop) satisfy −S=PSP^T for a
permutation P ⟹ (−1)^{n/2} Pf(S) = det(P) Pf(S): if Pf≠0 the converse-iso P is an ODD
permutation for n≡2 (mod 4). The fixed-point set spans from the CYCLOTOMIC pole
(round tournaments on roots of unity = the LRC worry-set, HYP-2089) to the RANDOM pole
(the Rado/generic tournament). The repo mined the cyclotomic pole; this is the random
pole.

Session: claude-2026-06-06-S690 (rado-tournament-converse-swap-pfaffian)."""
import sys; sys.stdout.reconfigure(line_buffering=True)
from itertools import product, combinations, permutations
import random

def edges(n): return [(i, j) for i in range(n) for j in range(i+1, n)]

def build_adj(n, bits):
    """bits over edges(n): 1 => i->j, 0 => j->i. returns adjacency bitmask rows."""
    adj = [0]*n
    for (i, j), b in zip(edges(n), bits):
        if b: adj[i] |= 1 << j
        else: adj[j] |= 1 << i
    return adj

def H(n, adj):
    size = 1 << n; dp = [[0]*n for _ in range(size)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(size):
        for v in range(n):
            c = dp[mask][v]
            if c:
                av = adj[v]
                for w in range(n):
                    if not (mask >> w & 1) and av >> w & 1: dp[mask | 1 << w][w] += c
    return sum(dp[size-1])

def skew(n, adj):
    S = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                S[i][j] = 1 if (adj[i] >> j & 1) else -1
    return S

def pfaffian(S):
    """exact Pfaffian of an even-dim skew matrix via the perfect-matching sum."""
    n = len(S); assert n % 2 == 0
    def rec(verts):
        if not verts: return 1
        i = verts[0]; rest = verts[1:]; total = 0
        for k, j in enumerate(rest):
            sign = -1 if (k % 2) else 1   # moving j past k elements to pair with i
            total += sign * S[i][j] * rec(rest[:k] + rest[k+1:])
        return total
    return rec(tuple(range(n)))

# ============================================================
print("=== (II) THE CONVERSE LOOP SWAPS H±Pf  (verify on n=6, the n≡2 mod 4 case) ===")
n = 6; es = edges(n); m = len(es)
from collections import Counter
pf_dist = Counter(); ok_swap = True; ok_H = True; ok_Q = True
examples = []
cnt = 0
for bits in product([0, 1], repeat=m):
    adj = build_adj(n, bits); S = skew(n, adj)
    h = H(n, adj); pf = pfaffian(S)
    # converse: reverse every arc -> complement bits
    cbits = tuple(1 - b for b in bits); cadj = build_adj(n, cbits); cS = skew(n, cadj)
    hc = H(n, cadj); pfc = pfaffian(cS)
    if hc != h: ok_H = False
    if pfc != -pf: ok_swap = False         # n=6: n/2=3 odd => Pf negated
    Q8 = h*h - pf*pf
    if Q8 < 0 or Q8 % 8 != 0: ok_Q = False
    pf_dist[abs(pf)] += 1
    if cnt < 3:
        examples.append((h, pf, hc, pfc, Q8//8)); cnt += 1
print(f"  enumerated all {1<<m} tournaments on n=6")
print(f"  H(T*) = H(T) [converse-invariant, THM-203]:           {ok_H}")
print(f"  Pf(T*) = -Pf(T) [converse negates Pf at n≡2 mod 4]:    {ok_swap}")
print(f"  ⟹ converse swaps the two outputs (H+Pf) ⟷ (H−Pf)")
print(f"  H²−Pf² = 8Q with Q∈ℤ≥0 [THM-174]:                     {ok_Q}")
print(f"  |Pf| distribution over all n=6 tournaments: {dict(sorted(pf_dist.items()))}")
print(f"  examples (H, Pf, H*, Pf*, Q):")
for (h, pf, hc, pfc, Q) in examples:
    print(f"     T:(H={h:>2},Pf={pf:+d})  T*:(H={hc:>2},Pf={pfc:+d})   "
          f"outputs H±Pf = {h+pf},{h-pf} -> swapped to {hc+pfc},{hc-pfc};  product=8·{Q}")
print(f"  branch point (two outputs merge): Pf=0 ⟹ # tournaments with Pf=0: {pf_dist[0]}")

print("\n=== contrast: n=4 (n≡0 mod 4) — the loop FIXES Pf, no swap ===")
n4 = 4; es4 = edges(n4); m4 = len(es4); fix = True
for bits in product([0, 1], repeat=m4):
    adj = build_adj(n4, bits); S = skew(n4, adj)
    pf = pfaffian(S)
    cbits = tuple(1-b for b in bits); cpf = pfaffian(skew(n4, build_adj(n4, cbits)))
    if cpf != pf: fix = False
print(f"  Pf(T*) = +Pf(T) for all n=4 tournaments (n/2=2 even): {fix}  ⟹ no output-swap at n≡0 mod 4")

print("\n=== self-converse fixed points: the converse-iso parity constraint (n=6) ===")
# find self-converse tournaments and check det(P) = (-1)^{n/2} = -1 when Pf != 0
def perm_sign(p):
    p = list(p); s = 1; seen = [False]*len(p)
    for i in range(len(p)):
        if not seen[i]:
            j = i; L = 0
            while not seen[j]: seen[j] = True; j = p[j]; L += 1
            if L % 2 == 0: s = -s
    return s
random.seed(1); checked = 0; sc_found = 0; parity_ok = True; pf0 = 0
samp = [tuple(random.randint(0,1) for _ in range(m)) for _ in range(4000)]
for bits in samp:
    adj = build_adj(n, bits); S = skew(n, adj); pf = pfaffian(S)
    target = [[-S[i][j] for j in range(n)] for i in range(n)]   # -S
    # is there P with P S P^T = -S ?  search permutations (n=6 -> 720)
    found = None
    for p in permutations(range(n)):
        ok = True
        for i in range(n):
            pi = p[i]
            row = S[i]
            for j in range(n):
                if row[j] != target[pi][p[j]]:  # (P S P^T)[pi][pj] = S[i][j]
                    ok = False; break
            if not ok: break
        if ok: found = p; break
    if found is not None:
        sc_found += 1; checked += 1
        if pf == 0: pf0 += 1
        else:
            if perm_sign(found) != -1: parity_ok = False
    if sc_found >= 40: break
print(f"  among sampled tournaments, found {sc_found} self-converse; for those with Pf≠0,")
print(f"  the converse-iso P is an ODD permutation (det P = −1 = (−1)^3): {parity_ok}  "
      f"(Pf=0 cases: {pf0})")

# ============================================================
print("\n=== (I) THE RANDOM/RADO TOURNAMENT: extension property converges, self-converse ===")
def has_k_extension(n, adj, k):
    """for every disjoint (A,B) with |A|+|B|=k, ∃ v∉A∪B with v->all A, all B->v?
       (the k-point existential-closure / extension axiom of the generic tournament)."""
    verts = range(n)
    for size in range(1, k+1):
        for A in combinations(verts, size):
            restA = [v for v in verts if v not in A]
            # B disjoint from A, sizes summing to k handled by trying all B with |B|=k-size
            for B in combinations(restA, k - size):
                S_ = set(A) | set(B)
                witness = False
                for v in verts:
                    if v in S_: continue
                    if all(adj[v] >> a & 1 for a in A) and all(adj[b] >> v & 1 for b in B):
                        witness = True; break
                if not witness: return False
    return True

random.seed(2)
print("  P(random tournament satisfies the k-extension axiom), k=1 then k=2:")
for n in range(3, 11):
    T = 3000; c1 = 0; c2 = 0
    for _ in range(T):
        bits = [random.randint(0,1) for _ in range(len(edges(n)))]
        adj = build_adj(n, bits)
        if has_k_extension(n, adj, 1): c1 += 1
        if has_k_extension(n, adj, 2): c2 += 1
    print(f"    n={n:>2}:  k=1 (no source/sink-pattern): {c1/T:.3f}   k=2: {c2/T:.3f}")
print("  ⟹ →1: the iid-arc random tournament a.s. satisfies ALL extension axioms")
print("     = the Fraïssé limit = the GENERIC countable tournament (the Rado graph's odd twin).")
print("  Self-converse: reversing iid-fair arcs is iid-fair ⟹ the generic tournament ≅ its")
print("  converse — it is a FIXED POINT of the loop, the RANDOM pole opposite the cyclotomic")
print("  worry-set pole (round tournaments on roots of unity, HYP-2089).")
