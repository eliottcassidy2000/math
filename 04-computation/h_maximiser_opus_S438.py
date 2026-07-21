#!/usr/bin/env python3
"""
H-MAXIMISER: the tournament with the MOST Hamiltonian paths (H(T), definitions.md:15).
Answers THM-1820 open Q1 ("what maximises H?") and tests the WOWII inflation/decoupling
motif on H (opus-2026-07-20-S438).

  A. exact max H + maximiser structure (score seq, skew spectrum, regular?) n=3..6 exhaustive;
     n=7 over the regular class (score all = 3; forced by max-min<=1) + random control.
  B. REFUTE "H is score-determined": same-score-different-H witnesses => THM-1820's
     "H Schur-concave in scores" is imprecise (H is NOT a function of the score sequence).
  C. VERIFY source/sink inflation is H-NEUTRAL: H(T+source)=H(T)=H(T+sink). The WOWII
     pendant (dominating/dominated vertex) does NOT pump H -- H-maximiser is
     inflation-unreachable, opposite to the alpha-graph where leaves pump alpha.
  D. Paley: is the quadratic-residue tournament the H-maximiser at n=3,7?
"""
import itertools
try:
    import numpy as np
    HAVE_NP = True
except Exception:
    HAVE_NP = False

# ---------- core ----------
def ham_paths(adj, n):
    """# directed Hamiltonian paths. dp[mask][v]=#paths covering mask ending at v."""
    size = 1 << n
    dp = [[0]*n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c:
                av = adj[v]
                base = dp
                for u in range(n):
                    if not (mask >> u) & 1 and av[u]:
                        base[mask | (1 << u)][u] += c
    return sum(dp[size-1])

def scores(adj, n):
    return tuple(sorted(sum(adj[i]) for i in range(n)))

def skew_abs_spectrum(adj, n):
    if not HAVE_NP: return None
    S = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                S[i][j] = 1.0 if adj[i][j] else -1.0
    ev = np.linalg.eigvals(np.array(S))
    return sorted(round(abs(e), 3) for e in ev)

def edges_iter(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(1 << len(pairs)):
        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: adj[i][j] = 1
            else:               adj[j][i] = 1
        yield adj

def canon(adj, n):
    """canonical relabelling key (brute force; fine for n<=7 on a handful of matrices)."""
    best = None
    for p in itertools.permutations(range(n)):
        key = 0; bit = 0
        for i in range(n):
            for j in range(n):
                if i != j:
                    key |= (adj[p[i]][p[j]]) << bit; bit += 1
        if best is None or key < best: best = key
    return best

def paley_adj(p):
    qr = set((x*x) % p for x in range(1, p))
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr: adj[i][j] = 1
    return adj

# ---------- A ----------
print("="*70); print("A. EXACT H-MAXIMISER, exhaustive n=3..6"); print("="*70)
maxdata = {}
for n in range(3, 7):
    best_H, best_adjs = -1, []
    by_score = {}
    for adj in edges_iter(n):
        h = ham_paths(adj, n)
        by_score.setdefault(scores(adj, n), set()).add(h)
        if h > best_H: best_H, best_adjs = h, [adj]
        elif h == best_H: best_adjs.append(adj)
    iso = {}
    for adj in best_adjs: iso[canon(adj, n)] = adj
    reps = list(iso.values())
    print(f"\n n={n}: max H={best_H}  ({len(best_adjs)} labelled, {len(reps)} iso class)")
    for adj in reps:
        s = scores(adj, n); reg = len(set(s)) == 1
        print(f"     score={s}{' REGULAR' if reg else ''}  |skew|={skew_abs_spectrum(adj,n)}")
    maxdata[n] = (best_H, reps, by_score)

# ---------- A (n=7 regular class) ----------
print("\n" + "="*70); print("A'. n=7: regular class (scores all=3, forced by max-min<=1) + control")
print("="*70)
n = 7
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
best_H, best_bits = -1, []
count_reg = 0
for bits in range(1 << len(pairs)):
    deg = [0]*n
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1: deg[i] += 1
        else:               deg[j] += 1
    if max(deg) - min(deg) != 0:   # regular only (all scores = 3)
        continue
    count_reg += 1
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1: adj[i][j] = 1
        else:               adj[j][i] = 1
    h = ham_paths(adj, n)
    if h > best_H: best_H, best_bits = h, [bits]
    elif h == best_H: best_bits.append(bits)

def bits_to_adj(bits):
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1: adj[i][j] = 1
        else:               adj[j][i] = 1
    return adj

iso = {}
for bits in best_bits: iso[canon(bits_to_adj(bits), n)] = bits_to_adj(bits)
reps7 = list(iso.values())
print(f" labelled regular tournaments on 7 vertices: {count_reg}")
print(f" max H over regular = {best_H}  ({len(best_bits)} labelled, {len(reps7)} iso class)")
pal7 = paley_adj(7); pal7_canon = canon(pal7, 7); pal7_H = ham_paths(pal7, 7)
print(f" H(Paley(7)) = {pal7_H};  maximiser iso-classes canon == Paley canon? ",
      [canon(a,7) == pal7_canon for a in reps7])
for adj in reps7:
    print(f"     score={scores(adj,7)}  |skew|={skew_abs_spectrum(adj,7)}")

# control: does ANY non-regular tournament beat best regular? random sample.
import random
random.seed(12345)
beat = 0; trials = 200000
for _ in range(trials):
    bits = random.getrandbits(len(pairs))
    adj = bits_to_adj(bits)
    if ham_paths(adj, 7) > best_H: beat += 1
print(f" random control ({trials} tournaments): #beating best-regular = {beat}")

# ---------- B. H is NOT score-determined ----------
print("\n" + "="*70); print("B. IS H SCORE-DETERMINED?  (refutes 'H Schur-concave in scores')")
print("="*70)
for n in range(4, 7):
    by_score = maxdata[n][2] if n <= 6 else None
    multis = {s: hs for s, hs in by_score.items() if len(hs) > 1}
    print(f" n={n}: {len(multis)}/{len(by_score)} score sequences carry MULTIPLE distinct H values")
    ex = sorted(multis.items(), key=lambda kv: -len(kv[1]))[:2]
    for s, hs in ex:
        print(f"     score={s}: H in {sorted(hs)}  <-- same scores, different H")

# ---------- C. source/sink inflation is H-neutral ----------
print("\n" + "="*70); print("C. INFLATION NEUTRALITY: H(T+source)=H(T)=H(T+sink)?")
print("="*70)
def add_source(adj, n):   # new vertex n beats everyone
    a = [row[:] + [0] for row in adj] + [[1]*n + [0]]
    return a
def add_sink(adj, n):     # new vertex n loses to everyone
    a = [row[:] + [1] for row in adj] + [[0]*n + [0]]
    return a
ok = True
for n in range(3, 6):
    for adj in itertools.islice(edges_iter(n), 0, 1 << (n*(n-1)//2)):
        h = ham_paths(adj, n)
        hs = ham_paths(add_source(adj, n), n+1)
        hk = ham_paths(add_sink(adj, n), n+1)
        if not (h == hs == hk):
            ok = False
            print(f"   COUNTEREXAMPLE n={n}: H={h} source={hs} sink={hk}"); break
print(f" source/sink inflation H-neutral on ALL tournaments n=3,4,5: {ok}")
print("   => the WOWII pendant does NOT pump H; the H-maximiser is inflation-UNREACHABLE.")

# ---------- D. summary sequence ----------
print("\n" + "="*70); print("D. MAX-H SEQUENCE and maximiser type")
print("="*70)
seq = [maxdata[n][0] for n in range(3,7)] + [best_H]
print(f" max H for n=3..7: {seq}")
print(f" Szele average n!/2^(n-1) for n=3..7: {[__import__('math').factorial(n)//2**(n-1) for n in range(3,8)]}")
