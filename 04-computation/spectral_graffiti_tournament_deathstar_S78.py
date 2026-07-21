#!/usr/bin/env python3
"""
spectral_graffiti_tournament_deathstar_S78.py   (HYP-8636)

Distinctive strand of the TournamentGraffiti program: the SPECTRAL / ALGEBRAIC
invariants that the fleet's classical zoos omitted.
  - kind-pasteur THM-1845 zoo: c3, H, beta(=tr), dom, kings, scc, scores, arb0.
  - klein THM-1850 directed zoo: tr, dichromatic, gamma, fas, #C3, H, diam.
  - Neither used: skew-determinant disc (THM-474), #distinct eigenvalues (my S76),
    spectral radius, or the {7,21} H-spectrum forbidden-value structure (my S70).

Builds those, mines inequalities against the classical ones, and checks the
Redei {7,21} H-spectrum gap. Reuses/credits the fleet's proved anchors
(klein gamma+tr<=n+1; kind-pasteur n-c3 <= beta <= smax+1).
"""
from math import comb
from itertools import combinations
try:
    import numpy as np
    HAVE_NP = True
except Exception:
    HAVE_NP = False

def sep(t): print("\n" + "=" * 72 + "\n" + t + "\n" + "=" * 72)

def all_tournaments(n):
    pairs = list(combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1: A[i][j] = 1
            else:             A[j][i] = 1
        yield A

def scores(A, n): return [sum(A[i]) for i in range(n)]

def c3(A, n):
    s = scores(A, n)
    return comb(n, 3) - sum(comb(x, 2) for x in s)

def ham_paths(A, n):
    full = (1 << n) - 1
    out = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            c = dp[mask][v]
            if c:
                for w in range(n):
                    if not (mask >> w & 1) and (out[v] >> w & 1):
                        dp[mask | 1 << w][w] += c
    return sum(dp[full][v] for v in range(n))

def largest_transitive(A, n):
    best = 0
    for mask in range(1 << n):
        S = [v for v in range(n) if mask >> v & 1]
        if not S: continue
        isc = sorted(sum(A[v][w] for w in S) for v in S)
        if isc == list(range(len(S))): best = max(best, len(S))
    return best

def num_kings(A, n):
    full = (1 << n) - 1
    out = [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]
    k = 0
    for v in range(n):
        r = out[v] | (1 << v)
        for w in range(n):
            if out[v] >> w & 1: r |= out[w]
        if r == full: k += 1
    return k

def spectral(A, n):
    """returns (spectral_radius, n_distinct_eigs, skewdet_int)."""
    M = np.array(A, dtype=float)
    ev = np.linalg.eigvals(M)
    rho = float(max(abs(ev)))
    # distinct eigenvalues by rounding real+imag to 4 dp
    seen = set()
    for z in ev: seen.add((round(z.real, 4), round(z.imag, 4)))
    ndev = len(seen)
    K = M - M.T                       # skew-adjacency, entries +-1 off diag
    skewdet = round(float(np.linalg.det(np.eye(n) + K)))   # det(I+K), integer
    return rho, ndev, skewdet

# ---------- build database ----------
print("building spectral invariant database n=3..6 ...", flush=True)
DATA = {}
for n in (3, 4, 5, 6):
    rows = []
    for A in all_tournaments(n):
        s = scores(A, n)
        row = dict(n=n, c3=c3(A, n), H=ham_paths(A, n), beta=largest_transitive(A, n),
                   kings=num_kings(A, n), smax=max(s), smin=min(s), srange=max(s)-min(s))
        if HAVE_NP:
            rho, ndev, skewdet = spectral(A, n)
            row.update(rho=rho, ndev=ndev, skewdet=skewdet)
        rows.append(row)
    DATA[n] = rows
    print(f"  n={n}: {len(rows)} tournaments done", flush=True)

# ---------- (1) the Redei {7,21} H-spectrum forbidden-value check ----------
sep("(1) Redei H-spectrum: which ODD values are achievable? the {7,21} gap (S70)")
for n in (3, 4, 5, 6):
    hs = sorted(set(d['H'] for d in DATA[n]))
    print(f"  n={n}: achievable H = {hs}")
allH = sorted(set(d['H'] for n in DATA for d in DATA[n]))
maxH = max(allH)
missing_odd = [k for k in range(1, maxH+1, 2) if k not in allH]
print(f"  union H (n<=6) = {allH}")
print(f"  ODD values in [1,{maxH}] NOT achieved at n<=6: {missing_odd}   (S70: forbidden = 7,21)")

# ---------- (2) spectral invariant ranges ----------
if HAVE_NP:
    sep("(2) spectral invariant ranges (distinctive: not in the classical zoos)")
    for n in (3, 4, 5, 6):
        rr = [d['rho'] for d in DATA[n]]; nd = [d['ndev'] for d in DATA[n]]
        sk = [d['skewdet'] for d in DATA[n]]
        print(f"  n={n}: spectral_radius [{min(rr):.3f},{max(rr):.3f}], "
              f"#distinct-eig {sorted(set(nd))}, det(I+K) {sorted(set(sk))}")

# ---------- (3) auto-mine comparabilities + tight bounds among all pairs ----------
sep("(3) MINER: which invariant inequalities hold on ALL n<=6? (survivors + tight)")
KEYS = ["c3", "H", "beta", "kings", "smax", "srange"] + (["ndev", "skewdet"] if HAVE_NP else [])
ALL = [d for n in DATA for d in DATA[n]]
def holds(f):
    w = [d for d in ALL if not f(d)]
    return (len(w) == 0, w[0] if w else None)

# fleet anchors (sanity: must hold)
anchors = [
    ("kind-pasteur: n - c3 <= beta", lambda d: d['n'] - d['c3'] <= d['beta']),
    ("kind-pasteur: beta <= smax + 1", lambda d: d['beta'] <= d['smax'] + 1),
]
print("  -- fleet anchors (must hold) --")
for name, f in anchors:
    ok, _ = holds(f); print(f"    [{'OK ' if ok else 'FAIL'}] {name}")

# my distinctive candidates involving spectral invariants
if HAVE_NP:
    cands = [
        ("ndev <= beta            (distinct eigs <= largest transitive sub)", lambda d: d['ndev'] <= d['beta']),
        ("ndev >= 1 + [c3>0]      (transitive has 1 eig; a cycle splits it)", lambda d: d['ndev'] >= 1 + (1 if d['c3']>0 else 0)),
        ("c3 <= C(n,3) - (ndev-1) (more distinct eigs bounds cycles?)", lambda d: d['c3'] <= comb(d['n'],3) - (d['ndev']-1)),
        ("H <= skewdet? ", lambda d: d['H'] <= d['skewdet']),
        ("H >= |det(I+K)|/2^(n-1) (H vs normalized skew-det disc)", lambda d: d['H'] >= abs(d['skewdet'])/2**(d['n']-1)),
        ("ndev-1 <= srange        (#distinct eigs vs score range)", lambda d: d['ndev']-1 <= d['srange']),
        ("rho >= (n-1)/2          (spectral radius >= avg score)", lambda d: d['rho'] >= (d['n']-1)/2 - 1e-9),
    ]
    print("  -- my spectral candidates --")
    for name, f in cands:
        ok, w = holds(f)
        msg = "HOLDS n<=6" if ok else f"FAILS (e.g. n={w['n']},c3={w['c3']},H={w['H']},beta={w['beta']},ndev={w['ndev']},rho={w['rho']:.2f})"
        print(f"    [{'++' if ok else '  '}] {name}: {msg}")

sep("done -- survivors are candidate theorems; failures locate the decoupling")
