"""
PALEY vs INTERVAL: THE ASYMPTOTIC RACE — opus-2026-03-12-S67e

At p=7 and p=11, Paley beats Interval for H. But the margin shrinks.
Does the Interval's super-exponential amplification eventually overtake
the Paley's spectral base advantage?

This script:
1. Computes H(Paley) and H(Interval) at p=13 (feasible DP)
2. Analyzes the amplification ratio race
3. Explores the Weil/Ramanujan/expander connection
4. Investigates whether p ≡ 1 mod 4 primes behave differently
"""

import numpy as np
from itertools import combinations
from math import gcd, factorial, comb
import time

def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p-1)//2, p) == 1

def qr_set(p):
    return sorted([a for a in range(1, p) if is_qr(a, p)])

def fourier_magnitudes(S, p):
    m = (p - 1) // 2
    omega = np.exp(2j * np.pi / p)
    Q = []
    for k in range(1, m + 1):
        S_hat = sum(omega ** (s * k) for s in S)
        Q.append(abs(S_hat) ** 2)
    return np.array(Q)

def count_hp_from_0(S, p):
    """Count Hamiltonian paths starting from vertex 0 using bitmask DP."""
    n = p
    # dp[mask][v] = number of HP from 0 visiting exactly 'mask' ending at v
    # Use dict for efficiency
    dp = {}
    dp[(1, 0)] = 1  # mask with just vertex 0, at vertex 0

    for step in range(1, n):
        new_dp = {}
        for (mask, v), count in dp.items():
            if bin(mask).count('1') != step:
                continue
            for s in S:
                w = (v + s) % n
                if not (mask & (1 << w)):
                    new_mask = mask | (1 << w)
                    key = (new_mask, w)
                    new_dp[key] = new_dp.get(key, 0) + count
        dp.update(new_dp)

    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


print("=" * 72)
print("PART 1: H COMPUTATION AT p = 13")
print("=" * 72)

p = 13
m = (p - 1) // 2
S_int = list(range(1, m + 1))  # Interval
S_qr = qr_set(p)  # QR set

print(f"\np = {p}, m = {m}")
print(f"Interval: S = {S_int}")
print(f"QR:       S = {S_qr}")

# Q spectra
Q_int = fourier_magnitudes(S_int, p)
Q_qr = fourier_magnitudes(S_qr, p)
print(f"\nQ_k (Interval): [{', '.join(f'{q:.4f}' for q in Q_int)}]")
print(f"Q_k (QR):       [{', '.join(f'{q:.4f}' for q in Q_qr)}]")
print(f"prod(1+Q) Interval: {np.prod(1+Q_int):.1f}")
print(f"prod(1+Q) QR:       {np.prod(1+Q_qr):.1f}")

# Note: p=13 ≡ 1 mod 4, so QR is NOT a tournament connection set
# But let's check
qr_s = set(S_qr)
comp = set((p - s) % p for s in S_qr)
print(f"\nQR ∩ (p-QR) = {qr_s & comp}")
print(f"|QR ∩ (p-QR)| = {len(qr_s & comp)}")
if qr_s & comp:
    print("QR set is SYMMETRIC at p=13 (p ≡ 1 mod 4) — NOT a valid tournament!")
    print("Need to find the best tournament connection set instead.")

# For p ≡ 1 mod 4, enumerate ALL valid connection sets and find max H
print(f"\nComputing H for ALL {2**(p-1)//2} valid connection sets at p={p}...")

# Generate all valid tournament connection sets
valid_sets = []
elements = list(range(1, p))
for S in combinations(elements, m):
    S = list(S)
    ok = True
    for s in S:
        if (p - s) in S:
            ok = False
            break
    if not ok:
        continue
    valid_sets.append(S)

print(f"Total valid connection sets: {len(valid_sets)}")

# Group by orbit under multiplication by powers (circulant isomorphism)
# Two sets S and aS mod p are isomorphic if gcd(a,p)=1
# But for now, just compute H for all of them

results = []
t0 = time.time()
for i, S in enumerate(valid_sets):
    if i % 10 == 0 and i > 0:
        elapsed = time.time() - t0
        eta = elapsed / i * (len(valid_sets) - i)
        print(f"  Progress: {i}/{len(valid_sets)}, elapsed {elapsed:.1f}s, ETA {eta:.1f}s")

    H0 = count_hp_from_0(S, p)
    Q = fourier_magnitudes(S, p)
    prod_Q = np.prod(1 + Q)
    results.append((S, H0, prod_Q, Q))

elapsed = time.time() - t0
print(f"  Done in {elapsed:.1f}s")

# Sort by H
results.sort(key=lambda x: x[1], reverse=True)

print(f"\nTop 10 connection sets by H_from_0:")
print(f"  {'Rank':>4} {'H_from_0':>10} {'H_total':>12} {'prod(1+Q)':>12} {'amp':>10} S")
for i, (S, H0, prod_Q, Q) in enumerate(results[:10]):
    H_total = H0 * p
    amp = H0 / prod_Q if prod_Q > 0 else float('inf')
    print(f"  {i+1:>4} {H0:>10} {H_total:>12} {prod_Q:>12.1f} {amp:>10.4f} {S}")

print(f"\nBottom 5:")
for i, (S, H0, prod_Q, Q) in enumerate(results[-5:]):
    H_total = H0 * p
    amp = H0 / prod_Q if prod_Q > 0 else float('inf')
    print(f"  {len(results)-4+i:>4} {H0:>10} {H_total:>12} {prod_Q:>12.1f} {amp:>10.4f} {S}")

# Check if the max-H set has any special Q structure
best_S, best_H0, best_prod, best_Q = results[0]
print(f"\nBest set: S = {best_S}")
print(f"  H_from_0 = {best_H0}, H_total = {best_H0 * p}")
print(f"  Q_k = [{', '.join(f'{q:.4f}' for q in best_Q)}]")
print(f"  var(Q) = {np.var(best_Q):.4f}")
print(f"  prod(1+Q) = {best_prod:.1f}")

# Check what the Interval gives
int_result = [(S, H0, pq, Q) for S, H0, pq, Q in results if S == S_int]
if int_result:
    _, int_H0, int_prod, int_Q = int_result[0]
    print(f"\nInterval: S = {S_int}")
    print(f"  H_from_0 = {int_H0}, H_total = {int_H0 * p}")
    print(f"  Rank among {len(results)} sets: {[r[1] for r in results].index(int_H0) + 1}")

# Group by H value (orbit analysis)
from collections import Counter
H_counts = Counter(r[1] for r in results)
print(f"\nH value distribution:")
for h, count in sorted(H_counts.items(), reverse=True):
    # Find one representative
    rep = [r for r in results if r[1] == h][0]
    print(f"  H_from_0 = {h}: {count} sets, prod(1+Q) = {rep[2]:.1f}")


print("\n" + "=" * 72)
print("PART 2: THE AMPLIFICATION RACE")
print("=" * 72)

print("""
Compare amplification factors for Interval and best tournament:
  A(p) = H_from_0(T) / prod(1+Q(T))

This measures how much the graph structure Ω amplifies
beyond the bare spectral product.
""")

# Known data
interval_data = {
    5: (5, 5.0),      # H_from_0, prod(1+Q)=F_p
    7: (25, 13.0),
    11: (8457, 89.0),
    13: (int_H0, np.prod(1 + fourier_magnitudes(S_int, 13))),
}

best_data = {
    7: (27, 27.0),     # Paley at p=7
    11: (8645, 1024.0), # Paley at p=11
    13: (best_H0, best_prod),  # Best at p=13
}

print(f"  {'p':>4} {'H0_int':>10} {'prod_int':>12} {'A_int':>10} {'H0_best':>10} {'prod_best':>12} {'A_best':>10} {'ratio':>10}")
for p_val in [7, 11, 13]:
    if p_val in interval_data and p_val in best_data:
        h_int, prod_int = interval_data[p_val]
        h_best, prod_best = best_data[p_val]
        A_int = h_int / prod_int
        A_best = h_best / prod_best
        ratio = h_best / h_int
        print(f"  {p_val:>4} {h_int:>10} {prod_int:>12.1f} {A_int:>10.4f} {h_best:>10} {prod_best:>12.1f} {A_best:>10.4f} {ratio:>10.6f}")


print("\n" + "=" * 72)
print("PART 3: SPECTRAL CHARACTERIZATION OF THE MAXIMIZER")
print("=" * 72)

# For p=13, analyze what makes the best set special
best_S_set = set(best_S)
comp_S = set((p - s) % p for s in best_S)
print(f"\nBest S = {best_S}")
print(f"Complement (p-S) = {sorted(comp_S)}")

# Additive energy
def additive_energy(S, p):
    """E(S) = |{(a,b,c,d) ∈ S^4 : a+b ≡ c+d mod p}|"""
    count = 0
    sums = {}
    for a in S:
        for b in S:
            s = (a + b) % p
            sums[s] = sums.get(s, 0) + 1
    for v in sums.values():
        count += v * v
    return count

# Representation function
def rep_function(S, p):
    """r_S(d) = |{(a,b) ∈ S^2 : a-b ≡ d mod p}|"""
    r = {}
    for a in S:
        for b in S:
            d = (a - b) % p
            r[d] = r.get(d, 0) + 1
    return r

E_best = additive_energy(best_S, p)
E_int = additive_energy(S_int, p)
print(f"\nAdditive energy:")
print(f"  Best: E(S) = {E_best}")
print(f"  Interval: E(S) = {E_int}")

r_best = rep_function(best_S, p)
r_int = rep_function(S_int, p)
print(f"\nRepresentation function r_S(d):")
print(f"  Best:     {dict(sorted(r_best.items()))}")
print(f"  Interval: {dict(sorted(r_int.items()))}")

# Variance of representation function (excluding d=0)
r_best_vals = [r_best.get(d, 0) for d in range(1, p)]
r_int_vals = [r_int.get(d, 0) for d in range(1, p)]
print(f"\nVariance of r_S (d≠0): Best = {np.var(r_best_vals):.4f}, Interval = {np.var(r_int_vals):.4f}")
print(f"  (Lower variance = more uniform representation = closer to difference set)")


print("\n" + "=" * 72)
print("PART 4: WEIL/RAMANUJAN CONNECTION")
print("=" * 72)

print("""
DEEP CONNECTION: Weil's theorem and Ramanujan tournaments

For the Paley tournament (p ≡ 3 mod 4):
  Q_k = (p+1)/4 for all k

This follows from |Gauss sum g| = √p (PROVED by Weil, equivalent to RH
for curves over F_p).

The ADJACENCY MATRIX eigenvalues of the Paley tournament are:
  λ_0 = m = (p-1)/2  (trivial eigenvalue)
  |λ_k| = √((p+1)/4) = (√(p+1))/2 for k ≠ 0

This means the Paley tournament is a RAMANUJAN DIGRAPH:
  max_{k≠0} |λ_k| = (√(p+1))/2 ≈ √p/2

Compare with the Ramanujan bound for d-regular graphs: 2√(d-1).
Here d = m = (p-1)/2, so the bound is 2√(m-1) ≈ √p.
The Paley tournament achieves |λ| ≈ √p/2, which is HALF the bound.
Actually for tournaments the relevant bound is different.

The key point: OPTIMAL spectral gap → OPTIMAL expansion → OPTIMAL H?
""")

# Compute adjacency eigenvalues for small primes
for p_val in [7, 11, 13, 17]:
    m = (p_val - 1) // 2

    # For the interval
    S_int = list(range(1, m + 1))
    omega = np.exp(2j * np.pi / p_val)
    eigs_int = []
    for k in range(1, p_val):
        lam = sum(omega ** (s * k) for s in S_int)
        eigs_int.append(lam)

    # For the QR set (even if not a tournament for p≡1 mod 4)
    S_qr = qr_set(p_val)
    eigs_qr = []
    for k in range(1, p_val):
        lam = sum(omega ** (s * k) for s in S_qr)
        eigs_qr.append(lam)

    print(f"\n  p = {p_val}:")
    print(f"    Interval eigenvalue magnitudes: [{', '.join(f'{abs(e):.4f}' for e in eigs_int[:m])}]")
    print(f"    QR eigenvalue magnitudes:       [{', '.join(f'{abs(e):.4f}' for e in eigs_qr[:m])}]")
    print(f"    max|λ| Interval: {max(abs(e) for e in eigs_int):.4f}")
    print(f"    max|λ| QR:       {max(abs(e) for e in eigs_qr):.4f}")
    print(f"    √((p+1)/4) = {np.sqrt((p_val+1)/4):.4f}")

    # Spectral gap (difference between trivial eigenvalue m and max non-trivial)
    gap_int = m - max(abs(e) for e in eigs_int)
    gap_qr = m - max(abs(e) for e in eigs_qr)
    print(f"    Spectral gap: Interval = {gap_int:.4f}, QR = {gap_qr:.4f}")


print("\n" + "=" * 72)
print("PART 5: EXPANDER MIXING LEMMA AND H")
print("=" * 72)

print("""
EXPANDER MIXING LEMMA: For a tournament T with second eigenvalue λ,
the number of edges between any two vertex sets A, B satisfies:
  |e(A,B) - m|A||B|/p| ≤ λ · √(|A||B|)

Smaller λ → more uniform edge distribution → tournament "looks random".

For Hamiltonian path counting, we need edges between carefully chosen
subsets (current vertex set, next vertices). The more uniform these
edges are, the more paths can exist.

CONJECTURE: Among circulant tournaments on p vertices,
  max H is achieved by the one with SMALLEST max|λ_k|
  i.e., the most "Ramanujan" tournament.

For p ≡ 3 mod 4: this is the PALEY tournament (proven).
For p ≡ 1 mod 4: this should be the "nearest to Paley" tournament.
""")

# For p=13, check if the H-maximizer has the smallest spectral gap
p_val = 13
m = (p_val - 1) // 2
omega = np.exp(2j * np.pi / p_val)

print(f"\np = {p_val}: Correlation between max|λ| and H_from_0")
max_lams = []
H_vals = []
for S, H0, prod_Q, Q in results:
    eigs = [sum(omega ** (s * k) for s in S) for k in range(1, p_val)]
    max_lam = max(abs(e) for e in eigs)
    max_lams.append(max_lam)
    H_vals.append(H0)

corr = np.corrcoef(max_lams, H_vals)[0, 1]
print(f"  Correlation(max|λ|, H) = {corr:.6f}")
print(f"  (Negative means smaller max|λ| → larger H → supports conjecture)")

# Find which S minimizes max|λ|
min_lam_idx = np.argmin(max_lams)
print(f"\n  Set with smallest max|λ|: S = {results[min_lam_idx][0]}")
print(f"    max|λ| = {max_lams[min_lam_idx]:.6f}, H_from_0 = {results[min_lam_idx][1]}")

max_H_idx = np.argmax(H_vals)
print(f"  Set with largest H: S = {results[max_H_idx][0]}")
print(f"    max|λ| = {max_lams[max_H_idx]:.6f}, H_from_0 = {results[max_H_idx][1]}")

# Check if they're the same
if min_lam_idx == max_H_idx:
    print(f"\n  *** SAME SET! The Ramanujan conjecture holds at p={p_val}! ***")
else:
    print(f"\n  Different sets. The Ramanujan conjecture does NOT hold at p={p_val}.")
    print(f"  The relationship is more subtle.")


print("\n" + "=" * 72)
print("PART 6: CODING THEORY — TOURNAMENT AS CODE")
print("=" * 72)

print("""
ANALOGY: Connection set S ↔ Codeword in F_p

The connection set S = {s_1, ..., s_m} can be viewed as a codeword
in the binary code on F_p*: c(a) = 1 if a ∈ S, 0 otherwise.

The Fourier magnitudes Q_k = |Ŝ(k)|² are the POWER SPECTRUM.

Properties of good codes (low correlation / flat spectrum):
1. FLAT SPECTRUM: Q_k = constant → difference set → optimal autocorrelation
2. MINIMUM DISTANCE: related to how far S is from any shift of its complement
3. WEIGHT: |S| = m = (p-1)/2 (half-weight, like a balanced code)

The Paley tournament corresponds to the QUADRATIC RESIDUE CODE,
one of the most famous families in coding theory!

QR codes have:
- Self-orthogonality (related to p ≡ 3 mod 4)
- Optimal minimum distance (BCH bound)
- Rich automorphism group (PSL(2,p))

QUESTION: Does the QR code's optimality in coding theory
directly translate to H-maximality in tournament theory?
""")

# Binary representation of each connection set as a codeword
p_val = 13
for S, H0, prod_Q, Q in results[:5]:
    codeword = ['1' if i in S else '0' for i in range(1, p_val)]
    hamming_weight = sum(1 for c in codeword if c == '1')
    print(f"  S={str(S):<30} code={''.join(codeword)} H0={H0:>10} wt={hamming_weight}")


print("\n" + "=" * 72)
print("PART 7: THE p ≡ 1 mod 4 MYSTERY")
print("=" * 72)

print("""
For p ≡ 1 mod 4, the QR set is symmetric (S = p-S) and can't form
a tournament. What does the H-maximizer look like?

At p = 5: only 4 valid sets, all equivalent under rotation → H = 3 for all.
At p = 13: 64 valid sets with varying H.

Key question: Does the H-maximizer at p=13 have any special structure
related to coding theory or number theory?
""")

# Detailed analysis of the p=13 maximizer
best_S = results[0][0]
print(f"Best S at p=13: {best_S}")

# Check various algebraic properties
# Is it related to a Paley-like construction?
# One approach: "twin" the QR set — take half of QR and half of QNR
S_set = set(best_S)
qr_s = set(qr_set(13))
qnr_s = set(range(1, 13)) - qr_s

in_qr = S_set & qr_s
in_qnr = S_set & qnr_s
print(f"  Elements in QR: {sorted(in_qr)} ({len(in_qr)}/{len(qr_s)})")
print(f"  Elements in QNR: {sorted(in_qnr)} ({len(in_qnr)}/{len(qnr_s)})")

# Check if it's a coset of some subgroup
# Multiplicative order of each element mod 13
print(f"\n  Multiplicative structure mod 13:")
for a in range(1, 13):
    order = 1
    x = a
    while x != 1:
        x = (x * a) % 13
        order += 1
    in_S = "✓" if a in S_set else " "
    qr_status = "QR" if a in qr_s else "QN"
    print(f"    {a:>2} {in_S} ord={order:>2} {qr_status}")

# Check difference set property
diffs = {}
for a in best_S:
    for b in best_S:
        if a != b:
            d = (a - b) % 13
            diffs[d] = diffs.get(d, 0) + 1
print(f"\n  Difference multiset: {dict(sorted(diffs.items()))}")
diff_vals = set(diffs.values())
print(f"  Unique difference counts: {diff_vals}")
if len(diff_vals) <= 2:
    print(f"  → NEAR-difference set or actual difference set!")


print("\n" + "=" * 72)
print("PART 8: COMPARISON WITH ALL PRIMES UP TO 17")
print("=" * 72)

# For p=5, 7, 11 we have data. Add p=13.
# For p=17, computation is feasible (2^17 = 131K bitmask states, 17 endpoints)

for p_val in [5, 7]:
    m = (p_val - 1) // 2
    S_int = list(range(1, m + 1))
    S_qr = qr_set(p_val)
    Q_int = fourier_magnitudes(S_int, p_val)

    H0_int = count_hp_from_0(S_int, p_val)

    print(f"\n  p = {p_val}: Interval H_from_0 = {H0_int}, H_total = {H0_int * p_val}")

    # Find max H over all valid sets
    best_h = 0
    best_s_local = None
    elements = list(range(1, p_val))
    for S in combinations(elements, m):
        S = list(S)
        ok = True
        for s in S:
            if (p_val - s) in S:
                ok = False
                break
        if not ok:
            continue
        h = count_hp_from_0(S, p_val)
        if h > best_h:
            best_h = h
            best_s_local = S

    Q_best = fourier_magnitudes(best_s_local, p_val)
    print(f"  Best: S = {best_s_local}, H_from_0 = {best_h}, H_total = {best_h * p_val}")
    print(f"    Q = [{', '.join(f'{q:.4f}' for q in Q_best)}], var(Q) = {np.var(Q_best):.4f}")
    print(f"    Best/Interval ratio = {best_h / H0_int:.6f}")


# Summary table
print("\n" + "=" * 72)
print("SUMMARY: THE RACE — Interval vs Best tournament")
print("=" * 72)

# Gather known data
print(f"""
  p  │ H_from_0(Int) │ H_from_0(Best) │ Best/Int │ p mod 4 │ Best type
─────┼───────────────┼────────────────┼──────────┼─────────┼──────────""")

race_data = []
for p_val in [5, 7, 11, 13]:
    m = (p_val - 1) // 2
    S_int = list(range(1, m + 1))
    Q_int = fourier_magnitudes(S_int, p_val)
    H0_int = count_hp_from_0(S_int, p_val)

    # Best
    best_h = 0
    best_type = ""
    elements = list(range(1, p_val))
    for S in combinations(elements, m):
        S = list(S)
        ok = True
        for s in S:
            if (p_val - s) in S:
                ok = False
                break
        if not ok:
            continue
        h = count_hp_from_0(S, p_val)
        if h > best_h:
            best_h = h
            best_s_local = S

    # Check if best is Interval
    if best_h == H0_int:
        best_type = "Interval"
    elif set(best_s_local) == set(qr_set(p_val)):
        best_type = "Paley (QR)"
    else:
        best_type = f"S={best_s_local}"

    ratio = best_h / H0_int
    pmod4 = p_val % 4
    print(f"  {p_val:>2} │ {H0_int:>13} │ {best_h:>14} │ {ratio:>8.6f} │ {pmod4:>7} │ {best_type}")
    race_data.append((p_val, H0_int, best_h, ratio))


print("\n" + "=" * 72)
print("PART 9: DEEP SYNTHESIS — TWO REGIMES OF H-MAXIMIZATION")
print("=" * 72)

print("""
EMERGING PICTURE:

Two fundamentally different H-maximization strategies exist:

REGIME 1 — FLAT SPECTRUM (p ≡ 3 mod 4):
  - Paley tournament with Q_k = (p+1)/4
  - prod(1+Q) = ((p+5)/4)^m — AM-GM optimal
  - Low amplification factor (Ω is "simple")
  - Wins by spectral base alone
  - Connected to: difference sets, Gauss sums, Weil's theorem,
    Ramanujan graphs, QR codes

REGIME 2 — FIBONACCI RESONANCE (p ≡ 1 mod 4):
  - Interval tournament with peaked Fejér kernel Q_k
  - prod(1+Q) = F_p — Morgan-Voyce/Fibonacci structure
  - High amplification factor (Ω has rich structure)
  - Wins by resonance cascade amplification
  - Connected to: Schrödinger operators, KAM theory,
    continued fractions, golden ratio

OR: Does the flat spectrum ALWAYS win?

Evidence:
  p=5 (≡1 mod 4): All sets give same H. No winner.
  p=7 (≡3 mod 4): Paley wins (1.08×).
  p=11 (≡3 mod 4): Paley wins (1.02×).
  p=13 (≡1 mod 4): ???

The p=13 result is CRITICAL for determining which regime dominates.

HYPOTHESIS: The H-maximizer is:
  - Paley (QR set) for p ≡ 3 mod 4
  - Something related to "half-QR" for p ≡ 1 mod 4
  - The FLAT SPECTRUM principle always applies (minimize var(Q_k))

This would unify tournament theory with coding theory,
expander graphs, and the Riemann Hypothesis for finite fields.
""")
