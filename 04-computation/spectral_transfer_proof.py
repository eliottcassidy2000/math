#!/usr/bin/env python3
"""
spectral_transfer_proof.py — opus-2026-03-12-S67c

DEEP PROOF DIRECTION: Using the transfer matrix / spectral theory approach.

Key insight from the cross-field analysis:
- The Interval tournament's Q_k are eigenvalues of a Jacobi matrix
- The Morgan-Voyce polynomial is the characteristic polynomial
- prod(1+Q_k) = F_p via the transfer matrix T = [[3,-1],[1,0]]
- H(T) is determined by the representation profile r_S

NEW IDEA: Can we express H(T) as a SPECTRAL DETERMINANT?

If H = det(something involving Q_k), then maximizing H reduces
to a spectral optimization problem with known tools.

Also explore: Is H related to a REGULATOR or L-value of
the real cyclotomic field Q(cos 2π/p)?
"""

import numpy as np
from math import comb
from collections import Counter

def legendre(a, p):
    a = a % p
    if a == 0: return 0
    ls = pow(a, (p-1)//2, p)
    return -1 if ls == p-1 else ls

def interval_Q(p):
    m = (p-1)//2
    return np.array([
        (np.sin(m*np.pi*k/p)/np.sin(np.pi*k/p))**2
        for k in range(1, m+1)
    ])

def all_circulant_H(p):
    """Compute H for all 2^m circulant tournaments at prime p."""
    m = (p-1)//2
    pairs = []
    seen = set()
    for s in range(1, p):
        if s not in seen:
            pairs.append((s, p-s))
            seen.add(s)
            seen.add(p-s)

    results = []
    for bits in range(2**m):
        S = set()
        for i in range(m):
            if bits & (1 << i):
                S.add(pairs[i][0])
            else:
                S.add(pairs[i][1])

        # Build tournament
        A = [[0]*p for _ in range(p)]
        for i in range(p):
            for j in range(p):
                if i != j and (j-i) % p in S:
                    A[i][j] = 1

        # DP for H
        n = p
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)): continue
                if dp[mask][v] == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    if A[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full_mask = (1 << n) - 1
        H = sum(dp[full_mask][v] for v in range(n))

        # Q_k values
        Qs = np.array([abs(sum(np.exp(2j*np.pi*s*k/p) for s in S))**2 for k in range(1, m+1)])

        # E(S)
        sums = Counter()
        for a in S:
            for b in S:
                sums[(a+b) % p] += 1
        E = sum(v*v for v in sums.values())

        interval = set(range(1, m+1))
        paley_set = set(s for s in range(1, p) if legendre(s, p) == 1)
        name = ""
        if S == interval: name = "INT"
        if S == paley_set: name = "PAL"

        results.append({
            'S': S, 'H': H, 'Q': Qs, 'E': E, 'name': name
        })

    return results

print("=" * 70)
print("SPECTRAL TRANSFER PROOF DIRECTION — opus-2026-03-12-S67c")
print("=" * 70)
print()

# ============================================================
# PART 1: Can H be expressed as a spectral function of Q?
# ============================================================

print("PART 1: H AS A FUNCTION OF THE SPECTRUM {Q_k}")
print("=" * 70)
print()
print("Testing: H = F(Q_1,...,Q_m) for various functional forms")
print()

for p in [7, 11, 13]:
    m = (p-1)//2
    results = all_circulant_H(p)

    print(f"p={p}, m={m}:")

    # Group by Q spectrum (sorted)
    from collections import defaultdict
    by_spectrum = defaultdict(list)
    for r in results:
        key = tuple(np.round(np.sort(r['Q']), 6))
        by_spectrum[key].append(r)

    # Check if H is determined by sorted Q spectrum
    determined = True
    for key, group in by_spectrum.items():
        Hs = set(r['H'] for r in group)
        if len(Hs) > 1:
            determined = False
            print(f"  Q-spectrum {key[:3]}... has multiple H values: {Hs}")

    if determined:
        print(f"  H IS determined by sorted Q spectrum")

    # Test specific functional forms:
    for r in results:
        if r['name'] in ('INT', 'PAL'):
            Q = r['Q']
            H = r['H']

            # Various spectral functions
            f1 = np.prod(1 + Q)  # prod(1+Q) — Fibonacci for Interval
            f2 = np.prod(Q + 1) / np.prod(Q)  # = prod(1+1/Q)
            f3 = np.sum(np.log(1 + Q))  # log prod
            f4 = np.sum(Q**2)  # IPR numerator
            f5 = np.prod(Q + 2)  # prod(Q+2)
            f6 = np.exp(np.sum(Q / (1 + Q)))  # exponential form

            print(f"  {r['name']:>3}: H={H}, prod(1+Q)={f1:.1f}, prod(Q+2)={f5:.1f}")
            print(f"       sum(Q²)={f4:.1f}, sum(log(1+Q))={f3:.4f}")

    print()

    # Try to find the functional relationship
    # At p=7: only 2 distinct H values, 2 distinct Q spectra
    # At p=11: 4 distinct H values
    # Fit H = polynomial in spectral invariants

    # Spectral invariants: e_j (elementary symmetric polys of Q)
    # For Interval: e_j = C(m+j, 2j)
    # For others: different values

    print(f"  Spectral invariants for each H-class at p={p}:")
    for r in results:
        Q = np.sort(r['Q'])[::-1]
        # e_1 = sum(Q), e_2 = sum pairs, etc.
        e1 = np.sum(Q)
        e2 = sum(Q[i]*Q[j] for i in range(m) for j in range(i+1, m))

        if r['name'] or r == results[0]:
            print(f"    {r['name'] or 'other':>5}: H={r['H']}, e1={e1:.2f}, e2={e2:.2f}, E={r['E']}")

    print()

# ============================================================
# PART 2: The P_m(t) and H connection
# ============================================================

print("=" * 70)
print("PART 2: CHARACTERISTIC POLYNOMIAL VALUES AND H")
print("=" * 70)
print()
print("P_m(t) = prod(t - Q_k). Does P_m evaluated at special points give H?")
print()

for p in [7, 11, 13]:
    m = (p-1)//2
    results = all_circulant_H(p)

    for r in results:
        if r['name'] in ('INT', 'PAL') or r == results[0]:
            Q = r['Q']
            H = r['H']

            # P(t) = prod(t - Q_k)
            def P(t):
                return np.prod(t - Q)

            # Try various t values
            vals = {}
            for t_try in [0, 1, 2, -1, -2, m, m+1, p, p-1, p-2]:
                vals[t_try] = P(t_try)

            print(f"  {r['name'] or 'other':>5} (H={H}):")
            for t, v in vals.items():
                ratio = H / v if abs(v) > 0.01 else float('inf')
                print(f"    P({t:3d}) = {v:20.2f}  H/P = {ratio:12.4f}")
            print()

# ============================================================
# PART 3: H modular structure
# ============================================================

print("=" * 70)
print("PART 3: H MODULAR ARITHMETIC (H mod small primes)")
print("=" * 70)
print()

for p in [7, 11, 13]:
    m = (p-1)//2
    results = all_circulant_H(p)

    print(f"p={p}:")
    for r in results:
        H = r['H']
        E = r['E']
        name = r['name']
        if name or r == results[0]:
            mods = [H % q for q in [2, 3, 4, 5, 7, 8, p]]
            print(f"  {name or 'other':>5}: H={H:>14d}  mod(2,3,4,5,7,8,p) = {mods}")

    # All H values mod p
    H_mod_p = set(r['H'] % p for r in results)
    print(f"  All H mod {p}: {sorted(H_mod_p)}")
    print()

# ============================================================
# PART 4: The representation profile → H formula
# ============================================================

print("=" * 70)
print("PART 4: REPRESENTATION PROFILE AND SPECTRAL DATA")
print("=" * 70)
print()
print("kind-pasteur showed: sorted (r_S(0),...,r_S(p-1)) determines H.")
print("r_S(t) = |{(a,b) ∈ S×S : a+b ≡ t mod p}|")
print()
print("Connection: r_S is the AUTOCORRELATION of the indicator of S.")
print("By convolution theorem: r_hat(k) = |S_hat(k)|² = Q_k")
print("So r_S and {Q_k} carry the same information!")
print()
print("But the SORTED r_S and SORTED Q may carry different info")
print("because the sorting permutations may differ.")
print()

for p in [7, 11, 13]:
    m = (p-1)//2
    results = all_circulant_H(p)

    for r in results:
        if r['name'] in ('INT', 'PAL'):
            S = r['S']
            Q = r['Q']

            # Compute r_S(t) = sum of 1_{a in S} * 1_{t-a in S}
            r_profile = []
            for t in range(p):
                count = sum(1 for a in S if (t - a) % p in S)
                r_profile.append(count)

            # r_hat(k) should equal Q_k
            r_hat = np.fft.fft(r_profile)
            Q_from_r = np.abs(r_hat[1:m+1])**2 / (p**2) * p  # Normalization

            # Actually: r_hat(k) = |S_hat(k)|² where S_hat(k) = sum_{s in S} omega^{sk}
            # And Q_k = |S_hat(k)|²
            # So r_hat(k) = Q_k (no extra normalization needed for the discrete case)

            # Direct DFT of r_profile
            r_dft = np.array([
                sum(r_profile[t] * np.exp(-2j*np.pi*k*t/p) for t in range(p))
                for k in range(1, m+1)
            ])

            print(f"  {r['name']:>3}, p={p}:")
            print(f"    r_S = {r_profile}")
            print(f"    sorted r_S = {sorted(r_profile)}")
            print(f"    Q_k = {np.round(Q, 4)}")
            print(f"    r_hat = {np.round(np.abs(r_dft), 4)}")
            print(f"    Match? {np.allclose(np.abs(r_dft), Q)}")
            print()

# ============================================================
# PART 5: The big question — analytic H formula
# ============================================================

print("=" * 70)
print("PART 5: IS H EXPRESSIBLE AS A DETERMINANT?")
print("=" * 70)
print()
print("H(T) counts Hamiltonian paths. For circulant tournaments,")
print("H = p · h where h counts paths starting at vertex 0.")
print()
print("Determinantal formula: h can sometimes be expressed as")
print("a permanent or pfaffian of a matrix built from the Q_k.")
print()

for p in [7, 11]:
    m = (p-1)//2
    results = all_circulant_H(p)

    print(f"p={p}:")
    for r in results:
        if r['name'] == 'INT':
            H = r['H']
            Q = r['Q']

            # H/p
            h = H // p
            print(f"  H(Int) = {H} = {p} * {h}")
            print(f"  h = {h}")

            # Try: det(I + Q-diagonal matrix)
            # This gives prod(1 + Q_k) = F_p, which is NOT H
            det_IQ = np.prod(1 + Q)
            print(f"  det(I + diag(Q)) = prod(1+Q) = {det_IQ:.1f} ≠ {H}")

            # Try: permanent of circulant matrix?
            # The permanent of a circulant matrix is related to H...
            # Actually, H(T) for a tournament T with adjacency matrix A
            # equals perm(A restricted to Hamiltonian path submatrices)

            # Simpler: H relates to the IMMANANT of the tournament matrix
            print(f"  Note: H is related to immanants of the adjacency matrix")

            # Key spectral relationship:
            # For circulant tournaments, the adjacency eigenvalues are
            # lambda_k = S_hat(k) = |S_hat(k)| * exp(i*angle_k)
            # and |S_hat(k)|^2 = Q_k

            S = r['S']
            lambdas = [sum(np.exp(2j*np.pi*s*k/p) for s in S) for k in range(p)]

            # Check: is there a simple relationship between H and lambdas?
            # For permanent: perm(A) relates to sum over all permutations
            # of prod a_{i, sigma(i)}
            # For H: H = sum over all Hamiltonian paths (special permutations)

            # The Fourier spectrum of a circulant matrix:
            # eigenvalues = lambda_k = sum_{s in S} omega^{sk}
            # So the permanent = sum_sigma prod lambda_{sigma(k)}? No...

            # Actually for circulant matrices of size n:
            # perm(C) = sum over all subsets T of eigenvalues...
            # This is complex. Let me just check numerically.

            print(f"  Eigenvalues |λ|: {np.round(np.abs(lambdas[:m+1]), 4)}")
            print(f"  sum |λ|^2 = {np.sum(np.abs(lambdas)**2):.1f} = m*(m+1)/2*p = {m*(m+1)/2}")
            print()

print()
print("=" * 70)
print("PART 6: SYNTHESIS — THE PROOF ROADMAP")
print("=" * 70)
print()
print("Based on all discoveries, the proof should go:")
print()
print("STEP 1: For circulant tournaments, H is determined by")
print("the representation profile r_S (kind-pasteur, verified)")
print()
print("STEP 2: The representation profile is the inverse DFT of Q_k.")
print("So H = F(Q_1,...,Q_m) for some specific function F.")
print()
print("STEP 3: The Walsh decomposition H = sum h_hat[S] chi_S(sigma)")
print("gives H as a polynomial in the sign variables sigma_k.")
print("The degree-4 term dominates for p >= 13.")
print()
print("STEP 4: The degree-4 coefficient h_hat[{i,j,k,l}] depends on")
print("the 4-point additive correlation of S, which is measured by")
print("the zero-sum index W = |{±i±j±k±l ≡ 0 mod p}|.")
print()
print("STEP 5: For the Interval S = {1,...,m}:")
print("  - E(S) = maximal additive energy (Freiman bound achieved)")
print("  - The representation profile r_S is a TRIANGLE function")
print("  - The degree-4 Walsh surplus is positive because the")
print("    triangle profile concentrates 4-point correlations")
print()
print("STEP 6: The hyperplane condition (sum h_hat >= 0 for all")
print("flip sets) can be verified by showing the positive h_hat")
print("terms (from zero-sum quadruples) outweigh the negative ones.")
print()
print("STEP 7: For p >= 13, the degree-4 energy ratio E4/E2 > 1,")
print("so the degree-4 terms control the H landscape.")
print("Combined with the hyperplane condition, this proves")
print("Interval maximizes H.")
print()
print("MISSING PIECES:")
print("  A. Analytic formula for h_hat[{i,j,k,l}] in terms of")
print("     additive 4-point correlations (character sum identity)")
print("  B. Proof that zero-sum surplus is positive for all p >= 13")
print("  C. Proof that E4/E2 > 1 for all p >= 13")
print("  D. These reduce to character sum estimates (Weil-type bounds)")
print()
