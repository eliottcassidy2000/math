#!/usr/bin/env python3
"""
gauss_magnitude_derivation.py -- Derive Walsh magnitudes analytically from Gauss sums

The goal: PROVE the Walsh magnitude formulas rather than just observing them.

Setup: For a circulant tournament on Z_p with connection set S ⊂ {1,...,p-1},
|S| = m = (p-1)/2, the orientation is encoded as sigma ∈ {±1}^m where
sigma_j = +1 means gap j is in S, sigma_j = -1 means gap p-j is in S.

The cycle count c_k(sigma) is a function on the orientation hypercube.
Its Walsh-Hadamard coefficients are:
  h_hat_k[{a,b}] = (1/2^m) sum_sigma c_k(sigma) * sigma_a * sigma_b

Key insight: For a circulant, the adjacency eigenvalues are
  lambda_t(sigma) = sum_{s in S} omega^{st}
where omega = e^{2pi i/p}, and S depends on sigma.

For the Paley base sigma_0 = (1,...,1) (gaps {1,...,m}), flipping chord j means
replacing gap j with gap p-j in S. The eigenvalue changes:
  lambda_t(flip j) = lambda_t(sigma_0) + omega^{(p-j)t} - omega^{jt}
                    = lambda_t(sigma_0) - 2i sin(2pi j t / p)

The trace formula: tr(A^k) = sum_t lambda_t^k gives closed walks of length k.
For k=3: c_3 = tr(A^3)/3 (exact for p >= 7).
For k=p: c_p involves simple Hamiltonian cycles (complex).
For general k: c_k involves inclusion-exclusion corrections.

HOWEVER: for Walsh coefficients, the KEY observation is that:
  Walsh of c_k = Walsh of (1/k) sum_t lambda_t^k  (when c_k = tr(A^k)/k exactly)

This holds for k=3. Let's check if the Walsh of tr(A^k)/k matches the Walsh of c_k
for k=5, 7, etc. (It won't be exact due to non-simple walks, but the SIGN might match.)

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def resonance_level(a, b, p):
    for q in range(1, p, 2):
        if (q * a - b) % p == 0 or (q * a + b) % p == 0:
            return q
        if q > 1 and ((a - q * b) % p == 0 or (a + q * b) % p == 0):
            return q
    return p


def compute_eigenvalues(S, p):
    """Compute circulant adjacency eigenvalues lambda_t = sum_{s in S} omega^{st}."""
    omega = cmath.exp(2j * cmath.pi / p)
    evals = []
    for t in range(p):
        val = sum(omega**(s * t) for s in S)
        evals.append(val)
    return evals


def trace_k(evals, k):
    """Compute tr(A^k) = sum_t lambda_t^k."""
    return sum(lam**k for lam in evals).real


def main():
    p = 11
    m = (p - 1) // 2
    pairs = [(j, p - j) for j in range(1, m + 1)]
    n_orient = 1 << m

    print("=" * 70)
    print(f"GAUSS SUM DERIVATION OF WALSH MAGNITUDES at p={p}")
    print("=" * 70)

    # For each orientation sigma, compute eigenvalues and tr(A^k)
    # Then compute Walsh coefficients of tr(A^k)/k

    all_traces = {}  # bits -> {k: tr(A^k)}

    for bits in range(n_orient):
        S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                   for i in range(m))
        evals = compute_eigenvalues(S, p)
        traces = {}
        for k in range(3, p + 1, 2):
            traces[k] = trace_k(evals, k) / k
        all_traces[bits] = traces

    # Compute Walsh coefficients of trace_k/k
    sigma = {}
    for bits in range(n_orient):
        sigma[bits] = tuple(1 if bits & (1 << i) else -1 for i in range(m))

    chord_pairs = [(a, b) for a in range(m) for b in range(a + 1, m)]

    print("\nDEGREE-2 WALSH OF tr(A^k)/k vs EXACT c_k:")
    print("(If these match, the trace formula captures the Walsh structure)")

    # Known exact c_k Walsh coefficients from ck_walsh_signs.py output
    exact_data = {
        5: {(0, 2): -5.5, (0, 3): -5.5, (1, 2): 5.5, (1, 4): 5.5, (3, 4): -5.5,
            (0, 1): 0, (0, 4): 0, (1, 3): 0, (2, 3): 0, (2, 4): 0},
        7: {(0, 1): -20.625, (0, 2): -15.125, (0, 3): -15.125, (0, 4): 20.625,
            (1, 2): 15.125, (1, 3): -20.625, (1, 4): 15.125,
            (2, 3): 20.625, (2, 4): 20.625, (3, 4): -15.125},
        9: {(0, 1): -5.5, (0, 2): -13.75, (0, 3): -13.75, (0, 4): 5.5,
            (1, 2): 13.75, (1, 3): -5.5, (1, 4): 13.75,
            (2, 3): 5.5, (2, 4): 5.5, (3, 4): -13.75},
        11: {(0, 1): -9.625, (0, 2): 48.125, (0, 3): 48.125, (0, 4): 9.625,
             (1, 2): -48.125, (1, 3): -9.625, (1, 4): -48.125,
             (2, 3): 9.625, (2, 4): 9.625, (3, 4): 48.125}
    }

    for k in range(3, p + 1, 2):
        print(f"\n  c_{k}:")
        for a, b in chord_pairs:
            ga, gb = a + 1, b + 1
            q = resonance_level(ga, gb, p)

            # Walsh of tr(A^k)/k
            h_hat_trace = 0
            for bits in range(n_orient):
                sig = sigma[bits]
                h_hat_trace += all_traces[bits][k] * sig[a] * sig[b]
            h_hat_trace /= n_orient

            # Exact (from data)
            h_hat_exact = exact_data.get(k, {}).get((a, b), 0)

            if abs(h_hat_trace) < 0.001 and abs(h_hat_exact) < 0.001:
                continue

            match = abs(h_hat_trace - h_hat_exact) < 0.01
            diff = h_hat_trace - h_hat_exact

            print(f"    ({a},{b}) q={q}: trace={h_hat_trace:>10.4f}, "
                  f"exact={h_hat_exact:>10.4f}, diff={diff:>10.4f} "
                  f"{'MATCH' if match else 'DIFFERS'}")

    # Analytical derivation of Walsh[trace_k/k][{a,b}]
    print(f"\n\n{'='*70}")
    print("ANALYTICAL WALSH FORMULA FROM EIGENVALUE PERTURBATION")
    print("=" * 70)

    # For the Paley base S_0 = {1,2,...,m}, eigenvalue:
    #   lambda_t = sum_{s=1}^m omega^{st}
    # When we flip chord j (replace j with p-j in S):
    #   delta_j(t) = omega^{(p-j)t} - omega^{jt} = omega^{-jt} - omega^{jt} = -2i sin(2pi jt/p)
    #
    # For orientation sigma: S(sigma) = base set with some chords flipped.
    #   lambda_t(sigma) = lambda_t(base) + sum_{j: sigma_j=-1} delta_j(t)
    # But we want the GENERAL form. Let's define:
    #   lambda_t(sigma) = sum_{s in S(sigma)} omega^{st}
    # where S(sigma) has gap j if sigma_j=+1, gap p-j if sigma_j=-1.
    #
    # lambda_t(sigma) = sum_j omega^{sigma_j * j * t}   (where sigma_j=+1 gives omega^{jt},
    #                                                     sigma_j=-1 gives omega^{-jt})
    # Wait: gap j -> omega^{jt}, gap p-j -> omega^{(p-j)t} = omega^{-jt}.
    # So lambda_t(sigma) = sum_{j=1}^m omega^{sigma_j * j * t}
    # but sigma_j ∈ {+1, -1}.
    #
    # More carefully: if sigma_j = +1, the gap is j and contributes omega^{jt}.
    #                 if sigma_j = -1, the gap is p-j and contributes omega^{(p-j)t} = omega^{-jt}.
    # So: lambda_t(sigma) = sum_{j=1}^m omega^{j * t} if sigma_j=+1
    #                        + omega^{-j*t} if sigma_j=-1
    #   = sum_{j=1}^m [sigma_j=+1] * omega^{jt} + [sigma_j=-1] * omega^{-jt}
    #   = sum_j (1+sigma_j)/2 * omega^{jt} + (1-sigma_j)/2 * omega^{-jt}
    #   = sum_j (omega^{jt} + omega^{-jt})/2 + sigma_j * (omega^{jt} - omega^{-jt})/2
    #   = sum_j cos(2pi jt/p) + i * sigma_j * sin(2pi jt/p)
    #
    # So lambda_t(sigma) = C(t) + i * sum_j sigma_j * sin(2pi jt/p)
    # where C(t) = sum_j cos(2pi jt/p) and the imaginary part depends on sigma.
    #
    # For t=0: lambda_0 = m (constant, independent of sigma).
    # For t != 0:
    #   C(t) = sum_{j=1}^m cos(2pi jt/p)
    #   sigma-dependent: D(t,sigma) = sum_j sigma_j * sin(2pi jt/p)
    #   lambda_t = C(t) + i * D(t,sigma)

    omega_root = cmath.exp(2j * cmath.pi / p)

    # Precompute C(t), sin(2pi jt/p)
    C = [0.0] * p  # C(t)
    sin_jt = [[0.0]*p for _ in range(m)]  # sin_jt[j][t] = sin(2pi (j+1) t / p)

    for t in range(p):
        for j in range(m):
            angle = 2 * math.pi * (j + 1) * t / p
            C[t] += math.cos(angle)
            sin_jt[j][t] = math.sin(angle)

    print(f"\n  C(t) = sum_j cos(2pi jt/p) for t=0..{p-1}:")
    for t in range(p):
        print(f"    C({t}) = {C[t]:.6f}")

    # The Gauss sum connects to C(t):
    # For Paley S = QR: S_hat(t) = (-1 + chi(t)*g)/2 for t != 0
    # where g = sum chi(s) omega^s = Gauss sum, |g|^2 = p.
    # C(t) = Re(lambda_t at sigma_0) where sigma_0 = (+1,...,+1) = base S = {1,...,m}
    # Note: this is NOT the Paley tournament!

    # At base sigma_0 = {1,...,m}: lambda_t = sum_{j=1}^m omega^{jt}
    # = omega^t * (omega^{mt} - 1) / (omega^t - 1)  for t != 0
    # Since m = (p-1)/2: omega^{mt} = omega^{(p-1)t/2}
    # = omega^{-t/2} (since omega^p = 1)
    # This is a geometric sum.

    # Actually: sum_{j=1}^m omega^{jt} = omega^t (1 - omega^{mt}) / (1 - omega^t)
    # = omega^t (1 - omega^{(p-1)t/2}) / (1 - omega^t)

    # For t != 0: this equals (-1 + delta_t) / 2 where delta_t depends on chi(t)
    # via the Gauss sum relationship.

    # Let's compute and verify:
    print(f"\n  Eigenvalues at base sigma_0 = {{1,...,{m}}}:")
    S_base = list(range(1, m + 1))
    evals_base = compute_eigenvalues(S_base, p)
    for t in range(p):
        lam = evals_base[t]
        chi_t = legendre(t, p)
        g = sum(legendre(s, p) * omega_root**(s) for s in range(1, p))
        pred = (-1 + chi_t * g) / 2 if t > 0 else complex(m, 0)
        match = abs(lam - pred) < 0.001
        print(f"    t={t}: lambda={lam.real:>8.4f}+{lam.imag:>8.4f}i, "
              f"chi(t)={chi_t:+d}, pred={pred.real:>8.4f}+{pred.imag:>8.4f}i, "
              f"{'MATCH' if match else 'NO'}  |lam|^2={abs(lam)**2:.4f}")

    # The base S = {1,...,m} is NOT the QR set, so the Gauss sum formula
    # S_hat(t) = (-1 + chi(t)*g)/2 does NOT apply to the base.
    # The Gauss sum formula applies only to S = QR (Paley).

    # For the base, S_hat_base(t) = sum_{j=1}^m omega^{jt} is a geometric series:
    # = (omega^t - omega^{(m+1)t}) / (1 - omega^t)
    # = omega^t * (1 - omega^{mt}) / (1 - omega^t)

    # Now: the WALSH COEFFICIENT of lambda_t^k w.r.t. {a,b} involves:
    # lambda_t(sigma)^k = (C(t) + i * D(t,sigma))^k
    # where D(t,sigma) = sum_j sigma_j sin_j(t).
    #
    # Walsh[lambda_t^k][{a,b}] = (1/2^m) sum_sigma lambda_t(sigma)^k * sigma_a * sigma_b
    #
    # Since lambda_t = C(t) + i*D(t,sigma), we need to expand the k-th power
    # and extract the degree-2 Walsh in sigma.

    # Key simplification: D(t,sigma) = sum_j sigma_j * s_j(t) where s_j(t) = sin(2pi(j+1)t/p).
    # lambda_t^k depends on sigma ONLY through D(t,sigma).
    # lambda_t = C + iD where C = C(t) is constant and D depends on sigma.
    #
    # Walsh[lambda_t^k][{a,b}] = Walsh[(C + iD)^k][{a,b}]

    # For degree-2 Walsh of a function f(D) where D is linear in sigma:
    # Walsh_2[f(D)][{a,b}] = s_a(t) * s_b(t) * E[f''(D)]
    # This is because sigma_a sigma_b = partial^2/partial(s_a s_b) when D is linear.
    # More precisely: D = sum_j sigma_j * s_j, and:
    # Walsh[(C+iD)^k][{a,b}] = (i^2) * s_a * s_b * k(k-1) * <(C+iD)^{k-2}>_sigma
    # + correction terms from higher-order interactions.
    # Actually this is only exact at leading order.

    # Let me compute EXACTLY: Walsh[(C+iD)^k][{a,b}] for small m.
    # With m=5 (p=11), sigma has 2^5=32 values.

    print(f"\n\n{'='*70}")
    print("WALSH OF EIGENVALUE POWERS: Walsh[lambda_t^k][{{a,b}}]")
    print("=" * 70)

    # For each t, compute Walsh[lambda_t^k][{a,b}]
    for k in [3, 5, 7, 9, 11]:
        print(f"\n  k={k}:")
        for a, b in [(0, 2), (0, 1)]:  # q=3 and q=5 representatives
            ga, gb = a + 1, b + 1
            q = resonance_level(ga, gb, p)

            # Compute sum over t of Walsh[lambda_t^k][{a,b}]
            total_walsh = 0
            contributions = []

            for t in range(p):
                # Walsh[lambda_t^k][{a,b}] = (1/2^m) sum_sigma lambda_t(sigma)^k sigma_a sigma_b
                walsh_t = 0
                for bits in range(n_orient):
                    sig = sigma[bits]
                    S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                               for i in range(m))
                    lam = sum(omega_root**(s * t) for s in S)
                    walsh_t += (lam**k).real * sig[a] * sig[b]
                walsh_t /= n_orient

                total_walsh += walsh_t
                if abs(walsh_t) > 0.01:
                    contributions.append((t, walsh_t))

            # Should equal k * h_hat_ck[{a,b}] since c_k = tr(A^k)/k for simple cycles
            # (only exact for k=3)
            print(f"    ({a},{b}) q={q}: sum_t Walsh[lam_t^{k}] / {k} = {total_walsh/k:.4f}")
            if contributions:
                # Group contributions by chi(t)
                chi_pos = sum(w for t, w in contributions if legendre(t, p) == 1)
                chi_neg = sum(w for t, w in contributions if legendre(t, p) == -1)
                chi_zero = sum(w for t, w in contributions if legendre(t, p) == 0)
                print(f"      QR contrib: {chi_pos/k:.4f}, QNR contrib: {chi_neg/k:.4f}, "
                      f"t=0: {chi_zero/k:.4f}")

                # Check: does chi(t) determine the sign of the contribution?
                for t, w in contributions[:6]:
                    chi_t = legendre(t, p)
                    print(f"      t={t}: chi(t)={chi_t:+d}, walsh_contrib={w:.4f}")

    # Now derive the ANALYTICAL form.
    # For c_3 = (1/3) sum_t lambda_t^3, the Walsh coefficient is:
    # h_hat_c3[{a,b}] = (1/3) sum_t Walsh[lambda_t^3][{a,b}]
    #
    # lambda_t(sigma) = C(t) + i * sum_j sigma_j s_j(t)
    # lambda_t^3 = (C + i*D)^3 = C^3 + 3C^2(iD) + 3C(iD)^2 + (iD)^3
    #            = C^3 - 3CD^2 + i(3C^2D - D^3)
    # We need the real part (since tr(A^k) is real):
    # Re(lambda_t^3) = C^3 - 3CD^2
    #
    # Walsh[C^3 - 3CD^2][{a,b}] = -3C * Walsh[D^2][{a,b}]
    # Since C is constant w.r.t. sigma, and D = sum_j sigma_j s_j(t).
    # Walsh[D^2][{a,b}] = Walsh[(sum_j sigma_j s_j)^2][{a,b}]
    #                    = 2 s_a(t) s_b(t)
    # (by the Walsh-Hadamard identity: Walsh[sigma_i sigma_j][{a,b}] = delta_{ij,ab})
    #
    # So: Walsh[Re(lambda_t^3)][{a,b}] = -3C(t) * 2 * s_a(t) * s_b(t) = -6 C(t) s_a(t) s_b(t)
    #
    # Therefore: h_hat_c3[{a,b}] = (1/3) sum_t (-6 C(t) s_a(t) s_b(t))
    #                             = -2 sum_t C(t) sin(2pi(a+1)t/p) sin(2pi(b+1)t/p)

    print(f"\n\n{'='*70}")
    print("ANALYTICAL FORMULA FOR c_3 WALSH")
    print("=" * 70)

    for a, b in chord_pairs:
        ga, gb = a + 1, b + 1
        q = resonance_level(ga, gb, p)
        chi_ab = legendre(ga * gb, p)

        # Analytical formula: -2 sum_t C(t) sin(ga*t*2pi/p) sin(gb*t*2pi/p)
        analytical = 0
        for t in range(p):
            analytical += C[t] * sin_jt[a][t] * sin_jt[b][t]
        analytical *= -2

        print(f"  ({a},{b}) q={q}: analytical = {analytical:.6f} "
              f"(expected 0 since c_3 is constant)")

    # c_3 Walsh = 0 for all pairs (since c_3 = C(p,3)/3 * 2 = 55 is constant).
    # This means sum_t C(t) sin(ga*t) sin(gb*t) = 0 for all (a,b).
    # This is an orthogonality relation!

    print(f"\n  Result: sum_t C(t) sin(ga*t) sin(gb*t) = 0 for all (ga,gb)")
    print(f"  This is EXPECTED: c_3 = p*(p-1)*(p-2)/6 is orientation-independent")

    # For k=5, the expansion is more complex:
    # lambda^5 = (C+iD)^5
    # Re(lambda^5) = C^5 - 10C^3D^2 + 5CD^4
    # Walsh[D^2][{a,b}] = 2*s_a*s_b
    # Walsh[D^4][{a,b}] = ? (involves 4th-order Walsh)
    # D^4 = (sum_j sigma_j s_j)^4 = sum of products sigma_j sigma_k sigma_l sigma_m * s_j s_k s_l s_m
    # Walsh_2[D^4][{a,b}] = sum over (j,k,l,m) with exactly {a,b} appearing as the
    # "unpaired" indices... This is combinatorially complex.
    # Actually: Walsh_2[D^4][{a,b}] = 4! * (# ways to pair 4 indices with 2 left as a,b)
    # Hmm, let me think differently.

    # D^4 = (sum_j sigma_j s_j)^4
    # The degree-2 Walsh coefficient at {a,b} picks out terms where the product of
    # sigmas equals sigma_a sigma_b. This happens when:
    # - Two of the four indices are a and b, and the other two are the same (cancel)
    # This gives: C(4,2) * s_a * s_b * (sum_j s_j^2) * (ways) - let me compute exactly.

    # Walsh_2[D^4][{a,b}] = (1/2^m) sum_sigma D^4 sigma_a sigma_b
    # D = sum_j sigma_j s_j(t). D^4 = sum_{j1,j2,j3,j4} sigma_{j1}...sigma_{j4} prod s
    # sigma_{j1} sigma_{j2} sigma_{j3} sigma_{j4} sigma_a sigma_b has average = 1 iff
    # the multiset {j1,j2,j3,j4,a,b} has each index appearing an even number of times.
    # (Because E[sigma_i^{n_i}] = 1 if n_i even, 0 if odd.)

    # So we need partitions of {j1,j2,j3,j4} such that combined with {a,b},
    # every index has even multiplicity.
    # Case 1: j1=j2=a, j3=j4=b. Count: C(4,2) = 6. Term: s_a^2 s_b^2.
    # Case 2: j1=a, j2=b, j3=j4=c (c != a,b). Count: C(4,1)*C(3,1) * sum_c = 12 * sum_c s_c^2. Wait.
    # Actually j1,j2,j3,j4 can repeat. Let me enumerate:
    # Need {j1,j2,j3,j4} + {a,b} = even multiplicities everywhere.
    # So {j1,j2,j3,j4} must contain a odd times and b odd times (to make them even when combined).
    # Since we have 4 indices, a appears 1 or 3 times, b appears 1 or 3 times, and the rest even.
    # Possible:
    # a x 1, b x 1: remaining 2 indices = c,c for some c. -> a,b,c,c
    #   Ways: C(4,1)*C(3,1)*1 = 12 (choose position for a, position for b, rest=c,c)
    #   Hmm, actually more careful: # orderings of (a,b,c,c) = 4!/(1!1!2!) = 12.
    #   Term: s_a s_b s_c^2. Sum over c: 12 * s_a s_b * sum_{c != a, c != b} s_c^2
    #   But wait, c can equal a or b: if c=a: ordering is (a,a,a,b) which has a appearing 3 times.
    #   Hmm I need to be more careful.

    # Let n_j = # times index j appears in {j1,j2,j3,j4}. n_j + (1 if j=a) + (1 if j=b) must be even.
    # For j != a and j != b: n_j must be even.
    # For j = a (if a != b): n_a must be odd.
    # For j = b (if a != b): n_b must be odd.
    # Total: sum n_j = 4. n_a = odd (1 or 3), n_b = odd (1 or 3), n_c = even for others.
    #
    # Case A: n_a=1, n_b=1. Remaining 2 indices have even multiplicities.
    #   -> one index c with n_c=2. (Or 2 indices each with 0, but that's already the case.)
    #   Multinomial: 4! / (1! 1! 2!) = 12. Sum over c: s_c^2 for each c (including c=a or c=b?)
    #   If c=a: n_a would be 1+2=3? No: c is a separate index that can equal a.
    #   Actually I should just sum over all c: the term is sum_c 12 * s_a * s_b * s_c^2.
    #   But wait, c ranges over all m indices, including a and b.
    #   n_a = 1 (from the "a" appearance) + n_c(c=a) = 1 + 2 = 3 (odd, OK!).
    #   So c=a is allowed and gives n_a=3, n_b=1. This is covered in case A.
    #
    # Case B: n_a=3, n_b=1. Remaining 0 indices with even multiplicities.
    #   Multinomial: 4! / (3! 1!) = 4. Term: s_a^3 * s_b.
    #
    # Case C: n_a=1, n_b=3. Similar: 4 * s_a * s_b^3.
    #
    # Case D: n_a=3, n_b=3. sum = 6 > 4, impossible.

    # So: Walsh_2[D^4][{a,b}] = sum over multisets (j1,j2,j3,j4) satisfying the parity condition:
    #   product of s_{j_i}
    #
    # = 12 * s_a * s_b * (sum_c s_c^2)  [Case A with c != a, c != b, n_c=2]
    # + 12 * s_a * s_b * s_a^2          [Case A with c = a]
    # + 12 * s_a * s_b * s_b^2          [Case A with c = b]
    # + 4 * s_a^3 * s_b                 [Case B]
    # + 4 * s_a * s_b^3                 [Case C]
    #
    # Wait, I'm double-counting. Let me redo this.
    # Actually, case A with c=a gives (j1,j2,j3,j4) where a appears 3 times and b once.
    # That's the SAME as Case B. So I should NOT include c=a in Case A separately.
    #
    # Correct enumeration:
    # We need n_a odd, n_b odd, all others even, sum n = 4.
    # Options: (n_a, n_b, rest) =
    #   (1, 1, 2): one other index c with n_c=2. #multisets: for each c != a,b: 1. Total: m-2.
    #     Multinomial coeff: 4!/(1!1!2!) = 12 for each c.
    #   (3, 1, 0): just a^3 b. Multinomial: 4!/(3!1!) = 4.
    #   (1, 3, 0): just a b^3. Multinomial: 4!/(1!3!) = 4.
    #
    # Walsh_2[D^4][{a,b}] = 12 * s_a * s_b * sum_{c != a, c != b} s_c^2
    #                      + 4 * s_a^3 * s_b + 4 * s_a * s_b^3
    #                      = s_a * s_b * [12 * sum_{c != a,b} s_c^2 + 4*s_a^2 + 4*s_b^2]
    #                      = s_a * s_b * [12 * S2 - 8*s_a^2 - 8*s_b^2 + 4*s_a^2 + 4*s_b^2]
    #   where S2 = sum_j s_j^2
    # Hmm wait, let me just compute directly:
    # = s_a s_b [12 sum_{c != a, != b} s_c^2 + 4 s_a^2 + 4 s_b^2]
    # = s_a s_b [12 S2 - 12 s_a^2 - 12 s_b^2 + 4 s_a^2 + 4 s_b^2]
    # = s_a s_b [12 S2 - 8 s_a^2 - 8 s_b^2]

    # OK this is getting complicated but is solvable. Let me verify numerically.

    print(f"\n\n{'='*70}")
    print("VERIFYING ANALYTICAL WALSH FORMULA FOR tr(A^k)/k")
    print("=" * 70)

    # For k=3, we derived: h_hat = -2 sum_t C(t) s_a(t) s_b(t)
    # For k=5: h_hat = (1/5) sum_t W5_ab(t) where W5_ab involves C(t), s_a(t), s_b(t), S2(t)

    # Let me compute S2(t) = sum_j sin^2(2pi(j+1)t/p)
    S2 = [0.0] * p
    for t in range(p):
        S2[t] = sum(sin_jt[j][t]**2 for j in range(m))

    # Formula for Re(lambda^5):
    # Re(lambda^5) = Re((C+iD)^5) = C^5 - 10C^3D^2 + 5CD^4
    # Walsh_2[C^5][{a,b}] = 0 (constant)
    # Walsh_2[-10C^3 D^2][{a,b}] = -10 C^3 * Walsh_2[D^2] = -10 C^3 * 2 s_a s_b = -20 C^3 s_a s_b
    # Walsh_2[5C D^4][{a,b}] = 5C * Walsh_2[D^4][{a,b}]
    #   = 5C * s_a s_b * [12 S2 - 8 s_a^2 - 8 s_b^2]

    # Total: Walsh_2[Re(lam^5)][{a,b}] per t:
    # = -20 C(t)^3 s_a(t) s_b(t) + 5 C(t) s_a(t) s_b(t) [12 S2(t) - 8 s_a(t)^2 - 8 s_b(t)^2]
    # = s_a(t) s_b(t) [-20 C^3 + 60 C S2 - 40 C s_a^2 - 40 C s_b^2]
    # = s_a(t) s_b(t) * C(t) * [-20 C^2 + 60 S2 - 40 s_a^2 - 40 s_b^2]

    # h_hat_{tr(A^5)/5}[{a,b}] = (1/5) sum_t of the above

    for a, b in [(0, 2), (0, 1)]:
        ga, gb = a + 1, b + 1
        q = resonance_level(ga, gb, p)

        # Numerical check
        numerical = 0
        for bits in range(n_orient):
            sig = sigma[bits]
            S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                       for i in range(m))
            evals = compute_eigenvalues(S, p)
            numerical += sum(ev**5 for ev in evals).real / 5 * sig[a] * sig[b]
        numerical /= n_orient

        # Analytical k=5
        analytical_5 = 0
        for t in range(p):
            sa = sin_jt[a][t]
            sb = sin_jt[b][t]
            Ct = C[t]
            S2t = S2[t]
            term = sa * sb * Ct * (-20 * Ct**2 + 60 * S2t - 40 * sa**2 - 40 * sb**2)
            analytical_5 += term
        analytical_5 /= 5

        print(f"  ({a},{b}) q={q}: numerical={numerical:.4f}, analytical={analytical_5:.4f}, "
              f"match={abs(numerical - analytical_5) < 0.01}")

    # For k=3 (as sanity check)
    for a, b in [(0, 2), (0, 1)]:
        ga, gb = a + 1, b + 1
        q = resonance_level(ga, gb, p)
        analytical_3 = 0
        for t in range(p):
            analytical_3 += -2 * C[t] * sin_jt[a][t] * sin_jt[b][t]
        analytical_3 /= 3  # divide by 3 for c_3 = tr/3
        numerical = 0
        for bits in range(n_orient):
            sig = sigma[bits]
            S = sorted(pairs[i][0] if bits & (1 << i) else pairs[i][1]
                       for i in range(m))
            evals = compute_eigenvalues(S, p)
            numerical += sum(ev**3 for ev in evals).real / 3 * sig[a] * sig[b]
        numerical /= n_orient
        print(f"  k=3: ({a},{b}) q={q}: numerical={numerical:.6f}, analytical={analytical_3:.6f}")

    # Now: the key SIGN PATTERN question.
    # The sign of h_hat_{c_k}[{a,b}] depends on the sum:
    # sum_t f_k(C(t), S2(t), s_a(t), s_b(t)) * s_a(t) * s_b(t)
    # The factor s_a(t) * s_b(t) = sin(2pi ga t/p) sin(2pi gb t/p)
    # Using product-to-sum: = (1/2)[cos(2pi(ga-gb)t/p) - cos(2pi(ga+gb)t/p)]
    # The (ga-gb) and (ga+gb) determine the Fourier mode.
    # For resonance q: q*ga ≡ ±gb mod p, so ga-gb and ga+gb relate to multiples of ga/q.
    # This connects the Walsh sign to the Gauss sum structure.

    print(f"\n\n{'='*70}")
    print("SIN PRODUCT STRUCTURE: s_a(t) * s_b(t) BY CHI(t)")
    print("=" * 70)

    for a, b in chord_pairs:
        ga, gb = a + 1, b + 1
        q = resonance_level(ga, gb, p)
        chi_ab = legendre(ga * gb, p)

        # Compute s_a(t) * s_b(t) for each t, group by chi(t)
        qr_sum = 0
        qnr_sum = 0
        for t in range(1, p):
            prod = sin_jt[a][t] * sin_jt[b][t]
            if legendre(t, p) == 1:
                qr_sum += prod
            else:
                qnr_sum += prod

        print(f"  ({a},{b}) q={q}, chi(ab)={chi_ab:+d}: "
              f"QR_sum={qr_sum:.6f}, QNR_sum={qnr_sum:.6f}, "
              f"diff={qr_sum-qnr_sum:.6f}, ratio={'N/A' if abs(qnr_sum)<0.001 else f'{qr_sum/qnr_sum:.4f}'}")

    # Also check: C(t) * s_a(t) * s_b(t) grouped by chi(t)
    print(f"\n  C(t) * s_a(t) * s_b(t) by chi(t):")
    for a, b in [(0, 2), (0, 1), (3, 4)]:
        ga, gb = a + 1, b + 1
        q = resonance_level(ga, gb, p)

        qr_sum = 0
        qnr_sum = 0
        for t in range(1, p):
            prod = C[t] * sin_jt[a][t] * sin_jt[b][t]
            if legendre(t, p) == 1:
                qr_sum += prod
            else:
                qnr_sum += prod

        print(f"  ({a},{b}) q={q}: QR={qr_sum:.6f}, QNR={qnr_sum:.6f}, diff={qr_sum-qnr_sum:.6f}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
