#!/usr/bin/env python3
"""
qr_transitivity_proof.py — Is the QR action on chord types ALWAYS transitive?

The eigenvector theorem requires that the QR subgroup acts TRANSITIVELY
on {1,...,m} (chord types). This was verified for p = 7, 11, 19, 23, 31.

QUESTION: Is this always true for p = 3 mod 4?

PROOF ATTEMPT:
  The QR action on chord types is: a * k mod p, reduced to {1,...,m} via
  k -> p-k if k > m. Two chord types k, k' are in the same orbit iff
  there exists a QR element a with a*k = k' or a*k = p-k' (mod p).

  Equivalently: k' = a*k mod p  OR  k' = -a*k mod p = p-a*k mod p.
  Since -1 is NQR for p=3 mod 4, the second case means k' = (-a)*k mod p
  where -a is NQR. So we need: for every k, k' in {1,...,m}, there exists
  b in (Z/pZ)* with b*k = k' or b*k = p-k' mod p. That is:
  k and k' are in the same orbit iff b = k'/k or b = -k'/k for some b in QR or NQR.

  Actually: k' = a*k mod p (with a QR) or k' = -a*k = (NQR)*k mod p (with -a NQR).
  So k' can be reached from k by ANY nonzero multiplication!
  Since (Z/pZ)* is a group that acts transitively on {1,...,p-1},
  and our orbit on {1,...,m} identifies k with p-k,
  the quotient action of (Z/pZ)* / {+1,-1} on {1,...,m} is transitive
  REGARDLESS of the QR structure!

Wait, but the action is by QR elements only, not all of (Z/pZ)*.
Let me reconsider.

The signed permutation P_a maps chord type k to:
  - a*k mod p  if a*k mod p <= m  (with sign +1)
  - p - (a*k mod p)  if a*k mod p > m  (with sign -1)

The UNDERLYING permutation pi_a maps k to the chord index, ignoring sign.
pi_a(k) = min(a*k mod p, p - a*k mod p).

For the orbit of k under {pi_a : a QR}: this is {min(a*k mod p, p-a*k mod p) : a QR}.

For this to be all of {1,...,m}, we need: for every target t in {1,...,m},
there exists a QR element a with a*k = t or a*k = p-t (mod p).

a = t/k or a = -t/k = (p-t)/k mod p.

a is QR iff chi(a) = 1.
chi(t/k) = chi(t)*chi(k)^{-1} = chi(t)*chi(k)  (since chi^2 = 1)
chi(-t/k) = chi(-1)*chi(t/k) = -chi(t/k)  (since chi(-1) = -1 for p=3 mod 4)

So: at least one of t/k and -t/k is QR (since they have opposite chi values).
Therefore: for every k and t, there exists a QR element a with pi_a(k) = t.
QED: The QR action on chord types is ALWAYS transitive for p = 3 mod 4!

Author: opus-2026-03-12-S62
"""

import math

def is_qr(a, p):
    if a % p == 0: return False
    return pow(a, (p - 1) // 2, p) == 1


print("=" * 70)
print("PROOF: QR ACTION ON CHORD TYPES IS ALWAYS TRANSITIVE (p = 3 mod 4)")
print("=" * 70)

print("""
THEOREM: For every prime p = 3 mod 4, the QR subgroup of (Z/pZ)* acts
TRANSITIVELY on the set of chord types {1, 2, ..., (p-1)/2}.

PROOF:
  Let m = (p-1)/2. We need to show: for every k, t in {1,...,m},
  there exists a QR element a such that a*k = t (mod p) or a*k = -t (mod p).

  Consider the two candidates:
    a_1 = t * k^{-1} mod p    (satisfies a_1 * k = t mod p)
    a_2 = -t * k^{-1} mod p   (satisfies a_2 * k = -t mod p)

  We have a_2 = -a_1 mod p, so chi(a_2) = chi(-1) * chi(a_1) = -chi(a_1)
  (using chi(-1) = -1 for p = 3 mod 4).

  Therefore: exactly one of a_1, a_2 is a QR.

  If a_1 is QR: take a = a_1, then a*k = t mod p, and since t <= m,
    chord type k maps to t.
  If a_2 is QR: take a = a_2, then a*k = -t = p-t mod p.
    Since p-t > m (because t <= m < p/2 < p-t), we have pi_a(k) = p-(p-t) = t.

  In both cases: pi_a(k) = t. QED.

COROLLARY: The fixed subspace of the QR signed permutation action
on R^m is EXACTLY 1-dimensional, spanned by the Paley orientation sigma_P.

COROLLARY: The interaction matrix J has exactly (m+1)/2 distinct eigenvalues,
with the Paley eigenvalue being simple (multiplicity 1) and all others
having multiplicity 2.
""")

# Verify the proof for many primes
print("VERIFICATION:")
for p in [3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 97, 103, 107]:
    if p % 4 != 3:
        continue
    m = (p - 1) // 2

    # Check: for every pair (k, t), at least one of t/k, -t/k is QR
    all_ok = True
    for k in range(1, m + 1):
        for t in range(1, m + 1):
            k_inv = pow(k, p - 2, p)  # Fermat's little theorem
            a1 = (t * k_inv) % p
            a2 = ((-t) * k_inv) % p  # = (p - t) * k_inv mod p

            a1_qr = is_qr(a1, p)
            a2_qr = is_qr(a2, p)

            if not (a1_qr or a2_qr):
                print(f"  FAILURE at p={p}, k={k}, t={t}!")
                all_ok = False
                break

            # Verify exactly one is QR
            if a1_qr and a2_qr:
                print(f"  BOTH QR at p={p}, k={k}, t={t}! (a1={a1}, a2={a2})")
                all_ok = False
                break
        if not all_ok:
            break

    status = "OK" if all_ok else "FAIL"
    print(f"  p={p:3d} (m={m:3d}): {status}")


# Now prove the STRONGER claim: the eigenvalue at sigma_P is maximal

print(f"\n{'=' * 70}")
print("IS THE PALEY EIGENVALUE ALWAYS MAXIMAL?")
print("=" * 70)

print("""
We've proved sigma_P is an eigenvector of J. But is its eigenvalue
lambda_0 always the LARGEST?

At p=7: lambda_0 = 7.0, others = -3.5. YES.
At p=11: lambda_0 = 561.0, others = 154.9, -435.4. YES.

QUESTION: Is lambda_0 > lambda_j for all j > 0, for ALL p = 3 mod 4?

From the eigenvalue formula:
  lambda_0 = (1/m) * [E[H * A^2] - m * E[H]]

Since A(sigma_P) = m (maximal) and H is monotone increasing in |A|
(at least at p=7,11), we expect lambda_0 to be large.

But we CANNOT prove lambda_0 is always maximal without understanding
the higher eigenvalues. The maximality depends on the CURVATURE of
H as a function of sigma, which involves all Walsh degrees.

CONJECTURE: lambda_0 is the unique maximal eigenvalue of J for ALL
p = 3 mod 4. This would mean Paley ALWAYS maximizes the degree-2
quadratic form, even when it doesn't maximize H overall.
""")

# Additional verification: check that the QR transitivity proof works
# by explicitly constructing the QR element for each (k,t) pair at p=19

print("\np=19: Explicit QR elements mapping chord types:")
p = 19
m = 9
print(f"  {'k->t':>6} {'a (QR)':>8} {'which':>6} {'a*k mod p':>10} {'maps to':>8}")
for k in [1, 5, 9]:
    for t in [1, 5, 9]:
        k_inv = pow(k, p - 2, p)
        a1 = (t * k_inv) % p
        a2 = ((p - t) * k_inv) % p

        if is_qr(a1, p):
            a = a1
            result = (a * k) % p
            mapped = result if result <= m else p - result
            which = "t/k"
        else:
            a = a2
            result = (a * k) % p
            mapped = result if result <= m else p - result
            which = "-t/k"

        print(f"  {k}->{t:>2} {a:>8} {which:>6} {result:>10} {mapped:>8}")
