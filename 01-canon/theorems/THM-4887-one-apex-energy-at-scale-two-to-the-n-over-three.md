---
id: THM-4887
title: "Smoothed energy of the one-apex counts at scales 2^(n/3) to 2^(n/2). For n >= 12 and 0 <= Delta <= 2^floor(n/3), the pairs (v, v') in [2^(n-1), 2^n)^2 with |g(v) - g(v')| <= Delta number at most (1.39 n + 19) 2^(4n/3), sharp up to O(n) at Delta = 2^floor(n/3); more generally E_n(Delta) <= 17 (n - 2m)(Delta + 1) 2^n for 2^(m-1) < Delta <= 2^m and (n-1)/3 <= m <= (n-4)/2, with an absolute constant at the top scale m = floor((n-4)/2). Hence the energy E_n is at most (1.39 n + 16) 2^(4n/3), the one-apex counts of length n take at least 2^(2n/3) / (5.56 n + 64) distinct values, and the number A(X) of Hamiltonian-path counts up to X is at least X^(2/3) / (13 ln X + 102) for every X >= 4096 (no earlier power bound beyond Session 150's sqrt(X)/4). Positive density by this route needs the inequality at Delta = 0 (linear energy), which is OPEN; density one and Chan-Pak OPEN"
status: >
  PROVED: the strip lemma, the smoothed-energy bound (D1) at every admissible split, (D2), the scale-invariant form
  (D5), the energy bound and the counting corollary (D3), and (D4). New: the method uses no arithmetic of g below the
  scale of the free block; it is an approximate-collision (thin-strip lattice point) count for the perturbed product set
  {g(h) + h c}, in the spirit of the multiplication-table and shifted-product-energy problems.
  INHERITED: the two-block identity is Session 193's concatenation law (1) (modular_sync/PROOF.md); same-shape identities
  are Session 206's for its tail function Q (session206/parameter-capacity/audit/REVIEW.md) and Session 16's one-block
  binary-product identity (BINARY_PRODUCT_DENSITY.md (4)). 1 <= g(l) <= l + 1 is in Session 158 binary_bank/NOTE.md.
  The counting lemma is the classical spacing lemma behind Weyl's inequality (Vaughan, The Hardy-Littlewood Method,
  Lemma 2.2), also THM-4886's step 3; the small-denominator count is a standard Dirichlet/Farey dissection.
  ELEMENTARY: the sharpness floor (Cauchy-Schwarz over intervals, valid for any multiset).
  The X^(2/3) exponent is new: there is no earlier power bound beyond Session 150's sqrt(X)/4 in the repository, and
  none was found in the literature (Chan-Pak, TCS 1015 (2024), section 6.5, Conjecture 6.2, records only the conjecture
  and Royle's verification to 80,557, MathOverflow 232751).
  The counting corollary is asymptotic. Order-join multiplicativity (THM-1862), the certified range (THM-4875 (D):
  every odd value up to 2^32 except 7, 21) and unique factorisation give more for every X < 2^14,751. Multisets of
  certified primes from all dyadic ranges suffice (the audit's rigorous f_smooth construction); 31 distinct primes from
  [2^31, 2^32] alone give C(98,182,656, 31) = 6.9e213 counts below 2^1000, against X^(2/3) = 10^200.7. The content is
  the exponent: a finite set of certified values has only polylogarithmically many products.
  NOT CLAIMED: linear energy E_n = O(2^n) (the statement at Delta = 0; larger-scale bounds do not imply it), positive
  density, density one, Chan-Pak (all OPEN).
  AUDIT: PASS WITH CORRECTIONS (3 majors on scope and wording, 10 minors, 4 upgrades; all applied),
  04-computation/experiments/one_apex_energy_20261106_audit_S83/AUDIT.md.
session: opus-2026-11-06-S83
source: 05-knowledge/results/one_apex_exponential_sums_bilinear_weyl_bound_20261106.md
related:
  - THM-4886
  - THM-4885
  - THM-4875
  - THM-1862
  - HYP-9479
  - 05-knowledge/results/session150_hamiltonian_spectrum_20261103/induction/RESULT.md
  - 05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md
  - 05-knowledge/results/feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md
  - 05-knowledge/results/session26_random_bridges_20261022/session26-tournament/CUT_ADDENDUM.md
  - 05-knowledge/results/session26_random_bridges_20261022/session26-tournament/ASYMPTOTIC_REVIEW.md
verification:
  - 04-computation/experiments/one_apex_fourier_20261106/d_identity.py
  - 04-computation/experiments/one_apex_fourier_20261106/d_lemma.py
  - 04-computation/experiments/one_apex_fourier_20261106/d_energy.py
  - 04-computation/experiments/one_apex_energy_20261106_audit_S83/b_strip.py
  - 04-computation/experiments/one_apex_energy_20261106_audit_S83/c_energy.py
  - 04-computation/experiments/one_apex_energy_20261106_audit_S83/d_constants.py
  - 04-computation/experiments/one_apex_energy_20261106_audit_S83/g_scales.py
---

# THM-4887 - smoothed energy of the one-apex counts at scales 2^(n/3) to 2^(n/2), and an X^(2/3) counting bound

**Notation.**
- g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2). g(v) is the Hamiltonian path count of the one-apex
  tournament of the binary word of v (Session 193 section 2).
- X_n = [2^(n-1), 2^n), the words of length n beginning with 1.
- For Delta >= 0, E_n(Delta) = #{(v, v') in X_n^2 : |g(v) - g(v')| <= Delta}. The energy is E_n = E_n(0).
- A(X) is the number of distinct integers at most X that are the Hamiltonian path count of some tournament.

## Statement

**(D1) General form.** Let m >= 1, C = 2^m, P = 2^(k-1) >= 8C, n = k + m and 0 <= Delta <= C. Then

    E_n(Delta) <= P C + 10 P^2 + P C^2 (35.74 + 12 ln(P/C)) + 3 P C (1 + ln(P/C)).

For Delta = 0 the constants improve to 8 P^2 and 33.74.

**(D2) Scale 2^(n/3).** For n >= 12 and 0 <= Delta <= 2^floor(n/3):

    E_n(Delta) <= (1.39 n + 19) 2^(4n/3).

- At Delta = 2^floor(n/3) this is sharp up to O(n): every multiset of 2^(n-1) values in [1, 2^n] has
  E_n(Delta) >= 2^(2n-2)(Delta + 1)/(2^n + Delta + 1), at least 2^(4n/3)/8 there.
- For Delta = o(2^(n/3)) it is not sharp. At Delta = 0 it exceeds the observed E_n = 3.25 |X_n| by about
  0.6 (1.39 n + 19) 2^(n/3).

**(D3) Energy and counting** (the Delta = 0 constants).
- E_n <= (1.39 n + 16) 2^(4n/3) for n >= 12.
  - The trivial bound is 2^(2n-2).
  - THM-4886's floor is E_n >= 2 |X_n| - 2^(n/2) - 1, from Session 16's converse involution.
  - The observed value is 3.25 · 2^(n-1) (HYP-9479 P2).
- The one-apex counts on X_n take at least 2^(2n/3)/(5.56 n + 64) distinct values, all at most 2^n.
- A(X) >= X^(2/3)/(13 ln X + 102) for every real X >= 4096.
- For comparison, the one-apex family alone has at most (ln 2/2) X + O(sqrt X) values up to X (Session 16
  BINARY_PRODUCT_DENSITY.md (2); Session 26 ASYMPTOTIC_REVIEW.md).

**(D4) The window near 0.** Let S_n(theta) = sum over v in X_n of e(theta g(v)). Then

    integral over |theta| <= 2^(-floor(n/3)-2) of |S_n(theta)|^2 dtheta <= 3.92 (1.39 n + 19) 2^n.

- The integral over [0, 1) is E_n. So this window adds only O(n 2^n) to Parseval's integral.
- The central peak |theta| <= 2^(-n-3) alone contributes at least 2^(n-5), because every phase there lies within pi/4
  of 0. So (D4) is sharp up to O(n): the window carries no L^2 mass beyond the central peak, up to O(n).
- This is not square-root size in mean square over the window: the window has length about 2^(-n/3), so its
  root-mean-square is at least about 2^(2n/3).
- Exactly, for n = 12..22, the window integral is 0.31 · 2^n and the peak 0.0618 · 2^n (audit `c_energy.out`).

**(D5) All scales from 2^(n/3) to 2^(n/2)/4** (the audit's U1). (D1) holds at every split with P >= 8C. Divided by
2^(n+m), for (n-1)/3 <= m <= (n-4)/2 and 0 <= Delta <= 2^m, it reads

    E_n(Delta) <= [17.87 + 6 ln 2 (n - 2m - 1) + 2.5 · 2^(n-3m) + 2^(-m-1) (4 + 3 (n - 2m - 1) ln 2)] 2^(n+m).

- For n >= 12 the bracket is at most 8.3 (n - 2m). So for 2^(m-1) < Delta <= 2^m:

      E_n(Delta) <= 17 (n - 2m)(Delta + 1) 2^n,

  i.e. linear smoothed energy, up to O(n - 2 log2 Delta), at every dyadic scale from 2^(n/3) to 2^(n/2)/4.
- At the top scale m = floor((n-4)/2) the bracket is an absolute constant:
  - at most 39.9 for n >= 12 (the maximum is at n = 13);
  - at most 34.9 for n >= 20;
  - tending to 30.35 for even n and 34.51 for odd n.
- The exact values E_n(2^m)/2^(n+m) are 0.61-0.63 for n = 18..24 and every admissible m (audit `g_scales.out`).

## Proof

1. **Two-block identity** (Session 193 eq. (1), with u the word of h and y the m-letter word of l, leading zeros
   allowed, H(y) = g(l)). For h >= 1, m >= 1 and 0 <= l < 2^m, write l = sum_(i<m) l_i 2^i with the least
   significant bit first. Then

       g(h 2^m + l) = g(h) + h C_m(l) + gamma_h(l),

   - C_m(l) = sum_(i<m) (1 - l_i) 2^(m-1-i) is the complement of l reversed, a bijection of [0, 2^m).
   - gamma_h(l) = g(l) - 1 + (h mod 2)[l < 2^(m-1)].
   - **Induction from the recursion** (audit). For m = 1 the identity is the recursion itself. To append a bit
     l = 2l' + b, use g(2V + 1) = g(V) and g(2V) = g(V) + V + (V mod 2) with V = h 2^m + l', together with
     C_(m+1)(2l' + 1) = C_m(l') and C_(m+1)(2l') = 2^m + C_m(l').
   - **Range of gamma.** 1 <= g(l) <= l + 1, by induction: g(2v + 1) = g(v); g(2v) <= 2v + 1 for v even; and for v odd,
     g(v) = g((v - 1)/2) <= (v + 1)/2, so g(2v) <= (v + 1)/2 + v + 1 <= 2v + 1.
     Equality g(l) = l + 1 holds only at l = 0 and at powers of 2, so 0 <= gamma_h(l) <= 2^m - 1 (audit U2). The proof
     below uses only gamma in [0, 2^m].

2. **Coordinates.** v = h 2^m + l, with h in [P, 2P) and l in [0, C), parametrises X_n. Write c = C_m(l).
   - **Same h, l != l'.** |g(v) - g(v')| >= h |c - c'| - C >= P - C >= 7C > Delta, so only the diagonal pairs count.
   - **h != h'.** If |g(v) - g(v')| <= Delta, then |h c - h' c' - D| <= W := C + Delta <= 2C, where D = g(h') - g(h).
     Hence such pairs number at most

         Z(h, h') = #{(c, c') in [0, C)^2 : |h c - h' c' - D| <= W}.

3. **Strip lemma.** Let B = P/C, a power of 2 and at least 8.
   - By Dirichlet there are integers a and 1 <= b <= B with |b h/h' - a| < 1/B.
   - Reduce a/b to lowest terms a'/beta. Then |h/h' - a'/beta| < 1/(beta B) <= 1/beta^2, and

         Z(h, h') <= (C/beta + 1)(2 beta W/h' + 3) <= 4 C^2/P + 3 C/beta + 7.

   *Proof.*
   - Since h' > 2W, each c has at most one c'. So Z(h, h') <= #{c in [0, C) : ||alpha c - delta|| <= t}, where
     alpha = h/h', delta = D/h' and t = W/h'. Dropping the constraint c' in [0, C) only enlarges the count.
   - **The classical spacing lemma.** Take a block c = c_0 + j with 0 <= j < beta.
     - alpha c - delta = x_0 + j a'/beta + j (alpha - a'/beta), and the last term has modulus below 1/beta.
     - The points x_0 + j a'/beta run over a translate of the grid (1/beta) Z mod 1.
     - So at most 2 beta (t + 1/beta) + 1 = 2 beta t + 3 values of j qualify.
   - At most C/beta + 1 blocks cover [0, C). Finally 2 beta W/h' <= 2B (2C)/P = 4, and 2 C W/h' <= 4 C^2/P.
   - The lemma holds for any admissible Dirichlet pair. The audit tested the largest reduced denominator, which is the
     strongest form.

4. **Few pairs have small denominators.** If beta(h, h') = beta, then |h beta - h' a'| < h'/B < 2C.
   - For fixed (h', a', beta), at most 4C/beta + 1 values of h qualify.
   - For fixed (h', beta), a' lies in an open interval of length (beta P + 4C)/h' <= beta + 1/2, so it takes at most
     beta + 1 values.
   - Summing over the P values of h', the number of pairs with beta(h, h') = beta is at most P (beta + 1)(4C/beta + 1).
   - Hence

         sum over h != h' of 1/beta(h, h') <= P sum_(beta <= B) (beta + 1)(4C/beta + 1)/beta
                                          <= P (4C ln B + 10.58 C + B + 1 + ln B).

   - So the lattice term 3C/beta + 7 is O(ln B) on average over pairs. That is the source of the O(n) factor.

5. **Sum.** Using C B = P, the term 3 C B P is 3 P^2, and

       E_n(Delta) <= P C + sum over h != h' of Z(h, h')
                  <= P C + P^2 (4 C^2/P + 7) + 3 C P (4 C ln B + 10.58 C + B + 1 + ln B),

   which is (D1). For Delta = 0, W = C and the lemma gives 2C^2/P + 3C/beta + 5, hence 8 P^2 and 33.74.

6. **(D2) and (D3).** Write n = 3t + r with 0 <= r <= 2, and take m = t and k = 2t + r. Then P/C = 2^(t+r-1) >= 8
   when n >= 12. Also P C = 2^(n-1), P^2 = 2^(4n/3) 2^(2r/3-2) and P C^2 = 2^(4n/3) 2^(-1-r/3).

   Dividing (D1) by 2^(4n/3) gives the table below. Each row includes the remaining terms; the remainder bound is taken
   at n = 12, 13, 14 respectively and decreases in n. Values are rounded up. The Delta <= C column was also computed at
   50 digits by the audit (`d_constants.out`).

   | r | Delta <= C (10 P^2, 35.74) | Delta = 0 (8 P^2, 33.74) |
   |---|---|---|
   | 0 | 16.54 + 1.3863 n | 15.04 + 1.3863 n |
   | 1 | 17.36 + 1.1004 n | 15.77 + 1.1004 n |
   | 2 | 18.72 + 0.8734 n | 16.83 + 0.8734 n |

   - Every entry is at most 1.39 n + 19 and 1.39 n + 16 respectively. The margins are at least 2.51 and 1.01.
   - **Sharpness at Delta = 2^t.** The values lie in [1, 2^n]. Cut [1, 2^n] into ceil(2^n/(Delta + 1)) intervals of
     length Delta + 1 and apply Cauchy-Schwarz to the interval counts.
   - **(D3).**
     - The number of distinct values is at least |X_n|^2/E_n >= 2^(2n/3)/(4(1.39 n + 16)).
     - For X >= 4096 take n = floor(log2 X) >= 12. Then A(X) >= (X/2)^(2/3)/(5.56 log2 X + 64), which is at least
       X^(2/3)/(13 ln X + 102).

7. **(D4).**
   - Let F(y) = max(0, 1 - |y|/(Delta + 1)) <= 1[|y| <= Delta].
   - Its transform is |sum_(j=0..Delta) e(j theta)|^2/(Delta + 1) >= 0, and it is at least 4(Delta + 1)/pi^2 for
     |theta| <= 1/(2(Delta + 1)).
   - So the window integral is at most (pi^2/4) E_n(Delta)/(Delta + 1).
   - Take Delta = 2^t. Then 1/(2(Delta + 1)) >= 2^(-t-2) and 2^(4n/3)/Delta <= 2^(n + 2/3).

8. **(D5).** Substitute P = 2^(n-m-1), C = 2^m and ln(P/C) = (n - 2m - 1) ln 2 into (D1) and divide by 2^(n+m). The
   condition P >= 8C is m <= (n-4)/2.
   - With m >= (n-1)/3 and n >= 12 we have n - 3m <= 1 and n - 2m >= 4, so the bracket is at most 8.3(n - 2m).
   - Since 2^m < 2(Delta + 1) for Delta > 2^(m-1), the stated form follows.

QED.

## What the method does and does not reach

- **Where the window comes from.** The proof uses about the low block only that gamma lies in [0, C], and no arithmetic
  of g there. A perturbation of size C forces the window W >= C. Per pair, the strip count then has an area term 4C^2/P
  and a lattice term 3C/beta + 7, which is O(ln B) on average.
- **Where it stops.** Balancing P^2 against P C^2 puts the smallest reachable scale at C = 2^(n/3). There the bound is
  sharp up to O(n). By (D5) every larger scale up to 2^(n/2)/4 is also linear, up to O(n - 2m).
- **What it does not give.** Linear energy E_n = O(2^n) is the statement at Delta = 0. Larger-scale bounds do not imply
  it, since E_n(0) <= E_n(Delta) runs the wrong way. Positive density by this route needs that inequality, which is
  OPEN. Heuristically, it would need arithmetic information on g below the scale 2^(n/3): where D = g(h') - g(h) sits
  relative to the lattice {h c - h' c'}.
- **Other routes.** Positive density of the spectrum could also come from other families: joins, banks, products.
  Linear energy of the one-apex family is sufficient, not necessary.

## Checks

**Session** (`04-computation/experiments/one_apex_fourier_20261106/`).
- `d_identity.py`:
  - the two-block identity on 200,000 random (h, m, l), with 0 violations;
  - C_m bijective for m <= 12;
  - g(v) <= v + 1 for v < 2^18.
- `d_lemma.py`: (n, m) = (10, 2), (11, 2), (12, 3), (13, 3), all ordered pairs, with the worst-case shift.
  - Strip lemma: 0 violations. The max ratio to the bound is 0.565 for W = C and 0.635 for W = 2C.
  - The small-denominator counts and the harmonic sum are within their bounds.
  - Same-h gaps exceed C.
  - E_n(0) and E_n(C) both lie below PC + sum Z* and below (D1).
- `d_energy.py`, n = 12..26:
  - E_n(2^floor(n/3))/2^(4n/3) lies in [0.39, 0.67], against the proved 35.7 to 55.1. The constants are loose by a
    factor of 53 to 141.
  - E_n/|X_n| rises from 2.93 to 3.2607 at n = 22, then falls slowly (3.2473 at n = 28; HYP-9479).

**Audit** (`04-computation/experiments/one_apex_energy_20261106_audit_S83/`).
- The identity, re-derived by induction and tested:
  - exhaustively for m <= 8 and h < 2^9;
  - on 100,000 random big-integer cases with m <= 64.
- The strip lemma with the strongest Dirichlet choice, on all pairs for (10,2), (12,4), (13,4), (14,4), (13,3), (15,5)
  and (16,5). (12,4) is the boundary case B = 8 of (D2). Random pairs were tested for (18,6) and (21,7). There were
  0 violations.
- Every constant at 50 digits.
- (D4) computed exactly by autocorrelation.
- (D5) against the exact E_n(2^m).
- The prime-product comparisons: `e_primes.out` and the rigorous `f_smooth.out`.
- The W = C totals of the session's earlier run are reproduced by `h_totals.py`:
  - lemma totals 120,050 and 2,027,274;
  - actual-shift totals 2,822 and 44,838, for (10, 2) and (13, 3).
