---
id: THM-4887
title: "Energy of the one-apex counts at scale 2^(n/3). For n >= 12 and 0 <= Delta <= 2^floor(n/3), the pairs (v, v') in [2^(n-1), 2^n)^2 with |g(v) - g(v')| <= Delta number at most (1.39 n + 19) 2^(4n/3), which is sharp up to a factor O(n). Hence the energy E_n is O(n 2^(4n/3)), the one-apex counts of length n take at least 2^(2n/3) / (5.56 n + 76) distinct values, and the number A(X) of Hamiltonian-path counts up to X is at least X^(2/3) / (13 ln X + 121) for every X >= 4096 (the previous all-X exponent was 1/2, Session 150). Positive density needs the same equidistribution below the scale 2^(n/3); density one and Chan-Pak OPEN"
status: >
  PROVED (new): the two-block identity (a special case of Session 193's concatenation law), the strip lemma, the
  smoothed-energy bound, its sharpness up to O(n), the energy bound and the counting corollary.
  The counting corollary is asymptotic. Order-join multiplicativity (THM-1862) and the certified range (every odd value
  up to 2^32 except 7, 21; THM-4875 (D)) already do better below 2^1000: products of 31 distinct primes from
  [2^31, 2^32] give C(98,182,656, 31) > 10^213 counts below 2^1000, against X^(2/3) = 10^200.7 there. The corollary's
  content is the exponent: finitely many certified values generate only polylogarithmically many products, and
  Session 150's family has exponent 1/2.
  NOT CLAIMED: linear energy, positive density, density one, Chan-Pak (all OPEN).
  AUDIT: pending.
session: opus-2026-11-06-S83
source: 05-knowledge/results/one_apex_exponential_sums_bilinear_weyl_bound_20261106.md
related:
  - THM-4886
  - THM-4885
  - THM-1862
  - THM-4875
  - HYP-9479
  - 05-knowledge/results/session150_hamiltonian_spectrum_20261103/induction/RESULT.md
  - 05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md
  - 05-knowledge/results/session26_random_bridges_20261022/session26-tournament/CUT_ADDENDUM.md
verification:
  - 04-computation/experiments/one_apex_fourier_20261106/d_identity.py
  - 04-computation/experiments/one_apex_fourier_20261106/d_lemma.py
  - 04-computation/experiments/one_apex_fourier_20261106/d_energy.py
---

# THM-4887 - energy of the one-apex counts at scale 2^(n/3), and an X^(2/3) counting bound

**Notation.**
- g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2); g(v) is the Hamiltonian path count of the one-apex
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

The bound is sharp up to a factor O(n): every multiset of 2^(n-1) values in [1, 2^n] has
E_n(Delta) >= 2^(2n-2) (Delta + 1) / (2^n + Delta + 1), which is at least 2^(4n/3) / 8 at Delta = 2^floor(n/3).

**(D3) Energy and counting.**
- E_n <= (1.39 n + 19) 2^(4n/3) for n >= 12. The trivial bound is 2^(2n-2), and the observed value is 3.25 · 2^(n-1)
  (HYP-9479 P2).
- The one-apex counts on X_n take at least 2^(2n/3) / (5.56 n + 76) distinct values, all at most 2^n.
- A(X) >= X^(2/3) / (13 ln X + 121) for every real X >= 4096.

**(D4) Mean-square cancellation near 0.** With S_n(theta) = sum over v in X_n of e(theta g(v)):

    integral over |theta| <= 2^(-floor(n/3)-2) of |S_n(theta)|^2 dtheta <= 3.92 (1.39 n + 19) 2^n.

The whole integral over [0, 1) is E_n >= 2^(n-1). So on this window, |S_n| has square-root size in mean square, up to
O(n). THM-4886 (A) gives only a pointwise power saving on the minor arcs.

## Proof

1. **Two-block identity.** For h >= 1, m >= 1 and 0 <= l < 2^m:

       g(h 2^m + l) = g(h) + h C_m(l) + gamma_h(l),

   where:
   - C_m(l) = sum_(i<m) (1 - l_i) 2^(m-1-i), with l_i the bits of l, i.e. the complement of l reversed. It is a
     bijection of [0, 2^m).
   - gamma_h(l) = g(l) - 1 + (h mod 2) [l < 2^(m-1)].

   This is Session 193's concatenation law H(uvw) with w empty: H(uv) = H(u) + H(v) - 1 + A(u) C(v) + p(u) eps(v).
   Here v is the m-letter word of l, leading zeros allowed, and H(v) = g(l).

   - 1 <= g(l) <= l + 1, by induction:
     - g(2v + 1) = g(v);
     - g(2v) <= 2v + 1 for v even;
     - for v odd, g(v) = g((v - 1)/2) <= (v + 1)/2, so g(2v) <= (v + 1)/2 + v + 1 <= 2v + 1.
   - Hence 0 <= gamma_h(l) <= 2^m.

2. **Coordinates.** v = h 2^m + l with h in [P, 2P) and l in [0, C) parametrises X_n. Write c = C_m(l).
   - **Same h, l != l':** |g(v) - g(v')| >= h |c - c'| - C >= P - C >= 7C > Delta, so only the diagonal pairs count.
   - **h != h':** if |g(v) - g(v')| <= Delta, then |h c - h' c' - D| <= W := C + Delta <= 2C, with D = g(h') - g(h).
     Hence such pairs number at most

         Z(h, h') = #{(c, c') in [0, C)^2 : |h c - h' c' - D| <= W}.

3. **Strip lemma.** Let B = P/C (a power of 2, at least 8). By Dirichlet there are integers a and 1 <= b <= B with
   |b h/h' - a| < 1/B. Reduce a/b to lowest terms a'/beta. Then |h/h' - a'/beta| < 1/(beta B) <= 1/beta^2, and

       Z(h, h') <= (C/beta + 1)(2 beta W/h' + 3) <= 4 C^2/P + 3 C/beta + 7.

   *Proof.*
   - Since h' > 2W, each c has at most one c'. So Z(h, h') <= #{c in [0, C) : ||alpha c - delta|| <= t}, with
     alpha = h/h', delta = D/h' and t = W/h'.
   - **Counting lemma.** In a block c = c_0 + j (0 <= j < beta):
     - alpha c - delta = x_0 + j a'/beta + j (alpha - a'/beta), and the last term has modulus below 1/beta;
     - the points x_0 + j a'/beta run over a translate of the grid (1/beta) Z mod 1;
     - so at most 2 beta (t + 1/beta) + 1 = 2 beta t + 3 values of j qualify.
   - At most C/beta + 1 blocks cover [0, C). Finally 2 beta W/h' <= 2B (2C)/P = 4 and 2 C W/h' <= 4 C^2/P.

4. **Few pairs have small denominators.** If beta(h, h') = beta, then |h beta - h' a'| < h'/B < 2C.
   - For fixed (h', a', beta), at most 4C/beta + 1 values of h qualify.
   - For fixed (h', beta), a' lies in an open interval of length (beta P + 4C)/h' <= beta + 1/2, so it takes at most
     beta + 1 values.
   - Summing over the P values of h', the number of pairs with beta(h, h') = beta is at most P (beta + 1)(4C/beta + 1).
   - Hence:

         sum over h != h' of 1/beta(h, h') <= P sum_(beta <= B) (beta + 1)(4C/beta + 1)/beta
                                           <= P (4C ln B + 10.58 C + B + 1 + ln B).

5. **Sum.** Use C B = P to absorb the term 3 C B P into 3 P^2:

       E_n(Delta) <= P C + sum over h != h' of Z(h, h')
                  <= P C + P^2 (4C^2/P + 7) + 3 C P (4C ln B + 10.58 C + B + 1 + ln B),

   which is (D1). For Delta = 0, W = C and the lemma gives 2C^2/P + 3C/beta + 5, hence 8 P^2 and 33.74.

6. **(D2).** Write n = 3t + r with 0 <= r <= 2, and take m = t, k = 2t + r. Then P/C = 2^(t+r-1) >= 8 when n >= 12, and:
   - P C = 2^(n-1);
   - P^2 = 2^(4n/3) 2^(2r/3-2);
   - P C^2 = 2^(4n/3) 2^(-1-r/3).

   Dividing (D1) by 2^(4n/3) gives the following, where each line includes the remaining terms, which total at most 0.32
   for n >= 12:

   | r | bound / 2^(4n/3) |
   |---|---|
   | 0 | 16.53 + 1.386 n |
   | 1 | 17.36 + 1.101 n |
   | 2 | 18.71 + 0.873 n |

   All three are at most 1.39 n + 19.

   **Sharpness.** The values lie in [1, 2^n], so cut [1, 2^n] into ceil(2^n/(Delta + 1)) intervals of length Delta + 1 and
   apply Cauchy-Schwarz to the interval counts.

7. **(D3).**
   - The first two items: E_n = E_n(0), and the number of distinct values is at least |X_n|^2 / E_n.
   - These values are Hamiltonian path counts of one-apex tournaments, and all are at most 2^n.
   - For X >= 4096 take n = floor(log2 X) >= 12. Then A(X) >= 2^(2n/3)/(5.56 n + 76) >= (X/2)^(2/3)/(5.56 log2 X + 76),
     which is at least X^(2/3)/(13 ln X + 121).

8. **(D4).**
   - Let F(y) = max(0, 1 - |y|/(Delta + 1)) <= 1[|y| <= Delta]. Its transform is
     |sum_(j=0..Delta) e(j theta)|^2 / (Delta + 1) >= 0, and it is at least 4(Delta + 1)/pi^2 for |theta| <= 1/(2(Delta + 1)).
   - So the integral of |S_n|^2 over that window is at most (pi^2/4) E_n(Delta)/(Delta + 1).
   - Take Delta = 2^floor(n/3). Then 1/(2(Delta + 1)) >= 2^(-floor(n/3)-2) and 2^(4n/3)/Delta <= 2^(n + 2/3).

QED.

## Why 2/3, and what positive density needs

- The proof uses about the low block only that gamma lies in [0, C], and no arithmetic of g there. A perturbation of size
  C forces the window W >= C.
- The strip count then has an area term 4C^2/P per pair and a lattice term O(1) per pair. Balancing P^2 against P C^2
  puts the free block at C = 2^(n/3), and the bound is sharp there up to O(n).
- Pushing E_n(Delta) to linear size at smaller Delta needs arithmetic information on g below the scale 2^(n/3):
  how D = g(h') - g(h) sits relative to the lattice {h c - h' c'}. This is exactly the input linear energy
  E_n = O(2^n) needs, and with it positive density (THM-4886, note section 5).

## Checks (`04-computation/experiments/one_apex_fourier_20261106/`)

- `d_identity.py`: the two-block identity on 200,000 random (h, m, l) (0 violations); C_m bijective for m <= 12;
  g(v) <= v + 1 for v < 2^18.
- `d_lemma.py`:
  - the strip lemma with the worst-case shift D over all ordered pairs h != h' for (n, m) = (10, 2), (11, 2), (12, 3),
    (13, 3): 0 violations, max ratio 0.565;
  - sum of Z over pairs 20,510 to 336,522 against lemma totals 120,050 to 2,027,274;
  - actual-shift totals 2,822 to 44,838.
- `d_energy.py`, n = 12..26:
  - E_n(2^floor(n/3)) / 2^(4n/3) lies in [0.39, 0.67] against the proved 35.7 to 55.1;
  - E_n / |X_n| rises from 2.93 to 3.26.
  The constants are loose by a factor of about 100; the exponent is the content.
