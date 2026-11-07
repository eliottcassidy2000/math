---
id: THM-4890
title: "Micro/macro balance of one-apex collisions. If two distinct reduced words (top bit 1, bottom bit 0) of length n have the same one-apex Hamiltonian path count, they differ only on a window [s, t]. With a = n-1-t the common-prefix length, carry depths dX, dXi of the window values and of their reversed complements, A = a - dXi and S = s - dX: |A - S| <= 2 or max(A, S) <= 3. So the window is centred, and the converse involution swaps A and S. Also: the energy of all words is at most 1 + 2 sum_j E^Y_j over reduced words; and with the bottom-half micro self-energies replaced by independent uniform random values, the expected energy is at most 109 |X_n| (n >= 12). So linear energy, which stays OPEN, is equivalent up to constants to a decorrelation between micro self-energies and macro strips."
status: >
  PROVED (new): (C) the balance lemma, (C') its converse symmetry, (D) the randomised-micro bound and the equivalence
  it gives, (E) the reduction to reduced words.
  INHERITED: (A) the two-block law is Session 193's concatenation law (1)
  (05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md), written symmetrically with the
  converse invariance g(C_m(l)) = g(l) (THM-4886's converse involution). (B) the three-block law is Session 193's (2).
  The strip lemma used in (D) is THM-4887 step 3 (the classical spacing lemma); 1 <= g(l) <= l + 1 is Session 158's
  (binary_bank/NOTE.md, quoted in THM-4887).
  FINITE-EXACT checks: (C) on every collision pair (converse pairs included) of reduced words for n = 10..20, 0
  violations. The sharper form |A - S| <= 1 or max(A, S) <= 3 also has 0 exceptions to n = 20; this is not proved.
  HYP-9482 P5 (pre-registered) confirmed |a - s| <= 1 + max(dX, dXi) at n = 22, 24, 26.
  NOT CLAIMED: linear energy, positive density, density one, Chan-Pak (all OPEN).
session: opus-2026-11-06-S84
source: 05-knowledge/results/one_apex_collision_balance_20261106/README.md
related:
  - THM-4887
  - THM-4886
  - THM-4885
  - HYP-9482
  - HYP-9480
  - 05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md
verification:
  - 05-knowledge/results/one_apex_collision_balance_20261106/check_balance.py
  - 05-knowledge/results/one_apex_collision_balance_20261106/s84_energy.py
---

# THM-4890 - micro/macro balance of one-apex collisions

## Objects

- **Counts.** g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2). The one-apex tournament of the word v has g(v)
  Hamiltonian paths (THM-4875, THM-4878).
- **Word sets.** X_n = [2^(n-1), 2^n). The reduced words are Y_n = {v in X_n : v even}, |Y_n| = 2^(n-2).
- **Reversed complement.** For a w-bit word x, C_w(x) = sum_(i<w) (1 - x_i) 2^(w-1-i). For v in X_n, v* = C_n(v).
  - C_w is an involution on [0, 2^w).
  - g(C_w(x)) = g(x): it is the converse tournament.
  - v -> v* maps Y_n onto Y_n.
- **Energies.** E_n = #{(v, v') in X_n^2 : g(v) = g(v')}, and E^Y_n is the same count on Y_n.

## Statement

**(A) Two-block law, symmetric form (INHERITED).** For h >= 1, m >= 1 and 0 <= l < 2^m, with c = C_m(l):

    g(h 2^m + l) = g(h) + g(c) + h c + [h odd][c odd] - 1.

Here h is the macro block and c is the micro block read backwards. The converse v -> v* exchanges the two blocks.

**(B) Window identity (INHERITED form).** Write v = h | x | sigma, with prefix h of a bits, window x of w >= 1 bits and
suffix sigma of s bits. Put tau = C_s(sigma). Then

    g(v) = g(h) + g(x) + g(tau) + h tau 2^w + h C_w(x) + tau x + eps - 2,
    eps = [x odd][tau odd] + [h odd][x_(w-1) = 0].

So for v' = h | x' | sigma with the same prefix and suffix:

    g(v) - g(v') = (g(x) - g(x')) + h Xi + tau X + (eps - eps'),   X = x - x', Xi = C_w(x) - C_w(x').

**(C) Balance lemma (PROVED, new).** Let v != v' in Y_n with g(v) = g(v').
- Let t and s be the highest and lowest bit positions where they differ. Put a = n - 1 - t and w = t - s + 1.
- Let x, x' be the window values (bits s..t), and X, Xi as in (B).
- Put dX = w - bitlength|X| and dXi = w - bitlength|Xi| (carry depths). Put A = a - dXi and S = s - dX.

Then a >= 1, s >= 1 and

    |A - S| <= 2   or   max(A, S) <= 3.

In particular |a - s| <= 2 + max(dX, dXi), unless a <= 3 + dXi and s <= 3 + dX. Without carries (dX = dXi = 0) the
window's midpoint (n - 1 + s - a)/2 lies within 1 of the word's midpoint (n - 1)/2, unless a, s <= 3.

**(C') Converse symmetry (PROVED).** (v, v') -> (v*, v'*) is again a collision pair of reduced words, and it maps
(a, s, dX, dXi) to (s, a, dXi, dX). So it maps (A, S) to (S, A): the micro and macro agreement lengths exchange roles.

**(D) Randomised micro self-energies (PROVED, new).** Let n >= 12, m = floor(n/2) - 2, C = 2^m and P = 2^(n-m-1).
So B = P/C is 8 or 16.
- Let F(l), 0 <= l < C, be independent random variables, each uniform on {0, 1, ..., C - 1}.
- Define g_F(h 2^m + l) = g(h) + h C_m(l) + F(l) + [h odd][l < 2^(m-1)] for h in [P, 2P).
- By (A), F(l) = g(l) - 1 recovers g, since 0 <= g(l) - 1 <= C - 1.

Then

    E_F #{(v, v') in X_n^2 : g_F(v) = g_F(v')} <= 109 |X_n|.

So the macro lattice {g(h) + h c} creates no excess collisions. Linear energy E_n = O(|X_n|) is equivalent to the
following decorrelation criterion (DC), up to absolute constants. The number of quadruples with h != h', l != l' and

    g(l') - g(l) = h c - h' c' - (g(h') - g(h)) + ([h odd][l < 2^(m-1)] - [h' odd][l' < 2^(m-1)])

is O(|X_n|). In the random model its expectation is at most 76 |X_n|.

**(E) Reduction to reduced words (PROVED, elementary).**
- X_n is the disjoint union of {1^n} and the layers Y_(n-k) 1^k, for 0 <= k <= n - 2, and g(y 1^k) = g(y).
- g(Y_j) lies in [2^(j-2) + 1, 2^j - 1], so every value occurs in at most two layers.
- Hence E_n <= 1 + 2 sum_(j=2..n) E^Y_j.
- If E^Y_j <= K |Y_j| for all j <= n, then E_n <= 2K |X_n| + 1.

## Proof

**(A), (B).** These are Session 193's laws (1) and (2), with g(C_m(l)) = g(l) substituted, and were checked on 200,000
random cases (`check_balance.py`).

**(C).**
- **Window values.** Both words have bit n-1 = 1 and bit 0 = 0, so 1 <= s <= t <= n - 2, giving a, s >= 1.
  The common prefix h has top bit 1, so h is in [2^(a-1), 2^a). The common suffix sigma has bit 0 = 0, so tau = C_s(sigma)
  has top bit 1 and lies in [2^(s-1), 2^s).
- **X and Xi.** The window values differ at their bit 0 and their bit w-1. So X is odd. Xi is also odd, because bit 0
  of C_w(x) is 1 - x_(w-1). In particular XXi != 0. By definition of the carry depths,
  2^(w-1-dX) <= |X| < 2^(w-dX), and likewise for Xi. Hence

      2^(A+w-2) <= h|Xi| < 2^(A+w),   2^(S+w-2) <= tau|X| < 2^(S+w).

- **Size of the right side.** By (B) the collision means h Xi + tau X = g(x') - g(x) + eps' - eps.
  - Since 1 <= g(y) <= y + 1 <= 2^w for y < 2^w, we have |g(x') - g(x)| <= 2^w - 1.
  - Also |eps' - eps| <= 2.
  - So |h Xi + tau X| <= 2^w + 1.
- **Case X Xi > 0.** Then 2^(A+w-2) + 2^(S+w-2) <= 2^w + 1, so 2^(A-2) + 2^(S-2) <= 3/2. Hence A <= 2 and S <= 2.
- **Case X Xi < 0.** Suppose A >= S + 3.
  - Then tau|X| < 2^(S+w) <= 2^(A+w-3), so h|Xi| - tau|X| > 2^(A+w-3).
  - This forces 2^(A-3) < 1 + 2^(-w) <= 3/2, so A <= 3 and max(A, S) <= 3.
  - The case S >= A + 3 is symmetric.
  - Otherwise |A - S| <= 2.
- **The corollary.** It follows from |a - s| <= |A - S| + |dXi - dX|.

QED.

**(C').**
- The map v -> v* reverses bit positions (i -> n-1-i) and complements, so it preserves reducedness. The differing
  positions of (v*, v'*) are the images of those of (v, v'). Hence a* = s and s* = a.
- The window of v* is C_w(x), so X* = C_w(x) - C_w(x') = Xi and Xi* = x - x' = X.
- g is invariant under the converse.

QED.

**(D).** Order the pairs ((h, l), (h', l')) by type.
- **Diagonal.** These number |X_n| = PC.
- **h = h', l != l'.** Here |g_F(v) - g_F(v')| >= h|c - c'| - |F(l) - F(l')| - 1 >= P - C > 0, so there are none.
- **h != h', l = l' (deterministic, since F cancels).** Put e = [l < 2^(m-1)]. The condition is
  (h - h') c = g(h') - g(h) + e([h' odd] - [h odd]).
  - For each ordered (h, h') and each e there is at most one c.
  - So there are at most 2P^2 = 2B PC <= 32 PC such pairs.
- **h != h', l != l'.**
  - Here F(l') - F(l) has P(= t) <= 1/C for every t, and the event needs |h c - h' c' - D| <= C, where
    D = g(h') - g(h). The parity terms are at most 1 and |F(l') - F(l)| <= C - 1.
  - THM-4887's strip lemma with W = C applies, because h' >= P >= 8C > 2W. Each ordered (h, h') has at most
    Z(h, h') <= 2C^2/P + 3C/beta + 5 such (c, c'), where beta = beta(h, h') is the reduced Dirichlet denominator.
  - THM-4887 step 4 gives sum over h != h' of 1/beta <= P(4C ln B + 10.58 C + B + 1 + ln B).
  - So the expected number is at most

        (1/C)[P^2(2C^2/P + 5) + 3CP(4C ln B + 10.58 C + B + 1 + ln B)]
            = PC [2 + 5B/C + 3(4 ln B + 10.58) + 3(B + 1 + ln B)/C] <= 76 PC,

    using B <= 16, ln 16 < 2.78 and C >= 16 for n >= 12.

Total: (1 + 32 + 76) PC = 109 |X_n|.

For the equivalence: with F(l) = g(l) - 1 the same three counts give E_n. The first two types contribute at most
33 |X_n|, and the third is exactly the count in (DC).

QED.

**(E).**
- **Layers.** Stripping the trailing ones of a word that is not 1^n leaves a word with top bit 1 and bottom bit 0, of
  length at least 2. Appending ones does not change g, since g(2v+1) = g(v).
- **Range.** For v = 2u in Y_j with u in X_(j-1): g(v) = g(u) + u + (u mod 2) >= 1 + 2^(j-2), and g(v) <= v + 1 <= 2^j - 1.
- **Overlap.** A value N lies in (2^(j-2), 2^j) for at most two j, and the value 1 comes only from 1^n.
- **Bound.** By Cauchy-Schwarz over the at most two layers, r_n(N)^2 <= 2 sum_j r^Y_j(N)^2.
- **Summation.** Summing, E_n <= 1 + 2 sum_j E^Y_j. Then sum_(j=2..n) |Y_j| = 2^(n-1) - 1 < |X_n|.

QED.

## Remarks

- **What (C) says.** It is a statement about where collisions can live. An exact coincidence of one-apex counts
  forces the two words to agree on a prefix and a suffix of nearly equal lengths. The only exceptions are carries,
  where the window values straddle a power of two. This is the macro (prefix) and micro (suffix, read backwards)
  balance that (A) and (C') make symmetric.
- **Data (HYP-9482, `outofsample_21_26.out`).**
  - At n = 26 the non-converse pairs with prefix length a = 1..10 number 0.347, 0.260, 0.153, 0.083, 0.043, 0.022,
    0.0115, 0.0058, 0.0029, 0.0014 per |Y|. This halving level law reflects the balance.
  - Deep, central windows are rare. The bulk of collisions are global (a, s <= 3).
- **What (D) says.** The macro bilinear lattice is not the obstruction to linear energy. Any proof must show that the
  micro self-energies g(l), on about n/2 bits, do not conspire with the strips {|h c - h' c' - D| <= C}; this is
  criterion (DC).
  - On the strips, c and c' are tied by a Beatty-type relation c' ~ (h/h') c.
  - The micro values enter through bit reversal: g(l) ~ l C_m(l)/2^m.
  - So (DC) is a decorrelation statement between digit reversal and multiplication. It is not provided here.
