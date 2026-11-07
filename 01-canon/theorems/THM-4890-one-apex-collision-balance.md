---
id: THM-4890
title: "Micro/macro balance of one-apex collisions. Let two distinct reduced words of length n (top bit 1, bottom bit 0) have the same one-apex Hamiltonian path count. They differ only on a window [s, t]. The window of one word starts and ends with 1, and the other's with 0. With A = a - dXi and S = s - dX, the carry-corrected common-prefix and common-suffix lengths, |A - S| <= 1, so |a - s| <= 1 + |dX - dXi|. The converse involution exchanges A and S. The energy of all words is at most 1 + 2 sum_j E^Y_j over reduced words. With the bottom-half micro self-energies replaced by independent uniform values, the expected energy is at most 109 |X_n| (n >= 12; the exact value is 1.64-1.99 |X_n|). Linear energy (OPEN) is equivalent, up to constants, to criterion (DC): the collisions with different macro and different micro blocks number O(|X_n|). In the random model their expectation is at most 76 |X_n|. The true count is 3.3-3.7 times the model value at n = 12..22, with about |X_n| of it forced by the converse involution. So (DC) is a bounded-correlation statement, not a decorrelation statement."
status: >
  PROVED (new):
  (C) the sharp balance lemma - the proof is the independent audit's upgrade U1, which supersedes the session's first
  version (|A - S| <= 2 or max(A, S) <= 3);
  (C') converse symmetry; (E) the reduction to reduced words;
  (D) the random-model bound - a new statement, with method and constants from THM-4887 steps 3-5.
  INHERITED:
  (A) is Session 193's concatenation law (1)
  (05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md), written symmetrically with the
  converse invariance g(C_m(l)) = g(l) (THM-4886's converse involution);
  (B) is Session 193's (2);
  the local formula (L0) refines THM-4885 (A) (1 <= g(v) - v v*/2^n <= 1 + #10(v));
  1 <= g(l) <= l + 1 is Session 158's (binary_bank/NOTE.md, quoted in THM-4887).
  FINITE-EXACT:
  (C), (C') and (E) on every collision pair of reduced words for n = 4..28 (audit, aud_pairs.c, 258.6 M pairs).
  A - S takes only the values -1, 0, 1.
  The same-sign case never occurs (proved in (C)).
  HYP-9482 P5 is a corollary of (C).
  NOT CLAIMED: linear energy, positive density, density one, Chan-Pak (all OPEN). In particular (D) is a statement
  about a random model only; the true micro values are positively correlated with the macro strips.
  AUDIT: PASS WITH CORRECTIONS (2 majors on reading, 13 minors, 4 upgrades; all applied),
  05-knowledge/results/one_apex_collision_balance_20261106/audit/AUDIT.md.
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
  - 05-knowledge/results/one_apex_collision_balance_20261106/dc_ratio.py
  - 05-knowledge/results/one_apex_collision_balance_20261106/audit/aud_pairs.c
  - 05-knowledge/results/one_apex_collision_balance_20261106/audit/aud_sharp.py
  - 05-knowledge/results/one_apex_collision_balance_20261106/audit/aud_randmodel.c
---

# THM-4890 - micro/macro balance of one-apex collisions

## Objects

- **The counts.** g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2). The one-apex tournament of the word v has
  g(v) Hamiltonian paths (THM-4875, THM-4878).
- **Words.** X_n = [2^(n-1), 2^n). The reduced words are Y_n = {v in X_n : v even}, with |Y_n| = 2^(n-2).
- **Reversed complement.** For a w-bit word x, C_w(x) = sum over i < w of (1 - x_i) 2^(w-1-i).
  - C_w is an involution on [0, 2^w).
  - g(C_w(x)) = g(x), because C_w(x) is the word of the converse tournament.
  - Write v* = C_n(v). The map v -> v* sends Y_n onto Y_n.
- **Energies.** E_n = #{(v, v') in X_n^2 : g(v) = g(v')}, and E^Y_n is the same count on Y_n.

## Statement

**(A) Two-block law, symmetric form (INHERITED).** For h >= 1, m >= 1 and 0 <= l < 2^m, with c = C_m(l):

    g(h 2^m + l) = g(h) + g(c) + h c + [h odd][c odd] - 1.

Here h is the macro block and c is the micro block read backwards. The converse v -> v* exchanges them.

**(B) Window identity (INHERITED form).** Write v = h | x | sigma, where the prefix h has a bits, the window x has
w >= 1 bits and the suffix sigma has s bits. Put tau = C_s(sigma). Then

    g(v) = g(h) + g(x) + g(tau) + h tau 2^w + h C_w(x) + tau x + eps - 2,
    eps = [x odd][tau odd] + [h odd][x_(w-1) = 0],

so for v' = h | x' | sigma, with X = x - x' and Xi = C_w(x) - C_w(x'):

    g(v) - g(v') = (g(x) - g(x')) + h Xi + tau X + (eps - eps').

**(L0) Local formula (refines THM-4885 (A)).** Let z be a w-bit string, leading zeros allowed, with bits b_1 ... b_w
read from the most significant end. Then

    g(z) = z C_w(z)/2^w + 1 + L(z),   L = #10(z) - R(z),   R = sum over i < j with b_i = 0, b_j = 1 of 2^(i-j-1),

and -1 < L(z) <= #10(z).

**(C) Sharp balance lemma (PROVED; the audit's U1).** Let v != v' in Y_n with g(v) = g(v'). Define:
- t and s, the highest and lowest bit positions where v and v' differ;
- a = n - 1 - t (the common-prefix length) and w = t - s + 1 (the window width);
- x and x', the window values (bits s..t), and X, Xi as in (B);
- the carry depths dX = w - bitlength|X| and dXi = w - bitlength|Xi|;
- A = a - dXi and S = s - dX.

Then a, s >= 1, and the following hold.
1. **Opposite signs.** X Xi < 0. Equivalently, the word with the larger window value has window bits 1 at both ends
   of the window, and the other word has 0 at both ends.
2. **Balance.** |A - S| <= 1.
3. **Width 1.** If w = 1 then ceil(h/2) = ceil(tau/2).
4. **Corollaries.**
   - |a - s| <= 1 + |dX - dXi|.
   - Without carries (dX = dXi = 0), the window's midpoint (n - 1 + s - a)/2 lies within 1/2 of the word's midpoint
     (n - 1)/2.
5. **Sign rule.** Whenever the windows are 1...1 against 0...0, A >= S + 2 implies g(v) < g(v'), and S >= A + 2
   implies g(v) > g(v').

**(C') Converse symmetry (PROVED).** (v, v') -> (v*, v'*) is again a collision pair of reduced words. It maps
(a, s, dX, dXi) to (s, a, dXi, dX), and so (A, S) to (S, A). The micro and macro agreement lengths exchange roles.

**(D) Random micro self-energies (PROVED; new statement, method and constants from THM-4887 steps 3-5).**
Setup:
- Let n >= 12, m = floor(n/2) - 2, C = 2^m and P = 2^(n-m-1). So B = P/C is 8 or 16.
- Let F(l), for 0 <= l < C, be independent random variables, each uniform on {0, 1, ..., C - 1}.
- Define g_F(h 2^m + l) = g(h) + h C_m(l) + F(l) + [h odd][l < 2^(m-1)] for h in [P, 2P).
- By (A), the choice F(l) = g(l) - 1 recovers g, since 0 <= g(l) - 1 <= C - 1.

Then

    E_F #{(v, v') in X_n^2 : g_F(v) = g_F(v')} <= 109 |X_n|.

The exact expectation is 1.64-1.99 |X_n| for n = 12..22 (audit U3), so 109 is far from sharp.

In the random model, therefore, the macro lattice {g(h) + h c} creates no excess collisions. Linear energy is
equivalent, up to absolute constants, to the following criterion.

**(DC)** The number of quadruples with h != h', l != l' and

    g(l') - g(l) = h c - h' c' - (g(h') - g(h)) + ([h odd][l < 2^(m-1)] - [h' odd][l' < 2^(m-1)])

is O(|X_n|).

The random model's expectation of this count is at most 76 |X_n| (exactly 0.44-0.60 |X_n| at n = 12..22). For the true
micro values the count is 3.3-3.7 times the uniform model and 1.6-1.8 times a model with the empirical law of g(l) - 1.
About |X_n| of it is forced by the converse involution. The ratios are still slowly rising at n = 22. So (DC) asserts
bounded correlation, not decorrelation (audit MAJOR-1; `dc_ratio_*.out`, `audit/aud_randmodel.out`).

**(E) Reduction to reduced words (PROVED, elementary).**
- X_n is the disjoint union of {1^n} and the layers Y_(n-k) 1^k, 0 <= k <= n - 2, and g(y 1^k) = g(y).
- g(Y_j) is contained in [2^(j-2) + 1, 2^j - 1], so every value occurs in at most two layers.
- Hence E_n <= 1 + 2 sum over j = 2..n of E^Y_j.
- So if E^Y_j <= K |Y_j| for all j <= n, then E_n <= 2K |X_n| + 1.

## Proof

**(A), (B).** These are Session 193's laws (1) and (2), with g(C_m(l)) = g(l) substituted. Both were checked
exhaustively for small lengths and on random cases (`check_balance.py`, `audit/aud_blocks.out`).

**(L0).**
- Read z with Session 193's machine. Each 0 at position i adds A_i + b_(i-1), where A_i is the value of the prefix
  above it.
- The sum of the A_i over the zeros is the part of z C_w(z)/2^w that pairs a 1 above a 0. The sum of the b_(i-1) over
  the zeros is #10(z).
- Group the zeros into maximal runs. A run contributes [preceded by a 1] - R_run, with 0 <= R_run < 1. Only the
  leading run can lack a preceding 1.
- Hence -1 < L(z) <= #10(z).
- Checked exhaustively for w <= 16 (`audit/aud_sharp.out`).

**(C).**
- **Ranges.** Both words have bit n-1 = 1 and bit 0 = 0, so 1 <= s <= t <= n - 2.
  - The common prefix h has top bit 1, so h is in [2^(a-1), 2^a).
  - The common suffix ends in 0, so tau = C_s(sigma) is in [2^(s-1), 2^s).
  - The windows differ in their bits 0 and w-1. So X is odd, and Xi is odd because bit 0 of C_w(x) is 1 - x_(w-1).
  - Swap the words if necessary so that X > 0, i.e. x_(w-1) = 1 and x'_(w-1) = 0.
- **Step 1 (opposite signs).** Suppose Xi > 0, i.e. x_0 = 0 and x'_0 = 1.
  - Then x is even and at least 2^(w-1), so g(x) >= 1 + 2^(w-2).
  - x' is odd and below 2^(w-1), so g(x') = g((x' - 1)/2) <= 2^(w-2).
  - With eps = 0 and eps' = [tau odd] + [h odd], (B) gives h Xi + tau X = g(x') - g(x) + eps' - eps <= 1.
  - But the left side is at least h + tau >= 2. This contradiction proves 1.
- **Step 2 (w = 1).** Here x = 1, x' = 0 and g(1) = g(0) = 1.
  - (B) reduces to tau + [tau odd] = h + [h odd], which is 3.
  - So h and tau lie in a common pair {2k-1, 2k}, and their bit lengths a and s differ by at most 1.
  - Since dX = dXi = 0 here, this is 2.
- **Step 3 (w >= 2: A >= S + 2 is impossible).**
  - **Sizes.** By Step 1, x = 1...1 and x' = 0...0. So X >= 3, |Xi| >= 3 and both are odd, hence
    |Xi| >= 2^(w-1-dXi) + 1. Also X <= 2^k - 1 with k = w - dX. Put y = C_w(x), which is below 2^(w-1), and note
    x' < 2^(w-1).
  - **Identity.** Substituting (L0) into (B) gives

        g(v) - g(v') = T X - H |Xi| - delta,
        T = tau + y/2^w,   H = h + x'/2^w,   delta = L(x') - L(x) + [h odd] - [tau odd].

  - **Lower bound on H|Xi| - TX.** If A >= S + 2, then H|Xi| >= h|Xi| >= 2^(A+w-2) + 2^(a-1) >= 2^(s+k) + 2^(a-1),
    while T X < (2^s - 1/2)(2^k - 1). So H|Xi| - T X > 2^(a-1) + 2^s + 2^(k-1) - 1/2 >= 2^(k-1) + 5/2.
  - **Bound on delta.** -delta < #10(x) + 2 <= 3 + (k-1)/2. When dX >= 1, the window is x = 1 0^(dX-1) z, with z a
    k-bit word ending in 1.
  - **Conclusion.** g(v) - g(v') < 1/2 + (k-1)/2 - 2^(k-1) < 0. This contradicts g(v) = g(v'), and it also gives the
    sign rule 5.
- **Step 4 (S >= A + 2 is impossible).** Apply Step 3 to (v*, v'*). This is a collision pair with A and S exchanged,
  by (C'). Steps 2-4 give 2.
- **Corollaries.** |a - s| <= |A - S| + |dX - dXi|. Without carries, |s - a| <= 1.

QED (independent audit U1; every step checked on all 1,007,162 collision pairs for n <= 20 and on 300,000 random
configurations biased toward carries).

**(C').**
- v -> v* reverses positions (i -> n-1-i) and complements, so it preserves reducedness.
- The differing positions of (v*, v'*) are the images of those of (v, v'). So a* = s and s* = a.
- The window of v* is C_w(x), so X* = Xi and Xi* = X.
- g is converse-invariant.

QED.

**(D).** Sort the pairs ((h, l), (h', l')) into four types.
- **Diagonal.** These number |X_n| = PC.
- **h = h', l != l'.** Here |g_F(v) - g_F(v')| >= h|c - c'| - |F(l) - F(l')| - 1 >= P - C > 0, so there are none.
- **h != h', l = l'.** These are deterministic, since F cancels. Put e = [l < 2^(m-1)]. The condition is
  (h - h') c = g(h') - g(h) + e([h' odd] - [h odd]). For each ordered (h, h') and each e there is at most one c, so
  there are at most 2P^2 = 2B PC <= 32 PC such pairs.
- **h != h', l != l'.**
  - Here P(F(l') - F(l) = t) <= 1/C for every t, and the event needs |h c - h' c' - D| <= C, where D = g(h') - g(h).
  - THM-4887's strip lemma applies with W = C, since h' >= P >= 8C > 2W. It bounds the number of such (c, c') per
    ordered (h, h') by 2C^2/P + 3C/beta + 5.
  - THM-4887 step 4 bounds the sum over h != h' of 1/beta by P(4C ln B + 10.58 C + B + 1 + ln B).
  - So the expected count is at most

        PC [2 + 5B/C + 3(4 ln B + 10.58) + 3(B + 1 + ln B)/C] <= 76 PC,

    using B <= 16 and C >= 16. The bracket is at most 75.72, at n = 13 (audit).

Total: (1 + 32 + 76) PC = 109 |X_n|.

With F(l) = g(l) - 1, the same partition gives E_n exactly. The first three types contribute at most 33 |X_n|, and the
fourth is the count in (DC).

QED.

**(E).**
- Stripping the trailing ones of a word other than 1^n leaves a word with top bit 1 and bottom bit 0, of length at least
  2. Appending ones does not change g.
- For v = 2u in Y_j: g(v) = g(u) + u + (u mod 2) >= 1 + 2^(j-2), and g(v) <= v + 1 <= 2^j - 1.
- A value N lies in (2^(j-2), 2^j) for at most two j, and the value 1 comes only from 1^n.
- Hence r_n(N)^2 <= 2 sum over j of r^Y_j(N)^2. Summing over N, and using sum over j of |Y_j| < |X_n|, gives the
  claim. The audit checked it for n = 4..28.

QED.

## Remarks

- **What (C) says.** An exact coincidence of one-apex counts forces the two words to agree on a prefix and a suffix of
  nearly equal lengths. The two windows are mirror images in sign: 1...1 against 0...0. The macro block (prefix) and
  the micro block (suffix, read backwards) are balanced, and (A) and (C') make this symmetric.
- **Data (HYP-9482, `outofsample_21_26.out`).** At n = 26, the non-converse pairs with prefix length a = 1..10 number,
  per |Y|:

  | a | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | pairs | 0.347 | 0.259 | 0.153 | 0.083 | 0.043 | 0.022 | 0.0115 | 0.0058 | 0.0029 | 0.0014 |

  - This is a halving law.
  - Deep central windows are rare. The bulk of collisions is global (a, s <= 3), so restricting attention to deep
    windows cannot prove linear energy.
- **What (D) says.** In the random model the macro bilinear lattice adds no excess collisions. For the true micro values
  the (DC) count is several times larger: about |X_n| of it from the converse involution, plus a positive correlation
  of about 1.4 beyond it (README section 4b).
  - A proof of linear energy must therefore show that this correlation stays bounded.
  - On a strip, c' ~ (h/h') c is Beatty-type, and the micro values enter through bit reversal: g(l) ~ l C_m(l)/2^m.
  - So (DC) is a bounded-correlation statement about digit reversal and multiplication. It is not provided here.
