---
id: THM-5141
title: "Counting exponent 3/4 for the Hamiltonian-path spectrum. The one-apex words of n bits take at least c 2^(3n/4) / (n^(9/2) log n) distinct Hamiltonian path counts, so the number A(X) of Hamiltonian-path counts of tournaments up to X is at least c X^(3/4) / (log X)^5 (previous exponent: 2/3, THM-4887). Ingredients: the graph of the bit reversal {(x, rev x)} has all signed-digit differences distinct, so its additive energy is exactly 6^w; the gcd sum S_w = sum over window pairs of gcd(X, Xi)/max(|X|, |Xi|) is at most C w 2^(3w/2); a collision of two words with the same prefix and suffix is one lattice-line equation c X - P Xi = gamma in the outer blocks"
status: >
  DRAFT of opus-2026-11-11-S93, written before the independent verification of this session. Every statement below
  carries its own label; nothing here is to be cited as PROVED until the AUDIT line at the end of this field says so.
  INTENDED LABELS: (A) PROVED; (B1) PROVED (elementary); (B2) PROVED (uses the van der Corput second-derivative
  estimate and the Erdos-Turan inequality, both CITED); (C) PROVED, on the block law, which is INHERITED (Session 193's
  concatenation law, THM-4890 (A)), re-proved here in the q-normalisation; (D), (E) PROVED.
  NOT CLAIMED: linear energy, positive density of the spectrum, density one, Chan-Pak (all OPEN). The method cannot
  give an exponent above 3/4 (remark 3).
  AUDIT: pending.
session: opus-2026-11-11-S93
source: 05-knowledge/results/s93_chanpak_wave_20261111/README.md
related:
  - THM-4887
  - THM-4890
  - THM-4886
  - THM-4885
  - THM-4875
  - THM-4878
  - HYP-9662
verification:
  - 04-computation/experiments/s93_chanpak_wave_20261111/coordinator/reversal_energy_checks.py
  - 04-computation/experiments/s93_chanpak_wave_20261111/coordinator/exponent_lemma_tests.py
  - 04-computation/experiments/s93_chanpak_wave_20261111/coordinator/eval_w3_w4.py
---

# THM-5141 - counting exponent 3/4 for one-apex Hamiltonian path counts

## Objects

- **Counts.** q(0) = 0, q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2), and x(h) = 2h + 1 - 2 q(h).
  With g of THM-4890 (g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2)) one has g(v) = v + 1 - 2 q(v) - (v mod 2)
  and x(h) = g(2h) (both by induction on the recursions). So {x(h) : h >= 0} is the set of Hamiltonian path counts
  of the one-apex tournaments (THM-4875, THM-4878).
- **Sizes.** 0 <= q(h) <= (h - 1)/2 for h >= 1 (induction), so h + 2 <= x(h) <= 2h + 1.
- **Reversal.** rev_w(x) = sum over i < w of x_i 2^(w-1-i) for 0 <= x < 2^w (x_i = bit i of x).
- **A(X)** = number of integers N <= X that are the number of Hamiltonian paths of some tournament.
- **Signed strings.** For eps in {-1, 0, 1}^w: D(eps) = sum eps_i 2^i, D*(eps) = sum eps_i 2^(w-1-i), z(eps) = number
  of zero digits.
- **Difference multiplicities.** r_w(a, b) = #{(u, u') in [0, 2^w)^2 : u - u' = a, rev_w u - rev_w u' = b}.
- **Window pairs.** An ordered pair (x, x') of w-bit words is *odd* if X = x - x' and Xi = rev_w x - rev_w x' are both
  odd (for w >= 2: the two words differ in the top bit and in the bottom bit). Its *height* is
  M = max(|X|, |Xi|)/gcd(X, Xi).
- **The gcd sum.** S_w = sum over odd pairs of gcd(X, Xi)/max(|X|, |Xi|) = sum over odd pairs of 1/M.

## Statement

**(A) Reversal graph (PROVED).**
1. No nonzero delta in {-2, ..., 2}^w satisfies sum delta_i 2^i = 0 and sum delta_i 2^(w-1-i) = 0.
2. eps -> (D(eps), D*(eps)) is injective on {-1, 0, 1}^w; r_w(D(eps), D*(eps)) = 2^z(eps), and r_w = 0 off the image.
3. sum over (a, b) of r_w(a, b)^2 = 6^w: the additive energy of {(x, rev_w x)} is exactly 6^w.
4. sum over (a, b) with a != 0 of r_w(a, b)^2 |a|^(-1/2) <= 42 * 6^w 2^(-w/2).

**(B1) The gcd sum, elementary (PROVED).** S_w <= C w (6 sqrt 2)^(w/2) = C w 2^(beta0 w), beta0 = (log2 6)/2 + 1/4 =
1.54248....

**(B2) The gcd sum (PROVED; van der Corput and Erdos-Turan CITED).** S_w <= C w 2^(3w/2).
(Numerically S_w is between 0.82 and 1.06 times 6^(w/2) for 4 <= w <= 15.)

**(C) Collision equation (PROVED on the INHERITED block law).** Let h != h' be n-bit words (top bit 1) with
x(h) = x(h'). Let t' >= 1 and s >= 0 be the lengths of their longest common prefix and suffix, w = n - t' - s, and
h = P | x | Z, h' = P | x' | Z with windows x, x' of w bits. Then (x, x') is an odd pair and, with
c = 2^(s+1) - rev_s(Z) in (2^s, 2^(s+1)] and P in [2^(t'-1), 2^t'),

    c X - P Xi = gamma,   gamma = 2 (q(x) - q(x')) + (P mod 2)(x_(w-1) - x'_(w-1)) + (c mod 2)(x_0 - x'_0),

and |gamma| <= 2^w. For a fixed odd pair the number of such (P, c) is at most
4 (1 + 2^max(t'-1, s) gcd(X, Xi)/max(|X|, |Xi|)). If neither word has a run of more than R equal bits, then
|s - t'| <= R + 3.

**(D) Energy of prefix families (PROVED).** Let n >= 8, R = ceil(log2 n) + 1, W = the set of n-bit words without a
run of more than R equal bits (|W| >= 2^(n-2)), and t = ceil(n/4). The number of ordered pairs h != h' in W with
x(h) = x(h') and the same first t bits is at most C n^(9/2) log n 2^n.

**(E) Counting (PROVED).** For n >= 8 the words of n bits take at least c 2^(3n/4)/(n^(9/2) log n) distinct values of x.
Hence for X >= 3

    A(X) >= #{x(h) : x(h) <= X} >= c X^(3/4) / (log X)^5.

With (B1) in place of (B2) the same proof gives the exponent beta0/(2 beta0 - 1) = 0.7398... in place of 3/4.

C and c denote positive absolute constants (not the same at each occurrence).

## Proof

### (A)

**1.** Let delta(z) = sum delta_i z^i. The hypotheses say delta(2) = 0 and delta(1/2) = 0. Division by the monic
z - 2 gives delta = (z - 2) delta_1 with delta_1 in Z[z], and delta_1(1/2) = 0. Since 2z - 1 is primitive, Gauss's
lemma gives delta_1 = (2z - 1) p with p in Z[z]. So delta = (2z^2 - 5z + 2) p. Suppose p != 0, with lowest nonzero
coefficient p_m and degree d; put p_i = 0 outside [m, d] and u_i = 2 p_i - p_(i-1). The coefficient of z^i in delta is

    delta_i = 2 p_i - 5 p_(i-1) + 2 p_(i-2) = u_i - 2 u_(i-1),   so   u_i = 2 u_(i-1) + delta_i.

Now |u_m| = 2 |p_m| >= 2, and |u_(i-1)| >= 2 implies |u_i| >= 2 |u_(i-1)| - 2 >= 2. So |u_i| >= 2 for all i >= m.
But u_(d+2) = 2 p_(d+2) - p_(d+1) = 0. Contradiction; hence p = 0 and delta = 0.

**2.** For a pair (u, u') put eps_i = u_i - u'_i. Then u - u' = D(eps) and rev_w u - rev_w u' = D*(eps). The pairs
with a given eps are obtained by choosing (0,0) or (1,1) at each zero digit: 2^z(eps) pairs. If
(D, D*)(eps) = (D, D*)(eps') then delta = eps - eps' is in {-2, ..., 2}^w and satisfies the two equations of 1, so
eps = eps'.

**3.** By 2, sum r_w^2 = sum over eps of 4^z(eps) = (4 + 1 + 1)^w.

**4.** By 2 the sum is sum over eps != 0 of 4^z(eps) |D(eps)|^(-1/2) (D(eps) != 0 for eps != 0: look at the lowest
nonzero digit). For 0 <= j <= w let F_j = sum of 4^z(eps) over the eps with |D(eps)| < 2^(w-j).
Claim: F_j <= (5/3) (2/3)^j 6^w. Write D(eps) = 2^(w-j) T + L, where T is the value of the top j digits (as a
j-digit signed string) and |L| <= 2^(w-j) - 1. If |D(eps)| < 2^(w-j) then |T| <= 1. The j-digit signed strings of
value 0: only the zero string (weight 4^j). Of value 1: exactly the strings (-1, ..., -1, 1, 0, ..., 0) with i >= 0
digits -1 below a digit 1 (induction on the lowest digit: it is odd; if it is 1 the rest has value 0; if it is -1
the rest has value 1), of total weight sum over i < j of 4^(j-1-i) = (4^j - 1)/3. Value -1: the same. So the top
digits have total weight < (5/3) 4^j, the lower digits at most 6^(w-j), which proves the claim. Since 1 <= |D| < 2^w,

    sum 4^z |D|^(-1/2) <= sum over j = 0..w-1 of 2^(-(w-j-1)/2) F_j
                       <= (5 sqrt 2 / 3) 6^w 2^(-w/2) sum over j >= 0 of (2 sqrt 2 / 3)^j < 42 * 6^w 2^(-w/2).

### The half split (used in (B1), (B2))

Let w >= 4, s1 = ceil(w/2), s2 = floor(w/2), ds = s1 - s2 in {0, 1}. Write x = v 2^s1 + y (0 <= y < 2^s1,
0 <= v < 2^s2) and u = rev_s1(y). Then rev_w(x) = u 2^s2 + rev_s2(v), and x -> (u, v) is a bijection. For a pair
(x, x') put

    a = u - u',  b = y - y' = rev_s1 u - rev_s1 u',  a' = v - v',  b' = rev_s2 v - rev_s2 v'.

Then X = a' 2^s1 + b and Xi = a 2^s2 + b'; the pair is odd iff b and b' are odd; the number of pairs with given
(a, b, a', b') is r_s1(a, b) r_s2(a', b'). If b is odd then a != 0, and if b' is odd then a' != 0.

For an odd pair let (k, k') = (X, Xi)/j with j = +-gcd(X, Xi) and k' >= 1; k and k' are odd and coprime, and
M = max(|k|, k'). For K >= 1 let N(K) be the number of odd pairs with M <= K. For P_u = (a, b) with b odd let
Deg_K(P_u) be the number of (a', b') with |a'|, |b'| < 2^s2, b' odd, and M <= K; define Deg'_K(P_v) symmetrically
(the number of (a, b) with |a|, |b| < 2^s1, b odd, M <= K). Since r r' <= (r^2 + r'^2)/2,

    N(K) <= (1/2) sum over P_u of r_s1(P_u)^2 Deg_K(P_u) + (1/2) sum over P_v of r_s2(P_v)^2 Deg'_K(P_v).      (1)

Since a pair with M in [2^i, 2^(i+1)) contributes at most 2^(-i),

    S_w <= sum over i >= 0 of 2^(-i) N(2^(i+1)),   and trivially N(K) <= 4^w.                                    (2)

**(D1)** Deg_K <= 4 K^2 and Deg'_K <= 8 K^2. Fix (k, k') (at most 2K^2 choices). Then (X, Xi) = j (k, k') with
j k = b (mod 2^s1), so j lies in one class mod 2^s1; and j k' = Xi lies in the open interval of length 2^(s2+1)
around a 2^s2, so j lies in an open interval of length at most 2^(s1+1): at most 2 values of j, each determining
(a', b'). For Deg': j k' = b' (mod 2^s2) and j k lies in an open interval of length 2^(s1+1) <= 4 * 2^s2: at most
4 values.

**(D2)** For 1 <= K <= 2^s2: Deg_K and Deg'_K are at most 3 * 2^s1 (3 + 2 ln K). Every (a', b') counted by Deg_K
arises from a triple (b', k', k): b' with |b'| < 2^s2, a divisor k' <= K of Xi = a 2^s2 + b', j = Xi/k', and k with
|k| <= K and j k = b (mod 2^s1). Given (b', k') there are at most 3 such k (a class mod 2^s1 in an interval of length
2K <= 2^(s1+1)); given k' there are at most 2^(s2+1)/k' + 1 values of b'. So
Deg_K <= 3 sum over k' <= K of (2^(s2+1)/k' + 1) <= 3 (2^(s2+1) (1 + ln K) + K). For Deg' exchange the roles: b with
|b| < 2^s1, a divisor k <= K of |X|, then at most 3 values of k' (a class mod 2^s2 in an interval of length
2K <= 2^(s2+1)): Deg'_K <= 3 (2^(s1+1) (1 + ln K) + K).

### (B1)

By (1), (A)3, (D1), (D2): N(K) <= 6^s1 min(8 K^2, 3 * 2^s1 (3 + 2 ln K)) for 1 <= K <= 2^s2. In (2) use this for
i + 1 <= s2 and N <= 4^w for i >= s2. With i0 = floor(s1/2):

    S_w <= 6^s1 [sum over i <= i0 of 32 * 2^i + sum over i > i0 of 2^(-i) 3 * 2^s1 (3 + 2 s2)] + 2 * 4^w 2^(-s2)
        <= (6 sqrt 2)^s1 (82 + 12 s2) + 2^(3w/2 + 2)  <=  C w (6 sqrt 2)^(w/2).

### (B2)

**The CRT form.** Fix P_u = (a, b), b odd, and a direction (k, k'). If (a', b') has this direction then
k' X = k Xi, that is 2^s2 (2^ds k' a' - k a) = k b' - k' b. Put e = 2^ds k' a' - k a. Then k b' - k' b = 2^s2 e and

    |e| < |k| + 2^ds k' <= 3K,   e = -k a (mod 2^ds k'),   2^s2 e = -k' b (mod |k|).

The moduli 2^ds k' and |k| are coprime, so e lies in one class e0 modulo m = 2^ds k' |k|, and e determines
(a', b'). Hence the number of (a', b') with direction (k, k') is at most 6K/m + 1, and it is 0 unless
||e0/m|| < 3K/m (|| . || = distance to the nearest integer). By the Chinese remainder theorem

    e0/m = -sgn(k) a/(2^ds k') - b inv(2^s1; |k|)/|k|   (mod 1),

where inv(2^s1; |k|) is the inverse of 2^s1 modulo |k|. So with theta_k = sgn(k) b inv(2^s1; |k|)/|k|:

    Deg_K(P_u) <= sum over (k, k') of 6K/(|k| k') + #{(k, k') : ||a/(2^ds k') + theta_k|| < 3K/(2^ds |k| k')}
               <= 12 K (1 + ln K)^2 + sum over odd k, |k| <= K, of #{odd k' <= K : ||a/(2^ds k') + theta_k|| < 3K/(|k| k')}.

**Discrepancy lemma.** Let a != 0 be an integer, ds in {0, 1}, L >= 1, theta real and 0 < dl. Then

    #{odd k' in (L, 2L] : ||a/(2^ds k') + theta|| <= dl} <= 2 dl L + C (|a|^(1/3) + L^(3/2) |a|^(-1/2) + 1 + ln L).

Proof. If L < |a|^(1/3) the count is at most L. Otherwise let F_h(l) = h a/(2^ds (2l+1)) for the l with 2l + 1 in
(L, 2L] (at most L values). F_h'' = 8 h a/(2^ds (2l+1)^3) has constant sign and lam <= |F_h''| < 8 lam with
lam = h |a|/(2^ds L^3). The van der Corput estimate (Kuipers-Niederreiter, Uniform Distribution of Sequences,
chapter 1, Theorem 2.7: |sum e(f(n))| <= (|f'(b) - f'(a)| + 2)(4/sqrt(rho) + 3) when |f''| >= rho) gives, with
|F'(b) - F'(a)| <= 8 lam L and the trivial bound when lam >= 1,

    |sum over l of e(F_h(l))| <= 56 L lam^(1/2) + 8 lam^(-1/2) + 6.

The Erdos-Turan inequality (the same book, chapter 2, Theorem 2.5) bounds the number of points of a sequence of N
reals in an interval of length 2 dl modulo 1 by 2 dl N + C N/H + C sum over h <= H of (1/h) |sum e(h u_n)| for
every integer H >= 1. With H = ceil(L |a|^(-1/3)):

    count <= 2 dl L + C [L/H + L (H |a|/L^3)^(1/2) + (L^3/|a|)^(1/2) + 1 + ln H]
          <= 2 dl L + C [|a|^(1/3) + L^(3/2) |a|^(-1/2) + 1 + ln L].

**The degree bound (D3).** Split the k' <= K into k' = 1 and the blocks (L, 2L], L = 1, 2, 4, ... < K, and apply
the lemma with dl = 3K/(|k| L) for each of the at most 2K values of k:

    Deg_K(P_u) <= C [K (1 + ln K)^2 + K (1 + ln K) |a|^(1/3) + K^(5/2) |a|^(-1/2)].

For Deg'_K(P_v), P_v = (a', b') with b' odd: put e = k a - 2^ds k' a'; then k' b = k b' + 2^s2 e,
|e| < |k| + 2^ds k', e = -2^ds k' a' (mod |k|), 2^s2 e = -k b' (mod k'), and the same computation gives the condition
||2^ds a'/k + theta'_(k')|| < 3K/(|k| k') with theta'_(k') = b' inv(2^s2; k')/k'. For fixed k' and each sign of k
the lemma applies to the sequence (+-2^ds a')/|k| in |k| (the lemma with the integer +-2^ds a' in place of a and
ds = 0), so Deg'_K satisfies the same bound with |a'| in place of |a| (and another absolute constant).

**The sum.** By (1), (A)3, (A)4 (applied at lengths s1 and s2, using |a| < 2^s1, |a'| < 2^s2) and (D3):

    N(K) <= C 6^s1 [K (1 + ln K)^2 + K (1 + ln K) 2^(s1/3) + K^(5/2) 2^(-s2/2)]            for all K >= 1,
    N(K) <= 3 * 6^s1 2^s1 (3 + 2 ln K)                                                    for 1 <= K <= 2^s2   (D2).

Since min(A K^(3/2), B/K) <= A^(2/5) B^(3/5), for 1 <= K <= 2^s2:

    N(K)/K <= C 6^s1 [(1 + s1)^2 + (1 + s1) 2^(s1/3) + (1 + s1)^(3/5) 2^(2 s1/5)] <= C (1 + s1)^2 (6 * 2^(2/5))^s1 <= C 8^s1,

because 6 * 2^(2/5) = 7.917... < 8. By (2),

    S_w <= sum over i < s2 of 2 C 8^s1 + sum over i >= s2 of 2^(-i) 4^w <= C w 2^(3w/2).

### (C)

**Block law.** For A >= 0, m >= 1, 0 <= B < 2^m:
2 q(A 2^m + B) = 2 q(A) + 2 q(B) + A rev_m(B) + (A mod 2) B_(m-1).
Proof: unrolling the recursion, q(h) = sum over i of h_i ceil(floor(h/2^(i+1))/2). For h = A 2^m + B the bits
i >= m contribute q(A). For i < m - 1 the inner floor is A 2^(m-i-1) + floor(B/2^(i+1)) with A 2^(m-i-1) even, and
for i = m - 1 it is A; so the bits i < m contribute q(B) + (A/2) sum over i < m of B_i 2^(m-1-i) + (A mod 2) B_(m-1)/2.

**Three blocks.** Applying the law twice to h = P 2^(w+s) + x 2^s + Z (0 <= x < 2^w, 0 <= Z < 2^s) gives

    x(h) = K(P, Z) + c x - P rev_w(x) - 2 q(x) - (P mod 2) x_(w-1) - (c mod 2) x_0,   c = 2^(s+1) - rev_s(Z),

where K(P, Z) does not depend on x (for s >= 1 the last bit Z_(s-1) equals c mod 2; for s = 0, c = 2). Z -> c is a
bijection from [0, 2^s) onto the integers of (2^s, 2^(s+1)].

**The equation.** For h, h' as in the statement the windows differ in their top bits and in their bottom bits (by
maximality of the common prefix and suffix; t' >= 1 because both top bits are 1, and w >= 1 because h != h'), so X
and Xi are odd. Subtracting the two three-block formulas gives c X - P Xi = gamma. |gamma| <= 2 (2^(w-1) - 1) + 2.

**Lattice line.** Fix the odd pair and the parities of P and c (4 choices). Then gamma is constant, and two
solutions (c0, P0), (c1, P1) satisfy (c1 - c0, P1 - P0) = j (Xi, X)/g, g = gcd(X, Xi), j an integer. Two values of P
in [2^(t'-1), 2^t') differ by at most 2^(t'-1) - 1, two values of c in (2^s, 2^(s+1)] by at most 2^s - 1. So the number
of solutions is at most 1 + g min(2^(t'-1)/|X|, 2^s/|Xi|) <= 1 + 2^max(t'-1, s) g/max(|X|, |Xi|).

**Carries.** Suppose 2^(w-rho) <= |X|, |Xi| < 2^w. Then |Xi| P <= |X| c + |gamma| < 2^(w+s+2) and
|Xi| P >= 2^(w-rho+t'-1), so t' <= s + rho + 2. And |X| c <= |Xi| P + |gamma| < 2^(w+t') while |X| c > 2^(w-rho+s),
so s <= t' + rho - 1. Hence |s - t'| <= rho + 2.

**Run-free words.** Suppose neither word has a run of more than R equal bits; we show |X|, |Xi| >= 2^(w-1-R), so
that rho = R + 1. If w <= R + 1 this holds because X and Xi are odd. Let w >= R + 2 and say x_(w-1) = 1,
x'_(w-1) = 0. The R + 1 bits of x below the top bit are not all 0, so x >= 2^(w-1) + 2^(w-2-R); the R + 1 bits of x'
below its top bit are not all 1, so x' <= 2^(w-1) - 1 - 2^(w-2-R). Hence X > 2^(w-1-R). The same argument applies to
rev_w x and rev_w x', whose top bits are x_0 != x'_0 and whose next R + 1 bits are x_1, ..., x_(R+1).

### (D)

The words of n bits with a run of at least R + 1 equal bits number at most (n - R) 2^(n-1-R) <= 2^(n-2) (2^R >= 2n),
so |W| >= 2^(n-2). Let E be the number of ordered pairs in the statement. Each has parameters (t', s) with t' >= t
and, by (C), |s - t'| <= R + 3; there are at most n (2R + 7) such (t', s). For fixed (t', s) the pairs inject into
the quadruples (x, x', P, c) of (C), so by the lattice-line bound their number is at most

    4 * 4^w + 4 * 2^max(t'-1, s) S_w,    w = n - t' - s <= n - 2t' + R + 3,    max(t' - 1, s) <= t' + R + 3.

Since n <= 4t <= 4t': 4^w <= 4^(R+3) 4^(n-2t') <= 4^(R+3) 2^n, and by (B2)
2^(t'+R+3) S_w <= C n 2^(t' + R + 3 + (3/2)(n - 2t' + R + 3)) <= C n 2^((5/2)(R+3)) 2^n. With 2^R < 4n this gives
E <= n (2R + 7) [4^(R+4) + C n (32 n)^(5/2)] 2^n <= C n^(9/2) log n 2^n.

### (E)

For a t-bit prefix p let W_p be the words of W with prefix p, E_p the number of ordered pairs h != h' in W_p with
x(h) = x(h'), and d_p the number of distinct values of x on W_p. By Cauchy-Schwarz |W_p|^2 <= d_p (|W_p| + E_p), and
again by Cauchy-Schwarz sum over p of d_p >= |W|^2/(|W| + E), since sum E_p = E. A value is taken in at most 2^(t-1)
prefix classes, so the number of distinct values of x on W is at least

    |W|^2 / (2^(t-1) (|W| + E)) >= 2^(2n-4) / (2^(t-1) C n^(9/2) log n 2^n) >= c 2^(3n/4) / (n^(9/2) log n)

(t <= n/4 + 1). These values are at most 2^(n+1) - 1. For X >= 2^9 take n with 2^(n+1) <= X < 2^(n+2); for
3 <= X < 2^9 adjust c.

With (B1): S_w <= C w 2^(beta0 w). Take t = ceil(n (beta0 - 1)/(2 beta0 - 1)) (>= n/4); then
t' + beta0 (n - 2t') <= n for t' >= t, the bound E <= C n^(9/2) log n 2^n holds as before, and the number of values
is at least c 2^(n-t)/(n^(9/2) log n), n - t >= n beta0/(2 beta0 - 1) - 1.

## Remarks

1. **Where the gain is.** With the trivial bound S_w <= 4^w the proof of (D), (E) gives the exponent 2/3 of
   THM-4887. The gain comes from the gcd sum: the progression term of the lattice line is large only when the
   difference vector (X, Xi) is a large multiple of a short vector, and (A)3 says that pairs of points of the
   reversal graph are rarely parallel to a given short vector (at most C 6^(w/2) pairs for each direction).
2. **Numerics.** S_w/6^(w/2) = 0.83, 1.02, 0.90, 1.06, 0.92, 1.06, 0.91, 1.05, 0.87, 1.01, 0.84, 0.97 for
   w = 4, ..., 15 (`exponent_lemma_tests_t3_15.out`; the same values for w <= 11 in `reversal_energy_checks_11.out`).
   The pairs of height 1 (X = +-Xi) carry (2/3) 6^(w/2) of it for even w and 0.816 * 6^(w/2) for odd w. The
   pre-registered W4 of this session (the same sum over pairs with one odd coordinate, tau = 13, 14) was confirmed.
3. **The cap.** The term 4 * 4^w in (D) counts one possible solution for every odd pair; it is at most 2^n only if
   n <= 4t. Beating 3/4 needs a proof that most window pairs have no solution (P, c) at all, which is the linear
   energy problem (THM-4890 (DC)), OPEN.
4. **Balance.** The inequality |s - t'| <= rho + 2 is the carry-corrected balance of THM-4890 (C) in the present
   normalisation, with weaker constants; it is re-proved here from the collision equation to keep the proof
   self-contained.

## Receipts

- `reversal_energy_checks_11.out`: (A)2 for w <= 14 (all 3^w strings), (A)3 exactly for w <= 11, the largest number
  of pairs in one direction (at most 1.53 * 6^(w/2)).
- `exponent_lemma_tests_13_22.out`: the block law; the collision equation, |gamma| <= 2^w and the oddness of X, Xi on
  all 2,026,664 unordered colliding pairs at n = 22 (and n = 8, ..., 20); |s - t'| <= R + 3 for run-free pairs
  (observed maximum R + 1); the half-split identities; (A)4 (observed constant 4.4 at w = 12, bound 42).
- `exponent_lemma_tests_t3_15.out`: S_w and the histogram of heights for w <= 15.
- `eval_w3_prefix_28.out`, `eval_w3_prefix_30.out`: colliding pairs per word with a common prefix of length at least
  t (0.147, 0.019, 0.0048 at t = 4, 7, 9), far below the bound of (D).
