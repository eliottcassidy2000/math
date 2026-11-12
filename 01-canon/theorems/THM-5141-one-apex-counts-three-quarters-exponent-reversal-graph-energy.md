---
id: THM-5141
title: "Counting exponent 3/4 for the Hamiltonian-path spectrum. The one-apex words of n bits take at least c 2^(3n/4) / n^4 distinct Hamiltonian path counts, so the number A(Y) of Hamiltonian-path counts of tournaments up to Y is at least c Y^(3/4) / (log Y)^4 (previous exponent: 2/3, THM-4887 and Session 216). Ingredients: digits -> (x, rev x) is a Freiman isomorphism of order 2 from the Boolean cube onto the graph of the bit reversal (all signed-digit difference vectors are distinct), so that graph has additive energy 6^w; the gcd sum S_w = sum over window pairs of gcd(X, Xi)/max(|X|, |Xi|) is at most C w 2^(3w/2); a collision of two words with the same prefix and suffix is one lattice-line equation c X - P Xi = gamma in the outer blocks"
status: >
  PROVED: (A)1, (A)2, (A)4, (B1), (B2), (C), (D), (E). (A)3 is (A)2 plus the KNOWN energy of the Boolean cube.
  CITED inside (B2): the van der Corput second-derivative estimate and the Erdos-Turan inequality
  (Kuipers-Niederreiter, Uniform Distribution of Sequences, chapter 1 Theorem 2.7 and chapter 2 Theorem 2.5; both
  statements checked against the book by the audit).
  KNOWN / INHERITED, re-proved or quoted here: the block law (Session 193's concatenation law (1) = THM-4890 (A) in the
  g-normalisation; in the present normalisation it is Session 206's identity for Q = 2q,
  session206/parameter-capacity/audit/REVIEW.md, also THM-5093 and Session 231 source/PROOF.md (1)); the collision
  equation is THM-4890 (B) (Session 193's law (2)) in the q-normalisation; the balance |s - t'| <= r + 2 is weaker
  than THM-4890 (C); the identification of g(v) with the Hamiltonian path count of the one-apex tournament
  (Session 193 section 2; THM-4875, THM-4878), on which the bound for A(Y) rests; the value 6^w is the additive
  energy of the cube {0,1}^w, the extremal value |A|^(log2 6) of Kane-Tao (Electron. J. Combin. 24(2) (2017) P2.31);
  the reversal graph is the Hammersley (0, w, 2)-net; the number of signed-digit strings of a given value is a
  Stern diatomic number (Monroe 2021); the lattice-line count is the elementary count of integer points of a line
  in a box; the CRT step is the reciprocity of Kloosterman fractions.
  NEW (no prior art found in the repository or in the literature searched; a search is not a clearance): the
  exponent 3/4 (earlier bounds: THM-4887, X^(2/3)/(13 ln X + 102); Session 216 density/audit/REVIEWED_RESULT.md,
  X^(2/3)/(60 (ln X)^(2/3)), which lists "exact gcd counting in the balanced-window equation" as explored but
  unproved; Chan-Pak, arXiv:2308.10214, section 6.5, state only Conjecture 6.2 and that almost-completeness is not
  known); the 2-isomorphism (A)1-(A)2; the bounds (B1), (B2) for the gcd sum.
  NOT CLAIMED: linear energy, positive density of the spectrum, density one, Chan-Pak (all OPEN). As written the
  proof does not give an exponent above 3/4 (remark 3).
  AUDIT: PASS WITH CORRECTIONS. Three independent proof-readers (one per part; each re-derived every line of its
  part and the rest at referee depth) found no fatal or major issue; 14 minor corrections and 20 notes, all applied
  (among them: the case w <= 3, a wrong power of n in the variant with (B1), a receipt that said R + 1 for an
  observed R + 2, missing credits, notation with double meanings). Eight earlier numerical verifiers were
  interrupted by a usage limit before reporting; their scripts and outputs are complete and every finite check
  passed (04-computation/experiments/s93_chanpak_wave_20261111/exponent34_verify/).
session: opus-2026-11-11-S93
source: 05-knowledge/results/s93_chanpak_wave_20261111/README.md
related:
  - THM-4887
  - THM-4890
  - THM-4886
  - THM-4885
  - THM-4875
  - THM-4878
  - THM-5093
  - HYP-9662
verification:
  - 04-computation/experiments/s93_chanpak_wave_20261111/coordinator/reversal_energy_checks.py
  - 04-computation/experiments/s93_chanpak_wave_20261111/coordinator/exponent_lemma_tests.py
  - 04-computation/experiments/s93_chanpak_wave_20261111/coordinator/one_apex_bruteforce.py
  - 04-computation/experiments/s93_chanpak_wave_20261111/coordinator/degree_statistics.py
  - 04-computation/experiments/s93_chanpak_wave_20261111/exponent34_verify/
---

# THM-5141 - counting exponent 3/4 for one-apex Hamiltonian path counts

## Objects

- **Counts.** q(0) = 0, q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2), and eta(h) = 2h + 1 - 2 q(h) (called x(h) in the
  scripts). With g of THM-4890 (g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2)) one has
  g(v) = v + 1 - 2 q(v) - (v mod 2) and eta(h) = g(2h) (induction on the recursions). So {eta(h) : h >= 0} = {g(v)}
  is the set of Hamiltonian path counts of the one-apex tournaments (INHERITED: Session 193 section 2, THM-4875,
  THM-4878; brute force for spines of length at most 10 in `one_apex_bruteforce_10.out`).
- **Sizes.** 0 <= q(h) <= (h - 1)/2 for h >= 1 (induction), so h + 2 <= eta(h) <= 2h + 1.
- **Reversal.** rev_w(x) = sum over i < w of x_i 2^(w-1-i) for 0 <= x < 2^w (x_i = bit i of x).
- **A(Y)** = number of integers N <= Y that are the number of Hamiltonian paths of some tournament.
- **Signed strings.** For eps in {-1, 0, 1}^w: D(eps) = sum eps_i 2^i, D*(eps) = sum eps_i 2^(w-1-i), z(eps) = number
  of zero digits.
- **Difference multiplicities.** r_w(a, b) = #{(u, u') in [0, 2^w)^2 : u - u' = a, rev_w u - rev_w u' = b}.
- **Window pairs.** An ordered pair (x, x') of w-bit words (leading zeros allowed) is *odd* if X = x - x' and
  Xi = rev_w x - rev_w x' are both odd (for w >= 2: the words differ in the top bit and in the bottom bit). Its
  *height* is M = max(|X|, |Xi|)/gcd(X, Xi).
- **The gcd sum.** S_w = sum over odd pairs of gcd(X, Xi)/max(|X|, |Xi|) = sum over odd pairs of 1/M.
- C_0 and c_0 denote positive absolute constants (not the same at each occurrence); ex(t) = exp(2 pi i t);
  ||t|| is the distance from t to the nearest integer.

## Statement

**(A) Reversal graph.**
1. (PROVED) No nonzero delta in {-2, ..., 2}^w satisfies sum delta_i 2^i = 0 and sum delta_i 2^(w-1-i) = 0.
2. (PROVED) eps -> (D(eps), D*(eps)) is injective on {-1, 0, 1}^w; r_w(D(eps), D*(eps)) = 2^z(eps), and r_w = 0 off
   the image. Equivalently, digits -> (x, rev_w x) is a Freiman isomorphism of order 2 from {0,1}^w onto the
   reversal graph {(x, rev_w x)}: the vectors (2^i, 2^(w-1-i)) are 2-dissociated.
3. (PROVED from 2; the value is KNOWN for the cube) sum over (a, b) of r_w(a, b)^2 = 6^w.
4. (PROVED) sum over (a, b) with a != 0 of r_w(a, b)^2 |a|^(-1/2) <= 42 * 6^w 2^(-w/2).

**(B1) The gcd sum, elementary (PROVED).** For every w >= 1: S_w <= C_0 w (6 sqrt 2)^(w/2) = C_0 w 2^(beta0 w),
beta0 = (log2 6)/2 + 1/4 = 1.54248....

**(B2) The gcd sum (PROVED; van der Corput and Erdos-Turan CITED).** For every w >= 1: S_w <= C_0 w 2^(3w/2).
(Numerically S_w/6^(w/2) lies in [0.75, 1.06] for 4 <= w <= 21; remark 2.)

**(C) Collision equation (PROVED; the block law is KNOWN and re-proved below).** Let h != h' be n-bit words (top
bit 1) with eta(h) = eta(h'). Let t' >= 1 and s >= 0 be the lengths of their longest common prefix and suffix,
w = n - t' - s >= 1, and h = P | x | Z, h' = P | x' | Z with windows x, x' of w bits. Then (x, x') is an odd pair and,
with c = 2^(s+1) - rev_s(Z) in (2^s, 2^(s+1)] and P in [2^(t'-1), 2^t'),

    c X - P Xi = gamma,   gamma = 2 (q(x) - q(x')) + (P mod 2)(x_(w-1) - x'_(w-1)) + (c mod 2)(x_0 - x'_0),

and |gamma| <= 2^w. For fixed (t', s) and a fixed odd pair, the number of (P, c) in these ranges that satisfy the
equation (all four parity classes together) is at most 4 (1 + 2^max(t'-1, s) gcd(X, Xi)/max(|X|, |Xi|)). If neither
word has a run of more than R equal bits, then |s - t'| <= R + 3.

**(D) Energy of prefix families (PROVED).** Let n >= 8, R = ceil(log2 n) + 1, W = the set of n-bit words without a
run of more than R equal bits (|W| >= 2^(n-2)), and t = ceil(n/4). The number Coll of ordered pairs h != h' in W with
eta(h) = eta(h') and the same first t bits is at most C_0 n^4 2^n.

**(E) Counting (PROVED).** For n >= 8 the words of n bits take at least c_0 2^(3n/4)/n^4 distinct values of eta.
Hence for Y >= 3

    A(Y) >= #{eta(h) : eta(h) <= Y} >= c_0 Y^(3/4) / (log Y)^4.

With (B1) in place of (B2) the same proof gives the exponent beta0/(2 beta0 - 1) = 0.7398... in place of 3/4, with
the same powers of n and of log Y.

## Proof

### (A)

**1.** Let delta != 0 and let m be the lowest index with delta_m != 0. Reducing sum delta_i 2^i = 0 modulo 2^(m+1)
shows that delta_m is even, so delta_m = +-2. Put U_i = sum over l <= i of delta_l 2^(i-l). Then U_m = delta_m,
U_i = 2 U_(i-1) + delta_i, and |U_(i-1)| >= 2 implies |U_i| >= 2 |U_(i-1)| - 2 >= 2. So |U_(w-1)| >= 2. But
U_(w-1) = sum delta_l 2^(w-1-l) = 0 by the second equation. Contradiction.

**2.** For a pair (u, u') put eps_i = u_i - u'_i. Then u - u' = D(eps) and rev_w u - rev_w u' = D*(eps). The pairs
with a given eps are obtained by choosing (0,0) or (1,1) at each zero digit: 2^z(eps) pairs. If
(D, D*)(eps) = (D, D*)(eps') then delta = eps - eps' is in {-2, ..., 2}^w and satisfies the two equations of 1, so
eps = eps'.

**3.** By 2, sum r_w^2 = sum over eps of 4^z(eps) = (4 + 1 + 1)^w. (This is the energy of the cube {0,1}^w
transported by the 2-isomorphism.)

**4.** By 2 the sum is sum over eps != 0 of 4^z(eps) |D(eps)|^(-1/2) (D(eps) != 0 for eps != 0: look at the lowest
nonzero digit). For 0 <= l <= w let F_l = sum of 4^z(eps) over the eps with |D(eps)| < 2^(w-l).
Claim: F_l <= (5/3) (2/3)^l 6^w. Write D(eps) = 2^(w-l) T + Lo, where T is the value of the top l digits (as an
l-digit signed string) and |Lo| <= 2^(w-l) - 1. If |D(eps)| < 2^(w-l) then |T| <= 1. The l-digit signed strings of
value 0: only the zero string (weight 4^l). Of value 1: exactly the strings (-1, ..., -1, 1, 0, ..., 0) with i >= 0
digits -1 below a digit 1 (induction on the lowest digit: it is odd; if it is 1 the rest has value 0; if it is -1
the rest has value 1), of total weight sum over i < l of 4^(l-1-i) = (4^l - 1)/3. Value -1: the same. So the top
digits have total weight < (5/3) 4^l and the lower digits at most 6^(w-l), which proves the claim. Since
1 <= |D| < 2^w,

    sum 4^z |D|^(-1/2) <= sum over l = 0..w-1 of 2^(-(w-l-1)/2) F_l
                       <= (5 sqrt 2 / 3) 6^w 2^(-w/2) sum over l >= 0 of (2 sqrt 2 / 3)^l < 42 * 6^w 2^(-w/2).

### The half split (used in (B1), (B2))

For 1 <= w <= 3 there are at most 4^(w-1) <= 16 odd pairs and each contributes at most 1, so S_w <= 16 and (B1), (B2)
hold. Let w >= 4, s1 = ceil(w/2), s2 = floor(w/2), ds = s1 - s2 in {0, 1}. Write x = v 2^s1 + y (0 <= y < 2^s1,
0 <= v < 2^s2) and u = rev_s1(y). Then rev_w(x) = u 2^s2 + rev_s2(v), and x -> (u, v) is a bijection. For a pair
(x, x') put

    a = u - u',  b = y - y' = rev_s1 u - rev_s1 u',  a' = v - v',  b' = rev_s2 v - rev_s2 v'.

Then X = a' 2^s1 + b and Xi = a 2^s2 + b'; the pair is odd iff b and b' are odd; the number of pairs with given
(a, b, a', b') is r_s1(a, b) r_s2(a', b'). If b is odd and r_s1(a, b) > 0 then a != 0 (a = 0 forces u = u', b = 0);
likewise a' != 0 if b' is odd.

For an odd pair let (k, k') = (X, Xi)/tau with tau = +-gcd(X, Xi) and k' >= 1; k and k' are odd and coprime, and
M = max(|k|, k'). For K >= 1 let N(K) be the number of odd pairs with M <= K. For a point p_u = (a, b) with
|a|, |b| < 2^s1 and b odd let Deg_K(p_u) be the number of (a', b') with |a'|, |b'| < 2^s2, b' odd, and M <= K; for
p_v = (a', b') with |a'|, |b'| < 2^s2 and b' odd let Deg'_K(p_v) be the number of (a, b) with |a|, |b| < 2^s1, b odd,
and M <= K. Since r r' <= (r^2 + r'^2)/2,

    N(K) <= (1/2) sum over p_u of r_s1(p_u)^2 Deg_K(p_u) + (1/2) sum over p_v of r_s2(p_v)^2 Deg'_K(p_v).      (1)

Since a pair with M in [2^i, 2^(i+1)) contributes at most 2^(-i),

    S_w <= sum over i >= 0 of 2^(-i) N(2^(i+1)),   and trivially N(K) <= 4^w.                                    (2)

**(Deg1)** Deg_K <= 4 K^2 and Deg'_K <= 8 K^2 for every K >= 1. Fix (k, k') (at most 2K^2 choices). Then
(X, Xi) = tau (k, k') with tau k = b (mod 2^s1), so tau lies in one class mod 2^s1 (k is odd); and tau k' = Xi lies in
the open interval of length 2^(s2+1) around a 2^s2, so tau lies in an open interval of length at most 2^(s1+1): at
most 2 values of tau, each determining (a', b'). For Deg': tau k' = b' (mod 2^s2) and tau k lies in an open interval
of length 2^(s1+1) <= 4 * 2^s2: at most 4 values. Hence, by (1) and (A)3, N(K) <= 8 K^2 6^s1 for every K >= 1.

**(Deg2)** For 1 <= K <= 2^s2: Deg_K and Deg'_K are at most 3 * 2^s1 (3 + 2 ln K). Every (a', b') counted by Deg_K
arises from a triple (b', k', k): b' odd with |b'| < 2^s2, a divisor k' <= K of Xi = a 2^s2 + b', tau = Xi/k' (odd,
since Xi is odd), and k with |k| <= K and tau k = b (mod 2^s1). Given (b', k') there are at most 3 such k (one class
mod 2^s1 in an interval of length 2K <= 2^(s1+1)); given k' there are at most 2^(s2+1)/k' + 1 values of b'. So
Deg_K <= 3 sum over k' <= K of (2^(s2+1)/k' + 1) <= 3 (2^(s2+1) (1 + ln K) + K). For Deg' exchange the roles: b odd
with |b| < 2^s1, a positive divisor kk <= K of |X| (X = a' 2^s1 + b), |tau| = |X|/kk, and the signed number
k'' = sgn(tau) k' in [-K, K], which lies in one class mod 2^s2 (|tau| k'' = Xi = b' (mod 2^s2), and |tau| is odd):
at most 3 values, and (b, kk, k'') determines (a, b). Given kk there are at most 2^(s1+1)/kk + 1 values of b. So
Deg'_K <= 3 (2^(s1+1) (1 + ln K) + K). Since K <= 2^s2 <= 2^s1, both bounds are at most 3 * 2^s1 (3 + 2 ln K), and by
(1) and (A)3, N(K) <= 3 * 6^s1 2^s1 (3 + 2 ln K) for 1 <= K <= 2^s2.

### (B1)

By (Deg1), (Deg2): N(K) <= 6^s1 min(8 K^2, 3 * 2^s1 (3 + 2 ln K)) for 1 <= K <= 2^s2. In (2) use this for
i + 1 <= s2 (so ln K <= s2) and N <= 4^w for i >= s2. With i0 = floor(s1/2) (i0 + 1 <= s2 because w >= 4):

    S_w <= 6^s1 [sum over i <= i0 of 32 * 2^i + sum over i > i0 of 2^(-i) 3 * 2^s1 (3 + 2 s2)] + 2 * 4^w 2^(-s2)
        <= (6 sqrt 2)^s1 (82 + 12 s2) + 2^(3w/2 + 2)  <=  C_0 w (6 sqrt 2)^(w/2).

### (B2)

**The CRT form.** Fix p_u = (a, b) (|a|, |b| < 2^s1, b odd, a != 0) and a direction (k, k'). If (a', b') has this
direction then k' X = k Xi, that is 2^s2 (2^ds k' a' - k a) = k b' - k' b. Put e = 2^ds k' a' - k a. Then
k b' - k' b = 2^s2 e and, because |b'| < 2^s2 and |b| < 2^s1,

    |e| < |k| + 2^ds k' <= 3K,   e = -k a (mod 2^ds k'),   2^s2 e = -k' b (mod |k|).

The moduli 2^ds k' and |k| are coprime, so e lies in one class e0 modulo m = 2^ds k' |k|, and e determines
(a', b'). Hence the number of (a', b') with direction (k, k') is at most 6K/m + 1, and it is 0 unless
||e0/m|| < 3K/m. By the Chinese remainder theorem

    e0/m = -sgn(k) a/(2^ds k') - b inv(2^s1; |k|)/|k|   (mod 1),

where inv(2^s1; |k|) is the inverse of 2^s1 modulo |k|. So with theta_k = sgn(k) b inv(2^s1; |k|)/|k| (theta_k = 0
for |k| = 1), since sum over (k, k') of 1/(|k| k') <= 2 (1 + ln K)^2:

    Deg_K(p_u) <= 12 K (1 + ln K)^2 + sum over odd k, |k| <= K, of #{odd k' <= K : ||a/(2^ds k') + theta_k|| < 3K/(|k| k')}.

(Dropping the coprimality of k and k' and enlarging the tolerance only increases the count.)

**Discrepancy lemma.** Let a != 0 be an integer, ds in {0, 1}, L >= 1, theta real and dl > 0. Then

    #{odd k' in (L, 2L] : ||a/(2^ds k') + theta|| <= dl} <= 2 dl L + C_0 (|a|^(1/3) + L^(3/2) |a|^(-1/2) + 1 + ln L).

Proof. Let l_1 <= l_2 be the least and the largest integer l with 2l + 1 in (L, 2L], and n_L = l_2 - l_1 + 1 <= L the
number of points u_l = a/(2^ds (2l+1)) + theta. If dl >= 1/2 or L < |a|^(1/3) the count is at most n_L <= L, which is
at most the right side. Otherwise the set {t in [0,1) : ||t|| <= dl} is the union of [0, dl] and [1 - dl, 1), so the
count is at most 2 dl n_L + 2 n_L D, where D is the discrepancy of the points (a closed interval is a decreasing
limit of half-open ones). The Erdos-Turan inequality (Kuipers-Niederreiter, chapter 2, Theorem 2.5: for every
positive integer H, D <= 6/(H+1) + (4/pi) sum over f <= H of (1/f - 1/(H+1)) |(1/n_L) sum over l of ex(f u_l)|) gives

    n_L D <= 6 L/(H + 1) + (4/pi) sum over f <= H of (1/f) |sum over l of ex(F_f(l))|,   F_f(l) = f a/(2^ds (2l+1)).

F_f'' = 8 f a/(2^ds (2l+1)^3) has the sign of a on [l_1, l_2] and lam <= |F_f''| < 8 lam with
lam = f |a|/(2^ds L^3), hence |F_f'(l_2) - F_f'(l_1)| <= 8 lam (l_2 - l_1) <= 4 lam L. The van der Corput estimate
(the same book, chapter 1, Theorem 2.7: for integers l_1 < l_2 and F twice differentiable on [l_1, l_2] with
F'' >= rho > 0 or F'' <= -rho < 0, |sum over l_1 <= l <= l_2 of ex(F(l))| <= (|F'(l_2) - F'(l_1)| + 2)(4/sqrt(rho) + 3))
gives for lam < 1 the bound (4 lam L + 2)(4 lam^(-1/2) + 3) <= 56 L lam^(1/2) + 8 lam^(-1/2) + 6; for lam >= 1, and
when there is a single point, the trivial bound n_L <= L is smaller than this expression. With
H = ceil(L |a|^(-1/3)) (so 1 <= H <= 2 L |a|^(-1/3)), sum over f <= H of f^(-1/2) <= 2 H^(1/2) and
sum of f^(-3/2) < 2.62:

    count <= 2 dl L + C_0 [L/H + L (H |a|/L^3)^(1/2) + (L^3/|a|)^(1/2) + 1 + ln H]
          <= 2 dl L + C_0 [|a|^(1/3) + L^(3/2) |a|^(-1/2) + 1 + ln L].

**The degree bound (Deg3).** For fixed k the term k' = 1 contributes at most 1 and a block (L, 2L],
L = 1, 2, 4, ... < K, contributes at most 6K/|k| + C_0 (|a|^(1/3) + L^(3/2) |a|^(-1/2) + 1 + ln L) by the lemma with
dl = 3K/(|k| L) (the tolerance 3K/(|k| k') is smaller than dl on the block). There are at most 1 + log2 K blocks,
sum over L of L^(3/2) <= 2 K^(3/2), at most 2K values of k, and sum over odd |k| <= K of 1/|k| <= 2 (1 + ln K). So

    Deg_K(p_u) <= C_0 [K (1 + ln K)^2 + K (1 + ln K) |a|^(1/3) + K^(5/2) |a|^(-1/2)].

For Deg'_K(p_v), p_v = (a', b') with b' odd, a' != 0: put e = k a - 2^ds k' a'; then k' b = k b' + 2^s2 e,
|e| < |k| + 2^ds k', e = -2^ds k' a' (mod |k|), 2^s2 e = -k b' (mod k'), e determines (a, b), and with m' = |k| k'

    e0/m' = -2^ds a'/|k| - sgn(k) b' inv(2^s2; k')/k'   (mod 1),

so the condition is ||2^ds a'/k + theta'_(k')|| < 3K/(|k| k') with theta'_(k') = b' inv(2^s2; k')/k'. For fixed k'
and each sign of k the lemma applies to the sequence in |k| (the lemma with the integer sgn(k) 2^ds a' in place of
a and ds = 0), so Deg'_K satisfies the same bound with |a'| in place of |a|.

**The sum.** By (1), (A)3, (A)4 (at lengths s1 and s2, with |a| < 2^s1, |a'| < 2^s2) and (Deg3), for all K >= 1:

    N(K) <= C_0 6^s1 [K (1 + ln K)^2 + K (1 + ln K) 2^(s1/3) + K^(5/2) 2^(-s2/2)].

For 1 <= K <= 2^s2 put T = C_0 6^s1 [(1 + ln K)^2 + (1 + ln K) 2^(s1/3)], A_1 = C_0 6^s1 2^(-s2/2) and
B_1 = 3 * 6^s1 2^s1 (3 + 2 ln K). This bound and (Deg2) give

    N(K)/K <= min(T + A_1 K^(3/2), B_1/K) <= T + min(A_1 K^(3/2), B_1/K) <= T + A_1^(2/5) B_1^(3/5).

Since ln K <= s2 ln 2 < s1 and 3 s1 - s2 <= 2 s1 + 1, A_1^(2/5) B_1^(3/5) <= C_0 6^s1 (1 + s1)^(3/5) 2^(2 s1/5); and
2^(s1/3) <= 2^(2 s1/5). Hence N(K)/K <= C_0 (1 + s1)^2 (6 * 2^(2/5))^s1 <= C_0 8^s1, because 6 * 2^(2/5) = 7.917... < 8.
By (2),

    S_w <= sum over i < s2 of 2 C_0 8^s1 + sum over i >= s2 of 2^(-i) 4^w <= C_0 w 2^(3w/2).

### (C)

**Block law (KNOWN).** For A >= 0, m >= 1, 0 <= B < 2^m:
2 q(A 2^m + B) = 2 q(A) + 2 q(B) + A rev_m(B) + (A mod 2) B_(m-1).
Proof: unrolling the recursion, q(h) = sum over i of h_i ceil(floor(h/2^(i+1))/2). For h = A 2^m + B the bits
i >= m contribute q(A). For i < m - 1 the inner floor is A 2^(m-i-1) + floor(B/2^(i+1)) with A 2^(m-i-1) even, and
for i = m - 1 it is A (and floor(B/2^m) = 0); so the bits i < m contribute
q(B) + (A/2) sum over i < m of B_i 2^(m-1-i) + (A mod 2) B_(m-1)/2.

**Three blocks.** Applying the law twice to h = P 2^(w+s) + x 2^s + Z (0 <= x < 2^w, 0 <= Z < 2^s) gives

    eta(h) = kappa(P, Z) + c x - P rev_w(x) - 2 q(x) - (P mod 2) x_(w-1) - (c mod 2) x_0,   c = 2^(s+1) - rev_s(Z),

with kappa(P, Z) = 2 P 2^(w+s) + 2 Z + 1 - 2 q(P) - 2 q(Z) - P 2^w rev_s(Z), which does not depend on x (for s >= 1
the top bit Z_(s-1) of Z is the lowest bit of rev_s(Z), so it equals c mod 2; for s = 0, c = 2). Z -> c is a
bijection from [0, 2^s) onto the integers of (2^s, 2^(s+1)].

**The equation.** For h, h' as in the statement, t' >= 1 because both top bits are 1, and w >= 1 because the
highest differing position is not below the lowest one. By maximality of the common prefix and suffix the windows
differ in their top bits and in their bottom bits, so X and Xi are odd. Subtracting the two three-block formulas
gives c X - P Xi = gamma. Since 0 <= q(y) <= 2^(w-1) - 1 for 0 <= y < 2^w and each parity term is -1, 0 or 1,
|gamma| <= 2 (2^(w-1) - 1) + 2 = 2^w.

**Lattice line.** Fix (t', s), the odd pair and the parities of P and c (4 choices). Then gamma is constant, and two
solutions (c0, P0), (c1, P1) satisfy (c1 - c0, P1 - P0) = mu (Xi, X)/g, g = gcd(X, Xi), mu an integer. Two values of
P in [2^(t'-1), 2^t') differ by at most 2^(t'-1) - 1, two values of c in (2^s, 2^(s+1)] by at most 2^s - 1. So the
number of solutions in one parity class is at most 1 + g min(2^(t'-1)/|X|, 2^s/|Xi|)
<= 1 + 2^max(t'-1, s) g/max(|X|, |Xi|).

**Carries.** Suppose 2^(w-r) <= |X|, |Xi| < 2^w. Then |Xi| P <= |X| c + |gamma| <= (2^w - 1) 2^(s+1) + 2^w < 2^(w+s+2)
and |Xi| P >= 2^(w-r+t'-1), so t' <= s + r + 2. And |X| c <= |Xi| P + |gamma| <= (2^w - 1)(2^t' - 1) + 2^w < 2^(w+t')
while |X| c > 2^(w-r+s), so s <= t' + r - 1. Hence |s - t'| <= r + 2.

**Run-free words.** Suppose neither word has a run of more than R equal bits; we show |X|, |Xi| >= 2^(w-1-R), so
that r = R + 1. If w <= R + 1 this holds because X and Xi are odd. Let w >= R + 2 and say x_(w-1) = 1,
x'_(w-1) = 0. The R + 1 bits of x below its top bit lie in the window and are not all 0, so
x >= 2^(w-1) + 2^(w-2-R); the R + 1 bits of x' below its top bit are not all 1, so x' <= 2^(w-1) - 1 - 2^(w-2-R).
Hence X > 2^(w-1-R). The same argument applies to rev_w x and rev_w x', whose top bits are x_0 != x'_0 and whose
next R + 1 bits are x_1, ..., x_(R+1).

### (D), (E) for a bound S_w <= C_0 w 2^(beta w), beta in {3/2, beta0}

Let t = ceil(n (beta - 1)/(2 beta - 1)) (t = ceil(n/4) for beta = 3/2; t >= n/4 in both cases).

**(D).** The words of n bits with a run of at least R + 1 equal bits number at most (n - R) 2^(n-1-R) <= 2^(n-2)
(2^R >= 2n), so |W| >= 2^(n-2). Each pair counted by Coll has parameters (t', s) with t' >= t and, by (C),
|s - t'| <= R + 3; there are at most n (2R + 7) such cells. In a cell the pairs inject into the quadruples
(x, x', P, c) of (C), so by the lattice-line bound, summed over the odd pairs, their number is at most

    4 * 4^w + 4 * 2^max(t'-1, s) S_w,    w = n - t' - s.

Since w <= n - 2t' + R + 3 and n <= 4t': 4 * 4^w <= 4^(R+4) 2^n. If s >= t' then
max(t'-1, s) + beta w = beta (n - t') - (beta - 1) s <= beta n - (2 beta - 1) t' <= n. If s <= t' - 1 then
max(t'-1, s) + beta w <= t' - 1 + beta (n - 2t' + R + 3) <= n + beta (R + 3) - 1. So
2^max(t'-1,s) S_w <= C_0 n (32 n)^beta 2^n in every cell (2^R < 4n), and

    Coll <= n (2R + 7) [4^(R+4) + C_0 n (32 n)^beta] 2^n <= C_0 n^4 2^n.

**(E).** For a t-bit prefix p let W_p be the words of W with prefix p, Coll_p the number of ordered pairs h != h' in
W_p with eta(h) = eta(h'), and d_p the number of distinct values of eta on W_p. By Cauchy-Schwarz
|W_p|^2 <= d_p (|W_p| + Coll_p), and again by Cauchy-Schwarz (over the nonempty classes)
sum over p of d_p >= |W|^2/(|W| + Coll), since sum Coll_p = Coll. A value is taken in at most 2^(t-1) prefix classes, so
the number of distinct values of eta on W is at least

    |W|^2 / (2^(t-1) (|W| + Coll)) >= 2^(2n-4) / (2^(t-1) C_0 n^4 2^n) >= c_0 2^(n-t)/n^4,

and n - t >= n beta/(2 beta - 1) - 1 (= 3n/4 - 1 for beta = 3/2). These values are at most 2^(n+1) - 1. For Y >= 2^9
take n with 2^(n+1) <= Y < 2^(n+2); for 3 <= Y < 2^9 adjust c_0.

## Remarks

1. **Where the gain is.** With the trivial bound S_w <= 4^w (beta = 2) the proof gives the exponent 2/3 of
   THM-4887. The gain comes from the gcd sum: the progression term of the lattice line is large only when the
   difference vector (X, Xi) is a large multiple of a short vector, and (A)3, through (1) and the count of (Deg1),
   shows that the odd pairs parallel to a given direction number at most 6^s1 + 2 * 6^s2 <= 3.27 * 6^(w/2).
2. **Numerics, and the pairs of height 1.** S_w/6^(w/2) = 0.83, 1.02, 0.90, 1.06, 0.92, 1.06, 0.91, 1.05, 0.87, 1.01,
   0.84, 0.97 for w = 4, ..., 15 (`exponent_lemma_tests_t3_15.out`) and 0.81, 0.94, 0.78, 0.91, 0.75, 0.89 for
   w = 16, ..., 21 (`exponent34_verify/de/de_gcdsum_1_21.out`): the ratio drifts down towards the share of the pairs
   of height 1. That share is exact (PROVED from (A)2): (D, D*)(reversed eps) = (D*, D)(eps), so X = Xi iff eps is a
   palindrome and X = -Xi iff eps is an antipalindrome; the odd pairs of height 1 number 4 * 6^(w/2 - 1) for even w
   and 2 * 6^((w-1)/2) for odd w. The pre-registered W4 of this session (the same sum over pairs with one odd
   coordinate, tau = 13, 14) was confirmed. The factor w in (B2) can be removed at the cost of a larger constant.
3. **The cap.** The term 4 * 4^w in (D) allows one solution (P, c) for every odd pair in each parity class; summed
   over the cells it is at most C_0 n^2 2^n only if t >= n/4 - O(log n). An exponent above 3/4 by this route needs,
   for widths w > n/2, a bound on the number of odd pairs that admit a solution at all. This is a restricted form of
   the linear-energy problem (THM-4890 (DC), OPEN), which would give exponent 1.
4. **Balance.** The inequality |s - t'| <= r + 2 is the carry-corrected balance of THM-4890 (C) in the present
   normalisation, with a weaker constant (THM-4890 (C) gives R + 2 for run-free pairs, which is attained); it is
   re-proved here from the collision equation to keep the proof self-contained.
5. **Subsets.** By (A)2 the Kane-Tao inequality transfers: every subset B of the reversal graph has additive energy
   at most |B|^(log2 6) (CITED plus (A)2; not used above).

## Receipts

- `coordinator/reversal_energy_checks_11.out`: (A)2 for w <= 14 (all 3^w strings), (A)3 exactly for w <= 11, the
  largest number of pairs in one direction (all pairs, direction (1, 1), 1.53 * 6^(w/2) at w = 11; limit 4/sqrt 6).
- `coordinator/exponent_lemma_tests_13_22.out`: the block law; the collision equation, |gamma| <= 2^w and the
  oddness of X, Xi on all 2,026,664 unordered colliding pairs at n = 22 (and n = 8, ..., 20); |s - t'| <= R + 3 for
  run-free pairs (observed maximum R + 2, at n = 22, R = 4); the half-split identities; (A)4 (observed constant 4.4
  at w = 12, limit 8.75, bound 42).
- `coordinator/exponent_lemma_tests_t3_15.out`: S_w and the histogram of heights for w <= 15.
- `coordinator/degree_statistics_6.out`, `_8.out`: exact degrees (weighted mean about 1.4 K, maximum about 5 K).
- `coordinator/eval_w3_prefix_28.out`, `_30.out`: colliding pairs per word with a common prefix of length at least
  t (0.147, 0.019, 0.0048 at t = 4, 7, 9), far below the bound of (D).
- `exponent34_verify/` (independent programs of the verifiers): a/ ((A)1 by a lattice search to length 26, (A)2 for
  w <= 15, (A)4 for w <= 24); b1/ ((Deg1), (Deg2) exhaustively for w <= 14, sampled at w = 18, 20); b2crt/ (the CRT
  form exhaustively for w <= 11; (Deg3) with constant 1); b2analytic/ (the discrepancy lemma on random instances);
  c/ and skeptic/ (equation, lattice-line bound, carries, run-free lemma on all colliding pairs for n <= 22; the
  chain of (E) for n <= 24); de/ (cell sums, distinct values for n <= 26, S_w for w <= 21); reader_*/ (the
  proof-readers' edge-case checks).
