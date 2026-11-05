# Exact 2-adic structure of the signed one-apex responses (THM-4883)

kind-pasteur-2026-11-05-S12. **PROVED:** T1-T5 below. **FINITE-EXACT:** the confirmatory checks P1-P6 of
`PREREGISTRATION.md` (outcomes in section 10, some predictions REFUTED). **Audit:** independent and adversarial, own
code, PASS WITH CORRECTIONS (`audit/AUDIT.md`; 2 majors, 14 minors, 6 upgrades, all applied; MISTAKE-724).
**Scope:** signed responses only. Nothing here is positive Hamilton-path coverage, C2/C3 range coverage or a step
toward the all-odd spectrum by itself; the C2 statement T4 prescribes an *offset* with a variable parent, exactly as
Session 198, section 5.

This note answers lane L2 of the dispatch ("determine the exact set {a n - b g(n)} for every rational slope"). It
builds on Session 198 `session198/rational-response-coverage/PROOF.md` (the case a even, b odd) and on THM-4882
(opus-2026-11-05-S80, audited PASS WITH CORRECTIONS at 15:08 with MISTAKE-722: both-odd slopes give every even
integer, even denominators are never complete, with a mod-4 confinement on one sign side and the sharp threshold
-3(b - a)). The new content is the exact 2-adic structure for even denominators (T1-T3), the C2 offsets of every left
source (T4) and clock uniqueness (T5). The base case (T1 for odd denominators, i.e. THM-4882 (3)(i)-(ii)) is re-proved
here by a different route (section 3), so this packet does not use THM-4882 as a proof dependency.

## 1. Objects and the physical dictionary

g(1) = 1, g(2n) = g(n), g(2n - 1) = g(n) + n - (n mod 2) (n >= 2). For integers alpha, beta with beta >= 1 put

    F_(alpha,beta)(n) = alpha n - beta g(n),
    S(alpha, beta)  = { F_(alpha,beta)(n) : n >= 3 odd },
    S*(alpha, beta) = S(alpha, beta) u { alpha - beta }          (the value at n = 1).

**Physical dictionary (KNOWN).** Session 198 (eq. (7)), Session 197 `path-polynomials/ONE_CORE_BINARY_INVERSE.md`,
Session 199 `hidden-structure/PROOF.md` section 4 (the zero-J class is exactly d = 1) and THM-4882 (1): the O/Z orbit
of the two-vertex seed
t = (D,E,U,d) = (1,0,2,1) is (g(n), n mod 2, n, 1); O sends n to 2n, Z sends n to 2n - 1. Its points are the one-apex
sources: a transitive chain x_1 -> x_2 -> ... and a root c with c -> x_1, and for each later x_i either c -> x_i (O)
or x_i -> c (Z). These are exactly the rooted tournaments with transitive root deletion (d = 1) and one root-start
Hamiltonian path: d = 1 forces V - c transitive, and s = 1 forces c -> x_1. The strong ones (E = 1) are the odd
n >= 3. Hence

    S(a, b) = { aU - bD : (D, E, U, d) a strong source with s = E = d = 1 },

and every element of S(a, b) is the signed response of an explicit literal tournament (the word of n: the binary
digits of n - 1 after the leading 1, digit 1 = O, digit 0 = Z). The checks in `verify.py` build these tournaments and
count their profiles by subset dynamic programming. The audit (upgrade 6) checked the converse exhaustively at small
order: of all tournaments of orders 2-6 rooted at a fixed vertex, the 2,141 with s = d = 1 all have profile
(g(U), U mod 2, U, 1), and the 1,070 with E = 1 are exactly the strongly connected ones, with odd U >= 3.

Write b = 2^s b' with b' odd, P_s for the 2^s residue classes modulo 2^(s+1) of the integers congruent to a - b
modulo 2, and a/b = 0.d_1 d_2 d_3 ... in binary. A class C in P_s is **full** if C is contained in S(a, b), and
**one-side bounded** if S(a, b) meets C in a set bounded above or bounded below (possibly empty).

## 2. Lemmas

**Lemma 1 (cone; THM-4882 (2), Session 198 (2)).** g(n) is odd, g(n) <= n, and 2 g(n) >= n + 3 for odd n >= 3.
*Proof.* Induction on n. g(2n - 1) - g(n) = n - (n mod 2) is even, so g is odd. g(2n) = g(n) <= n, and
g(2n - 1) <= g(n) + n <= 2n, hence <= 2n - 1 because it is odd. For the lower bound let m = 2n - 1 >= 3, so n >= 2:
if n is even, g(m) = g(n) + n >= n + 1; if n is odd (n >= 3), g(m) = g(n) + n - 1 >= (n + 3)/2 + n - 1 >= n + 1.
Thus g(m) >= n + 1 = (m + 3)/2. QED

**Lemma 2 (signs and extremes; the signs are Session 198 section 4 and THM-4882 (3), the attained extrema are added
here).** Let beta >= 1.
(a) If alpha >= beta: F(n) >= (alpha - beta) n, so S*(alpha, beta) is bounded below by alpha - beta (value at
n = 1) and S(alpha, beta) by 3(alpha - beta) = F(3) (g(3) = 3); both minima are attained.
(b) If 2 alpha <= beta: for odd n >= 3, F(n) <= (alpha - beta/2) n - 3 beta/2 <= 3(alpha - beta) = F(3), and
alpha - beta >= 3(alpha - beta); so S*(alpha, beta) has maximum alpha - beta < 0 and S(alpha, beta) has maximum
3(alpha - beta) < 0, both attained.
(c) Both kinds are unbounded on their open side: F(2^k - 1) = 2^(k-1)(2 alpha - beta) - alpha - beta tends to
+infinity in (a), and F(2^k + 1) = (alpha - beta)(2^k + 1) tends to -infinity in (b) (g(2^k + 1) = 2^k + 1,
g(2^k - 1) = 2^(k-1) + 1).
*Proof.* Lemma 1. QED

**Lemma 3 (branch decomposition; THM-4882 (2) iterated with g(2^(j-1) m') = g(m')).** For all integers alpha, beta,

    S(alpha, beta) = ( S(2 alpha - beta, beta) + beta - alpha )  u  U_(j >= 1) ( S*(2^j (2 alpha - beta), beta) - alpha ).

*Proof.* Write an odd n >= 3 as n = 2m - 1, m >= 2. If m is odd (m >= 3), g(n) = g(m) + m - 1 and
F_(alpha,beta)(n) = (2 alpha - beta) m - beta g(m) - alpha + beta. If m = 2^j m' with j >= 1 and m' odd (m' >= 1),
g(n) = g(m') + 2^j m' and F_(alpha,beta)(n) = 2^j (2 alpha - beta) m' - beta g(m') - alpha. Every pair (m odd >= 3),
(j >= 1, m' odd >= 1) occurs exactly once. QED

**Lemma 4 (homogeneity).** S(c alpha, c beta) = c S(alpha, beta) and likewise for S*.

**Lemma 5 (tail lemma; KNOWN in substance).** This is the block transfer of Session 198
`source-digit-blocks/DIGIT_BLOCK.md` (1) and Session 199 `hidden-structure/PROOF.md` Theorem 1, rewritten as a change
of slope. For an O/Z word T of length l with chronological letters eps_0, ..., eps_(l-1)
(eps = 1 for Z) put A_T = sum_i eps_i 2^i. For every strong head m (odd m >= 3),

    F_(alpha,beta)(m . T) = F_(2^l alpha - beta A_T, beta)(m) + kappa_T,

with an integer kappa_T depending only on alpha, beta, T. Here m . T is the cover count after appending T
(O: m -> 2m, Z: m -> 2m - 1).
*Proof.* One letter: F_(alpha,beta)(2m) = F_(2 alpha,beta)(m) since g(2m) = g(m); and
F_(alpha,beta)(2m - 1) = F_(2 alpha - beta,beta)(m) - alpha + beta (m mod 2) by the recurrence. Peel the letters of
T from the last one; the parity m mod 2 at each Z is the previous letter, or 1 for the first letter because m is
odd. The slope numerator becomes 2^l alpha - beta sum_i eps_i 2^i. QED

**Lemma 6 (clock; Session 198 (3)-(5), cyclic form).** Let alpha/beta in (1/2, 1) be reduced with beta odd,
k = ord_beta(2), B = 2^k, A = alpha (B - 1)/beta and W the word with eps_i = bit i of A. Then W ends with Z, is
mixed, and for every strong n, F(n . W) = F(n) - K with a constant K > 0, K even, and n . W strong.
*Proof.* B/2 <= A < B - 1 since 1/2 < alpha/beta < 1, so the last bit is 1 and W is mixed. Let
R = sum_i eps_i 2^(k-1-i), u* = R/(B - 1) and the periodic orbit u_0 = u*, u_(i+1) = 2 u_i - eps_i (so u_k = u_0 and
0 < u_i < 1, eps_i = floor(2 u_i)). With U_i the cover count before letter i, U_i = 2^i (U - u*) + u_i. On a strong
input E_0 = 1 = eps_(k-1), and E_i = eps_(i-1); so with cyclic indices D' = D + sum_i eps_i (U_i - eps_(i-1)) and
U' = BU - R. Substituting and using beta A = alpha (B - 1):

    F' = alpha U' - beta D' = F - beta Delta,     Delta = sum_i eps_i (u_i - eps_(i-1))   (cyclic).

Split the cyclic word into maximal runs of ones; a run of j ones followed by a 0 with orbit value v there has deficits
1 - u summing to (1 - v)(1 - 2^(-j)) (backwards, 1 - u_i = (1 - u_(i+1))/2 on ones) and adjacent-one corrections
j - 1, so it contributes v + (1 - v) 2^(-j) > 0 (Session 198's run argument, here with the cyclic term, which is
present because the input is strong). Hence K = beta Delta > 0; K is an integer (a difference of responses) and even
(both responses have the parity alpha - beta). The output ends with Z, so it is strong. QED

For a even the first letter is O and K is Session 198's K; for a odd the first letter is Z and K is THM-4882's
K_cyc = K_lin - b.

## 3. Base case: odd denominators (THM-4882 (3)(i)-(ii), re-proved)

**Theorem A0.** If alpha/beta in (1/2, 1) is reduced and beta is odd, then S(alpha, beta) is the set of all integers
congruent to alpha - beta modulo 2.
*Proof.* Inclusion: n and g(n) are odd. Closure: by Lemma 6, S(alpha, beta) - K is contained in S(alpha, beta).
Residues: fix j >= 1 with sigma = 2^j (2 alpha - beta) > beta. By Lemma 3, S(alpha, beta) contains
sigma L - beta g(L) - alpha for every odd L >= 1 (n = 2^(j+1) L - 1). By the reviewed JOINT_MODULAR_GADGET
theorem (Session 197, `session197/modular-gadget-universality/PROOF.md`, Theorem 1) with the even modulus K, for all
odd residues u, z modulo K there are arbitrarily large L with L = u and g(L) = z modulo K. Writing u = 1 + 2y,
z = 1 + 2x, the value is congruent to sigma - beta - alpha + 2(sigma y - beta x); since gcd(sigma, beta) =
gcd(2 alpha, beta) = 1, every residue modulo K of the parity alpha - beta occurs, with values >=
(sigma - beta) L - alpha arbitrarily large. Given r of that parity choose such a value v >= r with v = r (mod K) and
its strong n_0; then n_0 . W^((v - r)/K) has response r. QED

This is a route different from Session 198's (which uses the STRONG_SOURCE_INTERFACE seed followed by one O) and
from THM-4882's (the same seed with K_lin for the first period); only the JOINT_MODULAR_GADGET theorem is used.

## 4. T1: the 2-adic structure theorem

**Theorem T1.** Let a/b be a reduced rational slope, b = 2^s b' with b' odd.
1. If 1/2 < a/b < 1 and b' > 1, exactly one class C* in P_s is full, and every other class of P_s is one-side
   bounded.
2. In every other case (a/b <= 1/2, a/b >= 1, or b a power of two) no class is full and every class is one-side
   bounded.

Consequently S(a, b) is the set of all integers of its parity iff 1/2 < a/b < 1 and b is odd (Session 198 and
THM-4882); for an even b that is not a power of two it is one full class modulo 2^(s+1) plus one-side-bounded
remainders; for b a power of two, and for exterior slopes, it is one-side bounded in every class.

**Corollary (two-sided progressions; audit upgrade 4).** Every two-sided infinite arithmetic progression
{c + k e : k in Z} contained in S(a, b) lies in C*; in particular none exists when b is a power of two or the slope is
exterior. (Its intersection with each class of P_s is empty or again a two-sided progression, which a one-side
bounded class cannot contain.)

**Dyadic example (audit upgrade 5).** For 3/4 (s = 2) the classes modulo 8 are: 1 and 5 bounded below by 1 and -3,
3 and 7 bounded above by -5 (attained in class 3; class 7 has maximum -9). This refines THM-4882's two-sided
confinement modulo 4 at (3, 4).

*Proof.* Case 2, exterior slopes: Lemma 2 makes S one-signed, and a full class is unbounded on both sides. For
1/2 < a/b < 1 we prove 1 and the dyadic part of 2 together, by induction on s. The case s = 0 is Theorem A0 (P_0 has
one class). Let s >= 1; then a is odd. Put a_1 = a - b/2, b_1 = b/2, x_1 = a_1/b_1 = 2(a/b) - 1 in (0, 1); gcd(a_1,
b_1) = 1. By Lemmas 3 and 4,

    S(a, b) = Q_0 u U_(j>=1) Q_j,   Q_0 = 2 S(a_1, b_1) + b - a,   Q_j = 2 S*(2^j a_1, b_1) - a.

*Valuation routing.* Every element r of a piece has r + a even, and v_2(r + a) is determined by the piece:

| s      | piece                | reduced pair inside            | parity of its values | v_2(r + a) |
|--------|----------------------|--------------------------------|----------------------|------------|
| >= 2   | Q_0                  | (a_1, 2^(s-1) b'), odd/even    | odd                  | = 1        |
| >= 2   | Q_j, 1 <= j <= s - 2 | (a_1, 2^(s-1-j) b'), odd/even  | odd (scaled 2^j)     | = j + 1    |
| >= 2   | Q_(s-1)              | (a_1, b'), odd/odd             | even (scaled 2^(s-1))| >= s + 1   |
| >= 2   | Q_j, j >= s          | (2^(j-s+1) a_1, b'), even/odd  | odd (scaled 2^(s-1)) | = s        |
| = 1    | Q_0                  | (a_1, b'), even/odd            | odd                  | >= 2       |
| = 1    | Q_j, j >= 1          | (2^j a_1, b'), even/odd        | odd                  | = 1        |

(For Q_0 with s >= 2, r + a = 2(odd) + b with 4 | b; with s = 1, r + a = 2(odd) + 2b'.) The sets
G_nu = {v_2(r + a) = nu} (1 <= nu <= s) and G_(>s) = {v_2(r + a) >= s + 1} partition the odd integers; G_nu is a
union of 2^(s - nu) classes of P_s and G_(>s) is one class. So each class of P_s receives pieces from one group only:
G_1 gets Q_0 (s >= 2) or the family Q_(j>=1) (s = 1), G_(j+1) gets Q_j alone (1 <= j <= s - 2), G_s gets the family
Q_(j>=s), G_(>s) gets Q_(s-1) (s >= 2) or Q_0 (s = 1).

*Signs.* The piece Q_j has slope 2^j x_1 (j >= 0). Let j* be the unique j >= 0 with 2^j x_1 in (1/2, 1), if it
exists. Pieces with j < j* have slope < 1/2 (Lemma 2(b): bounded above), pieces with j > j* slope > 1 (Lemma 2(a):
bounded below). If no j* exists, some 2^(j_0) x_1 equals 1/2; this happens iff b' = 1 and a_1 = 1 (x_1 = 2^(1-s)),
and then j_0 = s - 2, the pieces j <= j_0 are bounded above, j > j_0 bounded below.

*Groups without the interior piece.* Each receives a single piece or the family Q_(j>=s). The family has one sign:
it misses j* in the cases below, and j >= s > j* (or j >= s > j_0 = s - 2), so it is bounded below. Every class in
such a group is one-side bounded.

*The interior group.* (i) j* >= s (s >= 2), or j* >= 1 (s = 1): Q_(j*) contains 2^s times the image of the reduced
pair (2^(j*-s+1) a_1, b') (resp. (2^(j*) a_1, b')), an even/odd interior slope with odd denominator b' > 1; by
Theorem A0 that image is all odd integers, so Q_(j*) contains 2^s (2Z + 1) - a = G_s (resp. 2(2Z + 1) - a = G_1):
this group is one full class. (ii) j* = s - 1 (s >= 2): the reduced pair (a_1, b') is odd/odd with b' > 1, its image
is all even integers, and Q_(s-1) contains 2^(s+1) Z - a = G_(>s), one full class. (iii) s = 1 and j* = 0:
Q_0 = 2 S(a_1, b') + b - a with an even/odd interior pair, so Q_0 = 2(2Z + 1) + 2b' - a = G_(>1), full.
(iv) s >= 2 and 0 <= j* <= s - 2: the interior piece is alone in G_(j*+1) and equals phi(S^(*)(a_1, 2^(s'') b'))
with s'' = s - 1 - j* in [1, s - 1] and phi(r) = 2^(j*+1) r + tau, tau = b - a (j* = 0) or -a (j* >= 1). phi is
increasing and maps the 2^(s'') classes of P_(s'') (odd integers modulo 2^(s''+1)) bijectively onto the 2^(s'')
classes of P_s inside G_(j*+1) (phi(r) + a = 2^(j*+1) r + b or 2^(j*+1) r). By induction the inner image has exactly
one full class iff b' > 1 and all other classes one-side bounded; phi transports both properties, and the single
extra value of S* changes neither. 

Every class not handled as full is one-side bounded, and a one-side bounded class is not full. In (i)-(iii) b' > 1
is automatic (an interior slope with odd denominator), and the no-j* case needs b' = 1; in (iv) the recursion
preserves b'. Hence exactly one full class iff b' > 1, none iff b' = 1. QED

**The classification algorithm** (called "the automaton" in `PREREGISTRATION.md`; it is the proof's recursion, not a
finite automaton, and has nothing to do with THM-4878's cut-state automaton). The proof is an algorithm
(`gen_pieces`/`classify` in `verify.py`): expand S(a, b) by Lemma 3 along the interior piece, reduce each pair, stop at
odd interior denominators (FULL), record every other piece with its sign and its attained extremum (Lemma 2), and cut
the family Q_(j>=s) once its pieces are positive (they share one class and one sign). Each class of P_s gets FULL, or
the sign and the extremum over the pieces of its group. That group extremum is attained in the group (Lemma 2), but
**not necessarily in each class of the group**: a piece of scale 2^nu with nu < s spreads over 2^(s - nu) classes and
its extreme value lies in one of them. The pre-registered P1 clause "attained in every class" is therefore false in
general (section 10).

## 5. T2: the explicit full class

**Theorem T2.** Let 1/2 < a/b < 1 with b' > 1, and let L be the least position L >= s + 1 with d_L = 1 (it exists
because the binary expansion of a/b is not finite). Let T be the chronological word d_(L-1), d_(L-2), ..., d_1
(1 = Z, 0 = O). Then for every strong head m, a (m . T) - b g(m . T) lies in C*, and

    { a n - b g(n) : n = m . T, m odd >= 3 } = C*.

Equivalently: C* is the class of a n - b g(n) for any odd n >= 3 whose O/Z word ends with the L letters
d_L, d_(L-1), ..., d_1 (d_L = 1, a Z), i.e. n - 1 = sum_(i<L) (1 - d_(i+1)) 2^i modulo 2^L and n > 2^L.
*Proof.* Lemma 5 with l = L - 1: A_T = sum_(l'=1)^(L-1) d_(l') 2^(L-1-l') = floor(2^(L-1) a/b), so
2^(L-1) a - b A_T = b frac(2^(L-1) a/b) = 2^s a'' where a''/b' = T^(L-1)(a/b) = 0.d_L d_(L+1)... is the doubling-map
iterate, reduced with denominator exactly b' (L - 1 >= s) and in (1/2, 1) (d_L = 1, and 1/2 has even denominator).
So F_(a,b)(m . T) = 2^s F_(a'',b')(m) + kappa_T, and Theorem A0 makes {F_(a'',b')(m)} the full parity class; the image
is 2^s (2Z + p) + kappa_T, one class of P_s, full; by T1 it is C*. The congruence form reads the word of n off the
binary digits of n - 1. QED

Examples (also tabulated in the receipt for b <= 24): 5/6 = 0.1101..., s = 1, L = 2, the word ends with ZZ, C* = 3
(mod 4); 7/12 = 0.10010101..., s = 2, L = 4, the word ends with ZOOZ, C* = 5 (mod 8); 11/12 = 0.1110..., L = 3, ZZZ,
C* = 7 (mod 8).

The full class is realized constructively: tail T after a strong head with the right residue modulo the head slope's
clock drop K'', followed by clock repetitions. Theorem A0 proves such heads exist with arbitrarily large values (its
branch heads 2^(j+1) L - 1 with gadget-built L); `verify.py` finds them in practice by direct search over odd heads
below 2^20, and builds witnesses this way for every target of C* in [-2000, 2000] that a direct search below 2^24
misses (P2). A0's own branch heads are explicit but can be very long when K'' is large (the gadget word has O(K''^2)
bits): the audit's implementation stayed under a 40,000-letter cap for 923 of 974 targets and 467 of 506 T2 heads, the
rest (K = 50, 62, 74, 98) being witnessed by direct heads.

## 6. T3: exact half-lines

**Theorem T3.** Let 1/2 < a/b < 1, b' > 1, s >= 1.
1. S(a, b) and C* agree on some half-line [c, infinity) iff d_2 = ... = d_(s+1) = 0 (i.e. a/b < 1/2 + 2^(-s-1)).
   Then S(a, b) n (c_+, infinity) = C* n (c_+, infinity) with
   c_+ = 5(a - b) for s = 1 and c_+ = F(2^s - 1) = 2^(s-1)(2a - b) - a - b for s >= 2,
   and c_+ is an element of S(a, b) outside C*. In particular S(a, b) n (0, infinity) = C* n (0, infinity).
   (PREREGISTRATION.md wrote max(5(a - b), F(2^s - 1)) for s >= 2; the maximum is always F(2^s - 1), since
   F(2^(j+1) - 1) = 2^j (2a - b) - a - b increases with j and already F(3) = 3(a - b) > 5(a - b); audit minor 2.)
2. S(a, b) and C* agree on some half-line (-infinity, c] iff d_1 = ... = d_(s+1) = 1 (a/b > 1 - 2^(-s-1)). Then
   S(a, b) n (-infinity, c_-) = C* n (-infinity, c_-) with c_- = (2^s + 1)(a - b) = F(2^s + 1), an element of
   S(a, b) outside C*.
3. Otherwise S(a, b) minus C* is unbounded above and below.

*Proof.* With the notation of section 4, digit d_(i+2) = 1 iff the interior index at the i-th level of the
recursion is 0. (1) If d_2 = ... = d_(s+1) = 0 then 2^s x_1 < 1, so j* >= s: case (i), C* = G_s, and the other
pieces Q_0, Q_1, ..., Q_(s-1) all have j < j*; their maxima are 5(a - b) (Q_0 at n = 5, Lemma 2(b) for S) and
F(2^(j+1) - 1) (Q_j at m' = 1), increasing in j. Conversely if some d_(j+2) = 1 with j <= s - 1, then j* <= s - 1
and the family Q_(j>=s) is a non-full piece bounded below and unbounded above (Lemma 2(c)). (2) If d_2 = ... =
d_(s+1) = 1 the recursion takes j* = 0 at the levels i = 0, ..., s - 1 (the slope at level i is
x_i = 2^i (a/b) - (2^i - 1)) and every non-full piece is a family Q_(j>=1) of positive slope at some level. Mapped to
the top level by r -> 2^i r + (b - a)(2^i - 1), the minimum of the level-i family is at j = 1, m' = 1, with value
2^i 3 (b/2^i)(x_i - 1) + (b - a)(2^i - 1) = (2^(i+1) + 1)(a - b), attained at n = 2^(i+1) + 1; the least is at
i = s - 1. Conversely, if some level has j* >= 1, its Q_0 is a non-full piece of slope < 1/2, unbounded below.
(3) Neither condition: both kinds of non-full pieces exist, each unbounded on its open side (Lemma 2(c)). QED

For b = 2 (mod 4) (s = 1) exactly one of 1, 2 holds: a/b < 3/4 gives S n (5(a - b), infinity) = {r : r = 2 - a
(mod 4)} on that half-line, and a/b > 3/4 gives S n (-infinity, 3(a - b)) = {r : r = -a (mod 4)} there. These are
the equality versions of THM-4882 (3)(iii) (its classes -(a + b) and -a modulo 4). The threshold 3(a - b) =
-3(b - a), attained at n = 3, was adopted into THM-4882 from its independent audit at 15:08 the same day
(concurrent; priority there); what is new here is the equality with a full class, the threshold 5(a - b) on the
positive side, and the general-s thresholds c_+ and c_-. Example: S(5, 6) n (-infinity, 0) = {r < 0 : r = 3
(mod 4)} u {-3}; S(7, 12) n (0, infinity) = {r > 0 : r = 5 (mod 8)}.

## 7. T4: C2 offsets of every left source

Let v = (D_v, E_v, U_v, d_v) be a nontrivial physical source and q = (D, 1, U, 1) a strong unit source. The
reviewed complete C2 theorem (`dependencies` of Session 198; BILINEAR_C2 (11)-(13)) realizes the parent
n = N(v, Oq) and the child m = M(v, Jq) as actual strong tournaments with exactly two Hamiltonian cycles
(Oq = (D, 0, 2U, 1) and Jq = (D - 1 + U, D - 1, D + U, D) are actual nontrivial sources).

**Theorem T4.** m - 2n = alpha_v D - beta_v U + c_v with

    alpha_v = D_v + 4 U_v - 2 d_v,  beta_v = alpha_v + E_v - d_v,  c_v = -D_v + U_v + d_v + 2 E_v.

alpha_v is odd, and E_v <= d_v for every source. If E_v < d_v, put gamma = gcd(alpha_v, d_v - E_v): over all strong
unit sources q the offsets m - 2n form exactly the residue class
U_v - D_v + E_v + 2 d_v (mod 2 gamma) (all odd offsets when gamma = 1), each offset with infinitely many parents.
If E_v = d_v the offsets are bounded above by c_v, and c_v is attained.
*Proof.* Substituting w = Oq into N and w = Jq into M (BILINEAR_C2 (11)-(12)) and subtracting gives the displayed
coefficients (direct expansion, also checked on 78,000 pairs, and symbolically by the audit). D_v is odd, so alpha_v is
odd. E_v <= d_v: deleting the root from a Hamiltonian path that ends at the root leaves a Hamiltonian path of V - c, an
injection (POSITIVE_SOURCE_GRAMMAR (1)). If E_v < d_v, the cone d_v <= D_v <= U_v (BILINEAR_C2 section 1) gives
alpha_v >= 3 D_v >= 3 d_v > 2(d_v - E_v), so beta_v/alpha_v in (1/2, 1); reduced by gamma it has odd denominator, and
Theorem A0 makes {(beta_v/gamma) U - (alpha_v/gamma) D} the set of all integers of parity (E_v - d_v)/gamma; the
offsets are c_v minus gamma times it, the class c_v + d_v - E_v = U_v - D_v + E_v + 2 d_v modulo 2 gamma. Unbounded D
(Theorem A0's witnesses grow) gives infinitely many parents. If E_v = d_v, beta_v = alpha_v and
m - 2n = c_v - alpha_v (U - D) <= c_v, with equality for the strong unit sources with U = D (cover n = 2^k + 1). QED

The class is odd, as it must be (m and n are odd): U - D + E is odd for every physical source (1 at t; O makes it
2U - D, Z leaves it unchanged, J makes it D; D is odd) (audit upgrade 1). For v = t this is Session 198's (9):
alpha = 7, beta = 6, c = 2. Every left source with d_v - E_v = 1 yields all odd offsets; examples are t, O t, every
O-born source over a source with d = 1, and every J-born source over a source with E = 1 (for instance over a strong
unit source), since J gives d' - E' = E. (A J-born source over an arbitrary strong source need not qualify:
t . JJZ = (9, 3, 9, 3) has the J-child (15, 6, 18, 9) with gamma = 3 and offsets 3 (mod 6) only; audit major 1.)
The parents are not prescribed: a fixed parent count with both children 2n - 1 and 2n + 1 is not obtained (dispatch
lane L3 stays open).

## 8. T5: clock uniqueness, and the witness-size question

**Theorem T5.** Fix a reduced slope a/b and an O/Z word w of length k with value A_w = sum_i eps_i 2^i. For every
strong unit source with cover count U, appending w changes aU - bD by lambda_w U + const_w with
lambda_w = a(2^k - 1) - b A_w. The change is independent of the source iff b A_w = a(2^k - 1). This forces
0 <= a/b <= 1 (since 0 <= A_w <= 2^k - 1), and for 0 < a/b < 1 it holds iff b is odd, ord_b(2) divides k and w is the
(k/ord_b(2))-th power of the ord_b(2)-letter word of a(2^(ord_b(2)) - 1)/b. For **1/2 < a/b < 1** that word is
Session 198's clock W, so the clock powers are the only source-independent moves; for even b there are none at any
slope. *Proof.* Lemma 5 with the head equal to the source: F_(a,b)(m . w) = F_(a,b)(m) + (a(2^k - 1) - b A_w) m +
kappa_w, and m ranges over infinitely many values. b A_w = a(2^k - 1) needs b | a(2^k - 1), i.e. b odd and
ord_b(2) | k, and then A_w = A (2^k - 1)/(2^(ord) - 1) with A = a(2^(ord) - 1)/b, whose bits repeat those of A. QED

Outside (1/2, 1) the source-independent words are not Session 198 clocks: 1/3 has ZO (change +1), 2/5 has OZZO, 1/5
ZZOO, 3/7 ZZO, and slopes above 1 (4/3, 6/5) have none (audit major 2; the first filing stated the interior "iff" for
every slope). Within the unit-deletion class, then, an interior signed functional has a single source-independent move
and a single drop K. This does not settle HYP-9475 (C): words of other slopes still change aU - bD by an amount
lambda_w U + const that depends on the source but can be used inside a target-dependent search. Whether the shortest
witness has O(log |r|) letters is HYP-9476 (C), with data in section 10.

## 9. Physical witnesses and what is not claimed

Every element of S(a, b) is the signed response of an explicit one-apex tournament (section 1). For b' > 1 every
element of C* has the recipe of T2 (tail word after a head with the right residue, then clock periods); for b' = 1
and for the non-full classes the elements are exactly the values of the listed one-signed pieces, whose exact
membership is a positive-coverage question for the one-apex sequence (not decided here; HYP-9476 (B) records their
densities). T4's pairs are actual
C2 graphs. No witness-size bound is claimed (the clock repetitions can be long), and no positive Hamilton-count,
C2/C3 range or all-odd spectrum consequence is drawn.

## 10. Verification (FINITE-EXACT) and pre-registered outcomes

See `README.md` for the run commands and `VALIDATION_NORMAL.json` / `VALIDATION_OPTIMIZED.json` (byte-identical) for
the receipts. Pre-registration: commit 2c11d95d56, before the confirmatory run.

- **P1 (structure census): structure CONFIRMED, per-class attainment REFUTED.** All 629 reduced slopes in (1/2, 1)
  with 3 <= b <= 64: one full class exactly when b is not a power of two (629/629); no class with mixed signs; zero
  bound violations over all odd n < 2^22; the digit description of T2 agrees with the automaton (629/629), and every
  odd n < 2^22 with the T2 suffix lands in C*. The clause "every bound equals the extremum of its class" is false:
  1,143 of the 1,902 non-full classes do not attain their group bound (first case: 3/4, class 7 mod 8, bound -5,
  maximum -9). Every one of the 759 groups attains its bound. This is the reason given at the end of section 4.
- **P2 (fullness): CONFIRMED.** 340 slopes with b <= 48, b' > 1: all 544,115 targets of C* with |r| <= 2000 have
  verified witnesses (320,748 by direct search below 2^24, 223,367 by the T2 recipe with at most 40 clock periods
  and at most 1,456 binary digits), each checked by the transfer and by the closed bit formula; 2,033 distinct
  witnesses below 2^12 counted literally (profile (g(n), 1, n, 1), one root-start path).
- **P3 (exact half-lines): CONFIRMED.** 186 slopes (69 positive-exact, 73 negative-exact, 44 neither), all as T3
  predicts, thresholds attained.
- **P4 (C2 offsets): CONFIRMED.** The offset identity on 78,000 pairs (1,950 bank sources x 40 strong unit sources);
  the slope/parity statement on all 1,950 sources; six literal parent/child pairs from five left profiles
  (1,0,4,1), (5,0,10,1), (11,0,24,5), (35,0,70,11), (23,0,48,11), offsets -1 and +1, graph orders 11-19, every graph
  with exactly two Hamiltonian cycles and Hamilton counts m = 2n + delta by independent subset DP.
- **P5 (densities): CONFIRMED.** 50 decided non-full classes of the 14 slopes have hit fractions between 0.097 and
  0.617 on 500 <= |r| <= 20000. The 9 further classes, which contain a slope-1 or slope-1/2 piece and were reported
  as lower bounds (0.526-0.546), are decided as well (audit upgrade 2): every value v of the slope-1 image {n - g(n)}
  is attained at some odd n <= 4|v| + 3 (by Lemma 3 with alpha = beta, an odd m gives F(2m - 1) = F(m) at a smaller
  argument, and m = 2^j m' gives F >= (2^j - 1) m' - 1), and every value v of the slope-1/2 image {n - 2 g(n)} at
  some n <= 4|v| + 3 (its values are 1 - 2g(m), m odd, and -1 - 2g(m'), the same for every j); after the scaling of
  each piece this is far below 2^22, and the audit found all 59 classes unchanged between n < 2^22 and n < 2^26. So
  P5 covers 59 classes, range 0.097-0.617.
- **P6 (witness size): partly REFUTED.** 3/5: all three clauses hold (every |r| <= 1024 found by length 34; maxlen
  9, 17, 22, 22, 32, 32 at j = 4..9; 184 of 522 clock-ending witnesses have t >= 2). 8/9: (a) 99.3% found and (c)
  160 of 433 hold, but (b) fails: maxlen(2^4) = 21 > 20 (target -15, decided). At j = 9 the depth-34 search left four
  targets undecided against the bound 3j + 8 = 35 (`verify.py` counted them as failures); the audit's depth-35 search
  shows -427, -425, -405 and 277 need at least 36 letters, so (b) fails at j = 9 too. The linear-in-j growth is not
  refuted; the additive constant 8 was too small for a six-letter clock. Clause (c) uses the first minimal witness in
  depth-first order (O before Z); that tie-break was fixed in `witness_length.c` after the pre-registration (margins
  184/522 and 160/433 against 1/2, so it cannot change the outcome).
- **T5 check:** all 1,253,070 words of length <= 12 over the 153 interior slopes with b <= 31 behave as T5 states
  (the audit: all 5,013,198 words of length <= 14).

## 11. Dependencies

- Session 198 `session198/rational-response-coverage/PROOF.md` (eq. (1), (7), the clock (3)-(5)); KNOWN, audited.
- Session 197 `session197/modular-gadget-universality/PROOF.md`, Theorem 1 (JOINT_MODULAR_GADGET); KNOWN, audited.
- BILINEAR_C2 (Session 197, frozen in Session 198's packet) for T4; `dependencies/bilinear.py`, `dependencies/core.py`
  are byte-identical copies (SHA256 in the receipt).
- THM-4882 (opus-2026-11-05-S80) for priority on the odd-denominator classification, the branch formulas (Lemma 3),
  the mod-4 confinement and the threshold -3(b - a); not used as a proof dependency (Theorem A0 re-proves the part
  that is used).
- POSITIVE_SOURCE_GRAMMAR (1) (Session 198) for E <= d in T4.

## 12. Audit

`audit/AUDIT.md` (independent agent, own code `audit/independent.py`, `audit/audit_census.c`, output
`audit/INDEPENDENT.json`): PASS WITH CORRECTIONS. Majors: (1) a false example in T4 (J-born sources over strong
sources need not have d - E = 1); (2) T5's "iff" was stated for every slope but holds as stated only inside (1/2, 1).
Minors: HYP-9476 sub-labels, the redundant maximum in c_+, the THM summary wording ("exactly one branch is interior"
is false in the dyadic case a = b/2 + 1), the unstated E <= d, two missing credits (Lemmas 2, 3), a sentence false for
b a power of two, a stale THM-4882 status, an over-reaching HYP sentence, the reason for P6 at j = 9, the T5 remark on
HYP-9475 (C), the description of the witness heads, the name "automaton", the post-registration P6 tie-break.
Upgrades adopted: the explicit T4 class and the parity lemma for U - D + E; the decided P5 classes; exterior density
evidence (HYP-9476 (B)); the corollary on two-sided progressions; the dyadic 3/4 example; the exhaustive dictionary
check. All applied; the receipts are unchanged (the corrections are textual). MISTAKE-724 records the errors.

**Concurrent work read at the merge (16:40).** Session 200 found the J-born example error (major 1) independently
from the pre-audit snapshot, with the same repair and the same counterexample profile (t . ZOJ = (9, 3, 9, 3), its
J-child (15, 6, 18, 9), offsets 3 (mod 6)): `session200/thm4883-example-review/REVIEW.md` (16:09). Session 200's
`session200/signed-block-transport/PROOF.md` (15:56, pending its audit) uses the slope-transfer identity of Lemma 5
with the fixed block O^l Z^j to move an odd-numerator slope with odd denominator to an even-numerator one, a third
proof of THM-4882 (3)(ii) with O(log b) extra vertices; T2 here (first pushed 15:24) is the even-denominator use of
the same identity, transferring to the odd-denominator descendant with the factor 2^s.
