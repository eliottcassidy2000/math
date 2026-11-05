---
id: THM-4883
title: "Signed one-apex responses S(a,b) = {a n - b g(n): n >= 3 odd}, the signed responses aU - bD of the strong unit-deletion sources: exact 2-adic structure. With b = 2^s b' (b' odd): (T1) S(a,b) contains a full residue class modulo 2^(s+1) iff 1/2 < a/b < 1 and b is not a power of two, then exactly one (C*), and every other class is met in a set bounded on one side; (T2) C* is read off the binary digits of a/b (the class of every n whose O/Z word ends with d_L ... d_1, L the first position > s with d_L = 1); (T3) S(a,b) equals C* on a positive half-line iff d_2..d_(s+1) are all 0, on a negative half-line iff d_1..d_(s+1) are all 1, with explicit attained thresholds, otherwise on neither; (T4) the C2 offset M(v,Jq) - 2N(v,Oq) of every nontrivial left source v is a signed response with odd denominator, filling the class U_v - D_v + E_v + 2d_v modulo 2 gcd(D_v + 4U_v - 2d_v, d_v - E_v) when E_v < d_v; (T5) a word moves aU - bD by a source-independent amount iff b A_w = a(2^k - 1): for 1/2 < a/b < 1 exactly the clock powers, none for even b"
status: >
  PROVED: T1-T5 (PROOF.md of the packet; the odd-denominator base case, THM-4882 (3)(i)-(ii) and Session 198, is
  re-proved there through the clock and the JOINT_MODULAR_GADGET theorem). KNOWN / CITED: the one-apex dictionary and
  the block transfer (Session 197 ONE_CORE_BINARY_INVERSE, Session 198 rational-response-coverage (7) and
  source-digit-blocks DIGIT_BLOCK (1), Session 199 hidden-structure, THM-4882 (1)); the branch formulas (THM-4882 (2));
  the sign statements outside (1/2, 1) (Session 198 section 4); the mod-4 confinement and the threshold -3(b - a) on
  one side for even b (THM-4882 (3)(iii), concurrent). FINITE-EXACT: pre-registered P1-P6 (P2-P5 confirmed; P1
  structure confirmed, its per-class attainment clause REFUTED; P6 partly REFUTED), HYP-9476. NOT CLAIMED: positive
  Hamilton-count coverage, C2/C3 range coverage, fixed-parent children, witness-size bounds; Chan-Pak remains OPEN.
  AUDIT: independent, adversarial, own code, PASS WITH CORRECTIONS (packet audit/AUDIT.md; 2 majors, 14 minors,
  6 upgrades, all applied; MISTAKE-724). Major 1 (the J-born example in T4) was also found independently by
  Session 200 (session200/thm4883-example-review/REVIEW.md). CONCURRENT: Session 200
  session200/signed-block-transport uses the same slope-transfer identity for odd denominators.
session: kind-pasteur-2026-11-05-S12
source: 05-knowledge/results/signed_response_2adic_structure_20261105/PROOF.md
related:
  - THM-4882
  - HYP-9475
  - HYP-9476
  - THM-4878
  - session198/rational-response-coverage/PROOF.md
  - session198/source-digit-blocks/DIGIT_BLOCK.md
  - session197/modular-gadget-universality/PROOF.md
  - session199/hidden-structure/PROOF.md
---

# THM-4883 - exact 2-adic structure of the signed one-apex responses

## Setting

g(1) = 1, g(2n) = g(n), g(2n - 1) = g(n) + n - (n mod 2). For coprime a, b (b >= 1) put
S(a, b) = { a n - b g(n) : n >= 3 odd }. The strong rooted tournaments with one root-start Hamiltonian path and
transitive root deletion (d = 1) are exactly the one-apex sources with odd cover count n >= 3, with profile
(D, E, U, d) = (g(n), 1, n, 1) (KNOWN; checked exhaustively at orders 2-6 by the audit), so S(a, b) is the image of
the signed response aU - bD on them and every element has a literal tournament. Write b = 2^s b' (b' odd),
a/b = 0.d_1 d_2 ... in binary, and P_s for the 2^s classes modulo 2^(s+1) of the integers of the parity of a - b.

## Statements (PROVED)

**T1 (2-adic structure).** S(a, b) contains a full class of P_s iff 1/2 < a/b < 1 and b' > 1; it then contains
exactly one, C*, and S(a, b) meets every other class of P_s in a set bounded above or below. In all other cases
(a/b outside (1/2, 1), or b a power of two) no class is full and every class is one-side bounded. For odd b this
is the statement that S(a, b) is all integers of its parity (Session 198 for a even, THM-4882 (3)(ii) for a odd).
Corollary: every two-sided infinite arithmetic progression contained in S(a, b) lies in C* (none for b a power of two
or exterior slopes). Proof: induction on s via the branch decomposition
S(a, b) = (S(2a - b, b) + b - a) u U_(j>=1) (S*(2^j(2a - b), b) - a); the 2-adic valuation v_2(r + a) routes each
branch to its own group of classes, at most one branch has an interior slope (none iff b' = 1 and a = b/2 + 1), the
others are one-signed by (n + 3)/2 <= g(n) <= n (odd n >= 3), and the recursion ends at an odd interior denominator
(full class) iff b' > 1. Example: 3/4 modulo 8 has classes 1, 5 bounded below and 3, 7 bounded above.

**T2 (explicit C*).** Let L >= s + 1 be least with d_L = 1. For every strong head m, appending the word
d_(L-1) ... d_1 (1 = Z) gives a n - b g(n) = 2^s (a'' m - b' g(m)) + kappa with a''/b' = frac(2^(L-1) a/b), so these
n alone fill C*: C* is the class of a n - b g(n) for any odd n with n - 1 = sum_(i<L) (1 - d_(i+1)) 2^i
(mod 2^L), n > 2^L. Examples: 5/6, C* = 3 (mod 4); 7/12, C* = 5 (mod 8); 11/12, C* = 7 (mod 8).

**T3 (exact half-lines, b' > 1, s >= 1).** S(a, b) agrees with C* on (c_+, infinity) iff d_2 = ... = d_(s+1) = 0,
with c_+ = 5(a - b) (s = 1) or c_+ = F(2^s - 1) = 2^(s-1)(2a - b) - a - b (s >= 2); on (-infinity, c_-) iff
d_1 = ... = d_(s+1) = 1, with c_- = (2^s + 1)(a - b) = F(2^s + 1); both thresholds are elements of S(a, b) outside
C*; otherwise S(a, b) minus C* is unbounded on both sides. For b = 2 (mod 4) exactly one side is exact: a/b < 3/4
gives S(a, b) n (5(a - b), infinity) = {r = 2 - a (mod 4)} there, a/b > 3/4 gives
S(a, b) n (-infinity, 3(a - b)) = {r = -a (mod 4)} there (the equality versions of THM-4882 (3)(iii), whose sharp
threshold -3(b - a) is concurrent).

**T4 (C2 offsets of every left source).** For a nontrivial physical source v and a strong unit source
q = (D, 1, U, 1), the actual C2 pair n = N(v, Oq), m = M(v, Jq) of the complete C2 theorem satisfies
m - 2n = alpha D - beta U + c, alpha = D_v + 4U_v - 2d_v (odd), beta = alpha + E_v - d_v, c = -D_v + U_v + d_v + 2E_v.
Since E_v <= d_v: if E_v < d_v the offsets fill exactly the class U_v - D_v + E_v + 2d_v modulo
2 gcd(alpha, d_v - E_v), each with infinitely many parents (every odd offset when d_v - E_v = 1, e.g. v = t, O t, an
O-born source over a source with d = 1, a J-born source over a source with E = 1); if E_v = d_v they are at most c,
with equality attained. Session 198's section 5 is the case v = t. No parent count is prescribed.

**T5 (clock uniqueness).** On strong unit sources, appending an O/Z word w of length k changes aU - bD by
(a(2^k - 1) - b A_w) U + const (A_w = sum_i eps_i 2^i); it is source-independent iff b A_w = a(2^k - 1), which forces
0 <= a/b <= 1, b odd and ord_b(2) | k. For 1/2 < a/b < 1 these words are exactly the powers of Session 198's clock, so
an interior signed functional has a single source-independent move in this class; for even b there is none at any
slope; for 0 < a/b <= 1/2 with b odd they are the powers of the word of a(2^ord - 1)/b (e.g. ZO for 1/3), which are not
Session 198 clocks.

## Verification

`05-knowledge/results/signed_response_2adic_structure_20261105/`: `verify.py` (normal and -O receipts identical),
`witness_length.c`; pre-registration commit 2c11d95d56. P1: 629 slopes b <= 64, full-class count and one-sidedness
exact, zero violations below 2^22; P2: all 544,115 targets of C* with |r| <= 2000 for 340 slopes b <= 48 witnessed
(direct search or the T2 recipe) and checked by two independent evaluations of g, 2,033 literal tournaments; P3: 186
half-line cases; P4: 78,000 offset identities and six literal C2 pairs (offsets +-1, five left profiles, exactly two
Hamiltonian cycles each); P5: 59 class densities; P6: witness lengths. The audit reproduced every receipt number with
its own code (census to 2^24, 66 further slopes up to b = 256 and 20 exterior slopes, 15,972,450 offset pairs, nine
literal C2 pairs). Outcomes and refuted clauses: HYP-9476.
