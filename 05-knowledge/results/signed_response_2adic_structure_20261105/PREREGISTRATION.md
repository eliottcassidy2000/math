# Pre-registration: exact 2-adic structure of signed one-apex responses (kind-pasteur-2026-11-05-S12)

Committed and pushed **before** the confirmatory run `verify.py` of this packet. Canon reservations: THM-4883,
HYP-9476 (commit 6cb451411f).

## Objects

g(1) = 1, g(2n) = g(n), g(2n - 1) = g(n) + n - (n mod 2). For coprime integers a, b with b >= 1,

    S(a, b) = { a n - b g(n) : n >= 3 odd }.

By Session 198 (rational-response-coverage, eq. (7)) and THM-4882 (1), S(a, b) is exactly the set of signed responses
aU - bD over the strong unit-deletion sources (the one-apex sources with odd cover count n >= 3, profile
(g(n), 1, n, 1)); every element has an actual O/Z source graph.

Write b = 2^s b' with b' odd and let P_s be the 2^s residue classes modulo 2^(s+1) of the integers congruent to a - b
modulo 2. A class C is **full** if C is contained in S(a, b); it is **one-side bounded** if S(a, b) meets C in a set
bounded above or bounded below. Write a/b = 0.d_1 d_2 d_3 ... in binary.

## Claims to be proved (proof in PROOF.md; not decided by computation)

- T1 (2-adic structure). For every rational slope, S(a, b) has a full class in P_s iff 1/2 < a/b < 1 and b' > 1 (b not
  a power of two). It then has exactly one, C*, and every other class of P_s is one-side bounded. For a/b outside
  (1/2, 1), or b a power of two, every class is one-side bounded.
- T2 (explicit C*). If b' > 1 and 1/2 < a/b < 1, let L be the least position >= s + 1 with d_L = 1. Then C* is the class
  of a n - b g(n) for every odd n >= 3 whose O/Z word (from the two-vertex seed) ends with the letters
  d_(L-1), ..., d_1 (chronologically; 1 = Z, 0 = O) preceded by a Z.
- T3 (exact half-lines, b' > 1, s >= 1). S(a, b) agrees with C* on a positive half-line iff d_2 = ... = d_(s+1) = 0; then
  S(a, b) on (c_+, infinity) equals C* there, with c_+ = 5(a - b) if s = 1 and c_+ = max(5(a - b), F(2^s - 1)) if
  s >= 2, attained. It agrees with C* on a negative half-line iff d_1 = ... = d_(s+1) = 1; then on (-infinity, c_-)
  with c_- = (2^s + 1)(a - b), attained at n = 2^s + 1. Otherwise on neither.
- T4 (C2 offsets for every left source). For every nontrivial physical source v and strong unit source q,
  M(v, Jq) - 2 N(v, Oq) = alpha D_q - beta U_q + c with alpha = D_v + 4U_v - 2d_v (always odd),
  beta = alpha + E_v - d_v, c = -D_v + U_v + d_v + 2E_v. If E_v < d_v the offsets over all such q fill exactly one
  residue class modulo 2 gcd(alpha, d_v - E_v) (every odd offset when the gcd is 1); if E_v = d_v they are bounded
  above by c.
- T5 (clock uniqueness, unit class). For a fixed reduced slope, an O/Z word changes aU - bD by an amount independent
  of the input unit source iff b is odd and the word is a power of Session 198's clock word (so b even has none).

## Confirmatory predictions (decided by `verify.py`, normal and `-O` runs must agree byte for byte)

Exploration before this file (scratch scripts, not committed) already ran a version of P1 on all 629 reduced slopes
with b <= 64 (n < 2^21, zero bound violations) and witness-length searches for 2/3, 4/5, 6/7. The predictions below
use new code, wider ranges or new statistics.

- **P1 (structure census).** For every reduced a/b in (1/2, 1) with 3 <= b <= 64, the automaton of PROOF.md section 4
  returns exactly one full class iff b is not a power of two, and a sign and bound for every other class. Over all
  odd n < 2^22: zero values violate a bound, and every bound equals the brute-force extremum of its class (attained).
- **P2 (fullness, exhaustive window).** For every reduced a/b in (1/2, 1) with b <= 48 and b' > 1, every target r in
  C* with |r| <= 2000 receives an explicit witness n with a n - b g(n) = r, checked by two independent evaluations of
  g (recurrence and closed bit formula); witnesses come from direct search (n < 2^24) or, if that fails, from the
  constructive recipe of PROOF.md (tail word + modular head + clock). Zero failures. Every witness with n < 2^12 is
  also checked on its literal tournament by subset dynamic programming (H = g(n), one root-start path, d = 1).
- **P3 (exact half-lines).** For every slope of P1 with b' > 1 and s >= 1, the digit criterion of T3 agrees with brute
  force over odd n < 2^22: in the predicted exact cases no element outside C* lies beyond c_+ (resp. c_-) and c_+
  (resp. c_-) is attained outside C*; in the "neither" cases elements outside C* above +200 and below -200 both occur.
- **P4 (C2 offsets).** The offset identity holds on a coordinate basis and on every pair (v, q) with v in the physical
  source bank (U_v <= 300) and q among 40 strong unit sources; the slope and gcd statements hold for every bank
  source; and at least six literal C2 parent/child graph pairs from three left sources other than t, with offsets
  -1 and +1, have subset-DP Hamilton counts n, m = 2n + delta and exactly two Hamiltonian cycles each.
- **P5 (density, HYP-9476 (A)).** For the 14 slopes 3/4, 5/6, 7/10, 9/10, 7/12, 11/12, 9/14, 13/14, 5/8, 7/8, 13/24,
  17/24, 19/24, 23/24, every non-full class whose pieces all have slopes different from 1 and 1/2 (so that odd
  n <= 2^22 decides the window exactly; the run checks the needed bound on n) and that has at least 50 elements of S
  in the window: the fraction of the class hit on its unbounded side within 500 <= |r| <= 20000 lies strictly between
  0.02 and 0.98. Classes containing a slope-1 or slope-1/2 piece are reported as lower bounds only. Prediction: no
  exception among the decided classes (no non-full class looks cofinite).
- **P6 (witness size, HYP-9476 (B)).** For the slopes 3/5 and 8/9 (not explored here), exhaustive search over all O/Z
  words of length <= 34: (a) at least 99% of the targets |r| <= 2^10 of the right parity are found; (b) the maximal
  minimal witness length over |r| <= 2^j is at most 3j + 8 for 4 <= j <= 9; (c) among targets whose minimal witness
  ends with t >= 1 full clock periods, at most one half end with t >= 2.

A failed prediction is recorded as REFUTED in HYP-9476, not repaired after the fact.
