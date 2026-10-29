---
id: THM-4847
title: "Fermat primes: what a proof of 'exactly five' must get past, and what tournaments can and cannot see. (A) Barriers. (1) Coprimality: a prime divides at most one Fermat number, so no finite set of primes covers {F_n : n >= 5}, and no covering-congruence proof exists. (2) Periodicity: for every modulus M the sequence F_n mod M is eventually periodic and is a non-unit for at most omega(M) indices, so no fixed-modulus certificate exists. (3) Algebra: X^(2^n)+1 = Phi_(2^(n+1))(X) is irreducible over Q, and the Aurifeuillean identity for 2^m+1 needs m = 2 (mod 4), i.e. n = 1. (4) Logic: 'exactly five' is the Pi^0_1 sentence 'for all n >= 5, F_n is composite', and it implies the open statement that infinitely many Fermat numbers are composite. (B) In the function-field analogue the statement is a classical theorem (Lidl-Niederreiter Thm 2.47; typed ANALOGY with the integer problem). Over F_q with q odd, f^(2^n)+1 = Phi_(2^(n+1))(f) is reducible for every nonconstant f in F_q[T] and every n >= 2, because (Z/2^(n+1))^* is not cyclic. In characteristic 2 this holds for n >= 1. The mechanism has no counterpart over Q. (C) Tournament counts are blind; three-state counts are Pepin's test. For n >= 2 and N = F_n, the number of labelled tournaments 2^C(N,2) and the number of circulant tournaments 2^((N-1)/2) are both 1 (mod N), whether F_n is prime or not. So every Fermat number passes THM-4746's tournament-count test. For n >= 1, the number 3^C(N,2) of labelled semicomplete digraphs (equivalently oriented graphs) on F_n vertices is -1 (mod F_n) iff F_n is prime. At a prime p the strong count is congruent to the total count (binomial vanishing in Wright's recurrence), so s_b(p) = (b/p): every Fermat prime p >= 5 has s_3(p) = -1, while s(p) = (2/p) is -1 at p = 5 and +1 at 17, 257, 65537 (audit correction). F_5, F_6, F_9 and F_11 fail both strong rotation tests at known prime factors (F_6, F_9, F_11 by the audit). (D) With THM-4644 (Fermat primes are exactly the primes all of whose circulant tournaments are rigid): tournament symmetry sees Fermat primality, tournament counting modulo N does not."
status: >
  PROVED (elementary; mac-mini-2026-10-29-growth): (A1)-(A4) and the counting forms in (C). Independently audited
  2026-10-29, corrections applied (04-computation/experiments/growth_gf_runs_fermat_20261029/audit/REPORT.md,
  MISTAKE-695). Two of those corrections: the strong-count sign in (C3) is (2/p), which is -1 at p = 5; and (B) is
  CLASSICAL (Lidl-Niederreiter, Finite Fields, Thm 2.47), typed ANALOGY.
  Classical inputs:
  - Goldbach's pairwise coprimality;
  - Euler-Lucas (every prime factor of F_n, n >= 2, is 1 mod 2^(n+2));
  - Pepin's test, with 3^((F_n - 1)/2) = (3/F_n) by quadratic reciprocity;
  - irreducibility of cyclotomic polynomials over Q, and the factorisation of Phi_m over F_q into factors of degree
    ord_m(q).
  FINITE-EXACT (fermat_camion.py; audit/fermat_primes_direct.py, fermat_extend.py):
  - the count test for F_2..F_11;
  - the F_5 strong tests at 641 (residue r = 257) and 6700417 (r = 641);
  - direct strong counts mod p (Wright's recurrence with exact binomials mod p^2):
    s_2(p) = -1, -1, +1, +1 and s_3(p) = 0, -1, -1, -1 at p = 3, 5, 17, 257;
  - F_6 at 274177 and 67280421310721, F_9 at 2424833 and F_11 at 319489 and 974849 all fail both strong tests
    (audit).
  CITED STATUS (checked 2026-10-29):
  - F_n is composite for 5 <= n <= 32 and F_33 is the smallest of unknown character (Wikipedia, Fermat number).
  - It is open whether infinitely many Fermat numbers are composite, and whether all of them are squarefree.
  - Boklan-Conway, arXiv:1605.01371 (Math. Intelligencer 2017), bound the expected number of new Fermat primes by
    about 2^-30. HYP-9410 has the repo's calibrated heuristic (about 7e-9 beyond F_32).
  INHERITED: THM-4644, THM-871 (circulant rigidity), THM-4746 (Camion tests, the semicomplete law), THM-4793/4800
  (Pepin as Collatz dynamics), THM-4801 (constructible polygons, initial-segment principle).
  - (A1) and (A3) restate S58's section 1.2 (05-knowledge/results/s58_fermat_pepin_collatz_20261026/README.md):
    covering congruences, algebraic and Aurifeuillean no-go, assembled there from classical facts.
  - The status facts and the Boklan-Conway figure are also in S58 section 1.1.
  CONCURRENT: the research packet 05-knowledge/results/four_vertex_glue_fermat_20261029/FERMAT.md (pushed 09:00 the
  same day, read after this checkpoint) has three parallel results:
  - (F6): the source 2 is blind, S^n(2) = -1 and S^(n+1)(2) = 1 mod F_n, prime or composite. This is the dynamical
    twin of our (C1).
  - The odd-order-unit certificate of compositeness, with an explicit F_5 witness x^3 = 1.
  - Exact squaring-digraph formulas.
  NEW HERE: (A2) as stated, (A4), the counting forms (C1)-(C3), the Camion computations and the routes table.
  (B) is classical.
  NOT PROVED, and believed far beyond current methods: that there are exactly five Fermat primes.
session: mac-mini-2026-10-29-growth
source: 05-knowledge/results/tournament_growth_gaussian_gf_runs_fermat_20261029.md
depends_on:
  - 01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md
  - 01-canon/theorems/THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md
related:
  - 01-canon/theorems/THM-871-fermat-rung-rigidity.md
  - 01-canon/theorems/THM-4793-the-fermat-line-under-collatz-valuation-two-descent-exact-merge-of-consecutive-lines-and-pepin-as-an-affine-cycle-type.md
  - 01-canon/theorems/THM-4800-pepin-as-the-valuation-two-collatz-branch-trivial-cycle-half-turn-fermat-line-arc-mersenne-mirror.md
  - 01-canon/theorems/THM-4801-collatz-meets-the-constructible-polygons-trunk-five-initial-segment-principle-collatz-closed-family.md
  - 05-knowledge/hypotheses/HYP-9410-collatz-line-calibration-of-the-fermat-prime-heuristic-algebraic-factorisations-and-immunity-of-genuine-fermat-lines.md
  - 05-knowledge/hypotheses/HYP-9445-run-observer-torsion-and-hamiltonian-cycles-star-law-fermat-routes-and-gaussian-gf-questions.md
scripts:
  - 04-computation/experiments/growth_gf_runs_fermat_20261029/fermat_camion.py
outputs:
  - 04-computation/experiments/growth_gf_runs_fermat_20261029/out/fermat_camion.out
---

# THM-4847 — Fermat primes: barriers, the function-field analogue, and the counting blindness of tournaments

`F_n = 2^(2^n) + 1`. The five known Fermat primes are `F_0, …, F_4 = 3, 5, 17, 257, 65537`. Their product is
`F_5 − 2 = 2^32 − 1`. "Exactly five" is equivalent to each of the following:
- every odd `n` with `φ(n)` a power of 2 divides `2^32 − 1`;
- there are exactly 31 odd constructible polygons (Gauss–Wantzel; THM-4801's initial-segment principle).

## 1. Statements

**(A) Barriers.**
1. *Coprimality.* `F_n − 2 = F_0 F_1 ⋯ F_{n−1}`, so a common prime factor of `F_m` and `F_n` (`m < n`) would divide 2.
   - Hence each prime divides at most one `F_n`, and no finite set of primes divides infinitely many of them.
   - So no Sierpiński-style covering proof can show that all `F_n` with `n ≥ 5` are composite.
   - By Euler–Lucas, a composite `F_n` (`n ≥ 2`) has a prime factor `p` with `2^(n+2) + 1 ≤ p ≤ √F_n` and
     `p ≡ 1 (mod 2^(n+2))`. Such witnesses are pairwise distinct and unbounded.
2. *Periodicity.* For any `M`, the sequence `x_{n+1} = x_n² mod M` is eventually periodic, so `F_n mod M` is
   eventually periodic. By (1), `gcd(F_n, M) > 1` for at most `ω(M)` indices. A certificate of the form "`F_n` shares
   a prime factor with the fixed modulus `M`" therefore decides compositeness for only finitely many `n`.
3. *Algebra.*
   - `X^(2^n) + 1 = Φ_(2^(n+1))(X)` is irreducible over `Q`, so the family `X^(2^n) + 1` has no algebraic
     factorisation. Any single composite number is of course the value at `X = 2` of some product of polynomials, so
     the barrier concerns uniform factorisations of the family.
   - The Aurifeuillean factorisation of `2^m + 1`, namely `4y⁴ + 1 = (2y² + 2y + 1)(2y² − 2y + 1)`, needs
     `m ≡ 2 (mod 4)`, and among Fermat exponents `m = 2^n` only `n = 1` qualifies (`F_1 = 5 = 1·5`).
   - HYP-9410 records that genuine Fermat lines carry no such factorisation.
4. *Logic.* Compositeness of `F_n` is decidable by a bounded search, so "exactly five" is the `Π⁰₁` sentence
   `∀n ≥ 5 ∃d < F_n (1 < d ∧ d | F_n)`.
   - If it is false, a single computation refutes it: Pépin's test at the new prime.
   - If it is true, it implies "infinitely many Fermat numbers are composite", which is open.

**(B) The function-field analogue is a classical theorem** (Lidl–Niederreiter, *Finite Fields*, Thm 2.47, on the
factorisation of cyclotomic polynomials over `F_q`; that `X⁴ + 1` is reducible modulo every prime is a textbook
example). It is typed ANALOGY with the integer problem. Let `q` be odd and `f ∈ F_q[T]` nonconstant.
- For every `n ≥ 2`, `f^(2^n) + 1` is reducible in `F_q[T]`.
- For `n = 1`, `T² + 1` is irreducible iff `q ≡ 3 (mod 4)`.
- In characteristic 2, `f^(2^n) + 1 = (f + 1)^(2^n)` is reducible for every `n ≥ 1` whenever `f + 1` is nonconstant.

So in `F_q[T]` there are no "Fermat primes" beyond `n ∈ {0, 1}`, at any base. Over `Q` the mechanism is absent,
because `Φ_(2^(n+1))` is irreducible for all `n`. The integer question concerns prime values of a different
irreducible polynomial for each `n`, at the fixed point `X = 2`. This is the opposite regime to
Bunyakovsky–Schinzel–Bateman–Horn, which fix the polynomial and vary the argument.

**(C) What tournament counts see.** Let `N = F_n`.
1. *Tournaments are blind (n ≥ 2).*
   - `2^(2^n) ≡ −1 (mod N)` gives `ord_N(2) = 2^(n+1)`, and `2^(n+1) | (N−1)/2 = 2^(2^n − 1)`.
   - Hence the circulant count `c(N) = 2^((N−1)/2) ≡ 1` and the labelled count `T(N) = 2^(C(N,2)) = c(N)^N ≡ 1 (mod N)`,
     prime or not.
   - Every Fermat number therefore passes THM-4746 §6's tournament-count test, and at a Fermat prime the Euler value
     `(2/N) = 1` carries no information.
2. *Three states are Pépin (n ≥ 1).* `3^(C(N,2)) ≡ −1 (mod N)` iff `N` is prime.
   - `3^(C(N,2))` counts labelled semicomplete digraphs on `N` vertices (each pair: `→`, `←` or `↔`), and equally
     labelled oriented graphs (each pair: `→`, `←` or nothing).
   - The circulant count is `3^((N−1)/2)`, the quantity in Pépin's test.
3. *Strong counts at a prime.* In Wright's recurrence `T(p) = Σ_{k=1}^{p} C(p,k) s(k) T(p−k)`, every `0 < k < p` term
   vanishes mod `p`. So `s(p) ≡ T(p) = b^(C(p,2)) ≡ (b/p) (mod p)` for any family with `b` states per pair and a
   linear condensation.
   - At a Fermat prime `p ≥ 5` this gives `s_3(p) ≡ −1` for strong semicomplete digraphs. For strong tournaments it
     gives `s_2(p) ≡ (2/p)`, which is `−1` at `p = 5` and `+1` at `17, 257, 65537`. At `p = 3`, `s_2(3) = 2 ≡ −1`,
     and `s_3(3) ≡ 0` (the semicomplete law needs `p ≥ 5`). Corrected in the audit; the draft said `s_2(p) ≡ 1` for
     every Fermat prime `p ≥ 5`.
   - The converse for the strong counts is not claimed. `F_5` fails both strong tests at both factors (§3).

**(D) Symmetry versus counting.** THM-4644: an odd prime `p` is a Fermat prime iff every circulant tournament on
`Z/p` is rigid (its automorphism group is `Z/p`), iff the multiplier group acts freely on the `2^((p−1)/2)`
connection sets. So the automorphism shadow of tournaments characterises Fermat primes exactly. By (C1), their
counting shadow modulo `N` is constant on all Fermat numbers.

## 2. Proofs

**(A1).** Goldbach's identity `Π_{i<n} F_i = F_n − 2` follows by induction from `(x − 1)(x + 1) = x² − 1` with
`x = 2^(2^(n−1))`. Euler–Lucas: if `p | F_n`, then `ord_p(2) = 2^(n+1)`. For `n ≥ 2`, `p ≡ 1 (mod 8)`, so 2 is a
quadratic residue, `2^((p−1)/2) ≡ 1`, and `2^(n+1) | (p−1)/2`.

**(A2).** A finite-state map is eventually periodic. The `ω(M)` bound is (A1).

**(A3).** `Φ_(2^(n+1))` is irreducible over `Q` (Gauss/Eisenstein at 2 after `X → X + 1`). For the Aurifeuillean
form, `2^m + 1 = 4y⁴ + 1` requires `m = 4j + 2`.

**(A4).** Trial division below `√F_n` is a bounded quantifier.

**(B).** `Φ_(2^(n+1))` factors over `F_q` into `φ(2^(n+1))/d` irreducible factors of degree `d = ord_(2^(n+1))(q)`. For
`n + 1 ≥ 3`, `(Z/2^(n+1))^* ≅ C_2 × C_(2^(n−1))` has exponent `2^(n−1) < 2^n = φ(2^(n+1))`, so there are at least two
factors `g`, `h`. Then `Φ(f) = g(f) h(f)`, and both factors have positive degree because `f` is nonconstant. The
case `n = 1` holds because `−1` is a square in `F_q` iff `q ≡ 1 (mod 4)`. In characteristic 2 the Frobenius gives
`X^(2^n) + 1 = (X + 1)^(2^n)`.

**(C1).** Shown in the statement. `ord_N(2)` divides `2^(n+1)` but not `2^n`, since `2^(2^n) ≡ −1 ≠ 1`.

**(C2).**
- (⇐) If `N` is prime, `N ≡ 1 (mod 4)` and `N ≡ 2 (mod 3)`, so `(3/N) = (N/3) = −1`. Then
  `3^(C(N,2)) = (3^((N−1)/2))^N ≡ (−1)^N = −1`.
- (⇒) Let `p | N` be prime (`p ≠ 3`, since `F_n ≡ 2 (mod 3)`). From `3^(N(N−1)/2) ≡ −1 (mod p)` we get
  `3^(N(N−1)) ≡ 1`, while the half power is not 1. So the 2-part of `ord_p(3)` equals the 2-part of `N(N−1)`, which
  is `2^(2^n)` because `N` is odd. Hence `2^(2^n) | p − 1`, so `p ≥ N`, so `p = N`.

  This is Pépin's argument applied to `N`-th powers.

**(C3).** `C(p,k) ≡ 0 (mod p)` for `0 < k < p`, and `T(0) = 1`. Euler's criterion gives `b^(C(p,2)) = (b^((p−1)/2))^p`.
For `b = 2` this recovers THM-4746's `s(p) ≡ (2/p)` in one line. The quasi-periodicity of THM-4746 Theorem 1 is the
stronger statement needed for composite orders.

## 3. Computation (FINITE-EXACT, `out/fermat_camion.out`)

**The count test.** It holds for `F_2, …, F_11`: `ord(2) = 2^(n+1)` divides `(N−1)²/2`.

**`F_5 = 641 · 6700417`.** The residues are `r_p ≡ F_5 (mod p − 1)`, which equals the cofactor reduced mod `p − 1`:
`r_641 = 257 = F_3`, and `r_6700417 = 641`.

| prime factor `p` | strong tournament test (`D(r) mod p`, passes iff 0) | strong semicomplete residue (passes iff `p − 1`) |
|---|---|---|
| 641 | `254` | `150` |
| 6700417 | `4217483` | `5459519` |

So `F_5` is neither a Camion pseudoprime nor a "semicomplete Camion pseudoprime".

**Direct strong counts at the Fermat primes** (audit, `audit/fermat_primes_direct.py`). These use Wright's
recurrence with exact binomials mod `p²`, without assuming the binomial vanishing, and reduce mod `p`:
`s_2(p) ≡ −1, −1, +1, +1` and `s_3(p) ≡ 0, −1, −1, −1` at `p = 3, 5, 17, 257`. An earlier draft listed
"`s_3(p) ≡ −1` at 5, 17, 257" from `3^(C(p,2))` alone, which is circular. `fermat_camion.py` now computes the
recurrence directly.

**`F_6`, `F_9` and `F_11`** (audit upgrade, `audit/fermat_extend.py`) also fail both strong rotation tests:
- `F_6` at 274177 (`r = 65537`: `D ≡ 215736`, `s_3 ≡ 123651`);
- `F_6` at 67280421310721 (`r = 274177`: `D ≡ 32409983229131`, `s_3 ≡ 35318451187950`);
- `F_9` at 2424833 (`r = 786433`: `D ≡ 522107`, `s_3 ≡ 623467`);
- `F_11` at 319489 (`r = 90113`: `D ≡ 61029`, `s_3 ≡ 273957`) and at 974849 (`r = 761857`: `D ≡ 350016`,
  `s_3 ≡ 751119`).

A single failing prime factor already decides each test. These are compositeness certificates; they prove nothing
about the remaining Fermat numbers.

## 4. Routes to "exactly five", and where each stands

| route | status |
|---|---|
| covering congruences / finite sieve | impossible (A1, A2) |
| algebraic or Aurifeuillean factorisation | impossible for `n ≥ 2` (A3) |
| a uniform failure of Pépin, `3^((F_n − 1)/2) ≢ −1` for all `n ≥ 5` | equivalent to the statement; no approach known |
| reformulations: circulant rigidity (THM-4644), Pépin dynamics (THM-4800), constructible polygons (THM-4801), (C2) | equivalent restatements; they move the problem, they do not shrink it |
| analytic: sieve and prime-number methods | no current method reaches doubly exponential sequences; even infinitely many composite `F_n` is open |
| heuristic | finitely many expected; at most about `2^-30` new Fermat primes (Boklan–Conway), about `7·10^-9` beyond `F_32` (HYP-9410) |
| function-field analogue | a classical theorem (B), by a mechanism that is absent over `Q` (ANALOGY) |
| tournaments | counting is blind (C1); symmetry is exact (D); three-state counting is Pépin (C2) |
| odd-order units | `F_n` is composite iff some unit `x ≠ 1` has odd order (concurrent FERMAT.md); a uniform construction of such a witness for all `n ≥ 5` would suffice and is not known |

The honest verdict is that "exactly five Fermat primes" is a `Π⁰₁` statement that no known method can approach. Any
proof would have to produce, for each `n ≥ 5`, a structure not visible to congruences of bounded modulus, to
polynomial algebra, or to tournament counting modulo `F_n`.
