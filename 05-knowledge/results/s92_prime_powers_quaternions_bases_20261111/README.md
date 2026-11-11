# S92 (2026-11-11): F1 at every odd modulus, the four-square series, and the algebras, bases and Collatz tables around it

Session `opus-2026-11-11-S92`. Canon: THM-5134, THM-5135. Hypotheses: HYP-9665, HYP-9666. Scripts and outputs:
`04-computation/experiments/s92_prime_powers_20261111/`. Audited 2026-11-11, PASS WITH CORRECTIONS (section 11).
The Chan-Pak conjecture (every odd number except 7 and 21 is a Hamiltonian path count) is NOT proved here and
remains OPEN: every local density below is positive, at 7 and 21 as elsewhere, and decides nothing about it. The
Collatz conjecture is OPEN as well, and nothing below bears on it: sections 6 and 7 relabel the arcs of `3n + 1` and
tabulate the known cycles, no more.

The owner's directive had seven parts. Short answers, with the section that holds the detail:

| asked | answer | section |
|---|---|---|
| the odd prime-power case of THM-5128, the series `sigma_1(n)/n` | PROVED at every odd modulus (THM-5134); the pre-registered tests are confirmed as registered, with partial blindness and uneven power | 1, 2 |
| quaternions versus bicomplex numbers, dihedral groups | prior work listed; the eight extensions of `F_2^2` put the four algebras in one table (THM-5135 (5); classical) | 3 |
| `(8/pi^2) sigma_1(n)/n` and the octonions / Fano plane | the series is the quaternion (rank 2) member of a ladder whose rank-4 member is E8, `240 sigma_3(n)`; the repository's octonion work concerns sign patterns, and no link to the series is known on the tournament side; no tournament family is known to carry the rank-4 law | 4 |
| Tao's toy example, all bases | taken to be the money-plant puzzle set by Terence Tao, as in Session266, which recorded and generalised it the same morning; the F1 joins are such quadratics; their local law holds in every base, and their second digit is equidistributed in every odd base (THM-5134 (G)) | 5 |
| multiplication tables modulo odd values and arcs of `3n + 1`, `3n - 1` | the dictionary exists (credited: Session255, the codex atlases, Session266); added: the F1 density is the additive square of the table, two relabellings and one cross-reference | 6 |
| three negative cycles against one positive, and the 4-vertex tournaments | a census of the known cycles; exact pieces (most of them KNOWN from Session266) and one coincidence, typed; no residue-level statement distinguishes the two Collatz sides, so no local density can see the difference | 7 |
| the pasted note on strong components as primes | profiles and generating function confirmed; the six strong classes of order 5 are all self-converse (KNOWN, Session224) | 8 |

A sweep of the repository by eight independent readers preceded the write-up; what it found is credited in
section 10. Section 9 types every bridge. Notation: `sigma_1(n)` is the sum of the divisors of `n`; `sigma` alone is
the sheet index (`c_sigma`); "sheet" alone means the F1 sheet (`s = 0`: sheet 0, the head sheet; `s >= 1`: sheet 1,
the tail sheet); the sign `eps` of a table label `m = 4x + eps` is its Collatz side (positive `eps = +1`, negative
`eps = -1`); `x_b(h) = 1 + b psi(h)` is the base-`b` state (`x_2 = x`); `X` is always the bound `N <= X`.

## 1. The law at every odd modulus (THM-5134)

F1: `N = x1 x2 + 2^s y1 y2 - [s = 0]` over pairs of one-apex states `(x, y) = (2 psi(h) + 1, 2h + 1)`. THM-5128
proved the class law by value at odd primes and squarefree moduli and left odd prime powers numerical ("not
written out", its (D3)). THM-5134 closes that step. For an odd modulus `m`:

- the type of a triple is `d = gcd(x1, y1, m)` for its smaller state; a triple of type `d` has `d | N - c_sigma`
  exactly (`c_1 = 0`, `c_0 = -1`), and its value is uniform on those classes with a power saving;
- the types have limiting shares `c_d = d^-2 prod_{p | m/d} (1 - p^-2)`;
- so the class factor tends to `L_m(e) = sum_{d | gcd(m, e)} d^-1 prod_{p | m/d} (1 - p^-2)`, `e = N - c_sigma`, which
  is multiplicative, with `L_(p^k) = (1 - p^-2) sigma_1(p^v)/p^v` at valuation `v < k`;
- the product over the odd primes is `S(N) = (8/pi^2) sigma_1(n)/n`, `n` the odd part of `N` (sheet 1) or of `N + 1`
  (sheet 0).

The factor for uniform states and the series were written down by the audit of THM-5128; the fixed-source law
without a rate is Session206's; the power saving and the sum over the smaller state are THM-5128 (B4)(iv), (C2),
(D2), and prime and squarefree `m` are THM-5128 (C1), (C2), (D2), (D3). New: the argument written out at the odd
moduli that are not squarefree (fibre count for a smaller state of content `d`, shares `c_d`, closed form `L_m`).

Data (exact censuses to `2^28`; 854,415,478 and 1,083,261,684 triples on the two sheets). The remainder of the
first-order law times `m`, at `2^28`: `0.02` at 9, `0.03` at 27 and 81, `0.004` at 25, `0.005` at 125, `0.10` at 49,
`0.04` at 45, `0.03` at 75, far above the sampling noise. It decreases at 5 and 7 (49: `0.163 / 0.141` at `2^24`,
`0.093 / 0.103` at `2^28`; 7 is the slow prime, as in THM-5128); at the powers of 3 it oscillates between `0.018` and
`0.042` from `2^24` to `2^28`, with no net decrease on sheet 0 (`0.021, 0.022, 0.024` at `2^24`, `0.024, 0.029, 0.031`
at `2^28`); the audit's census at `2^30` has `0.012 / 0.021` at 9 and `0.058 / 0.079` at 49
(audit/proof-5134-AE/summary_tables.out, audit/numerics-prereg/octaves.out). The factors are far from
their limits (the distinguished class modulo 9: `1.45` on sheet 0, `1.46` on sheet 1, against `35/27 = 1.296`), as
the `1/ln X` term predicts (NUMERICAL, audit, proof-5134-AE/a_constants.out: limit plus that term `1.441`,
`1.454`; at 49 the constants are not determined, as `A_7` in THM-5128).

Pre-registered (PREREGISTRATION.md; the registration commit 73c8437aac precedes the commit of the outputs):
window `[2^27, 2^28)`, eight moduli; a group is the set of the classes with the same `g0 = gcd(m, N - c_sigma)`. P1
(a-priori factors from the smaller states alone, within tolerance), P2 (flatness inside each group), P3 (closer to
the first-order value than to the limit at 9, 27, 25: 14 of 14, 12 of them distinct) CONFIRMED; P4 (remainder
smaller at `2^28` than at `2^26` in at least 8 of 12 pairs) CONFIRMED with exactly 8. The verdicts stand as
registered; the audit qualifies them (details in THM-5134, Numerical support).
Power of P1: the tolerances were set from the in-sample maxima; they are decisive at 9, 27, 81, 45; at 25, 125, 75
a law without the valuation split would also have passed; at 49 there is none (the tolerance `0.115` is of the
size of the predicted excess, `0.097` to `0.119`); so the valuation-resolved factor is tested at 3 only. P3 was
registered for 9, 27, 25 only and fails at 49 (0 of 4 groups in sample, on the window and cumulatively up to
`2^30`): at 49 the remainder still exceeds the first-order correction in every group with `g0 > 1` (at 63 and 21,
at `2^28`, in the groups `g0 = 7`, and `g0 = 21` at 63). P4 is unstable: two of its eight decreases are of
`0.00013`, below the noise scale, and the count is 4 of 12 from `2^26` to `2^27`, 10 of 12 from `2^27` to `2^28`.
(The tolerance at 49 is below three of the four in-sample values, not two as the registration says.)

*Blindness* (audit, numerics-prereg/prior_exposure.out). The registration (11:29) precedes the commit of the
outputs of the registered window and of `2^28` (11:43); the in-sample censuses at these moduli belong to the
registration commit. It was not blind at the prime level: the audit of THM-5128 had enumerated the same window at
prime moduli (its report prints the excess at 5), and THM-5128 (C) gives the remainders at `2^26` and `2^28` for 3,
5, 7. These numbers fix the `g0 = 1` group means of P1 at the six prime powers and the mean over the favoured
classes, and the prime-level direction agrees with the prime-power one in 9 of the 12 pairs of P4. Blind: the
favoured class split by valuation, flatness in the groups, 45 and 75.

What it does not give: an asymptotic for the multiplicity `r` of one `N`. On the odd `N` of `[2^24, 2^25)`, `r`
correlates with `sigma_1(n)/n` (`0.40`, `0.43`) and doubles across its quintiles (exploratory; HYP-9665 (B)). The
series varies by `1.61` over these quintiles and the slope of `r/mean(r)` against `S` is `1.44`, `1.43`: steeper than
the limit, as the positive `1/ln X` terms at 3 and 5 predict; the series explains 16 and 18 % of the variance of
`r`, and 8.07 % and 4.86 % of these `N` have `r = 0` (audit, proof-5134-AE/series_vs_counts_audit_25_24.out).

## 2. Three classical faces of `sigma_1(n)/n`

1. *The multiplication table.* The number of `(a, b)` modulo `m` with `a b = t` is
   `N_m(t) = sum_{d | gcd(t, m)} d phi(m/d)` (in the modular-centres packet, P4). The F1 seam is a sum of two table
   entries (`2^s` is a unit), and the additive square of the table is `(N_m * N_m)(c) = m^3 L_m(c)` (eleven odd
   `m <= 105`; audit: every `m <= 210`, even `m` included): `L_m` is the density of a sum of two table entries. The
   box version (2 x 2 matrices of determinant `t` in a box) is overnight_research_20261029/TWO_PORTS.md P5.
2. *Ingham (CITED).* A. E. Ingham, Some asymptotic formulae in the theory of numbers, J. London Math. Soc. 2
   (1927), 202-208: the number `D(n) = sum_{k < n} d(k) d(n - k)` of solutions of `a1 a2 + b1 b2 = n` in positive
   integers is `(6/pi^2) sigma_1(n) ln^2 n + O(sigma_1(n) ln n ln ln n)`; the secondary terms are Estermann's
   (1930). The Euler-product reading of the leading coefficient,
   `sigma_1(n)/(zeta(2) n) = prod_p (1 - p^-2) sigma_1(p^v)/p^v`, is THM-5135 (1) with `r = 2` at every prime, not
   Ingham's; for odd `n`, without the factor `3/4` at 2 it is `(8/pi^2) sigma_1(n)/n`. F1 is a thin analogue; the
   theorem is not used (on `[20000, 20400)` `D(n)/n` and `sigma_1(n)/n` correlate at `0.979`).
3. *Jacobi (CITED).* `r_4(n) = 8 sigma_1(n) = pi^2 n * S` for odd `n`: `S` is the product of the densities of the
   four-square form at the odd primes (real density `pi^2 n`, 2-adic density 1 at odd `n`, computed by the audit;
   that their product is the count is Siegel's formula for a genus of one class, CITED). The link: `a1 a2 + u b1 b2`
   and the sum of four squares are unimodular of the same discriminant over `Z_p`, `p` odd (classification by rank
   and discriminant, CITED; class counts equal modulo 9, 27, 25, 49, 121 with `u = 2`; audit: every unit `u` at
   fourteen moduli up to 49).

Divisor by divisor: `sigma_1(n)/n = sum_{d | n} 1/d`. The divisor `d` is the content of the state of the smaller
tournament (its count and its address share the factor `d`): the states with `d | g` have density `d^-2` (THM-5128
(B2)) and put their triples on the share `1/d` of the classes (`d | N - c_sigma`), `d` times as densely: `d^-2 * d = 1/d`.

## 3. Quaternions, split quaternions, bicomplex numbers, dihedral groups

Prior work found (all credited in section 10):

- THM-2887: 16 normalized cocycles on `F_2^2`, 8 with determinant commutator, three Arf-zero dihedral classes (the
  dihedral group of order 8, `D_4`; written D8 in THM-2887) and one Arf-one `Q_8` class; the four commutative
  classes are implied (16 - 8), not named; its eq. (4), a polarization identity, is the rule of squares used below.
- THM-5133 (D1) (kind-pasteur S23; audited the same day, PASS WITH CORRECTIONS, MISTAKE-762; (A), (B), (D1)
  unchanged by that audit): an antisymmetric sign table on `Z_2^n` is a tournament; on a line, cyclic = Hamilton's
  quaternions, transitive = split quaternions with the unit of square `-1` in the middle.
- THM-024, HYP-8205, THM-4818, THM-4743: for a self-converse `T`, `Aut+-(T)` (automorphisms and anti-automorphisms)
  is a split extension of an odd group by `Z_2` (otherwise `Aut+- = Aut`; "no quaternionic option", census
  HYP-8205); anti-automorphisms have order `2 (mod 4)`; it is the dihedral group of order 10 for the regular
  5-tournament and of order 6 (`S_3`) for the 3-cycle. The dihedral group of order `2p` (written `D_{2p}` in THM-127,
  THM-131, THM-586) acts on the Hamiltonian paths of the Paley tournament.
- THM-4826, Session116 (and an old theta-series output): Jacobi's `8 sigma_1(n)`, Hurwitz's `24 sigma_1(n)`, for grids
  of the Hurwitz order and sedenion zero divisors; never for path counts. MISTAKE-506: Hamilton's algebra ramifies
  only at 2 and infinity.
- The bicomplex numbers (tessarines, `C` tensor `C`) appear nowhere: in tracked files "bicomplex" is a double complex.

Added (THM-5135 (5), all by enumeration; classical mathematics):

| units of square `-1` | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| algebra | `R^4` | split quaternions `M_2(R)` | bicomplex `C + C` | Hamilton |
| group of signed units | `C_2^3` | `D_4` (order 8) | `C_4 x C_2` | `Q_8` |
| classes | 1 | 3 | 3 | 1 |
| sign table | signed triangle | transitive tournament | signed triangle | cyclic tournament |
| norm | - | determinant = F1 seam | product of two sums of two squares | sum of four squares |

- The class is the set of units of square `-1`; the algebra is commutative iff that set is even. So quaternions
  and bicomplex numbers share two imaginary units and differ in the third; the dihedral classes have one: one
  step below the bicomplex classes, two below the quaternion class, and the only other non-commutative type.
- The non-trivial coboundary reverses the three arcs: a class is a tournament up to converse, or a signed
  triangle up to a global sign. Commutative algebras have no tournament.
- *Norms, and on which orders.* On the Lipschitz order `Z[i, j, k]` (the span of the units) the number of
  quaternions of norm `n` is `r_2` additively convolved with itself, `8 sigma_1(n)` at odd `n` (Jacobi, CITED). On
  the maximal order `Z[i] x Z[i]` of `C + C` (16 units) the number of elements of norm `n` is by definition the
  multiplicative convolution of `r_2` with itself; on the lattice `Z[1, i, j, ij]` spanned by the sign table
  (`z1 = z2 (mod 2)`, index 4, signed units `C_4 x C_2`) it is half of that at odd `n` and 0 at `n = 2 (mod 4)`; both
  vanish at `n = 3 (mod 4)`. The norm of `C + C` is the algebra norm of degree 4, that of the quaternions the
  reduced norm; the determinant is the norm on `M_2(Z)` (on the span of the units it is `x0^2 + x1^2 - x2^2 - x3^2`,
  with the seam law at odd primes only) (audit, proof-5135/a2_squares.out [B], a7_lemmas.out [O]).
- In the join calculus the ordered join multiplies counts (rank 1; for uniform factors the divisor-function law,
  a comparison: section 8) and the one-arc join adds two products (rank 2; the law `sigma_1(n)/n`, THM-5134).
- The split-quaternion norm (F1, transitive triangle, `H = 1`) and the Hamilton norm (cyclic triangle, `H = 3`)
  have the same law at every odd prime (THM-5135 (2)); they differ at 2 and at infinity. Orientation, i.e. what
  THM-5133 calls flip against spin, is invisible to the odd local densities; so is the signature of a two-port
  kernel at the primes not dividing its determinant (THM-5135 (4)(a); modulo a prime this is THM-5128 (A), KNOWN).
- Two different "D versus Q" statements should not be confused: as symmetry groups of a tournament both `D_4`
  and `Q_8` are impossible (`|Aut|` is odd); as sign tables `D_4` is `TT_3` and `Q_8` is `C_3`.

## 4. The octonions and the Fano plane

Prior work: the octonion sign table is the Paley tournament `P_7` on the Fano plane (2026-09-17 / 09-22 notes: 16
alternative orientations with `H = 189`, 112 others with `H = 171`; THM-448 and THM-450: the doubling tower;
THM-5133 (E)-(G) and Session265 the same day: Fano and `PG(3,2)` parities, the line cube). All of it concerns
sign patterns (`F_2`-valued cocycles, parities); the series concerns residues of path counts at odd primes, and at
2 F1 has the flat law (THM-5128 (E)). No link between the two is known on the tournament side. A ladder of norms:

- rank 1 (`a b`): the divisor function, no series;
- rank 2 (`a1 a2 + b1 b2`, split quaternions): `sigma_1(n)/(zeta(2) n)`; definite model Hamilton / Hurwitz,
  ramified at 2, whence `8 sigma_1(n)` and `24 sigma_1(n)` at odd `n` only (recounted for odd `n <= 400` and `n < 60`);
- rank 4 (split octonions): `sigma_3(n)/(zeta(4) n^3)`; its definite model is E8, the integral octonions (Coxeter's
  order, one of seven isomorphic maximal orders), unramified at every prime:
  `(8 pi^4/3) n^3 * sigma_3(n)/(zeta(4) n^3) = 240 sigma_3(n)` for every `n` (THM-5135 (3)). The E8 form and the rank-4
  seam have equal class counts at the odd primes by the classification (CITED) and modulo every `2^k` by lifting
  from the 136 zeros modulo 2 (PROVED, THM-5135 (3)(b), proof supplied by the audit; enumerated modulo 2, 4, 8
  here and modulo 16 by the audit); shells `240, 2160, 6720` recounted. `Z^8` is not that model (MISTAKE-484).

That the real density times the local densities is the count is Siegel's formula for a genus of one class
(CITED). For a positive definite even unimodular form in `2r` variables the real density is
`(2 pi)^r n^(r-1)/(r-1)!`, and for even `r` its product with `sigma_(r-1)(n)/(zeta(r) n^(r-1))` is the Eisenstein
coefficient `(2r/|B_r|) sigma_(r-1)(n)`: `240 sigma_3(n)` at `r = 4`. At `r = 2`, where no such form exists, the same
expression is `24 sigma_1(n)`, Hurwitz's count at odd `n` (audit, citations/c6_followups.out F5, F6).

Ranks 1, 2, 4 are those of the split composition algebras, of dimensions 2, 4, 8 (Hurwitz 1898 for sums of
squares; Albert 1942, Jacobson 1958 for arbitrary non-degenerate forms; CITED); the laws hold for every rank and
the join calculus also has rank 3 (explicit kernels with three visits, found by the audit: THM-5135 Remark 2), so
the ladder is a reading, not a constraint. Two rank-4 forms of determinant `-1` occur in the join calculus (THM-5135
(4)(c)): the K-form (Session199; recorded in THM-4916 (1d)) and the corner of Session264's fusion law. THM-4916
computes the class law of the values of one word with at most two root changes (letter-count model; one free
state per value): factor `1 + (3p - 2)(p - 1)/(2p^3)` at the class 0 modulo an odd prime `p`, neither the value
`1 + p^-3 - p^-4` of two independent uniform states nor the value 1 of one uniform state under a non-zero linear
functional: its states are confined.
The fusion corner is unimodular in its displayed coordinates, which its states do not fill: on the lattice of the
states it has determinant `-4`, so it is unimodular at the odd primes only, and at 2 no family has the uniform law
(path counts are odd); at the odd primes nothing is claimed. Whether a tournament family carries the rank-4 law
at the odd primes is HYP-9666 (A), OPEN; the two-apex grammar `G(2,0)` of HYP-9660 (D) has class factors within
`0.0004` of 1 at `p = 5` (all words, `n = 16`), which excludes the value `1 + 5^-3 - 5^-4 = 1.0064` there. (The 240
labelled Paley tournaments against the 240 roots: section 9.)

## 5. Tao's toy example, and every base

The example was not attached to the message. It is taken to be, as in Session266, the puzzle set by Terence Tao
for Alex Bellos's column ("Can you solve it? A magical melody from the Mozart of math", The Guardian, 26 October
2026, section "The magic money plant"; confirmed by the audit on syndicated copies), recorded in
`session266_fiber_radix_20261111/radix/PROOF.md` (committed while this session was running): `C_t = a t^2 + b t + c`,
and `100 a + 10 b + c = 427` with digits gives `(4, 2, 7)`. Session266 proves the all-base version (injective iff the
day is at least the base), the nilpotent algebra `Q[nu]/(nu^3)` behind the recursion and "carries are the kernel
`(x - B)`". What this session adds:

- *The F1 joins are such quadratics.* The odd-cycle packing polynomial of `J(h1, h2, s)` is
  `1 + alpha_1 z + alpha_2 z^2`, `alpha_1 = psi1 + psi2 + (2^s y1 y2 - [s = 0])/2`, `alpha_2 = psi1 psi2` (THM-5070 (2) in
  these coordinates; 314 joins recounted, every_base.out (2)): an instance of Tao's quadratic, with constant term
  1. Its value at 2 is the path count, its value at `b` is `N_b = x_b(h1) x_b(h2) + b (2^s y1 y2 - [s = 0])/2`.
- *The digit bound.* For `b > max(alpha_1, alpha_2)` the base-`b` numeral of `N_b` is `(alpha_2)(alpha_1)(1)`: the
  count displays the packing numbers (the block lemma of THM-4940 is the case `alpha_2 = 0`). In base 2 this
  almost never holds and carries collide: 17 is the count of joins with `(alpha_1, alpha_2) = (6, 1)` and `(8, 0)`.
  Of the 10,400 joins with `h1 <= h2 < 64`, `s <= 4` the share with digits at most `b - 1` is
  `0.0002, 0.0004, 0.0014, 0.022` for `b = 2, 3, 10, 100` (bridges_exact.out S11).
- *The local law in every base* (THM-5134 (G), a PROVED limit law). At every odd modulus prime to `b` the base-`b`
  counts of F1 have the same class law, with the distinguished classes `0` (sheet 1) and `-b/2` (sheet 0);
  `N_b = 1 (mod b)`; the series runs over the primes not dividing `2b` (`(9/pi^2) sigma_1(n)/n` for `b = 3, 6`). The
  rates are those of the binary digit chain and do not depend on `b`: 7 is the slow prime in every base prime to 7.
  NUMERICAL, complete censuses (the draft's were truncated, section 11): sheet 1 (every_base.out), `N_b <= 2^24`,
  `b = 2, 3, 4, 5, 6, 10` (57,818,155 triples at `b = 2`): remainder of the first-order law at most `0.133`, at 7;
  sheet 0, `N_b <= 2^22`, `b = 2, 3, 5, 10` (10,532,579 triples at `b = 2`): every triple with a `p`-divisible smaller
  state has `2 N_b + b = 0 (mod p)` (0 violations), remainder at most `0.167`, at 7 (bridges_exact.out S10).
- *The second digit* (THM-5134 (G); found by the audit). The second base-`b` digit of `N_b` is `alpha_1 (mod b)`. For
  odd `b` it is equidistributed by value on each sheet, with a power saving (PROVED, a corollary of THM-5128
  (B4)(iv): modulo `b` it is an affine form in the residue vector of the larger state with a unit coefficient);
  for `b = 2` this is THM-5128 (E) with `j = 2`; for even `b > 2` it is HYP-9665 (C), OPEN. NUMERICAL (audit, sheet
  0 / sheet 1), `b = 3`: deviation `0.131 / 0.127`, `0.071 / 0.073`, `0.041 / 0.040`, `0.024 / 0.023` at
  `2^16, 2^20, 2^24, 2^28`, a factor `0.867` per binary digit, the rate `sqrt(3)/2` of the digit chain at 3.
- *Balanced states and early boosts.* Apart from `h = 0`, balanced states `x_b(h) = y(h)` exist only for `b = 2`
  (the fans) and `b = 3` (`h = 3, 45, 117, 309, ...`; with `h = 0`, 12 below `2^28` and 13 below `2^30` by the audit;
  finiteness OPEN). An early boost needs only a small divisible state, not a balanced one: `h = 3` (`psi = 2`,
  `y = 7`) is 7-divisible for every `b = 3 (mod 7)` (`x_b = y = 7` in base 3, `x_b = 21` in base 10; with `F` the
  share of the triples whose smaller state is `p`-divisible, `F p^2 = 3.37`, `3.28` on sheet 1 at `2^24`), and
  `h = 6` (`x_b = 26 = 2y`, `y = 13`; not balanced) boosts 13 in base 5 (`F p^2 = 3.38`).
- Two radices are in play and are independent: the base `b` of the numeral and the binary radix of the address
  `h` (visible in THM-5085 (2); Session266 separates the evaluation base from the digit bound).

## 6. Multiplication tables modulo odd `m` and the arcs of `3n + 1`, `3n - 1`

Prior work (credited): THM-4792 (the centre labels `+-1/2` of the table are the fixed points of `3x + 1`, `3x - 1`;
alternation iff `m | 2^t + 1`, Moree); Session255 (every arc of either map is the centre block of one table,
`m = 4x + eps`); the central-modular packet (`m -> 3m + 2 eps`, `m -> (m + eps)/2`); modular centres (cycle labels
`{5, 17, 9}`, `{3, 7}`, `{19, 55, 27, 79, 39}`); codex WEIGHTED_ATLAS and PUNCTURED_ATLAS (the product law, the triangle
`J(x, eps) = (3x + eps, -eps)`, 7 and 21 as the arcs `2 -> 1` of `3n - 1` and `5 -> 16` of `3n + 1`); Session266
collatz/PROOF.md sections 2 and 8 (codex, 11:33 the same day, before this session's scripts for this part, 11:43
and 12:04): in its notation `h = 4x + s` (there `h` is the count and `s = +-1` the sign called `eps` here), the
count `h` sent to the signed raw edge `x -> C_+(x)` (`s = +1`) or `-x -> -C_-(x)` (`s = -1`), `H = 1` a boundary
symbol; the code product; the code dynamics by `h mod 8`; the cycles `5 -> 17 -> 9 -> 5` and `3 -> 7 -> 3`, the latter
cut at 7; oldest, the star law `a * b = a + b + 2ab` on `(H - 1)/2`, `1 * n = 3n + 1` (paired_prime_collatz_20261003).

One coordinate puts the two Collatz sides together (a restatement of Session266's signed source, `n = eps x`;
bridges_exact.out S7). Every odd `m >= 1` is `|4n + 1|` for exactly one integer `n = (chi_4(m) m - 1)/4`; `m` labels
the arc of `3n + 1` on `Z` leaving `n`. Then: an odd step is `m* -> 3 m* + 2` and an even step `m* -> (m* + 1)/2`
for `m* = 4n + 1`; multiplying labels is `n * k = 4nk + n + k`; tripling a count is `n -> -3n - 1`. The Collatz
side is the sign of `n`, that is `chi_4(m)`, not the F1 sheet (`N mod 4` is flat on each F1 sheet, THM-5128 (E)).
Added:

- *The density is the additive square of the table* (section 2, item 1).
- *The favoured classes in arc coordinates.* For a count `N = |4n + 1|` and an odd prime `p`: `p | N` iff
  `n = -1/4 (mod p)`; `N = -1 (mod p)` iff `n = -1/2` (for `n >= 0`) or `n = 0` (for `n < 0`). So sheet 1 of F1
  favours the arcs leaving the class `-1/4`, and sheet 0 the arcs leaving the class of the fixed point `-1/2` of
  `3n + 1` (positive side) or the class 0 (negative side). The argument of the series of sheet 0 is `2n + 1` for
  `n >= 0` and the odd part of `n` for `n < 0`. A relabelling; no dynamics.
- *Fan primes are the alternation primes.* KNOWN (THM-5128 (F), THM-5126 (A); S91 note, section 3): a prime has a
  balanced `p`-divisible state iff `p | 2^t + 1` for some `t` iff `ord_p(2)` is even, the first such state being the
  fan `h = 2^(t-1)`, `t = ord_p(2)/2`. Added: for a prime modulus this is the alternation condition of THM-4792 (c)
  (Moree); among `p < 200` the fan is the first `p`-divisible state exactly for `p = 3, 5, 11, 13, 17, 43, 113` (the
  S91 audit had `3, 5, 11, 13, 17, 43` and `257`, `65537`; among `200 < p < 400` they are 241, 257, 331, by the
  audit).
- *Lines* (small_tournaments_collatz_lines.out). `c_-(N)` = odd part of `N - 1` (line `1 + c 2^s`, THM-5126),
  `c_+(N)` = odd part of `N + 1` (anti-line, THM-5130). By definition `c_-(3N) = T_-(N)`, `c_+(3N) = T_+(N)`, the
  Syracuse maps of `3n - 1`, `3n + 1`: the ordered join with a triangle (`N -> 3N`) puts a count on the line
  labelled by its own image under the Syracuse map of `3n - 1` and on the anti-line labelled by its image under
  that of `3n + 1` (the remark `1 * n = 3n + 1` of 2026-10-03 with odd parts; cf. codex's `J`). It is not a map on
  line labels: `c_-(3N) = T_+(c_-(N))` only for `N = 3 (mod 4)`, `c_+(3N) = T_-(c_+(N))` only for `N = 1 (mod 4)`.

## 7. Negative against positive cycles, and the four tournaments on four vertices

Prior work: THM-4472 (the two converse classes with `H = 3` are the forward map and the inverse tree of one pair:
converse is time reversal composed with negation; their reading as the two Collatz sides, called sheets there,
is REFUTED there); the conjugacy `x -> -x` between `3n + 1` on the negatives and `3n - 1` on the positives, under
which every residue-level datum is invariant (THM-4474 (E), THM-4521, the functional-uniqueness note);
collatz_five_mirrors_20260929, Proposition 5 (word reversal fixes the cycles of `-1` and `{-5, -7}`, a rotation,
and sends the cycle of `-17` to the rational point `-13801/139`) and THM-4849 (circulant tournaments are rational
cycles, on which reversal is the converse); THM-5131 (the 48 labelled self-converse 4-tournaments have `H = 1` or
`5`, the 16 others `H = 3`); Session266 collatz/PROOF.md sections 3, 4, 8 (`D = 2^L - 3^r = 1, -1, -1, -139` for the
four known non-zero cycles, the positive / negative distinction being the sign of `D`, which residues forget; the
four 4-vertex classes with `H = 1 + 2 c_3` and their converse types; the weights at order 4).

Exact statements (bridges_exact.out S7, small_tournaments_collatz_lines.out):

- The five known cycles of `3n + 1` on `Z` (`{0}`, `{1, 4, 2}`, `{-1, -2}`, the cycles of `-5` and of `-17`) have 29 arcs
  with 29 distinct labels; exactly one label, 7 (the arc `-2 -> -1`), is not a path count. The positive trivial
  cycle has labels `5, 17, 9`, all attained; the negative trivial cycle has labels `3, 7` and is cut at 7 (KNOWN:
  Session266 section 8, PUNCTURED_ATLAS, as is the arc `5 -> 16` of the other known non-value, 21); the trunk arcs
  `(4^k - 1)/3 -> 4^k` have labels `(4^(k+1) - 1)/3 = 5, 21, 85, ...`. New: the count over all five known cycles,
  including the 18 labels of the cycle of `-17`. Other cycles and other non-values are not excluded (both OPEN).
- KNOWN (Session266 sections 4, 8; THM-791 (a), quoted by THM-4930 and THM-5131), restated in the coordinate `n`:
  the path counts `3, 5` of the 4-vertex tournaments are the labels of the arcs `-1 -> -2`, `1 -> 4`; at order 4,
  self-converse iff `H = 1 (mod 4)` iff `n >= 0` iff the number of 3-cycles is even, the converse pair (`H = 3`)
  being the only class on the negative side (tested in small_tournaments_collatz_lines.out); for even order in
  general only one direction holds: self-converse implies `H = 1 (mod 4)`.
- Added: `H = 1` read as the arc `0 -> 0` (so `1, 3, 5` are the arcs leaving `0, -1, 1`, one in each of the three
  cycles through them), and the label sets at orders 5 and 6: the integers covered are `[-4, 3]` without `-2` at
  order 5 and `[-11, 11]` without `-10, -9, -2, 5` at order 6.

What is not there: a map from cycles to tournaments. "Three known negative cycles against the one known positive
cycle" and "four tournaments on four vertices" agree only as counts: COINCIDENCE (also Session266's verdict). No
residue-level statement distinguishes the two Collatz sides (THM-4474 (E)); the numbers of known cycles on the
two sides are not local data, and the local densities of this session are of that blind kind: a reading, not a
theorem about cycles. The denominators `D = 2^L - 3^r` of the three small cycles are units (`4 - 3`, `2 - 3`,
`8 - 9`; Session266 section 3), of the Catalan / Pillai type under which kind-pasteur files the sign-flip criterion
`2^delta - 3^(delta-1) = 1` of THM-5130 (B) (solutions `delta = 1, 2`: `2 - 1`, `4 - 3`; kps_s22 note, sections 3, 10).
Only `4 - 3 = 1` is common to both: ANALOGY; no map between cycles and gadgets is claimed.

## 8. The pasted note: strong components as primes

- Confirmed by enumeration: at order 4 the arc-usage profiles (number of Hamiltonian paths through each arc) are
  `{1,1,1}`, `{2,2,2,1,1,1}` twice (a converse pair), `{4,3,3,3,1,1}`; the classes satisfy `B(x) = 1/(1 - S(x))`,
  `S(x) = x + x^3 + x^4 + 6x^5 + 35x^6 + ...` the series of strong classes (classical); enumeration confirms
  `1, 1, 2, 4, 12` at orders 1 to 5 here and 56 at order 6 by the audit (`S` truncated at `6x^5`, as printed in
  small_tournaments_collatz_lines.out, is exact to order 5 only: it would give 21 at order 6).
- The closing question is KNOWN (Session224 birth ledger, census/birth_fibers.csv rows `H = 9, 11, 13, 15`, and
  THM-4922 (3); 8 self-converse classes of 12 is A002785, also in Session224 and HYP-8205): the six strong
  classes of order 5 have `H = 9, 9, 11, 13, 15, 15` and are all self-converse; recounted here with the
  automorphism orders `1, 1, 1, 1, 3, 5`. The other 4 of the 12 classes are two converse pairs (non-strong).
- Prior work: the free monoid on strong tournaments with `H` multiplicative is THM-1862 / THM-1936 (classical:
  Moon; used in THM-1370 (1)), its series `T = 1/(1 - S)` is in THM-4816 / THM-4817; arc usage is the
  consecutive-position count of THM-050 and Session265's `W(u, v)` (it determines the strong components and their
  order up to reversal); Session265 also corrects two statements of the note (a general vertex inflation does not
  multiply `H`; the total weight of a boundary, not each boundary arc, is `H`) and shows the limit of the analogy
  (63 is irreducible among path counts); the product law on table codes is codex's. The mac-mini note
  character_of_two_records_20261111 (same day, pre-audit) types the same pasted text with these credits.
- In the language of section 4, "strong components as primes" is the rank-1 seam, with the divisor-function law
  `(1 - 1/p)(v + 1)` for uniform factors; no equidistribution of strong counts is known: a comparison, not a theorem.

## 9. The bridges, typed

EXACT = a proved or enumerated identity; MECHANISM = a common cause with a stated loss; ANALOGY = same shape, no
map; COINCIDENCE = equal numbers only. S7, S8, S10, S11 are sections of bridges_exact.out.

| bridge | type | evidence |
|---|---|---|
| F1 seam = four-square form at every odd prime | EXACT | THM-5134 (F); class counts modulo 9, 27, 25, 49, 121 |
| `S(N) = r_4(n)/(pi^2 n)` | EXACT (Jacobi and Siegel CITED) | `r_4 = 8 sigma_1` for odd `n <= 400`; the densities computed by the audit |
| `S` = Ingham's leading coefficient without the prime 2 | EXACT (CITED); F1 as a thin Ingham problem: ANALOGY | section 2 |
| `L_m` = additive square of the multiplication table | EXACT | eleven moduli; audit: every `m <= 210` |
| rank-4 law over all primes = E8 = integral octonions | EXACT for uniform states (odd primes: classification CITED; at 2: PROVED by lifting, audit) | THM-5135 (3) |
| a tournament family with the rank-4 law at the odd primes | OPEN (at 2 the law fails in every family) | HYP-9666 (A); THM-4916 (one free state, confined) is not one; `G(2,0)` is excluded at 5 |
| sign tables against norm forms | EXACT for the forms (classical): the norm form sees only the units of square `-1`; on one line this is the orientation, on the Fano plane it is not (128 orientations, one form); definite and split forms agree at odd primes, for E8 also at 2 | THM-5135 Remark 1; no statement about path counts |
| `240 = 7!/21` labelled Paley tournaments, 240 roots of E8 | COINCIDENCE (both `16 x 15`) | no equivariant bijection known (fano_code, section 4) |
| eight extensions = subsets of imaginary units; quaternion / bicomplex / dihedral | EXACT (classical) | THM-5135 (5) |
| quaternion norm additive, norm of `Z[i]^2` multiplicative in `r_2` | EXACT (definitions on `Z[i, j, k]` and `Z[i]^2`; `8 sigma_1(n)` is Jacobi, CITED, checked for odd `n <= 400`) | section 3 |
| ordered join = rank 1, one-arc join = rank 2 | MECHANISM (laws for uniform states; rank 1 is a comparison only) | THM-5135 (1) |
| `Aut+-(T)` never contains `D_4` or `Q_8` | KNOWN | THM-024, THM-4818, HYP-8205 |
| `Out(Q_8) = S_3 = Aut+-(C_3)` | EXACT (classical; not recomputed) | THM-2887 |
| the Collatz side of a table label is `chi_4` of the count | EXACT (restatement); KNOWN | Session266 section 8; S7 |
| Tao's quadratic: the packing polynomial of a join is an instance; digit bound = no carry | EXACT | THM-5070 (2), every_base.out (2); S11; THM-4940 |
| base-`b` law | PROVED limit law (THM-5134 (G)); EXACT parts: the formula for `N_b` (brute force) and type `d` implies `d` divides the class; NUMERICAL: complete censuses, remainders up to `0.17` | every_base.out, S10 |
| second digit in odd bases | PROVED (THM-5134 (G), found by the audit); even `b > 2` OPEN (HYP-9665 (C)) | audit censuses to `2^28` |
| tripling = `n -> -3n - 1`, products = `4nk + n + k` | EXACT, KNOWN: star law `a * b = a + b + 2ab` on `(H - 1)/2` with `1 * n = 3n + 1` (paired_prime_collatz_20261003, sections 2, 3); `(x, eps)` form Session255, codex, Session266 | S7 |
| one known non-value among the 29 arc labels of the known cycles | EXACT (uses THM-1370 and THM-4875); the cut at 7 KNOWN | S7; Session266 section 8 |
| `1, 3, 5` = arcs at `0, -1, 1`; self-converse iff `n >= 0` at order 4 | EXACT (small); KNOWN except the arc `0 -> 0` | S7, small_tournaments_collatz_lines.out; Session266 sections 4, 8 |
| three known negative cycles and one known positive against four 4-tournaments | COINCIDENCE | THM-4472; THM-4474 (E); Session266 sections 3, 4 |
| cycle denominators `2^L - 3^r` against the sign-flip criterion of THM-5130 (B) | ANALOGY (Catalan / Pillai; one common equation) | Session266 section 3; kps_s22 note |
| fan primes = alternation primes of THM-4792 (c) | EXACT (arithmetic); KNOWN except the cross-reference | S8; THM-5128 (F) |
| `sigma_1(n)/n` against the tournament abundancy `Pi_2/H` of kind-pasteur S22 | ANALOGY | `a(C_3) = 4`, `sigma_1(3)/3 = 4/3` |
| Session265's prime-power lifts, exceptions 5 and 761 | unrelated (a name clash) | THM-5134 Remark 6 |
| 7 at the minimum of the series of sheet 0 | COINCIDENCE (true of every `2^k - 1`) | THM-5134 (E) |

## 10. Credits and prior work

- THM-5128 with its audit (the factor for uniform states, the series, the enumeration at 9, 25, 27, 49; (A), (B2),
  (B4)(iv), (C1), (C2), (D2), (D3), (E), (F) as cited in sections 1, 2, 3, 5, 6), the S91 note, MISTAKE-761.
- Session206 (weighted parameter capacity), Session193 modular_sync, Session 25: the fixed-source law and the
  digit chain. Session223, THM-4916: valuation-resolved factors for root changes; Session199: the K-form (recorded
  in THM-4916 (1d)).
- kind-pasteur: THM-5133 (audited 2026-11-11, PASS WITH CORRECTIONS, MISTAKE-762) and HYP-9664 ((C) REFUTED for
  `n >= 5` by that audit); THM-5130, HYP-9662 and the kps_s22 note; THM-5070 (2), THM-5085 (2), THM-5126 (fans,
  lines), the modular-centres packet; from June, THM-447, THM-448, THM-450 (the doubling tower).
- Codex lanes: Session255, Session260 (the norm form `x x' + 2 y y'`), Session261 (the kernel with determinant 329),
  Session263, Session264 (fusion law), Session265 (`W(u, v)`, the Fano cube, 63), Session266 (Tao's puzzle, radix;
  collatz/PROOF.md sections 2, 3, 4, 8: signed-edge coordinate `h = 4x + s` in its notation, code product and
  dynamics, the cut at 7, `D = 2^L - 3^r` for the four known non-zero cycles, the four 4-vertex classes), the
  central-modular, odd-tables and signed-tables packets (WEIGHTED_ATLAS, PUNCTURED_ATLAS, PHASE_CLOCKS).
- Other packets: Session224 (birth ledger: the six strong classes of order 5, the 8 self-converse classes of 12;
  with THM-4922 (3), A002785); paired_prime_collatz_20261003 (star law); overnight_research_20261029/TWO_PORTS.md
  P5; collatz_five_mirrors_20260929 (Proposition 5). mac-mini: THM-5131, THM-5132, HYP-9663 (the mod-4 race,
  self-converse tournaments), the note character_of_two_records_20261111. claudebox: THM-480 (read; not used).
- Older canon: THM-2887 (with its eq. (4) and `Out(Q_8) = S_3`), THM-024, HYP-8205, THM-4818, THM-4743, THM-127,
  THM-131, THM-586, THM-4472, THM-4474, THM-4521, THM-4792, THM-4826, THM-4849, THM-791 (a), THM-4930, THM-4940,
  THM-4756, THM-1370, THM-1862, THM-1936, THM-4816, THM-4817, THM-050, THM-4875, Session112 (norm signs),
  Session116, the 2026-09-17 and 2026-09-22 Paley / Fano / octonion notes, MISTAKE-484, MISTAKE-506.
- External (CITED): Jacobi (1829/1834; `r_4(n) = 8 sum_{d | n, 4 not | d} d`, `r_8(n) = 16 sum_{d | n} (-1)^(n+d) d^3`);
  Hurwitz, Vorlesungen ueber die Zahlentheorie der Quaternionen (1919) (24 times the sum of the odd divisors) and
  Nachr. Goettingen 1898 (`1, 2, 4, 8`; Albert 1942, Jacobson 1958); Ingham 1927 (section 2); Estermann, Proc.
  London Math. Soc. (2) 31 (1930), 123-133; Siegel's formula for a genus of one class; Coxeter, Integral Cayley
  numbers, Duke Math. J. 13 (1946), 561-578; OEIS A004009 (`Theta_E8 = 1 + 240 sum sigma_3(n) q^(2n)`), A002785;
  Gronwall, Trans. AMS 14 (1913), 113-122; Moree, Acta Arith. 80 (1997), 197-212; O'Meara sec. 92-93,
  Conway-Sloane ch. 15, Cassels ch. 8 (unimodular forms over `Z_p`, `p` odd and `p = 2`); Serre, A Course in
  Arithmetic, ch. V (even unimodular definite lattices have rank divisible by 8); Moon (strong components).

THM-5133 (G), as audited, identifies `CD_4` with `T_15` of THM-448 (isomorphism found by its audit) and records that
`CD_5` is not `T_31`. Recomputed here independently (cd_tower_iso.py: explicit isomorphisms for `n = 2, 3, 4`; at
`n = 5` neither `T_31` nor its converse, both with 21 automorphisms and the same census of 4-subsets). The split-octonion
and split-sedenion values 45 and 24,453,597 also occur in THM-450 (7) (cross-reference only, not tested).

## 11. Audit

Independent audit, 2026-11-11, by eight auditors (proof-5134-AE, proof-5134-FH, proof-5135, citations,
numerics-prereg, bridges-note, prior-art, hypotheses), each with its own code: every verdict PASS WITH
CORRECTIONS; 111 findings (6 MAJOR, 64 MINOR, 41 NOTE), all applied; MISTAKE-763. Report:
`05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/AUDIT_REPORT.md`. The majors:

- the base-`b` census of sheet 0 (bridges_exact.py S10) covered 27 to 36 % of the triples and was described as
  complete (two findings; that of sheet 1 in every_base.py was incomplete too, a minor); fixed and regenerated;
- no credit to the Collatz lane of Session266, which answered the same part of the directive first (sections 6, 7);
- the six strong classes of order 5 being self-converse was given as new; it is KNOWN from Session224 (section 8);
- the pre-registration was only partly blind, which neither the registration nor the draft said (section 1);
- the falsifier of HYP-9665 (A) was already met by exact data: it tested the oscillating remainder at the odd
  part of the modulus, not the conditional flatness; it is replaced.

What the audit added: the second digit in odd bases (PROVED, THM-5134 (G)); the Smith-form law of two-port
kernels at every prime (PROVED, THM-5135 (4)(b); it was HYP-9666 (C)); E8 at 2 modulo every `2^k` by lifting
(PROVED, THM-5135 (3)(b)); conditional-flatness data for even moduli (HYP-9665 (A): flat to `0.0002` at 12 at
`2^28`); the first-order constants `A_(m,d)` and the mass constants in base `b`; a census at `2^30`; explicit kernels
of rank 3; the series of THM-5134 (G) restricted to the primes not dividing `2b`.

The corrections to this note were then checked against the outputs and recomputed with independent code
(audit/verify-README/: the exact statements of sections 2 to 8; the censuses of sections 1 and 5, at `2^28` and
`2^30` included; the multiplicities of section 1).

## 12. Reproduce; next obligations

    cd 04-computation/experiments/s92_prime_powers_20261111
    python3 odd_modulus_census.py cum 24            # exact remainder of the first-order law (2^28: 35 s, 8 threads)
    python3 odd_modulus_census.py eval 28 27 prereg_window_28_27.json
    python3 seam_rank_exact.py; python3 bridges_exact.py; python3 every_base.py; python3 cd_tower_iso.py
    python3 small_tournaments_collatz_lines.py; python3 series_vs_counts.py 25 24

every_base.py is the receipt of THM-5134 (G), cd_tower_iso.py of the last paragraph of section 10; `cum 28` takes
35 s in the stored output and under two minutes on a loaded machine (audit). Next: (1) write out even moduli and
the second digit in even bases `b > 2` (HYP-9665 (A), (C): every ingredient is available, audit); (2) uniformity
in the modulus (it needs at least HYP-9660 (B) with a norm bound), which is necessary, not sufficient, for a
statement about single values: a second-moment input is needed as well; (3) a rooted family with rank-4 states
equidistributed at the odd primes, or a proof that none of the known ones is (HYP-9666 (A)).
