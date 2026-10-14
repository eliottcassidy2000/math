# Collatz continued fractions and irrationality exponents

2026-10-14. **PROVED scoped encodings and equivalences; FINITE-EXACT controls;
CITED supplied analytic inputs; universal Collatz remains OPEN.**

The positive integers can be embedded faithfully into continued fractions,
with one Collatz step becoming one continued-fraction shift. Every number in
this natural encoding has irrationality exponent2. A different, explicitly
engineered encoding has exponent2 exactly when its source reaches1. The
difference between these constructions identifies the missing proof task:
an exponent must retain the desired terminal predicate before an exponent
bound can prove termination.

## Inheritance and the objects retained

The closest mechanisms are the actual-source certificates in
[native tail recursion](collatz_native_tail_recursion_20261013c.md), the
exact cycle/continued-fraction gates in
[the Pillai note](collatz_mod6_20260921_pillai_convergents_cycle_gates.md),
and the two-place distinction in
[the earlier zeta5 exploration](collatz_posets_dags_zeta5_20260927.md).
The canonical hostile is the signed fixed point−1: arbitrarily long
positive Mersenne prefixes shadow it. The corrected near miss is a
coefficient or frequency condition promoted to an actual integer cycle;
the affine carry and source guard are still required. The least-used
sidecar explored here is the complete score of a verifier-specific
statistical model, supplied in
[the information note](information_observers_collatz_20261014a.md).

The concept board consists of the actual source, its parity address,
continued-fraction denominator growth, terminal status, and retained
statistical information. The anchor is source-specific Collatz closure;
the niche is exact exponent engineering; the wildcard is the quadratic
field of a parity continued fraction. The applied research cards are
"Separate local support from bounded-height coverage" and
"Keep observer, recurrence class, and finite head" from META-PATTERNS.

Throughout use the shortcut map

    T(n)=n/2                         if n is even,
    T(n)=(3n+1)/2                   if n is odd.

This is the version with1↔2. The odd-only Syracuse map instead has1→1.
Connectivity of the entire positive graph is part of the conjecture, not
an assumption of the construction.

## 1 The unique two-cycle is elementary

If n=2m and T²(n)=n, then an even m gives n=0; an odd m gives
(3m+1)/2=2m, hence n=2. If n is odd, put m=(3n+1)/2. An even m gives
(3n+1)/4=n, hence n=1; an odd m gives (9n+5)/4=n, hence n=−1.
Thus the only positive solutions of T²(n)=n are1,2. There are no positive
fixed points. This proves uniqueness of the bidirectional edge, not
exclusion of longer cycles or unbounded forward paths.

## 2 A lossless encoding on which every exponent equals two

Write b_j(n)=T^j(n) mod2 and define

    beta_n=[0;1+b_0(n),1+b_1(n),1+b_2(n),...].

All digits are1 or2, so beta_n is an irrational real in(0,1). Let
G(x)=1/x−floor(1/x) be the Gauss map. Deleting the first digit proves

    G(beta_n)=beta_(T(n)).

This is a genuine dynamical semiconjugacy, and it is injective on positive
integers. To prove injectivity without a probabilistic assumption, fix a
length-L parity word w. On that word

    T^L(n)=(P_w n+B_w)/Q_w,
    P_w=3^(number of ones), Q_w=2^L.

Starting at(P,B,Q)=(1,0,1), an even letter sends this to(P,B,2Q); an odd
letter sends it to(3P,3B+Q,2Q). The unique legal source residue is

    n = -B_w P_w^(-1) mod Q_w.

The residue is sufficient as well as necessary: reduce the final
integrality congruence to each prefix; the odd remaining multiplier is
invertible modulo that prefix's power of2. Equivalently, each legal
prefix residue has exactly two lifts modulo its next power of2, giving
the two possible next parities. Equal infinite words therefore force
equal residues modulo every2^L, hence equal positive integers.

Decoding a bounded source needs its height bound: L bits recover n once
0≤n<2^L is known. Injectivity of an infinite code does not supply a uniform
finite stopping rule for reconstructing an unbounded source.

For any infinite simple continued fraction alpha=[a0;a1,...], with
convergent denominators q_j, the classical formula is

    mu(alpha)=2+limsup log(a_(j+1))/log(q_j).

Since beta's digits are bounded, **mu(beta_n)=2 for every positive n**,
without using Collatz termination. The full address distinguishes every
source; the single scalar exponent distinguishes none. Moreover, any
exact Gauss-shift encoding has mu constant along each orbit, because
deleting a finite prefix is an invertible rational fractional-linear
change. An assignment mu(beta_n)=n cannot obey this intertwining law.

The direct-state code also has exponent2, even if an orbit hypothetically
escapes. The stronger local-growth theorem and its explicit finite
denominator bound are proved in
[continued-fraction state codes](continued_fraction_collatz_codes_20261014a.md).
This is not restricted to bounded digits.

## 3 Periodicity retains information that the exponent loses

The root codes are explicit:

    beta_1=[0;overline(2,1)]=(sqrt(3)-1)/2,
    beta_2=[0;overline(1,2)]=sqrt(3)-1.

By the classical Lagrange theorem and injectivity above,

    beta_n is quadratic irrational
      iff its continued fraction is eventually periodic
      iff the actual T-orbit of n is eventually periodic.

The last equivalence uses G^(r+s)(beta_n)=G^r(beta_n), then injectivity
at T^(r+s)(n) and T^r(n). Quadraticity alone does not identify which cycle.
The exact root criterion is stronger:

    n reaches1 iff beta_n has a finite binary-digit prefix
                    followed by overline(2,1).

An odd shift of the alternating tail can be incorporated into the prefix.
Thus the correct target is membership in one specified eventual-tail
class, not merely mu=2 or even arbitrary eventual periodicity.

The inverse branches make the arithmetic guard visible. If h_a(x)=1/(a+x),
then h_1(beta_m)=beta_(2m) always, while h_2(beta_m) is the code of the
positive integer(2m−1)/3 only when m=2 mod3. The symbolic Cantor set permits
both branches everywhere; the integer graph retains the congruence guard.

The same prefix residues define a unique2-adic source for every infinite
binary word. For example all ones give the signed source−1 and
[0;overline(2)]=sqrt(2)−1. For every L, the positive integer2^L−1 has the
same first L parities, since

    T^j(2^L−1)=3^j 2^(L−j)−1                 (0≤j≤L).

This supplies an exact hostile to any argument that identifies all lawful
symbolic paths with positive-integer trajectories. A growing supply of
finite positive witnesses does not make their2-adic limit positive.

**Wildcard, decisive field hostile.** Among all digit words over{1,2} of
lengths1..16, the periodic codes in Q(sqrt3) were exactly the alternating
words at even lengths. This pattern fails at length18. The primitive word
`111111121121211122` has CF matrix

    [[29872,12553],[18471,7762]],
    trace²−4det = 3*21728².

Its periodic real is therefore in Q(sqrt3), with polynomial
`12553 x²+22110 x−18471=0`, but its tail is nonalternating. The associated
parity word has(P,Q,B)=(243,262144,402560), giving the actual rational
fixed source `402560/261901`. Exact replay produces18 distinct positive
rational states of odd denominator261901. This is a primitive rational
cycle, not an integer counterexample. A quadratic field label still loses
the specified tail class; the source integrality guard cannot be discarded.
The separate implication beta_n∈Q(sqrt3)⇒root restricted to positive
integer sources is not decided by this rational hostile.

## 4 A computable exponent that detects termination exactly

There is also a positive answer to the proposed engineering program.
Fix an integer M≥3. Define z_n=[n;2,a2,a3,...] as follows. At stage j≥1,
compute the actual orbit through time2^j, retaining any earlier visit to1.
Using the already constructed denominator q_j, set

    a_(j+1)=1                    if1 has been visited,
    a_(j+1)=q_j^(M−2)            otherwise.

Each digit is computable from finite work; no termination oracle is used.
The Fibonacci lower bound for q_j gives a uniform effective convergence
bound for the convergents, so z_n is a uniformly computable real as well.
The doubling simulation budget is only a practical compression of the
checking schedule. Every finite hitting time is eventually inspected.

If n reaches1, the code has an all-ones tail and mu(z_n)=2. If n never
reaches1, every subsequent digit is q_j^(M−2), so the displayed limsup
formula gives mu(z_n)=M exactly. Consequently

    Collatz for every positive source
      iff mu(z_n)=2 for every n
      iff mu(z_n)<M for every n.

This is a proved reformulation, not an independently established exponent
bound. The gate stores the terminal predicate by design. It works just
as well for any computable process with a detectable terminal state;
therefore one cannot expect an unrestricted algorithm that determines
this exponent from every computable continued-fraction generator.

It also explains what a successful analytic bridge must provide: an
independent upper bound below M for these particular codes, or for a
similarly faithful code whose terminal predicate has been proved. Showing
mu=2 for the natural state/parity code cannot substitute for this step.
An explicit representation of every z_n as an invertible rational
fractional-linear transform of a logarithm covered by the supplied theorem
would suffice, but no such representation is established.

The proposed recipe a_(j+1)=a_j², a1≥2, gives mu=3, not4: the sum of prior
logarithms is asymptotic to2log(a_j). To prescribe mu=M, growth relative to
the entire denominator q_j is the correct control. The separate recipe
a_j=2^(j!) does give infinite exponent. Proofs and finite hostile controls
are in the state-code note. Finite initial peaks never determine mu;
changing finitely many continued-fraction digits preserves it.

## 5 The logarithmic cycle target and the supplied papers

The supplied [logarithm preprint](https://arxiv.org/html/2610.10192v2),
Theorem1 and Corollaries6–7, asserts mu(log r)=2 for positive rational
r≠1 and its stated nonzero rational combinations/rational-power extensions.
It is accepted as a cited input here. It gives neither a uniform constant
in r nor the conclusion for the quotient log3/log2. In particular,
log((k+1)/k) tends to0 although each such logarithm has exponent2 under the
theorem. A cycle uses a varying argument2^A/3^r, and proximity of its
logarithm to the integer0 has denominator1. An asymptotic rational-
approximation threshold for each fixed argument does not exclude it.

For an actual odd cycle x0,...,x_(r−1), with total valuation A, the exact
identity is

    A log2−r log3 = sum_j log(1+1/(3x_j)) >0.

The useful bridge is a lower bound for this linear form against a
source-specific upper bound retaining actual orbit states and distinctness.
[The cycle-prefix packet](collatz_log_gap_20261014a.md) develops this into
an exact exclusion gate, strengthening the old repeated-minimum estimate.
It still concerns cycles, not all possible escaping paths.

"All logarithms" needs an argument class: for a Liouville number L,
log(exp L)=L has infinite exponent. The hypothesis mu(zeta(2k+1))=2 for
k≥1 remains an additional hypothesis here; it includes irrationality of
every one of those values. A zero-free half-plane of a zeta or L-function
does not itself give rational-approximation estimates for these values.

The other supplied sources suggest specific operations:

* [Ramanujan's inequality repository](https://github.com/blueberryvertigo/ramanujan-prime-inequality-and-missing-factor)
  assembles overlapping certified ranges and a final analytic tail. The
  transferable design is a complete coverage ledger with overlap margins.
  For Collatz, neither the last finite census nor sparse grounded families
  currently supplies that final all-source tail.
* [Primitive roots for every admissible integer base](https://github.com/openai/math/blob/main/preprints/Primitive-roots-for-every-admissible-integer-base-October-4-2026/primitive-roots-all-integer-bases.pdf),
  Theorem1.1 and sections1.3–2, separates possible order-index divisors by
  constructing controlled predecessors and excluding each range. Its
  base2 consequence, accepted as supplied, gives infinitely many primes
  with one doubling orbit on nonzero residues. Through the independently
  scoped THM-4631 level factorization this gives one fast-mode orbit at
  prime level; lifting to p² still requires its own order condition.
  It does not transfer a prime-counting density into pointwise Collatz
  coverage. The information note gives exact finite orbit examples.
* The quoted [Weil-group extension](https://hexagonmath.org/2610.00200v1)
  is a supplied claim; the full page returned403 during this session.
  Its described descent/induction pattern suggests proving a primitive
  obstruction and transporting it only through predicate-preserving
  operations. No L-function has been constructed whose zero-free region
  certifies these Collatz addresses.
* The [boundary-function archive](https://github.com/micah541/kaczynski-mathematical-works)
  and S6 automorphism example are developed in the information note: keep
  the approach path or incidence interpretation when changing observers.

The logarithm paper's determinant method also suggests a precise task:
find a nonzero arithmetic determinant tied to actual residual sources,
then force an upper bound below its arithmetic minimum. Word counts or
dimension counts alone do not provide the necessary nonvanishing or
source-independent thresholds. This remains a proposed construction.

## 6 Fisher information can be made verifier-specific

Let p be a strictly positive law on a finite set of source records, and
V their Boolean certificate verdict. Tilt by

    p_theta(x)=p(x) exp(theta V(x))/Z(theta).

For a deterministic compression S, the exact Fisher-information loss is

    I_full−I_S = E_theta Var_theta(V | S).

It vanishes iff the verdict is constant on every S-fibre, precisely when
the verifier can be evaluated from the compressed record. This supplies
a rigorous target-dependent test for the user's information perspective.
The full derivation, an ancillary-information hostile, and equal-cost
Collatz words with different source cells are in
[the information packet](information_observers_collatz_20261014a.md).

A point mass describes an exactly known value, not a theorem about every
other source. For Gaussian location measurements with variance sigma²,
Fisher information is1/sigma² and grows as the variance shrinks. The
limiting Dirac family is singular for the usual density/score setup;
zero variance should not be equated with zero Fisher information.

## 7 Reproduction and the next proof obligations

Run

    python 04-computation/experiments/collatz_parity_cf_boundary_20261014a.py

The saved JSON records112874 exact checks. Normal and python−O runs agree
byte for byte; normalized-LF SHA256:
`4ca6504b09511429637d3834aaad072ea1b262bd8f4e6fc0bb82c78a9658b9ed`.
The universe is all8190 parity words of lengths1..12, each with three
positive lifts; Mersenne shadows through64; all sources1..2048 for a
14-digit gated-code prefix; and every binary digit word of lengths1..16
for the separate field exploration, followed by the explicit length18
field hostile and its complete rational orbit. Positive controls include
1,2,3,27; hostiles include signed−1, the nonintegral cycle, and an explicitly
abstract never-hit process. Fourteen malformed interface inputs are rejected.
No finite check is used to infer an infinite exponent or termination.

The next obligations are different and should remain separate:

1. For cycle exclusion, enlarge the source-owned prefix gate while
   retaining carry, distinctness and an exact integer clock test.
2. For universal convergence, prove membership in the specified root-tail
   class, or independently bound the stopping-gated exponent below M.
3. For certificate compression, preserve the verifier on every source
   fibre; test this with its own tilted statistical family, not an
   unrelated average or score.

These are exact connections with testable interfaces. None supplies the
still-missing all-source Collatz closure.
