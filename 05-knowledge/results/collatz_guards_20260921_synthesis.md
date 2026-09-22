# Collatz guard counterarguments: reset density, bounded discrepancy, and squarefree word obstructions

**Date:** 2026-09-21. **Status:** PROVED scoped results; CITED restricted
one-cycle theorem; FINITE-EXACT computations; VERIFIED Lean reductions.
**Global Collatz convergence remains OPEN.** None of the results below
classifies all positive or negative cycles. No priority claim is made.

The new [user-supplied blueprint](../reference/COLLATZ-GUARDS-2026-09-21-SOURCE.md)
correctly moves toward ordered affine carries, but its proposed global
proof still lacks a descent argument for every starting integer. Exploring
the retained coordinates produces more than counterexamples: an exact
irrational reset density, a proof excluding bounded scaling discrepancy,
and sharp squarefree restrictions on finite halving words. The Lean work
separately repairs the actual trajectory counters and formalizes the
remaining global obligation without assuming it.

## 1. What was inherited, and what changed

| Live concept | Closest earlier mechanism | Retained coordinate / new result |
|---|---|---|
| Guarded affine words | [Previous affine audit](collatz_blueprint_20260921_affine.md) | The exact binary cylinder replaces a mod-nine valve claim |
| Arbitrary finite paths | [Inverse completion](arithmetic_braids2_20260917_inverse_completion.md) | Positive realizations of every finite prefix need not give a positive realization of an infinite word |
| Mechanical clocks | [HYP-2456 Beatty-Pell address](../hypotheses/HYP-2456-beatty-pell-crossover-word.md), still an OPEN synthesis | Separate the proven address-clock operation from carry and ordinary height |
| Squarefree paths | [Joint linear-form sieve](arithmetic_braids2_20260917_squarefree_symmetry.md) | Prime-square zero positions decide finite-word admissibility |
| Signed dynamics | [Arithmetic braids](arithmetic_braids_20260917_collatz.md) | The sign of the affine carry changes a normalized positive sum into a convergent subtractive sum |
| Global proof target | [Existing Lean audit](../../04-computation/lean/CollatzBlueprintAudit/README.md) | Count odd operations separately from total steps; retain the universal quantifier |

The anchor is the descent target, the niche is arithmetic capacity, and
the wildcard is squarefree symbolic restrictions. The canonical hostile
remains a very long all-one growth word. The corrected near miss is
inferring each orbit's behavior from source density. The underused
coordinates recovered here are the odd cofactor of n+1, the count of
ordinary integers below a height, and prime-square phase positions.

The comparison with Beatty-Pell work transports a **method**, not a
theorem about unrelated objects: both have an irrational floor clock
and extra state needed to recover the visible arithmetic. Here the
clock is derived directly from powers of two and three. No map from
tournament counts, E8 fibers, or the isolated prime sum 196 to every
Collatz orbit has been constructed.

## 2. The pasted counterarguments: strongest survivors

| Claim in the attachment | Status and decisive correction |
|---|---|
| Ordered carry `2^K T^L(n)=3^L n+B_L` | **PROVED**, with actual accelerated steps and ordered cumulative exponents |
| `5 mod9` forces the 23 reset | **REFUTED**: 95 has `95->143->215->323`; the exact reset word `(1,1,5)` requires `23 mod256` |
| The carry alone forces descent | **REFUTED**: the 23 and 95 three-step paths both have carry 19, but totals K=7 and K=3 |
| Baker's theorem leaves only one positive rational cycle | **REFUTED**: the inherited rational cycle is `5/7->11/7->5/7`; general integer-cycle exclusion has not been proved |
| Exactly three negative cycles | **UNESTABLISHED**: three known cycles do not establish exhaustiveness |
| The 196 prime-sum collision forces global scaling bounds | **UNSUPPORTED IMPLICATION**: the equality holds, but no map to a trajectory or inequality is supplied |
| The displayed Lean passes an axiom audit | **REFUTED as supplied**: it contains `sorry`, merges declarations, mixes clocks, and states a false all-natural conclusion at zero |

The legitimate restricted use of the cycle literature is explicit in
the [valve note](collatz_guards_20260921_valves.md): a first-reset return
would be a one-rise/one-fall cycle, and the cited Steiner result excludes
nontrivial cycles of that shape. This is far weaker than excluding all
cycles. The new density theorem below does not need that external result.

The prime sum also has a repaired exact explanation. If p_i is the ith
odd prime and c_i counts the odd composites at most p_i, then
`1+sum_(i<=k)p_i=(k+1)^2+2*sum_(i<=k)c_i`. At k=11 this reads
`196=12^2+52=14^2`. The missing coordinate is the cumulative count of
skipped odd composites. The extra 7 in the attachment's list would give
203. This identity supplies no Collatz descent estimate.

## 3. A reset law with a genuine irrational clock

For odd n write `n=2^(r+1)u-1`, u positive odd. There are exactly r initial
halving exponents one, followed by the first exponent
`k=1+v_2(3^(r+1)u-1)>=2`. If y is this block's endpoint, then

```text
2^(r+k)y=3^(r+1)n+3^(r+1)-2^(r+1).
```

Put `D=2^(r+k)-3^(r+1)`. If D<0 the block grows. If D>0, integrality
forces y<=n, with equality only at `D*u=2^(k-1)-1`. The spatial law of
the exact block is `Pr(r,k)=2^(-r-k)`. It follows, entirely elementarily,
that strict descent at the first reset has relative natural density

```text
P=sum_(L>=1)2^(-floor(L*log_2(3)))
 =0.713725497675891258032337198700... .
```

The binary digit positions of P form a Beatty sequence. In particular
P is irrational. The same density holds inside any fixed odd-modulus
source row, including `5 mod9`. These are source-counting statements;
the actual next block is determined by the previous cofactor and carry.

For example, the family `n_a=3*2^(2a+1)-1` extends 23. Its first-reset
endpoint descends for a=1,2,3 and **exceeds its start for every a>=4**.
Thus the motivating micro-example has an infinite regime with opposite
first-reset behavior. The [full proof and exact interval for P](collatz_guards_20260921_valves.md)
keep both the positive mechanism and its failure boundary.

## 4. Bounded discrepancy cannot be a positive integer trajectory

For an actual positive odd orbit define `K_j` as the first j halving
exponents summed, and put `q_j=2^K_j/3^j`. The ordered carry becomes

```text
n_j*q_j=n_0+(1/3)*sum_(i<j)q_i.                           (S1)
```

**PROVED:** no positive integer orbit has q_j bounded above and bounded
away from zero. Equivalently, `K_j-j*log_2(3)` cannot stay in any fixed
bounded interval, even eventually.

The proof first derives two-sided linear height bounds from hypothetical
bounded q. Repeated states are impossible under that hypothesis, since
each cycle would multiply q by a nonunit ratio of a power of two to a
power of three. A tail would therefore be a positive-density set of
distinct integers, each of which never subsequently falls below a fixed
fraction of itself. But that exceptional set has density zero: the
finite halving-cylinder formula and a Chebyshev bound give a direct
proof. This establishes the needed bridge from a source-counting result
to this particular hypothetical orbit; it does not assume random
sampling along trajectories.

**Corollary:** an infinite mechanical halving word with cumulative
`K_j=floor(j*log_2(3)+rho)-floor(rho)` cannot belong to a positive
integer orbit. Every finite prefix is nevertheless realized by infinitely
many positive integers. Their nested cylinders converge in the 2-adic
sense, and their least positive representatives tend to infinity.

The [discrepancy note](collatz_guards_20260921_discrepancy.md) also proves
capacity bounds for general `an+b`, sharpened by the invariant
`gcd(n,b)`. For b<0, positivity itself gives
`sum q_j<=a*n_0/abs(b)`, hence q_j tends to zero. This applies to positive
`3n-1` orbits and to `3n-5` orbits that remain positive. It does not
classify their cycles. The stronger bounded-strip exclusion extends
to every fixed positive odd b in `3n+b`; its density proof does not
extend unchanged to multiplier five.

## 5. Squarefree itineraries have a finite local test

Squarefreeness of the starting number alone is compatible with every
finite exponent word. Requiring **every displayed odd node** to be
squarefree is a different predicate with genuine forbidden words.

For a word of length L, intersect its source cylinder with `5 mod9`.
The L+1 nodes become integer linear forms. At a prime p>=5 each form
forbids one source residue modulo p^2; at two and three there is no
obstruction. Consequently only primes with `p^2<=L+1` can cover all
source residues. Checking these finitely many primes is both necessary
and sufficient for positive-density all-squarefree realizations, by
the fixed-word linear-form sieve.

**PROVED sharp boundary:** every exponent word of length at most 23
has all-squarefree realizations with source `5 mod9`. At length 24,
the word consisting of 24 exponents all equal to three is forbidden:
the map `(3n+1)/8` cycles through every residue modulo 25, and its
25 displayed nodes must include a multiple of 25.

For the motivating three-step reset, set `R(n)=(27n+19)/128`. Its
slope gap is 101 and its carry 19 is a unit modulo 101. Therefore R
acts as one full cycle modulo `101^2`. Tracking the three intermediate
phase roots, including the last endpoint, gives the exact theorem:

```text
r successive copies of (1,1,5), with all 3r+1 nodes squarefree,
have positive-density sources in5mod9  iff  r<=4985.
```

At r=4986 one node must be divisible by `101^2`. All other primes have
explicit avoiding residues. This is an arithmetic obstruction attached
to a specific word, not a universal trigger encountered by every orbit.
The same mechanism gives a whole constant-exponent family; see the
[squarefree proof and finite certificates](collatz_guards_20260921_squarefree.md).

This connection also explains why squarefree density cannot by itself
be a descent potential. Arbitrarily long growing words can remain
squarefree, while long enough repetitions of a contracting word cannot.
The relevant product is a joint linear-form sieve depending on the
word, not the ambient single-integer constant `6/pi^2`.

## 6. What is actually formalized

The extended [Lean package](../../04-computation/lean/CollatzBlueprintAudit/README.md)
has 35 root-imported audited theorems. The new
[ClockAudit module](../../04-computation/lean/CollatzBlueprintAudit/CollatzBlueprintAudit/ClockAudit.lean)
defines actual counts for the standard map C:

```text
K+O=t,                    2^K*C^t(n)=3^O*n+B.
```

Here O counts odd operations, K counts even divisions, and t counts
all standard steps. It proves the affine identity, the corresponding
local descent equivalence, and the equivalence between convergence of
all positive starts and a **universally quantified** actual-tally margin.
The latter margin remains unproved.

The distinction matters even at n=2. The attachment's free metadata
can choose fake K=3,B=2 to satisfy its local equation, despite the true
one-step counts being K=1,O=0,B=0. With the true K, its incorrect `3^t`
clock fails. A generic positive map with one descending point and a
different fixed point refutes the unrestricted local-to-global inference.
It does not, by itself, refute the still-open Collatz-specific positive
equivalence. The attachment's global statement on **all naturals** is
directly false at zero, and that refutation is also formalized.

The new analytic and sieve proofs in sections 3--5 are independently
audited written proofs with exact controls; they are **not claimed as
Lean-formalized** in this checkpoint. The Lean build contains no `sorry`
or custom axiom dependencies; allowed dependencies are `propext` and
`Quot.sound`. The actual audit records, rather than source appearance,
are the evidence for this statement.

## 7. Further hypotheses and decisive tests

1. **OPEN: control the cofactor dynamics of successive resets.** With
   `A=3^(r+1)u+2^(k-1)-1`, the next block has
   `r'=v_2(A)-k` and `u'=A/2^(k+r')`. This exact map keeps the dependence
   lost by the independent *spatial* law of r and k. A useful candidate
   certificate would depend on this cofactor and accumulated carry.
   It must survive the family n_a, every finite-word completion, and
   repeated contracting cylinders before claiming global drift.
2. **OPEN: quantitative excursion capacity.** The bounded-strip proof
   supplies an honest density-to-orbit bridge. For an unbounded
   discrepancy excursion, a proof must replace the lost linear height
   bound and the fixed fractional non-descent predicate by quantitative
   bounds that can still be compared. The basic conditional test is
   exact: if N distinct orbit values lie below H(N) in an exceptional
   set E, then `|E intersect [1,H(N)]|>=N`. Any proposed upper bound
   must contradict this inequality on the same domain and scale.
3. **OPEN: organize the squarefree finite-word language.** The new local
   test decides each finite word and detects the shortest forbidden
   length. Study how these prime-square exclusions interact with block
   concatenation and discrepancy excursions. A finite automaton for all
   word lengths is not supplied; the tested prime range grows with L.
   Long squarefree all-one words are a mandatory hostile to any claimed
   global contraction consequence.

These are defined research objects and tests, not premises smuggled into
a convergence proof. The immediate formalization path is the general
finite-word source-cylinder theorem, followed by the finite-prime
squarefree criterion and the density/capacity implication. Formalizing
the global conclusion still requires new mathematics controlling every
actual unbounded excursion.

## Reproduction and checkpoint

Run `python 04-computation/experiments/collatz_guards_20260921_verify.py`.
The [manifest](collatz_guards_20260921_manifest.json) records normal and
optimized exact replays, independent modular paths, the clean Lean
build/axiom audit, and hashes of the proof notes and source archive.
The [transcript](collatz_guards_20260921_verification.out) records commands.
The earlier blueprint manifest is refreshed against the extended formal
package; its original mathematical scope is unchanged.

Independent reviews checked the reset classification/density, all-prime
sieve quantifiers and sharp endpoint cutoff, and the strengthened
bounded-discrepancy proof. Computational universes and positive/hostile
controls are declared in each script. No finite census is used as a
proof that all trajectories converge or that all cycles are known.
