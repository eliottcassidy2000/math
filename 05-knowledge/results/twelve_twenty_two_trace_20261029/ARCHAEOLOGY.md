# Recovery atlas: earlier12 and22, accepted maps, and rejected shortcuts

**Status: source recovery and scoped synthesis.** The mathematical proofs
of this packet are in [POWER_CLOCKS.md](POWER_CLOCKS.md) and
[FROBENIUS_DUAL.md](FROBENIUS_DUAL.md). Historical prose is not promoted
by appearing in this atlas.

## Search scope and selection

`search_archive.py` searched digit-block occurrences of12 or22 in Markdown
under `01-canon/theorems`, `05-knowledge/results`, and `07-reflections`.
At cached HEAD `2cc2ce454a`,9462 candidate documents were readable;
3692 contained both integers. `archive_search.json` records the roots,
pattern, full both-number path inventory and zero read failures. Dates,
equation labels and identifier fragments can produce hits. This is not
a claim to have deeply read all candidates or exhausted every connection.

The selected chains below were inspected, then compared with current
corrections and incoming Session104 work. Their range is March to October,
including niche rejected models as well as successful constructions.
The decisive test is a common operator or a preserved predicate, not the
frequency with which a small integer occurs.

## March: divisor-prime selection and a necessary repair

The [March16 Bernoulli fixed-point reflection](../../../07-reflections/bernoulli-fixed-point.md)
uses F(n)=product of primes p with p-1 dividing n, for positive even n.
The classical denominator input is the
[von Staudt--Clausen theorem, DLMF24.10.1](https://dlmf.nist.gov/24.10).
The exact selector gives:

|n|Selected primes|F(n)|
|---|---|---:|
|6|2,3,7|42|
|8|2,3,5|30|
|12|2,3,5,7,13|2730|
|22|2,3,23|138|
|42|2,3,7,43|1806|
|1806|2,3,7,43|1806|

The old claim that {2,3} was self-selecting was false:7-1 divides6.
The suggested freeze whenever n+1 is composite fails already at n=8.
The valid chain6 to42 to1806 to1806 uses **every divisor**, not just the
largest one. The reflection is repaired in place and the mechanism is
recorded in `01-canon/MISTAKES.md` under October29.

This supplies a real arithmetic connection to the new clock calculation.
At a split prime p for the seventh-root-trace algebra, p is congruent to
plus or minus1 modulo7. If also p-1 divides12, every nonzero coordinate
of u has twelfth power1. Intersecting this split-prime condition with
the complete Bernoulli12 selector leaves exactly13. Thus the selector
detects a sufficient full-collapse prime for the12-clock. It does not
explain all partial losses at29 and43, identify Bernoulli iteration with
counter iteration, or classify every power-order index.

The October28 recovery
[Borwein ladders, owner numbers, and vSC twelve](../borwein_ladder_owner_numbers_vsc_twelve_20261028.md)
also explains the user's older set {4,2,5,3,7,13}: it is the unfiltered
set {d+1:d divides12}; its prime-filtered version omits4. This is the
same selector operation. No Monster or tournament classification follows
from the cardinality of the resulting set.

## June:12 is a trace,22 an odd-sector dimension

The [June29 signed-cycle/free-monoid synthesis](../signed-cycle-index-verified-SC-extended-and-the-free-monoid-bridge-opus-20260629.md)
records at order6 the exact counts A=56, SC=12, NS=22, and V_merged=34.
For the permutation involution R on the basis of tournament classes,

\[
\operatorname{tr}R=SC=12,\quad
\dim V_+=(A+SC)/2=34,\quad\dim V_-=(A-SC)/2=22.
\]

SC is the number of fixed basis objects, not the dimension of the fixed
linear subspace. A swapped pair also contributes one even vector. The
October24 repair already states this in
[THM-587, per-level signed cycle index](../../../01-canon/theorems/THM-587-per-level-signed-cycle-index-antipodal-euler-number.md)
and the [counting precision notes](../session52_saturation_residues_metagraph_20261024/counting-precision/PRECISION_NOTES.md).
This session propagates that existing correction to the June source.
No topological Euler or Betti interpretation is available without a
specified complex and differential.

The surviving free-monoid map is concrete: a tournament's strong
components form a totally ordered word; converse reverses the word and
converses each component. Fixed words count self-converse classes. The
map preserves this fixed-class predicate, not an unspecified LRC degree.
The new arithmetic distinction is analogous in a precise linear sense:
an observer trace and a carrier dimension are different invariants.

## August: common gauges, folded7, and K4

[THM-3497, Berggren T4 calibration languages and harmonic splitting](../../../01-canon/theorems/THM-3497-berggren-t4-calibration-languages-and-harmonic-splitting.md)
distinguishes the variable gauge language with12 states from the fixed
drift language with192 states. Its harmonic17/96 quantity retains a
specified clock and is not an unguarded physical LRC frequency. The
lesson relevant here is the quantifier: choosing a transporter for each
word differs from choosing one common source gauge. No identification of
that12-state machine with C^12 has been supplied.

[THM-3524, folded C7/K4 common-drift transporter rank obstruction](../../../01-canon/theorems/THM-3524-folded-c7-k4-common-drift-transporter-rank-obstruction.md)
folds negation orbits {0},{1,6},{2,5},{3,4}. Multiplication by2 cycles the
three nonzero paired orbits, precisely the seventh-root trace phases.
The output is a weighted looped symmetric K4 support, not a tournament.
Its common-source rank obstruction is3 versus4; separate Fourier
nonvanishing tests do not create one common transporter.

This is the closest older obstruction to P7: separate clocks can each
look adequate while their common-source span is deficient. The models
are not declared isomorphic. A decisive transfer would have to intertwine
their source maps and joint operators, which their different ranks prevent
in the proposed full-rank form.

## September: compensated prefixes, square-sum growth, and Paley choices

The [September6 planar Jacobian board](../planar_jc48_sep06_board.md)
records a compensated weight22 construction that changes earlier prefix
rows12 through15. It is compatible with the old least-weight claim only
because that claim froze the prefix. This is a source/operation change,
not an equality between two invariants. No planar Jacobian theorem is
obtained by replacing the frozen-prefix guard with a new clock observer.

The September22 square-sum chains, recovered with exact sources in
[Session104's three chains](../session104_descent_exterior_20261029/summands/THREE_CHAINS_RECOVERED.md),
include the three components at Q12 and their mergers at13 and14.
The distinct-parent closure from {2,3}, with holes {1,4,6}, is a different
operation. The consecutive-square identity
21^2+22^2+23^2+24^2=25^2+26^2+27^2 belongs to the all-k balance law
with starting point k(2k+1), as proved in
[the summand-square note](../session104_descent_exterior_20261029/summands/SUMMAND_SQUARES.md).
None of these growth operations is the matrix C until an intertwiner is
provided.

The [September22 Paley/Fano/octonion design](../collatz_mod6_20260922_paley_fano_octonion_design.md)
provides the seven-point aut21/H189 control. Its14 cyclic triples split
into two Fano planes; an octonion construction retains a chosen plane.
F3 gives an actual new map from the counter algebra modulo2 to this
tournament chart, preserving Frobenius and recording converse. This is
stronger than matching7,21 or189; it still carries no Collatz-orbit map.

## October: Pell collisions, ramification, and mixed conjugates

The [October19 parabolic cyclotomic cycle](../parabolic_cyclotomic_cycles_20261019/PARABOLIC_CYCLOTOMIC_THREE_CYCLE.md)
contains the very same cubic q(d). F1 spells out its finite-set affine
map and the derivative obstruction, retaining the independent audit's
boundary. This is an exact algebraic connection, not a global dynamical
conjugacy.

The [October22 Monster/levels/square-sum note](../monster_levels_square_sum_7_21_20261022.md)
distinguishes the group label Fi22 from15 Ogg primes and from field
degrees. Its square-sum graph adds22 with lower neighbors3 and14, because
22+3=25 and22+14=36. A separate Pell collision condition involves
x^2-2y^2=-1. It can be placed in the same quadratic order as the g22 unit,
but a graph-growth interpretation needs the endpoint/height predicates;
a common norm equation alone does not preserve adjacency.

The [October28 Ramanujan recovery](../cursor_cm_phase_20261028/RAMANUJAN.md)
and [trace-clock classification](../cursor_cm_phase_20261028/TRACE_CARRIER.md)
already establish the quadratic Smith formula and exceptional primes.
P6 uses these as a hostile control, not a rediscovery. The elementary
identity involving P(q), P(-q), P(q^2) changes the parameter by a factor4;
nome sign, nome squaring and Fricke inversion are distinct. Ramanujan's
original g22 and g130 factors are pinned in
[his1914 paper](https://ramanujan.sirinudi.org/Volumes/published/ram06.pdf).
The g130 full-conjugate warning now yields P7's exact mixed-readout
counterexample. A one-bit quadratic sign rule does not remove a second
independent quadratic factor.

## Incoming October29: native ports remain part of the source

The [Session104 common-rank audit](../session104_descent_exterior_20261029/audit_incoming/INCOMING_COMMON_RANK_AUDIT.md)
distinguishes12 ordered ports and22 native two-run covers at four vertices.
An aligned native matrix has B=I+N with N^3=0; re-marking can destroy the
common acyclic grading. A product with a specified re-marking has trace22,
while some cyclic unions retain trace12 and require longer words to
detect. Individual unipotence is insufficient for a common acyclic order.
These actual port objects are not the12/22 time samples in P4.

The [positive incidence note](../session104_descent_exterior_20261029/positive/POSITIVE_INCIDENCE_AND_ATOMS.md)
keeps12 owner tags at the tight1..13 LRC boundary, attached to six safe
locations among182 events, with zero interval measure. Tags, locations,
event counts and measure are distinct observations; no12-dimensional
phase space is inferred. The boundary birth term remains necessary.

The [overnight terminal-selector result](../overnight_rank_return_20261029/TERMINAL_SELECTOR.md)
raises rank18 to24 while leaving the clock polynomial unchanged. It
supplies the closest native warning for the new power-clock proofs:
recovering a generator from its powers says nothing about whether the
chosen readout sees all its amplitudes. P3's trace-versus-counter contrast
makes this warning explicit in dimension3.

Finally, the twelve-vertex81/27 hostile retains equal three-run mass
80614764 but unequal ranks430/435 modulo3; its two-run ranks24/22 have
the opposite comparison. Source, ordered ports, joint law and actual
complemented relation survive this synthesis. No scalar clock or field
degree is substituted for them.

## What unifies these accepted connections

The useful common object is a source module together with its allowed
operators, readouts and integer lattice. A missing coordinate can be:
an odd clock phase, a prime-power precision digit, a mixed amplitude,
a relative transporter, or a realization guard. Each has a different
cheapest exposing operation. P1--P5 show when clock algebra alone gives
recovery for free; P6--P7 prove why that inference fails outside its
hypotheses. This is a collection of typed theorems and tests, not an
assertion that all appearances of12 and22 share one universal structure.
