# Independent audit: fixed three-ear radix-16 obstruction

## Finding

**PASS. No mathematical error identified in the obstruction, occupancy identities, spectral filtration, or complementary-template formula.** The result excludes a fixed repeatable bank of two or more distinct balanced radix-16 digit actions at an actual seed H>1 in this three-ear/transitive-spine architecture. It does not exclude any tournament Hamiltonian-path count and does not address banks whose admissible templates depend on phase/state so that a physical template cannot necessarily repeat indefinitely.

The independent audit uses only standard-library integer/rational arithmetic and imports no producer module, LP solver, symbolic algebra package, or discovery code. It rebuilds the forward transition by assigning each new spine vertex to each admissible run, then obtains row actions by transposition. Thus the independently implemented transition direction is different from the producer's direct row-action code.

## Exact evidence

`audit.py`, normal and optimized, gives byte-identical receipts:

- **280,608 explicit checks**;
-324 states independently reconstructed;
-9,362 independent subset-DP comparisons: all4,681 spine words of lengths0 through4, for both ear-core types;
-48 prefix-normalization row invariances;
-248,832 occupancy quotient coordinate equalities and all80 initial cross-order occupancy relations;
-52 exact rational separating certificates, with16,848 nonnegative coordinate inequalities and strict negative normalization constants;
- exact universe coverage: all26 nonzero eligible two-letter templates for each core;
-16 complementary-template operator identities,5,184 coordinate equalities;
- eight corrupt separating certificates and two corrupt identity tests rejected;
- exact universal cubic coefficient4/27.

`spectral_and_example.py` adds175 checks, also byte-identical normal/optimized. It independently checks common local triangular flags for all32 letter/run configurations, spectral multiplicities, and the four warning-example values via full tournament subset DP.

Inputs are only `farkas_certificates.json` and `farkas_occupancy_certificates.json`; their SHA-256 values are in the receipts. All transfer, acceptance, occupancy, spectral, and example data are reconstructed independently.

## Transfer and normalization

For each of six possible ear orders, a path partitions the transitive spine into four increasing runs. Every ordered path has one such assignment and every accepted assignment yields one path. Prefix endpoint, two middle endpoints, and suffix-start status therefore suffice. A run can start only when its preceding ear beats its first vertex; every already-started run can receive a new larger spine vertex, overwriting its endpoint status. Acceptance additionally checks that each preceding endpoint beats its following ear, or that an empty middle run is replaced by the corresponding core arc. These are exactly the implemented324 states and acceptance conditions.

The prefix-only assignment is unique for every spine word, giving each order's normalization ell_o(v)=1. The independent checker verifies ell_o M_a=ell_o for every order and letter, and checks the initial values. Nonnegativity follows from the counting interpretation.

## Growth exclusion and fixed-bank quantifiers

If a two-letter template has a common set bit, the corresponding ear beats every appended vertex during any repetition. In each ear order, the increasing run immediately before that ear can contain no new vertex: if it did, its last vertex would be new and lose to the ear. At most three runs are therefore available to every appended vertex. For a seed spine of length n, counting all possible old assignments and six ear orders gives

`H(seed (ab)^k) <= 6*4^n*3^(2k)`.

A repeatable balanced recurrence gives `H_k=16^k*(H_0+d/15)-d/15`, which is at least `16^k*(H_0-1)+1` when H_0>1 and |d|<=15. The growth rates contradict each other. Exactly27 pairs remain because each of three bit-columns may be00,01,10; the other37 are excluded.

The certificate exclusions then leave at most the physical template00. At a fixed seed, one physical template has a single output count, so it cannot realize two distinct digits. This proves the stronger two-digit-bank obstruction, and hence the stated sixteen-digit obstruction. The argument genuinely uses arbitrary self-repetition within the fixed bank. It does not apply merely because a template appears once in a phase-dependent construction. H_0>1 is genuinely used in the growth step.

## Rational separators and occupancy

An affine16 recurrence implies `R M^j v=0` for all j>=0, with `R=c(M-I)(M-16I)`. Dividing each finite row by its positive gcd does not alter its zero set or introduce a sign change. The independent checker rebuilds these rows directly, clears denominators of all certificate multipliers using a positive common multiple, then checks every coordinate by integer arithmetic.

For each accepted separator, the combined row is coordinatewise nonnegative, while its value on a normalized state satisfying the recurrence and any appended occupancy relations is a strictly negative constant. This is a contradiction. No floating-point tolerance, integrality assumption on v, or undocumented rank assumption is used.

For occupancy, the prefix is labeled P and the other three runs are labeled by their preceding ears. A nonempty run remains nonempty; an empty ear-labeled run can start exactly when that letter's ear bit is1. This sixteen-state quotient is independent of ear order. The independently verified intertwining identities and common empty-spine initial values prove the80 cross-order relations for all reachable states, not merely enumerated examples.

The first bank covers precisely the14 eligible pairs with union-support size1 or3, for both cores. The occupancy bank covers the12 eligible pairs with support size2, for both cores. Their union is exactly the26 nonzero eligible pairs. No existence conclusion for00 is inferred from its survival.

## Spectral filtration

An independent common basis for the three-state prefix operator is

`(e_good-e_bad, e_empty-e_good, e_bad)`.

Its diagonal is(0,0,1). For either middle operator, use

`(e_good-e_bad, e_good, e_empty)`;

the diagonal is(0,1,0). For the suffix, use `(e_started,e_empty)`, with diagonal(1,0). Every letter's local matrix is upper triangular in these same bases; the audit checks this exactly.

The full per-order operator is the Kronecker sum of the four local operators. Its fixed diagonal multiplicities are the coefficients of `(2+t)^3(1+t)`, namely8,20,18,7,1. A product of L letters remains triangular with diagonal eigenvalues0,1,2^L,3^L,4^L. The six independent order blocks give the stated characteristic polynomial. This does not assert diagonalizability or commutation of the letter matrices.

## Complementary templates and4/27

For all eight complementary templates and both cores, the independently computed factored row equals3840 times the sum of the six normalization rows:

`c M(M-I)^3(M-4I)^3(M-9I)^2(M-16I)=3840*sum ell_o`.

Consequently, every normalized actual seed satisfies the claimed scalar inhomogeneous recurrence with right side23040. Applying one further(E-1) yields characteristic roots0 (multiplicity1),1 (multiplicity4),4 (multiplicity3),9 (multiplicity2),16 (multiplicity1). Thus for k>=1 the solution has exactly the allowed exponential-polynomial form; a zero-root transient can affect k=0. Rational initial counts and rational characteristic roots imply rational coefficients by an invertible rational confluent-Vandermonde system.

On its cubic polynomial term, the factor outside(E-1)^3 contributes

`1*(1-4)^3*(1-9)^2*(1-16)=25920`,

and the third forward difference of k^3 is6. Therefore its coefficient is necessarily `23040/(6*25920)=4/27`. Distinct characteristic modes cannot cancel this nonzero cubic contribution. A pure affine16 recurrence contains only a16^k mode and a constant, so it is incompatible. This complementary-template contradiction needs normalization but not positivity.

## Finite-search scope

The explicit cyclic example35751 with template70 was independently checked by subset DP:

`215,3443,55091,878159`.

The first two increments fit16H+3; the third falls short by3300. The producer's much larger enumeration of all words through length7 and its coverage maxima were **not rerun** in this audit. Those search statistics are expressly motivational and are not used in the independently verified obstruction proof. Likewise the producer's additional200 longer random transfer comparisons were not reused; this audit performed its own complete short-word comparisons and exact structural proofs.

## Reproduction

```
python audit/audit.py
python -O audit/audit.py
python audit/spectral_and_example.py
python -O audit/spectral_and_example.py
```

Portable receipts: `receipts/independent_audit.normal.out`, `receipts/independent_audit.optimized.out`, `receipts/spectral_example.normal.out`, and `receipts/spectral_example.optimized.out`. The only relocation edit to the independent audit is its source-relative certificate directory (`../code`).
