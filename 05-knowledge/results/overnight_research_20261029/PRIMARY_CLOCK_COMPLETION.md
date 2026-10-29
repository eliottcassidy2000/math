# Primary clock layers, adjoining roots, and compatible pairing completions

2026-10-29. **PROVED scoped module theorems / FINITE-EXACT native modules.**
This completes the C7 higher-primary obligation in
[CLOCK_MINIMAL_PAIRING](CLOCK_MINIMAL_PAIRING.md). It does not complete
the extra affine/run clocks at2/3 or provide an arithmetic Sha realization.

**Continuation:** [FULL_NATIVE_COMPLETION](FULL_NATIVE_COMPLETION.md)
now gives the full126-element vertex/run minimum at2, index2^40, and a
universal index3^10 lower bound at3, after the
[root-extension obstruction](TERNARY_EXTENSION_OBSTRUCTION.md) excludes3^8.
The formulas below still concern
their specified cyclic clocks. Root attachment and pure completion
remain different operation sets.

## Inheritance and concept board

Closest proved mechanism: reciprocal constituents in C4 of
CLOCK_MINIMAL_PAIRING, together with the actual complete Smith certificate
in [THREE_RUN_GLUE R4](THREE_RUN_GLUE.md). Canonical hostile: the labelled
cyclic-triple twins with automorphism orders81/27. Corrected near miss:
forgetting the index54 S3 lattice creates four spurious Z/6 factors at T5.
Least-used relevant sidecar: the original-generator recovery maps in the
actual [two-saturation certificate](../session97_overnight_arithmetic_carriers_20261029/halving/phase2/README.md).

The live concepts are native endpoint ports, reciprocal clock factors,
primary depth, root adjoining, and compatible reduction. Anchor: the
smallest paired module retaining C7 and every source element. Niche:
distinguish adding roots from adding independent generators. Wildcard:
ask whether minima at separate precisions assemble into an injective
ordinary quotient tower. The last question has an exact obstruction below.

## P1. The all-exponent minimum with an odd cyclic clock

Let V be a finite abelian p-group with a specified C_n action, where n is
odd and p does not divide n. A completion means a C_n-equivariant injection
V -> W, where W is a finite abelian p-group with a perfect alternating
pairing W x W -> Q/Z preserved by C_n. The ambient clock must also satisfy
C^n=1. No pairing on V is prescribed, and V need not be a direct summand.

Factor X^n-1 over Z_p into its coprime unramified irreducible factors f.
The resulting idempotents decompose V into modules over the unramified
discrete valuation rings O_f. Write each partition in descending order:

    V_f = direct sum_i O_f/(p^lambda_f,i).

Let f* be the monic reciprocal, d_f=deg(f), and pad partitions by zeros.
There are three cases, giving the exact minimum ambient partitions:

1. For f=X-1, duplicate the odd-position exponents:
   (lambda_1,lambda_1,lambda_3,lambda_3,...).
2. For an unordered reciprocal pair f!=f*, put
   mu_i=max(lambda_f,i,lambda_f*,i) on BOTH sides.
3. For a nontrivial self-reciprocal factor, retain the original partition.

The exact minimum logarithmic index is therefore

    Delta(V) = sum_j(lambda_0,2j-1-lambda_0,2j)
             + sum_{pairs f!=f*} d_f sum_i |lambda_f,i-lambda_f*,i|.  (P1)

Here zero denotes the trivial factor. This extends the elementary formula
C6, including p=2 and all primary exponents. These are elementary module
arguments, not a claim of new historical representation theory.

**Proof of necessity.** For modules over an unramified DVR, an injection
of type lambda into type mu requires lambda_i<=mu_i for every i.
Indeed the p-socle of p^(h-1)V injects into that of p^(h-1)W; their residue
field dimensions count the parts of height at least h. Those inequalities
for all h are equivalent to the componentwise inequalities. A perfect
invariant pairing pairs an f-isotype only with its reciprocal. Dual
isotypes therefore have identical partitions. The trivial isotype has
each cyclic multiplicity even: split off a primitive hyperbolic pair of
maximum order and iterate. These constraints give the stated lower bound.

**Attainment and source map.** An O/(p^e) source factor embeds into an
O/(p^E) target factor, E>=e, by x -> p^(E-e)x. Pair matched reciprocal
factors by the perfect trace pairing and pair the trivial factors in
hyperbolic planes. All embeddings are O-linear, so the original clock
and source elements are retained, even when roots were adjoined.

For a nontrivial self-reciprocal factor of degree2r, the unramified
involution bar sends the clock root to its inverse. Choose a unit
c=z-bar(z) with bar(c)=-c; such a unit exists by taking a residue outside
the fixed field. The pairing

    (x,y) -> Tr(c*x*bar(y))/p^e mod Z

on O/(p^e) is perfect and invariant. It is genuinely alternating also
at p=2: c*x*bar(x) has relative trace zero in characteristic zero before
reduction, rather than merely satisfying a skew-symmetry identity modulo2.
The unramified trace pairing is perfect, completing the construction.

## P2. An exact recurrence through all primary heights

Let r_f(h) count the parts of lambda_f of size at least h. Define I_H to
be the minimum index for the SINGLE source V/p^H V. Then

    I_0=1,
    I_H/I_(H-1) = p^[ (r_0(H) mod2)
                     + sum_{pairs f!=f*} d_f |r_f(H)-r_f*(H)| ].   (P2)

This follows by summing heights of the Ferrers diagrams in P1. Taking
componentwise maxima and duplicating odd positions commute with truncating
each exponent at H. Consequently P2 is an all-height exact recurrence,
not a fit to the four observed Paley layers.

The word SINGLE is essential: numerical minima at separate heights need
not assemble into an injective system under ordinary quotient reduction.

## P3. Actual Paley7 layers and minimum indices

Let M be the210-by210 three-run endpoint matrix for Paley7 with directed
steps{1,2,4}. The named ports are ordered triples of distinct vertices.
Translation tau:v->v+1 has exactly30 free orbits of length7 on them.

At p=2, write eta^2+eta+2=0 with eta even in Z_2 and eta_bar=-1-eta.
Then Phi7 has reciprocal cubic factors

    f+ = X^3-eta X^2+eta_bar X-1,
    f- = X^3-eta_bar X^2+eta X-1.                       (P3)

Modulo2 these are X^3+X+1 and X^3+X^2+1. The exact O-module partitions
of the actual cokernel are:

| Prime | Trivial factor | f+ | f- |
| ---: | --- | --- | --- |
| 2 | (2,2,1,1,1,1,1,1,1,1,1,1) | (4,4,3,1,1) | (4,1,1,1,1,1) |

At p=3, Phi7 is irreducible of degree6. Its entire source component is
zero; the full3-primary group has trivial C7 action and partition
(4,2,2,2,2,2,2,1,1). These recover the independently established group
orders2^80 and3^18 and ALL integer elementary-divisor multiplicities.

| h | 2-primary isotype ranks (trivial,+,-) | Added binary index exponent | Cumulative exponent |
| ---: | --- | ---: | ---: |
| 1 | (12,5,6) | 3 | 3 |
| 2 | (2,3,1) | 6 | 9 |
| 3 | (0,3,1) | 6 | 15 |
| 4 | (0,2,1) | 3 | 18 |

For p=3 the trivial ranks are9,7,1,1, so the added index exponents are
1,1,1,1. The minimum C7 completion orders are therefore2^98 and3^22.
At2 both reciprocal target partitions become(4,4,3,1,1,1); the trivial
partition is already alternating. All maps from named ports are retained.

Combining with the six already certified elementary odd-primary pieces
in CLOCK_MINIMAL_PAIRING gives the minimum C7-equivariant index for the
ENTIRE native cokernel:

    2^18 * 3^4 * 11^6 * 79^3 * 5861^6 * 392809^3.       (P4)

Primes5 and17 contribute1 for this C7 action and, by P6 below, for the
full vertex/run group. This is not yet a minimum for that larger group:
its modular actions still need compatibility tests at2 and3. The four one-sided
odd-primary components already have full-clock optimality by C4's
isotropic lower bound and functorial hyperbolic double.

## P4. Adjoining roots is cheaper than a direct-sum repair

The simplest hostile is V=Z/2+Z/4. Adding independent summands to make
each cyclic multiplicity even costs index8. Adjoining a square root of
the first generator gives the injection

    (a,b) -> (2a,b) in (Z/4)^2,

of index2. The source is retained; its embedding is not a direct summand.
Thus the DIRECT-SUM qualification in THREE_RUN_GLUE R6 is substantive.

For Paley7 the exact comparison is:

| Source | Minimum injective index | Minimum direct-sum index |
| --- | ---: | ---: |
| 2-primary, forgetting clock | 2^2 | 2^8 |
| 2-primary, retaining C7 | 2^18 | 2^30 |
| 3-primary, retaining C7 | 3^4 | 3^4 |

The unclocked2-primary partition is(4^9,3^3,2^2,1^31). Its minimum
paired ambient partition is(4^10,3^2,2^2,1^32), explaining the index4.
The C7 penalty is therefore a real constraint of the native operation,
not an inference from square order.

For a reciprocal pair put m_f(e)=r_f(e)-r_f(e+1). A direct-sum repair
adds d_f sum_e e*|m_f(e)-m_f*(e)| to the logarithmic index, whereas P1
adds d_f sum_h |r_f(h)-r_f*(h)|. In the native dyadic pair the height
imbalance is(-1,2,2,1): the two costs are30 and18. Opposite imbalances
at different depths can be absorbed by root adjoining but cannot cancel
when each exponent must be repaired by an independent summand.

## P5. Ordinary quotient compatibility forces a split injection

Fix an injection V -> W of finite p-primary O-modules. Its induced maps

    V/p^h V -> W/p^h W

are injective for EVERY h iff V is a direct summand of W as an O-module.
The forward condition is purity, V intersect p^h W = p^h V. Here is an
elementary splitting proof for this bounded finite setting. Choose an
element v of maximum order p^e in V. Purity implies p^(e-1)v is not in
p^e W. In a cyclic decomposition of W this forces a unit coordinate
of v in a summand of EXACT length e: greater length contradicts the
order of v, and shorter length kills p^(e-1)v. Projection to that
coordinate splits off O v from W and V. The remaining intersection is
again pure; induction splits V. The reverse direction is immediate.
Apply this proof separately in each unramified clock isotype to obtain
an equivariant splitting.

Consequently the minimum FIXED completion retaining injectivity through
all ordinary quotient levels is the direct-sum minimum, not necessarily
P1. In general its logarithmic index is

    sum_e e*(m_0(e) mod2)
      + sum_{pairs f!=f*} d_f sum_e e*|m_f(e)-m_f*(e)|.  (P5)

The example Z/2+Z/4 proves that the loss is unavoidable: its minimum
ambient group is(Z/4)^2, and EVERY order-two element of that ambient
group reduces to zero modulo2. Thus no choice of minimum embedding
can give an injective lower-level map. This is a universal obstruction,
not just failure of the displayed coordinate choice. A full direct-sum
completion gives a compatible tower, but its individual lower levels
need not be the smallest completions of the truncated sources.

Torsion subgroups W[p^h] and quotients W/p^h W have different transition
maps; switching them does not repair a claim about ordinary reductions.
Nor does numerical recurrence P2 by itself supply compatible pairings
or bound an inverse limit. These are the missing obligations when trying
to transport this finite-module construction toward a Selmer tower.

## P6. Full native clocks add no cost at5 or17

The separate checker `full_clock_forms.py` resolves the previously open
elementary cases. It uses the actual quotient projections from
three_run_glue.json and the four literal port permutations generated by
v->v+1, v->2v, a swap of two run slots and a three-cycle of run slots.
They generate Aut(Paley7) x S3 of order126; the vertex factor's order21
was independently enumerated in the earlier certificate.

For each generator R, the checker constructs C with BR=CB on ALL210
named port columns. It solves C^T J C=J on the alternating matrix space
simultaneously for all four generators, then verifies the resulting J
literally and checks full rank.

| Prime | Source dimension | Dimension of invariant alternating forms | Nondegenerate forms | Minimum added dimension |
| ---: | ---: | ---: | ---: | ---: |
| 5 | 12 | 1 | 4 | 0 |
| 17 | 4 | 1 | 16 | 0 |

Every nonzero scalar in each one-dimensional form space is nondegenerate.
The receipt retains both the full space basis and a named perfect form.
Thus ALL six elementary odd-primary pieces now have exact minima for
the full126-element native vertex/run clock group: zero at5/17 and
hyperbolic doubling at11/79/5861/392809. Only the higher-primary2/3
full-clock extensions remain open. Native complement is an additional
operation, not an element of this automorphism group.

## Connections, tests, and remaining obligations

The source-to-target map is explicit: named endpoint coefficients ->
30 orbit polynomials -> unramified factors -> O-linear quotient coordinates
-> the root-adjoining embeddings of P1. It preserves all source classes
and the supplied C7 clock. The factor chart and its inverse preserve
the original labels; no Frobenius orbit is collapsed to just its degree.
It does not automatically preserve native complement operations or
extend the other vertex/run clocks. Source names and the target relation
remain required sidecars for those operations.

The concrete arithmetic analogue is the Session97 two-saturation result:
its13 halving steps preserve original rational points via explicit
identities2Q=sum P_j and inverse recovery. Our root-adjoining module map
preserves elements by the same kind of inclusion, rather than equating
the old and new groups. This is a structural analogy only. No map from
the tournament cokernel into those Mordell--Weil groups, their Selmer
groups or Sha has been constructed. The arithmetic certificate contains
the rational points, local characters and torsion exclusions that the
finite-module formula lacks.

Incoming ef06067902 supplies a distinct native bridge in
[COMPLEMENT_KERNEL K1--K5](../overnight_rank_return_20261029/COMPLEMENT_KERNEL.md).
For every strong source through seven, a selected Hamilton-path complement
makes the paired weight observer integrally left-invertible. Coprime full
minors give a Bezout left inverse, so this weight inclusion remains
injective over EVERY coefficient ring. This is the free-lattice analogue
of the splitting requirement in P5. The objects differ: incoming work
embeds arbitrary port weights into two responses; here we embed a finite
response cokernel into a paired group. No cokernel-to-weight identification
or transfer of the pairing is implied.

The incoming B5 bad-path plane is fixed by vertex automorphisms but has
nontrivial standard S3 run action. Thus extra vertex phase cannot recover
its amplitudes, just as retaining C7 alone does not settle our extra
run-clock problem. Its selected successful path DOES recover them; its
bad fixed path does not improve by repetition. This sharpens the next
native test beyond the earlier T5 example, without claiming that all
physically realizable attachment laws fill the arbitrary-weight space.

Incoming K5 also proves the minimum unnamed weight count
(m^210+14m^70+6m^30)/21 for its successful Paley7 observer. Its precision
recurrence has characteristic roots p^210,p^70,p^30. That recurrence
counts orbit states of the retained paired path frames; P2 instead
counts the minimum index of a finite paired completion. In both cases
the recurrence is legitimate because the actual source and operation
were specified, but only the incoming construction supplies the stated
native complement reconstruction. Neither sequence is A000568.

Reproduce with `python -B primary_clock_layers.py` from this directory.
The checker independently enumerates all three-run covers, verifies
every900 native7-by7 translation block, and supplies the reversible7-by7
factor chart modulo2^6 and3^6. It computes O-linear Smith row maps U,
stores their inverses, checks that EVERY source relation maps to zero,
and independently checks expanded integer matrices with local elimination.
All pivots occur before the precision limit. The receipt preserves
all210 port labels,30 translation orbits, source-row maps and embeddings.
Literal reciprocal trace matrices verify the alternating dyadic pairing.

The finite minimum census covers126 abstract partitions and1,225 ordered
reciprocal partition pairs with parts at most4 and the stated length
bounds. It tests the root-adjoining hostile and replays the full81/27
singleton-slice rank decks. Normal and optimized Python produce identical
receipts; no proof gate depends on `assert`. This is a same-agent
alternative-path audit, not an independent human or agent review.

The dyadic injective problem is now closed in the continuation above.
Next test ternary index3^10 attainment and the pure full-clock problem. For an
arithmetic bridge, first specify whether its transition maps are ordinary
quotients, torsion inclusions or genuinely different Selmer maps, then
supply a commuting source diagram. P5 rules out silently substituting
the cheapest finite-level completion for that diagram. Sha finiteness,
LRC(14), Collatz and a general native attachment-minimality theorem remain
open.
