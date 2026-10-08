# Source-owned ternary entries in the binary parameter complement

**PROVED:** an all-reset inverse-entry family, its exact source/exponent
guards, CRT fusion with every finite binary guard, and the endpoint inverse
depth obstruction. **FINITE-EXACT:** the declared labelled bank and controls.
**CONDITIONAL:** a smaller child's ROOT proof is still an external obligation.
**OPEN:** universal coverage and grounding every parameter. No large source
is materialized and no production API searches for ROOT.

## 1. Inheritance and the two places a guard can be applied

The closest mechanism is the corrected ternary-refinement boundary in
[THM-4594, maximal class-decided sieve and sign barrier](../../01-canon/theorems/THM-4594-the-maximal-class-decided-collatz-sieve-and-the-sign-barrier.md).
Its sign obstruction concerns unrefined binary classes. The negative-cycle
inverse entries and their exact normal forms are inherited from
[guarded child normal forms](collatz_child_normal_forms_20261007.md).
The elementary inverse-word algebra is not new; the present useful output
is a reset-indexed bank with no omitted positive native head, its exact
pullback to the actual residual parameter, and a labelled fusion interface.

The canonical hostile is `31 -> 27` as a smaller inverse child, while
`27 --(1,2)-->31`: a paid dependency does not itself ground either source.
The corrected near miss is assuming a new ternary refinement is free at an
already expanded endpoint. The least-used sidecar is its full affine carry.
The anchor is complement coverage, the niche is the two-base guard product,
and the wildcard is a lower bound on inverse depth after expansion. The
live board is original source, ordered carry, native phase, endpoint,
payment and terminal obligation.

Keep the actual child family from
[marked periodicity](collatz_bott_marked_periodicity_20261007d.md):

\[
 E(t)=F+Tt,\quad F=924745897,\quad T=2^{32},\quad
 n(t)=2^{E(t)}-1,\qquad t\ge0.                         \tag{1}
\]

| Source | Map/target | Preserved predicate | Required sidecar / loss |
|---|---|---|---|
| Supplied original n | Exact inverse-word child h | Positive odd integer, h<n, actual h-to-n word | Child ROOT obligation remains |
| Source exponent E | Ternary parameter phase | Exact native inverse guard | Word and immutable source, not just residue |
| Binary and ternary cells | CRT intersection | Both congruences at the same t | Binary height cuts remain separate |
| Expanded endpoint Y | Its ternary observations | Exact affine carry residues | Source payment is not endpoint payment |

## 2. An all-height family indexed by its terminal reset

Let a>=2 be even and choose the **least** k>=0 such that

\[
 Q=2^{k+a}<P=3^{k+1}.
\]

Necessarily k>=1. The preceding failed inequality gives

\[
                         1<P/Q<3/2.                 \tag{2}
\]

For the positive valuation word `w=1^k,a`, direct composition gives

\[
 F_w(h)=\frac{Ph+B}{Q},\quad B=P-2^{k+1}>0,
 \qquad h=\frac{Qn-B}{P}.                            \tag{3}
\]

Its exact inverse guard on positive odd supplied n is

\[
                         n\equiv BQ^{-1}\pmod P.     \tag{4}
\]

**There is no omitted finite positive head.** A native n=1 would give an
odd integer h in the interval (-1,1): the lower bound follows from Q>0,
and the upper bound from Q<P and `2^(k+1)/P<1`. This is impossible. Every
native positive odd n is therefore at least3. By (2),
`Qn-B >=3Q-P+2^(k+1)>0`. Thus h is positive and odd. Since Q<P and B>0,
it is also strictly below n.

Final odd integrality in (3) forces all intermediate formal states to be
integers and odd. Indeed a dyadic noninteger cannot regain integrality
under `(3x+1)/2^a`, and an even intermediate would create such a noninteger
at the next positive valuation. Positivity persists from h>0. Consequently
w is the exact actual word from h to n. Since n>1 and an actual ROOT state
cannot leave ROOT, no earlier ROOT occurs in this prefix.

This is a strict common-future receipt `n --()--> n <--w-- h`, paying the
immutable n. It is not a forward descent from n. A supplied first-hit ROOT
word for h must begin with w, by deterministic dynamics; cutting that
prefix exports a checked ROOT suffix for n. Without that supplied proof,
the receipt ends at an open obligation.

The minimal-k condition is a useful guard, not decorative. It proves (2)
and all-positive native safety. We do not claim arbitrary larger k has the
same finite-head argument. Odd terminal a cannot give a Mersenne endpoint:
every actual endpoint of a word ending in odd a is2 modulo3, whereas a
positive Mersenne number is0 or1 modulo3.

## 3. Exact phases at the original source

Putting `n=2^E-1` in (3) gives the equivalent congruence

\[
                    2^{E+k+a}\equiv Q+B\pmod P.       \tag{5}
\]

Here `Q+B=P+2^(k+1)(2^(a-1)-1)` is a3-unit. The elementary order identity
`ord_(3^(k+1))(2)=2*3^k` shows that (5) has exactly one E class modulo
`2*3^k`. For example, `v3(4^(3^j)-1)=j+1` follows by repeated binomial
expansion, proving the full unit-group order. Modulo3, `Q+B=2^(k+1)`;
subtracting the parity k+a of Q's exponent shows that this E class is odd.

Because F is odd and `gcd(2^31,3^k)=1`, substitution of (1) gives **one**
parameter class modulo `3^k`, with no change to the supplied source. This
calculation uses modular powers and a digit-lift logarithm, not an expansion
of `2^E`. The script independently verifies the full base-two cycles through
3^9, and checks exact order and congruences for larger selected resets.

The following table records a compact concrete bank. G17 is the inherited
inverse `(2048n-2363)/2187`, with its actual forward word retained. G1 has
no positive Mersenne source, so it contributes no phase to this bank.

| Entry | k or word length | Actual word cost | Exact t class |
|---|---:|---:|---|
| G17 | length7 |11|696 modulo729|
| a=2, the inherited G5 |k=1|3|1 modulo3|
| a=4 |k=5|9|84 modulo243|
| a=6 |k=8|14|1007 modulo6561|
| a=8 |k=11|19|34786 modulo177147|
| a=10 |k=15|25|1025955 modulo14348907|
| a=12 |k=18|30|111363584 modulo387420489|

In particular a=4 gives `h=(512n-665)/729` and `E=391 mod486`.
Its least positive odd native source is91 and its child is63, with actual
word `(1,1,1,1,1,4)`. This is different from asserting either integer's
ROOT fate. The a=6 entry has `E=2283 mod13122` and
`h=(16384n-19171)/19683`.

The old G5/G17 phase union has parameter density244/729. The a=4 and a=6
cells are disjoint from it and each other. The resulting four-cell union is

\[
 \tau_6=\frac{2224}{6561},\qquad
 \tau_6-\frac{244}{729}=\frac{28}{6561}.               \tag{6}
\]

Through a=12, all seven **labelled alternatives** are retained. The a=8
cell is contained in G5's cell, so only six cells are needed to count the
union. The exact union mass is

\[
 \tau_{12}=\frac{131325004}{387420489},\qquad
 \tau_{12}-\frac{244}{729}=\frac{1653400}{387420489}.    \tag{7}
\]

The least uncovered parameter is t=0. Counting containment does not license
discarding a receipt with a different child: that alternative might have
the available ROOT proof. `coverage_core` is a counting-only operation;
`entries` returns all labelled records. This matches the retained-obligation
principle of [guard cover defects](collatz_guard_cover_defects_20261007e.md).

These are new entries relative to the explicitly declared G1/G5/G17 bank,
not a novelty claim about arbitrary inverse words or a claim of independence
from every prior Collatz rule. The quantities in (6)–(7) are natural/Haar
densities of the exponent parameter, not positive densities of Mersenne
integers among all integers.

## 4. Fusion improves every finite dyadic complement by an exact amount

Let D be any finite union of dyadic parameter cells and let T be the finite
ternary union above. Write their parameter densities as beta and tau.
At a common precision D is a subset modulo2^h; T is a subset modulo3^b.
CRT gives a bijection of their residue product with residues modulo
`2^h*3^b`. Therefore

\[
 \operatorname{dens}(D\cap T)=\beta\tau,\quad
 \boxed{\operatorname{dens}(D\cup T)=\beta+\tau-\beta\tau.} \tag{8}
\]

This is a counting theorem, not a stochastic independence assumption about
actual orbit steps. Each selected inverse entry is pointwise paid at its
supplied n. If D also consists of authenticated paid guards, (8) is the
exact union of those paid parameter sets. Finite source-height cutoffs in
D remain part of each receipt and can only be dropped for the asymptotic
density statement, not for pointwise acceptance.

There is a controlled countable extension. Suppose finite dyadic truncations
D_M exhaust a guard grammar, and the omitted rules lie in finitely many
dyadic parent cells of total mass epsilon_M tending to zero, uniformly after
their finite height cuts. Finite CRT computes each `dens(D_M intersect T)`.
The same parent cells bound the upper density of the omitted intersection
by epsilon_M, so squeezing proves (8) for the full grammar as well. The
run-surgery grammar in
[parameter cover](collatz_parameter_cover_20261007e.md) has such a tail:
indices k>=M lie in one parent of mass `2^(-M-3)`. This uniform tail, rather
than merely existence of the union's natural density, licenses its fusion.

For an individual binary cell `t=b mod2^h`, `fuse` solves its intersection
with a selected inverse entry, returning the unique class modulo
`2^h*3^k`. This operation does not choose a nearby source. At a supplied t,
both predicates must actually hold. A concrete complement class is
`t=327 mod486`: it is odd and satisfies the a=4 inverse entry, so every
member lies outside the old deep binary guard `t=0 mod2^47` while having
this paid inverse child. Relative to G5/G17 plus that single binary guard,
a=4/a=6 add exactly `(28/6561)(1-2^-47)` of the parameter space.
Coverage against a larger binary bank uses its own beta in (8).

The symbolic child API returns residues of the same actual child by
computing its numerator modulo `P*m` before division by P. Reducing only
modulo m would lose the ternary quotient information, especially when
3 divides m. Actual literal receipts use the supplied integer n; no
symbolic residue is passed off as its ROOT certificate.

## 5. Why applying inverse guards only after expansion is too late here

Let Y(t) be the endpoint after the old child's initial ones and its known
15-letter prefix. Retain its exact carrier
`P0=14348907, Q0=2^32, B0=4309583573`. Then

\[
 Y(t)=\frac{2\,3^{E(t)+14}+B_0-P_0}{2^{32}}.           \tag{9}
\]

For every b<=E(t)+14, its residue modulo3^b is the fixed rational residue
`(B0-P0)/2^32`. In particular **Y=5 modulo9 for every t**. The disjoint
G1/G5/G17 router therefore always takes exactly one G1 step; its child
`(2Y-1)/3` is divisible by3. That child has no positive odd predecessor
at all, because `2^a h-1` is then never divisible by3. Binary refinement
of t cannot alter this terminal fact. Other first inverse branches at Y
are not ruled out by this argument.

There is a broader size obstruction allowing every inverse valuation
letter a>=1. Put `g_a(x)=(2^a x-1)/3`; these maps are increasing and
`g_a(x)>=g_1(x)` for x>=0. For any r-step positive inverse path at Y,

\[
 h\ge g_1^r(Y)=(2/3)^r(Y+1)-1.
\]

If `0<=r<=E-16`, the right side is positive and at least its value at
r=E-16. From (9), `B0-P0+2^32>0`, and hence

\[
 g_1^{E-16}(Y)+1
 > 2^E\frac{3^{30}}{2^{47}}>2^E.                     \tag{10}
\]

Thus **no positive inverse word of length at most E-16 at this endpoint
can pay the original source M_E**, whatever its valuation costs or native
guards. A method following this expanded prefix first would need at least
E-15 inverse odd steps even to escape this size obstruction. This is a
necessary bound, not existence of a paying path at that depth. The constant
16 is the first integer that works in this elementary comparison:
`3^29<2^46`; the script gives a rational E=80 hostile to replacing16 by15
in the same lower-bound argument. It does not assert optimality among legal
inverse paths.

The positive source-level entries of Sections2–4 avoid the expansion and
retain the original payment comparison. This is the meaningful distinction
between reusing a source guard and merely enriching an endpoint observer.

## 6. Reproduction and boundaries

Run normally and optimized:

```text
python -B -X utf8 04-computation/experiments/collatz_complement_guard_fusion_20261007e.py
python -B -O -X utf8 04-computation/experiments/collatz_complement_guard_fusion_20261007e.py
```

The 27,274 exact controls include independent full base-two cycles through
3^9; all resets2,4,...,20 at twelve ordinary native sources each; literal
word verification and symbolic quotient checks; exhaustive independent
verification of the four-cell union modulo6561; all binary residues through
six bits intersected with the a=4 class; endpoint ternary residues through
twenty digits; rational size-bound controls; and malformed or forged inputs.

A single explicit supplied proof `11 --(1,2,3,4)-->1` tests discharge for
the inherited entry `13 ->11`; it exports `(3,4)` at13. The new entries'
all-height claims concern paid dependencies, not newly discovered ROOT
families. Unknown and wrong-source child proofs are not synthesized.
No observation, density or lower inverse-depth estimate closes their
remaining terminal obligations.

Normal, optimized and saved stdout agree. LF-normalized SHA256:
`72d74873360e5f60bbf87613b194ed4b7dbe3429686df64ca508528677dada79`.
