# Phase transport, six signed frames, and the different ABC maps

2026-10-26. **PROVED** elementary deductions F1-F4 and the stated arithmetic
continuation failures. **INHERITED PROVED/CITED** carousel decomposition,
normal-basis existence and the sharp 4p theorem. **FINITE-EXACT** receipts
in controls-normal.json. No claim of ABC, Collatz, or perfect-number closure.

## 1. Retain phase as marked ancestry when the line leaves the family

For `T_w=Q_m[A_(w_i)]`, odd m>=3, the inherited
[THM-4804 normal-basis dictionary](../../../01-canon/theorems/THM-4804-complement-lines-shatter-every-module-and-the-carousel-fibre-as-frobenius-orbits-of-a-normal-basis.md)
uses `beta_i=F^i(alpha)` in `F_(2^m)` and `x_w=sum_i w_i beta_i`.
The least word period d is the field degree and

    Aut(T_w)=C3^m semidirect C_(m/d).

The [previous marked-port law A1](../rooted_lines_frobenius_sheets_20261026/ARITHMETIC.md)
already shows that `(x_w,beta_i)` has orbit length m: marking a block kills
the residual cyclic factor. There are d block-rooted types, and d each of
singleton-rooted and triangle-vertex-rooted types. Rooting a triangle vertex
also removes its internal factor three. A full ordered path has trivial
stabilizer. These are inherited results, not newly inferred from the census.

**F1 (PROVED, phase transport with retained ancestry).** Name the ancestral
block map `c:V->Z/m` and let the active path start in block i. A rotation-free
coordinate is

    z=F^(-i)(x_w),       r_j=c(P_j)-i in Z/m.

Rotating decoration, owners and path together preserves this record. If a
legal re-marking Q starts in block i+delta, then

    z'=F^(-delta)(z),    r'_j=c(Q_j)-i-delta.                  (F1)

This follows by substitution. Along a sequence of such frame changes the
exponents add; the field element sees delta modulo d, while the normal-basis
port and owner labels retain its value modulo m. Dropping either the owner
map or the common frame can therefore lose relative phase.

Actual complement does not rotate w or apply field Frobenius. It changes
the tournament relation and retains its vertex/path identities. Once the
target leaves the decoration family, c is **ancestral provenance**, not a
claim that its fibres remain modules. The proved module-shattering and
uniform-source flux results forbid that interpretation. A target-only
classifier may drop the ancestry; a process asked to track original phase
must keep it. Re-marking then transports c on actual vertices by (F1).

For two separately classified cyclic decorations of periods d,e, the
incoming [Session68 aligned-phase theorem, section6](../session68_tournament_composition_classification_20261026/session68-arithmetic/REPORT.md)
gives gcd(d,e) joint phases of period lcm(d,e). A marked block is its e=m
case. This is the cyclic instance of the relative-frame double cosets in
README C1, not a replacement for their noncyclic alignment information.

**Incoming noncyclic extension (INHERITED PROVED).** The final
[Session69 field classifier, sections2-3](../session69_prime_ports_power_aliquot_20261026/session69-tournament/REPORT.md)
uses the full Paley7 action `GammaL1(8)=C7 semidirect C3`. Write its binary
decoration as `w(x)=a+Tr(bx+cx^3)`. For b nonzero the ordinary invariant is
`c/b^3` modulo Frobenius, with the b=0 cases treated separately. This gives
12 ordinary types; quotienting by rotations alone gives20 and loses the
larger symmetry. Marking slot z gives `(a,B,C)=(a,bz,cz^3)` modulo simultaneous
Frobenius, with point stabilizer order `3/[F2(B,C):F2]`. Thus **joint** field
degree, rather than the degree of either marginal alone, controls the
remaining symmetry. There are48 block-pointed types and96 vertex-rooted
types. These counts are inherited from the separately audited incoming
packet, not rerun by our checker.

Its marked modular flip `(a,B,C)->(a+1,B+1,C+1)` is a closed field-coordinate
operation; it is not actual path complement. The same packet also gives
two distinct Paley28 classes with equal H, full automorphism order and
complete raw flux histogram. This strengthens the 81/27 hostile: correcting
the normalization still does not recover the tournament or its line data.

## 2. The 81/27 witness and sharp 4p threshold

The controls construct `C3[A+,A+,A+]` and `C3[A-,A+,A+]`, where
`A+=C3->1` and A- is its converse. Both have exactly the same 67 labelled
cyclic triples and H=462915. Their automorphism orders are 81 and 27, so
their numbers of presentations are 5715 and 17145. The structured groups
are complete because the canonical four-vertex blocks and their quotient
are recoverable, as proved in the inherited carousel classification.

The decoration words 000 and 100 have degrees 1 and 3. Their block-rooted
type counts are 1 and 3; their full vertex-rooted type counts are 2 and 6.
The root stabilizer orders are respectively 27 for singleton roots and
9 for triangle roots, on both sides. What changes is the number of root
orbits. This explains how the same local rooted stabilizers coexist with
the different unrooted normalization.

The independently recomputed endpoint path matrix gives these nonzero
labelled-path flux counts for `Delta c3=score(first)-score(last)-1`:

| Delta c3 | order81 source | order27 source |
|---:|---:|---:|
|-4|0|4293|
|-3|73710|44064|
|-2|0|86535|
|-1|267624|132084|
|0|0|120258|
|1|121581|64017|
|2|0|11664|

Thus matching triangle fields, H and scalar cover counts cannot replace
actual line responses. Division by 81 or 27 gives orbit counts for these
flux classes, not automatically degrees of converse-merged nodes.

The [sharp fibre theorem F2](../tournament_composition_classification_20261026/fibre/RESULT.md)
keeps the p-part of automorphism order fixed in a labelled C3 fibre below
4p, for odd prime p. At m=p a degree jump 1->p changes one child chirality
and removes a factor p, realizing the sharp bound. This is an **unmarked
type-change cost**. Marking a vertex or retaining an alignment is a different
operation; its stabilizer can fall below 4p. C4's linear growth in required
alignment bits concerns continuation information, not the p-part threshold.

## 3. F2: a genuine S3 bridge between two six-state carriers

In `F8=F2[t]/(t^3+t+1)` choose normal alpha=1+t (bit value3). Its conjugates
are beta=(3,5,7), with sum1. Write e_i for the binary coordinate vector.
The six nonconstant words are `{e_i,1+e_i:i in Z/3}`. Define

    F(e_i)=e_(i+1),   R(e_i)=e_(-i),   S(w)=1+R(w).

Then `F^3=S^2=1` and `SFS=F^-1`. F has two three-cycles, and S exchanges
them, so this is the regular six-point S3 action. In field coordinates F
is Frobenius; R is linear reflection and S is affine reflection. S is not
a field automorphism (it sends 0 to1).

**F2 (PROVED).** A signed, path-marked triangle frame is `(epsilon,i)`,
epsilon=+/-1, with path `(i,i+epsilon,i+2epsilon)` in the corresponding
orientation of C3. The map

    (+1,i) -> e_i,             (-1,i) -> 1+e_i                 (F2)

is S3-equivariant for vertex rotation `i->i+1` and vertex reflection
`(epsilon,i)->(-epsilon,-i)`. This is immediate on the two generators.
It realizes, with Frobenius coordinates, the signed-frame S3 carrier
already proved in [Session25 section3](../session25_spectrum_wiring_20261021/session25-wiring/RESULT.md).

The distinction from six equally many states is testable: plain word
complement commutes with F and gives C3 x C2, a cyclic order-six action.
Affine reflection S instead inverts F and gives S3. Retaining the operation
distinguishes these two six-state structures.

Nor does (F2) transport the actual metagraph complement. Complementing a
path-marked triangle produces a transitive tournament; complementing a
uniform 12-vertex decoration leaves its triangle-count family. The map
preserves the specified relabelling action, not the complement operation.

**Normal versus degree-three boundary.** All six nonconstant field elements
have degree three, but only the three single-bit elements form normal
conjugate triples in this chosen basis. The three two-bit elements have
trace zero and are dependent. This matches the incoming
[marked 1-2-4 audit](../aliquot_power_124_20261026/tournaments/TOURNAMENTS_124.md):
an exponent Frobenius orbit can be a rank-two Fano line. Degree three alone
is not the normal-basis hypothesis. The affine S3 bridge may exchange
normal and nonnormal orbits.

**F3 (PROVED).** Even the complete distribution of Hamming weight can miss
Frobenius degree. At m=9, words `100100100` and `111000000` both have weight
three, but degrees3 and9 and residual cyclic orders3 and1. All moments of
their weight agree. Full finite-moment inversion recovers the weight law;
it does not recover the labelled word from that coarser random variable.

## 4. Incoming arithmetic resolves the two ABC maps

The repository advanced during this session. The audited
[Session69 power report, section6a](../session69_prime_ports_power_aliquot_20261026/session69-power/REPORT.md)
identifies the requested power/doubling map. For n=2^a u, u odd,
`a>=1` and `u>2^a`, its square has later odd doubling row u^2. Dividing
`u^2-n=u(u-2^a)` by u gives

    (u-2^a)+2^a=u.                                           (P)

This is a bijection to positive primitive additive triples with a
distinguished binary middle term at least2 and odd last term: n=BC.
It is lossless on its stated domain. The prime-factor radical of the
difference remains additional information for an ABC inequality.

By contrast the inherited aliquot normalization, for n>=2, is

    (n/g)+(s(n)/g)=sigma(n)/g,   g=gcd(n,s(n)).                (A)

It preserves abundance order but forgets scale. For **every** perfect n,
it is `1+1=2`, whether or not n is one of the known even examples. Its
missing gcd reconstructs n and its outgoing aliquot edge.

A third, legitimate but different comparison is to split off the largest
power of two dividing n and normalize:

    1+(u-1)=u,  n=2^a u, u>1.                               (D)

This remembers only the odd core and is invariant under doubling. We
initially considered it as a candidate interpretation; the incoming report
identifies (P) as the requested construction. Keep all three maps distinct:

| n | power/doubling (P) | aliquot (A) | odd-core split (D) |
|---:|---|---|---|
|6|1+2=3|1+1=2|1+2=3|
|12|outside its positive domain|3+4=7|1+2=3|
|28|3+4=7|1+1=2|1+6=7|
|72|1+8=9|24+41=65|1+8=9|
|120|7+8=15|1+2=3|1+14=15|

The last row gives a useful third origin of `1+2=3`: 120 is triply perfect.
An equation without its source map does not identify which phenomenon was
observed. The incoming reports also prove the special prime-square/perfect
intersection `25->6->6`; the equality of (A) at6 and28 forgets that incoming
structure. This result is credited to the incoming work, not reproved here.

**F4 (PROVED, arithmetic continuation hostiles).** The odd-core quotient is
not closed under aliquot continuation: 6 and12 both have odd core3, but
their aliquot successors 6 and16 have odd cores3 and1. The perfect quotient
under (A) is not closed under doubling followed by (A): 6 and28 first give
the same triple, while12 and56 give `3+4=7` and `7+8=15`.

For the specified even-perfect family `N=2^(r-1)q`, `q=2^r-1` prime,
the elementary divisor-sum product gives `s(2N)=2^(2r)` and hence

    (A)(2N) = q + 2^r = 2^(r+1)-1.

The erased exponent reappears in the continuation. This is the same
operation-kernel test as README C3, on a different carrier. It gives no
conjugacy between aliquot dynamics and tournament complement.

## 5. One further exact 6/28 bridge and its stopping point

The [recovered multiperfect dimension calculation](../session65_cryptid_tournament_interfaces_20261026/session65-root/multiperfect/REPORT.md)
uses `D_r=binom(2^r,2)=2^(r-1)(2^r-1)` and
`sigma(D_r)/D_r=sigma(2^r-1)/2^(r-1)`.
This yields6 and28 at r=2,3 and the triply perfect120 at r=4. It links
unordered pairs of a power-of-two set to a precise divisor-sum predicate.
An orientation of those pairs is additional data; their count does not
determine a tournament, a port response, or a relative phase.

Thus there are actual maps connecting several appearances of six: (F2)
is equivariant, (P) is a guarded bijection, (A) is a scale quotient, and
`binom(2^r,2)` is a counting identity. Their preserved predicates differ.
