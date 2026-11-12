# Stable port obstructions and the forbidden Hamiltonian spectrum

**PROVED explicit obstructions / FINITE-EXACT order-nine classification;
Chan–Pak OPEN. Research continues in this packet.**

Two different proposed composition principles fail on actual strong unit
seeds. Individual stable ports need not be jointly stable, first at order
nine in the complete strong bank. More decisively, even stable ports in
different contractible blocks need not be jointly stable: an order-eleven
seed contracts to two nonisomorphic terminal births. Thus global ordinary
stable-twin confluence is false. These are failures of a seed-generation
strategy, not counterexamples to the Hamiltonian spectrum conjecture.

## 1. Exact objects and conventions

Rows encode outgoing neighbours, vertices start at zero. H counts directed
Hamiltonian vertex sequences. For a selected set O, W is the number of
ordered two-path covers with the end of the first path and start of the
second both selected. The seed defect is delta=H-W. A unit seed has
delta=1. U is the inherited selected tail slope. A port is unselected;
its length m means replacement by a transitive module TT_m. Stability
means delta=1 for every positive integer length. Joint stability permits
all indicated lengths to vary independently.

An eligible contraction removes one vertex from an unselected transitive
twin block containing a stable port. A birth has a stable port and no
eligible contraction. We count ordinary marked isomorphism classes; converse
classes remain distinct. Canonical unselected twin quotients retain frozen
block lengths, as in the previous atlas.

## 2. J1 The first physical mixed obstruction

**PROVED, with a complete finite skeleton certificate.** Take

```
rows = [16,1,67,7,142,31,187,303,127]
O = {1,3}; ports = {5,8}.
```

This strong marked tournament has H=265, W=264, U=908. If the port lengths
are a,b, its exact defect is

```
delta(a,b) = 1 - (2^a-2)(2^b-2).                         (J1)
```

The signed Stirling skeleton has only coefficients a_(1,1)=1 and
a_(2,2)=-1. The lift factor for r visits to TT_m is r! S(m,r), so
F_1=1 and F_2=2^m-2 prove (J1) at every positive length. The complete
skeleton is computed by two different algorithms: signed words with one
selected cut, and separate path-end/path-start convolution. Direct
permutation counting at the base gives 265 and264; literal graph counts
on the entire 1..4 by1..4 grid agree. See [obstruction.py](obstruction.py)
and [its exact receipt](obstruction.json.gz).

Both axes are free. The exact positive unit locus is a=1 or b=1;
at (2,2), H=1177, W=1180 and delta=-3. Its stable-face complex has
two singleton facets and obstruction ideal <xy>. This realizes physically
the kind of mixed obstruction previously available only as a formal tensor.
The second marked type in the order-nine census is its converse.

The full slope for the displayed orientation is

```
U(a,b)=38+15*2^b+14*2^a+122*2^(a+b)+54*3^a*2^b.
U(a,1)=68+258*2^a+108*3^a;
U(1,b)=66+421*2^b.
```

These give valid seed-preserving rays; the two-dimensional grid is invalid.
The failure is the omitted mixed coefficient, not a parity or primality
obstruction.

## 3. J2 Contractible blocks can fail to commute

**PROVED explicit counterexample; minimality is still under audit.** Start
from the following marked quotient, which is not itself a unit seed:

```
Q = [2,100,49,7,75,409,173,31,223]
O = {1,4,5}; ports = {0,8}.
```

Its complete signed skeleton is

```
5 F_1(a)F_1(b)-2 F_1(a)F_2(b)-2 F_2(a)F_1(b)+F_2(a)F_2(b),
```

so its exact defect is

```
delta(a,b)=1+(2^a-4)(2^b-4).                             (J2)
```

Both directions are individually stable at a=b=2. The literal order-eleven
unit seed is

```
rows = [6,4,200,99,15,151,1843,347,63,1471,447]
selected mask = 100
stable ports = {0,1,9,10}
contractible blocks = {0,1}, {9,10}
H=2729; W=2728; U=22400.
```

Contracting the first or second block gives lengths (1,2) or (2,1).
Both are strong order-ten unit seeds with H=1347. In the first, only
the singleton port0 is stable; in the second, only the singleton port9
is stable. Neither has an eligible contraction. The two terminal births
are nonisomorphic and are converses of one another. Double contraction
would give delta(1,1)=5, so the local diamond cannot close. Joint growth
to (3,3) instead gives delta=17.

This identifies the missing coordinate: a stable direction can depend
on an unrelated block having the particular frozen length2. The positive
length identity along either line does not justify moving the other line.
The unit locus is exactly a=2 or b=2. It is two crossing rays, each with
a different terminal birth. Uniqueness of an ordinary birth root is false.
The previous conditional confluence theorem remains logically valid; its
universal pair-stability hypothesis is refuted.

## 4. Complete order-nine census

**FINITE-EXACT.** The universe is all191,536 ordinary tournament classes,
all nonempty selected sets, with actual marked automorphism orbits. The
source is the inherited full order-nine catalogue, checksum retained in
[unit_nine.json.gz](unit_nine.json.gz). No converse merging occurs.

There are178,133 strong tournament classes,15,427 strong unit-source
classes and16,443 strong marked unit seeds. Exactly324 have a stable port.
Every stable set is acyclic. Port counts are certified through the proved
run bound R=1+min(d+,d-); a full positive-length identity follows from R
successful lengths, rather than from the two-clone screen.

| Stable ports |0|1|2|3|4|5|7|
|---|---:|---:|---:|---:|---:|---:|---:|
| Marked seeds |16119|150|90|56|23|4|1|

There are202 births and122 descendants. The122 descendants agree by full
marked canonical keys with the independently generated previous layer.
The new birth ranks are150 with one stable port,44 with two, and8 with
three. Of the44 two-port births,42 are jointly stable and2 are the converse
pair (J1). All eight three-port births are jointly stable. Thus four
independent birth clocks first occur at order at least ten; Session288's
explicit order-ten four-clock core attains that minimum.

The complete through-nine birth sequence is
0,0,1,0,4,8,25,60,202, with total300.
Eighty port presentations pass the length-two screen but fail a later
length. They are retained as hostile controls, not counted as stable.

## 5. J3 Exact union of the 300 birth families

**PROVED all-length family count, based on the finite complete birth bank.**
The maximal stable faces give302 cylinders on276 marked twin quotients.
All quotient automorphisms are trivial. No cylinders from different births
intersect at any positive lengths. The two cylinders within each (J1)
birth intersect only at that birth. Fixed-length conflicts certify the
disjointness, without a length cutoff.

With G_8 the prior98-birth generating function, the exact new union is

```
G_9(z)=G_8(z)+154*z^9/(1-z)+42*z^9/(1-z)^2
                    +8*z^9/(1-z)^3-2*z^9.
```

Its order-nine coefficient is324, and for every n>=10 its coefficient is
4n^2+n-7. This is the union under the initially certified stable faces;
newly emergent ports are not included. It does not enumerate every larger
unit seed or every tournament. [birth_union.py](birth_union.py) checks the
canonical cylinder intersections and the exact rational-function identity.

## 6. Forbidden values and current research boundary

Let F_n be the positive odd integers not attained by a tournament of at
most n vertices. Padding makes F_(n+1) a subset of F_n. The target is
intersection_n F_n={7,21}; every finite F_n is still infinite above the
finite maximum path count. Convergence must therefore mean stabilization
on each finite value interval, not equality to a finite set at one order.

The prime and7-times-prime reduction is inherited from
[THM-4094, Hamiltonian matching deficit and two-prime-lane completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
The exact minimal-order strong-component recursion is inherited from
[Session224, first-birth fibers](../session224_first_occurrence_20261107/proof/PROOF.md).
Neither scalar factorization nor a prime-wheel survivor supplies the
missing physical endpoint selector. The99*2^k+23 wheel remains a search
policy with explicit divisibility certificates, not a completeness proof.

The spectrum conjecture is Conjecture6.2 of
[Chan and Pak](https://arxiv.org/html/2308.10214v2#S6.SS5), who credit
bof and Gordon Royle for the2016 formulation. It remains OPEN here.

## Reproduction and dependencies

From the repository root, run the Python files in this packet in order:
census_nine.py unit; census_nine.py stable; joint_nine.py;
obstruction.py; birth_union.py. Normal and optimized replays use separate
receipts. The order-ten and two-block probes are continuing computations;
their final scope and audits will be added before close-out.

Inherited mechanisms: the previous atlas's endpoint tables and contraction
lemma; Session288's selected-endpoint run bound and signed acyclic-port
tensor. The new unrestricted and contractible counterexamples are checked
with independent split-cover kernels and literal graph counts. Finite
results are not promoted into a universal coverage theorem.
