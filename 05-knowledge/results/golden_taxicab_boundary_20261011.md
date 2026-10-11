# Authenticated six-colour boundaries from golden and taxicab coordinates

2026-10-11. **PROVED, scoped and self-audited:** G1–G4 and the elementary
arithmetic comparisons. **FINITE-EXACT:** complete six-port table, exact
coordinate incidences and list-count controls. **INHERITED:** the stronger
two-cube divisor theorem and Moser rosette count. **OPEN:** a universal plane
six-colouring, a seven-chromatic unit-distance witness, universal Collatz,
and planar JC(2). No prize conclusion follows from these low-chromatic gadgets.

## Inheritance and an actual next object

The immediate suppliers are [weighted boundary relations](boundary_context_weights_20261011d.md)
and [unit rhombus interfaces](unit_rhombus_interfaces_20261011c.md). The former
already proves that a six-colour equality tester at unit-separated ports
is equivalent to a finite unit-distance obstruction to six-colourability.
The abstract K6 selector is not planar unit geometry. The latter already
gives a sufficient marked two-triangular-carrier six-colouring. Neither
should be replaced by a weaker new statement disguised as a solved frontier.

The next constructive target here is weaker and fully achievable:
**a nontrivial complete six-colour boundary relation, its exact extension
weights, and actual coordinates with every unit edge authenticated**.

The closest taxicab mechanisms are
[THM-463, two-cube divisor/split-axis theorem](../../01-canon/theorems/THM-463-two-cube-representations-are-a-divisor-property-on-the-split-axis.md)
and [THM-434, Moser rosette count](../../01-canon/theorems/THM-434-moser-ladder-unit-count-12-plus-rE-t.md).
We use their actual formulas, not their broad prose about the significance
of coincident records. The corrected hostile is the old Paley1729 modular
analogy: the next Paley ratios contain inert factors, so it is not a
splitting law. The [cubic-lens recovery](../../07-reflections/the-cubic-lens-where-1729-is-a-theorem-and-where-it-is-a-coincidence.md)
and [1729 growth-law recovery](../../07-reflections/mining-1729-from-taxicab-to-the-excess-constant-kps-S134.md)
preserve that history without making it a theorem premise.

Portfolio: six-colour boundary construction / ROOT height envelope /
elliptic pencil transport. Board: **hidden versus exposed vertex, available
colour set, exact unit incidence, trace plus norm, source height, relative
differential**. The reusable move is to change which coordinate is exposed,
then recalculate the complete interface before gluing.

## G1. A seven-vertex induced gadget with a genuine six-colour relation

Let c be a hidden centre and b_1,...,b_6 six distinct points on its unit
circle, no two at distance1. The induced graph is K_(1,6). All six leaves
are ordered boundary ports. For a fixed assignment of six colour names,
let j be the number of distinct colours on the ports. Then its extension
weight is exactly

    w(b_1,...,b_6)=6−j.                                    (1)

Its Boolean boundary relation accepts precisely j<=5. It rejects exactly
the all-distinct partition, even though no two boundary vertices are
adjacent. This is the required authenticated six-colour relation.

**Proof.** The centre may take any colour absent from the boundary, and
each such choice uniquely extends the assignment. There are6−j choices.

There are Bell(6)=203 boundary partitions; 202 extend. The labelled
colouring count is

    sum_pi (6)_(number of blocks) w(pi)=6*5^6=93750.

The full203-row exact table is saved in the output. The number of boundary
assignments rejected is6!=720 out of6^6. With at most five neighbours, a
single hidden vertex can never reject an otherwise proper six-colour
boundary assignment; thus six is the exact arity threshold for this
one-hidden-vertex mechanism. No global minimum-gadget claim is made.

This reverses the earlier hidden-leaf example. Hiding a leaf merely
multiplies the count by5; hiding the centre of six exposed leaves creates
a genuine forbidden boundary pattern. Context includes the direction of
forgetting, not just the unlabelled star.

### Golden coordinates

Set phi=(1+sqrt5)/2 and let zeta have Cartesian coordinates

    zeta=(phi/2, sqrt(10−2sqrt5)/4).

It is rotation by36 degrees. Take c=0 and b_i=zeta^(i−1), i=1,...,6,
using the usual complex-plane multiplication. All are at unit distance
from c. A chord among tenth roots of unity has length1 only if its angle
is60 or300 degrees, neither a multiple of36. Thus no additional unit
edges occur. The checker verifies all45 pairs of the *full* ten-point
decagon by exact algebraic-number arithmetic; any six-point subset works.

The golden ratio is load-bearing in this coordinate realization. The
boundary relation itself depends only on the resulting incidence graph;
it does not require golden arithmetic in every realization.

### Taxicab coordinates

The eight positive-sector solutions of a^2+ab+b^2=1729 are

    (3,40), (8,37), (15,32), (23,25),
    (25,23), (32,15), (37,8), (40,3).                        (2)

Map each to the unit vector

    b(a,b)=((a+b/2)/sqrt1729, b sqrt3/(2sqrt1729)).          (3)

All lie strictly within one sixty-degree sector, so their pairwise
angular differences are less than60 degrees: all pairwise distances are
strictly less than1. Alternatively the exact squared distance is

    [ (a−c)^2+(a−c)(b−d)+(b−d)^2 ]/1729,

which the checker verifies lies strictly between0 and1 for every pair.
The centre and any six of(2) therefore realize exactly the same G1 relation.

This is a concrete map from the inherited taxicab norm data to an actual
geometric boundary object. THM-434's fuller construction uses
omega_t=((2t−1)+i sqrt(4t−1))/(2t) and alpha(1−omega_t) for N(alpha)=t.
At t=1729 these are a common rotation of(3), because
|1−omega_t|=1/sqrt(t). Thus our eight chosen directions sit inside the
inherited60-unit rosette after a rotation. We do not reprove the full
rosette classification.

The loss is explicit: once we pass to the star relation, its norm labels,
eight-versus-six choice, and arithmetic angles disappear. Retain(2) and(3)
if future gluing depends on coordinates. High unit-vector count by itself
does not force high chromatic number; this star remains bipartite.

## G2. Joining hidden centres creates a joint available-colour constraint

Take hidden unit-separated centres u=(0,0), v=(1,0). Give each five
private boundary leaves, with no edges except the central edge and the
ten leaf-centre edges. The following *induced* realization is exact:

    left leaves:  u+((1-s^2)/(1+s^2), 2s/(1+s^2)), s=2,3,4,5,6;
    right leaves: v+((1-s^2)/(1+s^2), 2s/(1+s^2)), s=7,8,9,10,11.

The saved output supplies all12 rational points and all11 edges. Every
pair is tested for equality and unit distance; there are no accidental
unit edges. Distinct vertices are not allowed to collapse.

Let S and T be the sets of colours used on the two leaf groups, and
A=[6]\S, B=[6]\T the available centre colours. Its exact weight is

    w=|A||B|−|A intersection B|
     =(6−|S|)(6−|T|)−(6−|S union T|).                      (4)

**Proof.** Choose one colour in each available set, deleting the pairs
with equal colours because uv is an edge. The support is exactly

    A nonempty, B nonempty, |A union B|>=2.                 (5)

For this five-leaf geometry A and B are automatically nonempty. It rejects
precisely the assignments in which each group uses five distinct colours
and the two groups use the *same* five-colour set. If their omitted colours
differ, there is exactly one extension. Thus it implements a conditional
inequality of the missing colours, without pretending that its leaves
are forced to use five distinct colours.

The number of colours globally present on the ten ports is insufficient:
two assignments using five colours can have extension weight0 or positive.
The intersection data in(4) are the missing coordinate. There are120
forbidden canonical boundary partitions (the bijections between two
five-element groups), or86400 labelled assignments. The graph itself is
a tree and two-colourable. Its useful constraint concerns *precoloured
boundary extensions*, not its unconditioned chromatic number.

## G3. A geometric compiler for available-colour relations

Start with any finite unit-distance graph on distinct real-algebraic
points, whose actual unit incidences are known. At any vertex z, one may
adjoin any finite number of private unit leaves, with no other new unit
edges, preserving algebraicity of all coordinates.

**Proof.** Candidates are z+((1-s^2)/(1+s^2),2s/(1+s^2)) for rational s.
There are infinitely many. For each other existing vertex w, its unit
circle intersects z's circle in at most two points. Coincident existing
points are another finite set to avoid. Thus only finitely many candidates
are forbidden. Choose and check a fresh rational parameter, then repeat.
The argument also avoids edges to earlier leaves. It does not move or
identify old ports. QED.

Expose the leaves and hide the old core. A boundary assignment gives each
core vertex an available list A_v. Its extension weight is *exactly* the
number of proper list-colourings of the old core with these lists. This
is a bijection of extension sets, not a probabilistic independence claim.

For a unit equilateral triangle core, write the three lists as A,B,C.
The exact formula is

    w=|A||B||C|−|A intersect B||C|−|A intersect C||B|
      −|B intersect C||A|+2|A intersect B intersect C|.       (6)

It is positive exactly when every list is nonempty, each pair has union
size at least2, and the full union has size at least3. These are the
three-set distinct-representative conditions; they also follow directly
by inspecting the at most three selected colours. The checker compares
(6) with an independent literal choice count for all64^3 list triples.

Five private leaves realize any nonempty proper available subset of[6],
by repeating a colour as needed on the leaves. They do not realize the
empty forbidden set; zero private leaves would give the full available
palette. Six private leaves also permit an empty available list. These
scope distinctions matter when using the broader64^3 algebraic test bank.

This compiler gives actual geometric relations for every supplied core,
but it does not force the useful leaf assignments to occur in an arbitrary
external colouring. To compose geometric gadgets, one still needs rigid
motions matching their named ports and a fresh audit of all cross-edges.
An arbitrary abstract identification of ports need not be geometrically
realizable.

### A 25-vertex golden gadget detects a whole-cycle obstruction

The compiler was also executed on a regular *unit* pentagon. In the
golden coordinates above, its five hidden core points are

    c_i = 2 zeta^(2i)/sqrt(10−2sqrt5),        i=0,...,4.

Their induced unit graph is exactly C5; its diagonals have length phi.
At each c_i the deterministic rational-parameter search adjoins four
private leaves. All25 exact algebraic points, the20 chosen parameters,
and the25 actual unit edges are saved; every one of300 point pairs was
checked. Its ordered boundary is the20 leaves.

Let A_i be the available list at c_i. With J the6-by6 all-ones matrix
and D_A its diagonal membership matrix, the complete extension weight is

    trace[D_(A_0)(J−Id)D_(A_1)(J−Id)...D_(A_4)(J−Id)].      (6b)

Expanding the trace sums exactly the choices around the closed cycle;
each diagonal enforces its list and each J−Id its unit edge. This is a
formula for *every* boundary assignment, not just a sampled sufficient
condition. Four leaves at each centre realize any available list of
size2..5.

If all A_i={0,1}, every individual edge has two proper choices and every
proper path fragment is colourable, but the whole weight is0. If all
A_i={0,1,2}, the weight is30. Explicitly assign each group of four leaves
the colours (2,3,4,5) for the zero case and (3,4,5,3) for the positive case.
Thus the gadget retains a genuinely global
boundary constraint missed by checking each edge in isolation. The
1024-case list bank compares(6b), evaluated by a separate path recurrence,
with literal enumeration using lists {0,1},{0,2},{1,2},{0,1,2}.

This is an authenticated geometric instance of the prior
[odd-cycle/parity mechanism](plane_odd_cycle_channels_20261011b.md), not a
new proof of the general odd-cycle theorem. The induced graph has
chromatic number3. What is new to this continuation is the explicit
geometric20-port relation and its extension weights, ready for a context
that really supplies those lists.

The next constructive question is a *pair of
compatible geometric boundary layouts* whose complete relations conflict,
or whose weights combine to constrain an unresolved core.

## G4. What cubes, phi, and1729 actually share

For any two numbers alpha,beta, retain both trace s=alpha+beta and
product p=alpha beta. Then

    alpha^3+beta^3=s^3−3ps.                                 (7)

The pair(s,p) recovers the unordered roots of z^2−sz+p. The cube sum alone
does not. This exact loss-and-recovery map is the common structure.

For the golden conjugates phi and psi=1−phi, (s,p)=(1,−1), so their cubic
trace is4; specifically phi^3=2+sqrt5 and psi^3=2−sqrt5. For the two
taxicab representations,

    (1,12): (s,p)=(13,12),
    (9,10): (s,p)=(19,90),

and both values in(7) are1729. The coincidence is a collision after
forgetting(s,p), not an identity of the underlying quadratic fields.

In fact there is **no algebraic integer alpha in Q(sqrt5) with
Tr(alpha^3)=1729**. To see this, its trace s and norm p are integers and
Delta=s^2−4p=5b^2 for an integer b. Also

    1729=s(s^2+3Delta)/4,

so s is a positive divisor of1729 and s^3<=4*1729. Only s=1,7,13,19
remain. The required discriminants (4*1729−s^3)/(3s) are respectively

    2305, 313, 121, 1.

None is five times an integer square. The ordinary taxicab roots with
traces13 and19 are distinct rationals, so they are not exchanged by
golden conjugation. This is a useful hostile to an overly broad golden
explanation, while preserving the exact common trace identity(7).

The inherited arithmetic-to-geometry bridge is instead

    x^3+y^3=(x+y)(x^2−xy+y^2),

whose cofactor is the Eisenstein norm. At1729 the good divisors13,19 give
cofactors133=7*19 and91=7*13. THM-463 already proves the full primitive
splitting statement and its much larger exact audit; our gadget is a
consumer of that stronger theorem.

There is a second exact golden comparison. In Z[phi]/(1729), phi has order

    lcm(16,28,18)=1008=16*63,

from its orders modulo7,13,19. The first two primes are inert for the
golden polynomial, while19 splits (roots5 and15). All three split for the
Eisenstein polynomial. Thus the same prime triple has different roles in
the two rings. Cubing phi gives order336. These are verified CRT clocks;
no dynamical coverage, Fermat-primality test or Collatz conclusion is
inferred from the numerical appearance of63.

The earlier [golden solid note](collatz_coincidence_atlas_20260927.md)
already gives the cube-in-dodecahedron edge ratio phi and volume ratio
1−sqrt5/5, plus the antipodal graph quotients. Those are three-dimensional
geometric facts; projection of a solid into the plane need not preserve
its unit edges. Here the golden decagon supplies a direct plane realization
instead of using such a projection.

Finally, the two-cube curve has a precise elliptic encoding. For
x^3+y^3=N and s=x+y nonzero set

    U=12N/s,    V=36N(x−y)/s.

Then V^2=U^3−432N^2. At N=1729 the two positive representations yield
(1596,52668) and(1092,3276), with a sign change of V if the summands are
reversed. These are points on a constant j=0 cubic. The planar source
thread's Legendre j-invariant varies with its parameter; sharing genus1
does not provide a pencil isomorphism. The [new pencil test](legendre_pencil_normalizer_20261011.md)
retains both j and the relative differential and rules out every nontrivial
base permutation under symplecticity.

## What was gained and what must still be supplied

| Source and operation | Preserved target | Information lost | Required next sidecar |
|---|---|---|---|
| Golden or norm coordinates to a unit star | Exact six-colour extension relation | Angles, norm labels | Marked coordinates for geometric gluing |
| Hide one centre, then a joined pair | Complete extension count | Internal colour choices | Available-set intersections; explicit lift |
| Pair of cube roots to trace and norm | Unordered roots and cube sum | Ordering | A sign for the difference |
| ROOT word to height envelope | Valid unconditional upper bound | Actual ROOT support | Independent lower readout |
| Elliptic pencil to j | Generic isomorphism class over an algebraic closure | Differential scale and descent | Area-form compatibility |

The [Collatz continuation](collatz_height_envelope_20261011.md) proves an
exact height envelope and a source-only sibling discount, lowering its
remaining error to O(q^−7/(log q)^2). Its independent positive lower readout
is still missing. The analogous geometric progress is an authenticated
interface; its decisive compatible external context is still missing.
The analogy concerns retained constraints and actual witnesses, not an
asserted reduction between Collatz and plane colouring.

## Incoming work changes the construction frontier

Before filing, incoming commit c1e7e404ed was read in full across its five
mathematical notes. Its [quantized equality guards](unit_equality_quantized_guards_20261011e.md)
recover a stronger **necessary27-vertex floor** for the unit-separated
six-equality target, using two cited edge bounds. This supersedes the
earlier elementary eleven-vertex filter as a search floor. Our25-vertex
pentagon gadget does not approach or contradict that target: it has a
twenty-port relation and chromatic number3, not two ports forced equal
under every six-colouring.

The [amplifier theorem](existence_count_amplifiers_20261011e.md) closes
finite two-terminal serial/parallel repetition with private interiors:
positive separating support stays positive at palette at least3. The
new list gadgets therefore retain *all* ports and their shared core;
reducing them prematurely to a favourable two-port ratio would discard
their useful constraints. The exact zero in the pentagon's conditioned
relation is an available building block, not an unconditional two-port
equality gadget.

The [Moser-field certificate](unit_equality_candidates_20261011e.md)
implements an inherited all-size four-colouring obstruction to a large
rotation/sum bank. Merely enlarging that bank cannot produce a seven-colour
witness. We do not claim that leaving that particular field suffices:
the explicitly constructed golden and taxicab graphs above are still
two- or three-colourable. A future high-chromatic core needs a separate
field screen and chromatic audit before private-leaf compilation.

Finally, [conditional runner survival](lonely_runner_conditional_survival_20261011e.md)
and the [incoming synthesis](existence_count_equality_and_runner_20261011e.md)
keep exact conditional counts, the original clock, and closed boundary
points. They sharpen the common lesson without identifying the problems:
a fixed graph has a finite count quantum; a fixed runner source has a
complete clock; our present Collatz height estimate has no bound on the
unknown ROOT word length and therefore supplies no analogous finite
denominator quantum. This is the remaining missing coordinate, rather
than a license to turn a small residual into universal positivity.

## Reproduction and verification scope

[Script](../../04-computation/experiments/golden_taxicab_boundary_20261011.py),
[JSON including the full six-port table](golden_taxicab_boundary_20261011.json).

    python3 -B 04-computation/experiments/golden_taxicab_boundary_20261011.py
    python3 -B -O 04-computation/experiments/golden_taxicab_boundary_20261011.py

Controls authenticate all pairs of the taxicab sector, all golden decagon
pairs, and all pairs of the rational double-star. The203 star partitions
are complete; all4096 available-mask pairs and262144 triangle list triples
are compared with separate counting paths. The25-vertex pentagon gadget
has a complete incidence audit and1024 list tests plus two explicit
positive/hostile patterns. Empty/full lists are hostile
controls beyond the five-leaf realizability range. The all-source proofs
and finite Collatz universe are recorded in their companion note. These
are exact checks and self-audits, not an independent peer audit or Lean proof.
The joint run has **573,281 always-active checks**. Ordinary and optimized
Python produce byte-identical JSON. The incoming scripts' independent
audits are inherited evidence for those notes, not audits of this work.
