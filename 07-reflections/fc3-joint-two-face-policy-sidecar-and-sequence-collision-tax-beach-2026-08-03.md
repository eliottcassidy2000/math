# FC3 joint two-face policy: the face marker is one bit, not decoration

**Status:** FINITE-EXACT discovery package, not yet a canonical theorem.

**Scope:** only the promoted support-`(1,2)`, bank-`I2` and
support-`(1,3)`, bank-`I2` physical faces of THM-3249, using lawful response
rows 2 and 10.  This does not prove `FC(3)`, `SFC(3)`, another support/bank,
or that an adaptive row choice is one scalar Gaussian-moment functional.

Exact companion:

~~~text
04-computation/fc3_joint_two_face_policy_sidecar_beach_20260803.py
05-knowledge/results/fc3_joint_two_face_policy_sidecar_beach_20260803.out
~~~

## 1. Inheritance pass and concept board

The closest proved mechanism is
[THM-3256](../01-canon/theorems/THM-3256-optimal-two-row-threshold-policy-injective-signed-trace-and-factored-distance-enumerator.md):
axis-threshold dynamic programming compresses the row-2/10 selector on the
large face while the signed edit trace remains injective.  The canonical
hostile examples are the sharp two-state fixed-gauge clutches in
[THM-3249](../01-canon/theorems/THM-3249-cross-support-upset-atlas-local-sections-and-no-constant-gauge.md)
and
[THM-3254](../01-canon/theorems/THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go.md).
The corrected near miss is the thought that transporting the same two rows
also transports a face-blind classifier.  The least-used relevant sidecar is
the persistent face/reset marker.

The live concept board was:

| object | predicate | operation | lost coordinate | cheapest test |
|---|---|---|---|---|
| tagged states `(f,n)` | lawful row 2 or 10 | forget `f` | ambient support face | search identical exclusive vectors |
| axis trees | classifier size | add one face test | margins and directions | exact leaf/depth recurrence |
| signed routes | reconstruct start | Parikh abelianization | face/reset | intersect displacement boxes |
| distance sequence | count starts by reset distance | ordinary generating function | labels and order | multiply coordinate factors |
| previous-chart bit | future switch control | initialize/transport bit | face at first decision | hostile state `(5)` |

## 2. Nineteen literal collisions refute every face-blind policy

Pad the small face into eight multiplicity coordinates by putting
`n_6=n_7=n_8=0`.  On each face, discard overlap states because either row is
lawful there and retain only exclusive row-2 or row-10 requirements.  The
small and large faces contribute respectively `71=52+19` and
`407=304+103` exclusive records.

There are exactly **19** identical padded vectors with opposite exclusive
requirements.  Every one is oriented the same way:

~~~text
support (1,2): row 2 required;
support (1,3): row 10 required.                         (1)
~~~

The smallest hostile is the literal one-pole state

~~~text
(5),       n=(0,0,0,0,1,0,0,0).                       (2)
~~~

Consequently there is no function of the untagged count vector—regardless of
circuit depth, lookup size, or test language—which selects a lawful member of
`{2,10}` on both faces.  This is stronger than failure of an axis-threshold
tree.

It also settles the first-decision meaning of a one-bit memory proposal.  A
previous-chart bit with the same initialization on both faces sees identical
input at `(2)` and cannot help.  A bit initialized differently on the two
faces is exactly a persistent face bit; it is sufficient, and zero bits are
impossible.  This does **not** solve the separate THM-3244 problem of using
previous-chart memory to realize its globally optimal two-switch routes.

This is the precise phase-marker-style connection contract:

~~~text
source:             tagged physical state (f,n)
target:             padded multiplicity vector n in N^8
map:                forget f
preserved:          the literal pole submultiset
destroyed:          response-row availability
minimal sidecar:    one face/reset bit
hostile test:       state (5).                           (3)
~~~

The resemblance to a phase marker is an information-flow analogy, not a map
from FC to LRC.

## 3. Exact policy frontier with the bit restored

The primitive tests are `n_j<=h`: 11 on the small box, 16 on the large box,
and the same 16 after padding their union.  Add the single Boolean test
`f=(1,2)`.  For a constrained record set `S`, the exact recurrence is

~~~text
L_d(S)=min_T [L_(d-1)(S intersect T)+L_(d-1)(S\T)],      (4)
~~~

with cost one for a pure leaf and infinity for a mixed set at depth zero.
The unbounded version of `(4)` gives the global leaf minimum.

| policy class | minimum depth | leaves at that depth | global minimum leaves | first depth attaining it |
|---|---:|---:|---:|---:|
| small face, axes | 4 | 8 | 8 | 4 |
| large face, axes | 5 | 15 | 15 | 5 |
| joint, axes only | impossible | — | impossible | — |
| joint, axes plus face bit | 6 | 23 | 22 | 8 |

A gate-first composition of the two separate trees has depth 6 and 23 leaves,
so it is simultaneously optimal at minimum depth.  One recovered minimum-depth
tree is more economical operationally: it first tests `n_8=0`; only that
branch asks for the face bit.  It consults the bit on all `238` small
nonreset starts but only `2,159/4,318` large starts.

There is a genuine Pareto tradeoff.  Interleaving two conditional occurrences
of the same face test reduces the leaf count from 23 to the absolute minimum
22, but the first such tree has depth 8.  It consults the face marker on only
`198/238` small and `599/4,318` large starts.  Thus one bit is enough, while
one early query is not always leaf-optimal.

## 4. Three different collision taxes

The two tagged banks contain `4,558=239+4,319` starts.  Forgetting the face at
the raw-state level identifies every small vector with a large-face vector,
leaving 4,319 vectors and a tax of 239.

For signed Q-monotone edits, let

~~~text
q_12=(1,2,1,1,1,0,0,0),
q_13=(1,0,2,1,1,1,1,1).                               (5)
~~~

On either face the signed Parikh vector is `q_f-n`, hence reconstructs `n`
once `f` is known.  Equality across faces forces

~~~text
m=n+(q_13-q_12)=n+(0,-2,1,0,0,1,1,1).                 (6)
~~~

The admissible choices are

~~~text
n_1: 5 choices,  n_2: 2,  n_3: 2,  n_4: 2,  n_5: 2,
n_6=n_7=n_8=0,                                        (7)
~~~

and neither empty-state exclusion binds.  Therefore the cross-face signed
displacement intersection has exactly

~~~text
5*2*2*2*2=80                                          (8)
~~~

elements.  The face-forgotten signed Parikh support has `4,478` elements;
the face-tagged support has all `4,558`.  This is the finite cardinality
analogue of the support/multiplicity collision tax in THM-2000/2005: labels
give a multiset of traces, and forgetting the marker charges one for each of
80 doubled values.

For one deterministic minimum-depth controller with least-coordinate
tie-breaking, ordered signed words retain more than their Parikh vectors:
only 22 cross-face pairs have identical complete words (including the two
empty reset words).  Its untagged/tagged continuation-class counts are

| retained output | untagged | face-tagged |
|---|---:|---:|
| chart | 313 | 345 |
| unsigned coordinate | 4,125 | 4,208 |
| chart plus coordinate | 4,520 | 4,548 |
| signed edit | 4,536 | 4,558 |

The number 22 is controller-dependent; the 80 signed-Parikh collision tax is
policy-independent.  Unsigned Parikh vectors are much coarser: the small and
large faces have 96 and 1,536 values, with the entire small set inside the
large one.  Face tagging therefore gives 1,632 unsigned classes, not 4,558.

## 5. The face marker is also the catalytic sequence variable

Coordinatewise multiplication gives

~~~text
D_12(z)=(1+2z+z^2+z^3)(1+z)^4(1+2z)-z^6,             (9)

D_13(z)=(1+2z+z^2+z^3)(1+z)^4
        (1+z^2)(1+z+z^2)(1+2z)^2-z^8.                (10)
~~~

The subtracted monomials remove the forbidden empty physical state.  Their
sum has the compact shared-factor form

~~~text
D_joint(z)
=(1+2z+z^2+z^3)(1+z)^4(1+2z)
  [1+(1+z^2)(1+z+z^2)(1+2z)]-z^6-z^8,                (11)
~~~

with coefficients

~~~text
(2,19,82,220,426,648,803,821,687,467,250,101,28,4).  (12)
~~~

If a downstream consumer needs to retain the face, the correct sequence
object is simply the catalytic refinement

~~~text
u D_12(z)+v D_13(z).                                  (13)
~~~

Setting `u=v=1` is lawful for distance counting and unlawful for selecting a
response row.  The same one-bit marker is therefore a control sidecar before
abelianization and a catalytic variable after it.

## 6. Hostile controls and honest frontier

The companion reconstructs all 239 small-face response vectors from the
original pinned product-Gamma formula and replays the promoted large-face
direction bank.  It verifies both complete covers, reproduces THM-3256's
large-face `5/15` optimum as a positive control, exhausts all threshold-tree
subproblems, checks both Pareto trees on every nonreset state, and compares
the route histogram with an independently multiplied coordinate-factor
polynomial.  The script contains no assertion node, float, randomness, or
scratch cache; normal and optimized Python runs must match the frozen output.

Next exact problems exposed by this package are:

1. derive the 19-vector conflict locus symbolically from response-sign
   inequalities rather than a census;
2. optimize joint policies for the other 23 cross-support covering row pairs;
3. determine whether a previous-chart bit, used for future rather than face
   memory, can attain THM-3244's sharp global two-switch bound;
4. test whether face-aware threshold depth stays bounded over further
   product-Gamma faces; and
5. carry `u D_12+v D_13` rather than its specialization whenever a later
   response consumes the face label.
