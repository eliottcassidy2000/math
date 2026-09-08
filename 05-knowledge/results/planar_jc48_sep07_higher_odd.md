# Higher finite odd cusps in the birational (4,6) normalization class

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[classification and transport audit](planar_jc48_sep07_higher_odd_audit.md)
and the separate [exact representative braid audit](planar_jc48_sep07_higher_braid_audit.md)
both pass. Every curve in the declared class is excluded as whole Keller
nonproperness support. Planar JC remains OPEN.

## 1. Exact class, inheritance, and the three remaining types

Consider an irreducible complex affine curve with birational polynomial
normalization `(U(t),V(t))`, `deg U=4`, `deg V=6`. Its finite singularities
are exactly one unibranch cusp analytically equivalent to `v^2=u^m`,
where `m>=7` is odd, and `N>=2` ordinary nodes. The classification below
has precisely these possibilities:

| Finite cusp | Infinity cusp | Number of finite nodes |
|---|---|---:|
| `(2,7)` | `(2,7)` | 4 |
| `(2,7)` | `(2,9)` | 3 |
| `(2,9)` | `(2,7)` | 3 |

In particular **no finite cusp `(2,m)` with odd `m>=11` occurs in this
degree pair**. That exclusion follows from a nonvanishing local coefficient,
not only from the genus bound.

The closest mechanism is the audited
[(2,5) classification and connected-family transfer](planar_jc48_sep07_twofive_sextics.md),
with its [full incidence/topology audit](planar_jc48_sep07_twofive_sextics_audit.md).
The actual consumer is the
[two-positive-meridian retained-sheet obstruction](planar_jc48_sep07_next_braid.md),
together with the [odd-cusp actual node ledger](planar_jc48_sep06_odd_cusp.md).
The canonical hostile is loss of an actual chart by dividing by a vanishing
coefficient; the preceding `beta=1/4` good-curve control already demonstrates
that danger. Here the corrected near miss is to discard the exceptional
`U=t^4` curve when setting the cubic coefficient of `U` to one.

The needed sidecar is the **actual coefficient family before normalization
by that cubic coefficient**. Its irreducibility makes the exceptional
curves part of the same connected good locus. The five concepts are finite
characteristic coefficients, infinity type, birational normalization,
connected good parameter sets, and positively marked affine meridians.
No literature priority claim is made. The earlier (2,5) result is not being
repeated; the new objects are the higher-cusp exhaustion and the family
that retains its exceptional chart.

## 2. Complete local exhaustion

Move the cusp parameter and image to zero, make `U` monic, and remove
the fourth-degree term of `V` by subtracting a scalar multiple of `U`.
Exactly as in the inherited classification, the normal form is

```text
U=t^4+a t^3+b t^2,
V=c t^6+e t^5+f t^3+g t^2,             c!=0.           (1)
```

These are degree-preserving affine target operations and a parameter
translation. Both linear terms vanish at a cusp. If `b=0,a!=0`, either
the multiplicity exceeds two or the local type is `(2,3)`, so this case
does not belong to the present class.

Suppose first that `b!=0`. The absence of a finite third characteristic
term forces `f=ga/b`. Set `d=g/b`. The local polynomial
`V-dU+(d/b^2)U^2` has fifth coefficient `e+2ad/b`. To obtain a cusp of
odd type at least seven, one must therefore impose

```text
e=-2ad/b.                                             (2)
```

If `a=0`, (2) makes both `U,V` even, contradicting birationality.
Thus rescaling the parameter by `a` and the two target coordinates
normalizes `a=1,c=2`. Reusing `b,d` for the resulting parameters gives

```text
U=t^4+t^3+b t^2,
V=2t^6-(2d/b)t^5+d t^3+db t^2,          b!=0.           (3)
```

After removing the sixth even term by a multiple of `U^3`, the seventh
coefficient is

```text
C7=-(6b^2+(4b+3)d)/b^3.                              (4)
```

If it is nonzero, the finite cusp is `(2,7)`. If it vanishes, necessarily

```text
d=-6b^2/(4b+3),                     b!=0,-3/4.          (5)
```

The apparent excluded denominator is not another higher cusp: at `b=-3/4`
the coefficient (4) is the nonzero constant eight for every `d`.
On (5), the next even removal has coefficient
`6(b+2)/(b^5(4b+3))` multiplying `U^4`. After it, the ninth coefficient is

```text
C9=-44/(b^2(4b+3)) != 0.                              (6)
```

The constant numerator has a direct coefficient explanation. On (5),
the removed sixth coefficient multiplying `U^3` is
`-4/(b^2(4b+3))`. Before the `U^4` removal the ninth coefficient is
`4(1+6b)/(b^2(4b+3))`; the removal contributes
`24(b+2)/(b^2(4b+3))`. Their difference is exactly (6).
The preceding coefficients through order eight all vanish. Thus (6)
gives type `(2,9)` and proves that the finite eleventh boundary does not
exist in this chart. These local subtractions are singularity tests,
not additional global normal-form operations.

The infinity chart is `z=1/t`, `X=U/V`, `Z=1/V`. Its orders are two and
six, with contact six along the original infinity line. In (3), the
seventh characteristic coefficient is

```text
I7=-(2d+3b)/(2b).                                     (7)
```

If it vanishes, `d=-3b/2`, and direct expansion gives

```text
C7=9/(2b^2) != 0,
k4=-6-24b,
[z]((Z/X^3-4)/X-k4)=7,
[z^9](Z-4X^3-k4X^4)=7/16 != 0.                       (8)
```

Hence this is exactly finite seven/infinity nine, with no infinity-eleven
specialization. On the finite-nine family (5), (7) becomes
`-9/(2(4b+3))`, which is never zero: finite nine always has infinity seven.

Finally suppose `a=b=0`. Then `U=t^4`; multiplicity two requires `g!=0`.
The fifth coefficient of `V^2-g^2U` is `2gf`, so `m>=7` forces `f=0`.
Birationality forces `e!=0`, because otherwise both coordinates are even.
After scaling `c` to two, this exceptional family is

```text
U=t^4,       V=2t^6+e t^5+g t^2,       eg!=0.           (9)
```

Here `U-V^2/g^2` has first nonzero term `-(2e/g)t^7`; the infinity
seventh coefficient is `e/2`. These curves are therefore finite seven /
infinity seven, and cannot yield a higher finite type.

The projective closure is a rational sextic with its unique infinity
preimage at `t=infinity`. Its arithmetic genus ten is exactly

```text
10=(m_finite-1)/2+(m_infinity-1)/2+N,
```

under the declared singularity inventory. This gives the node counts
in the table and completes the exhaustion.

## 3. A single irreducible family retains the exceptional chart

For finite seven/infinity seven, normalized parameters `(b,d)` alone
do not include (9). Instead keep the actual coefficient family

```text
M={2ag+e b^2=0} subset A4_(a,b,e,g),
U=t^4+a t^3+b t^2,
V=2t^6+e t^5-(eb/2)t^3+g t^2.                        (10)
```

The equation is irreducible: it is primitive and linear in `e` over
`C[a,b,g]`, with coprime coefficient `b^2` and constant term `2ag`.
On the multiplicity-two locus `b!=0` or `g!=0`, the variety is smooth,
as its derivatives with respect to `e` and `a` are respectively `b^2`
and `2g`. The two charts are explicit affine-coordinate charts with one
coordinate inverted.

On `b!=0`, (10) gives exactly (1)–(2) with `d=g/b`. The finite seventh
coefficient before normalizing `a` is

```text
C7_b=-a(6b^2+(3a^2+4b)d)/b^3.                        (11)
```

On `g!=0`, solve `a=-eb^2/(2g)`. Using `V` as the order-two local
coordinate and removing its fourth and sixth even terms gives

```text
C7_g=-e(3b^3e^2+24b^2g+16g^2)/(8g^3),
I7=(e-3a)/2=e(3b^2+2g)/(4g).                         (12)
```

On the overlap the exact relation is `C7_g=-(b/g)C7_b`; their
nonvanishing conditions agree. Hence the multiplicity-two, finite-seven,
infinity-seven locus is a Zariski open subset of this one irreducible
smooth coefficient family. At `b=0,g!=0`, (10) forces `a=0`, and (12)
reduces to `C7_g=-2e/g`, `I7=e/2`. Every exceptional curve (9) belongs
to this same open locus, as an interior parameter point.

The explicit bridge holds fixed nonzero `e,g` and varies `b`, with
`a=-eb^2/(2g)`. All coefficients of the actual polynomial pair remain
regular at `b=0`. If an exceptional curve is good, intrinsic openness
below ensures nearby nonzero `b` are good too. This does not infer good
fibres from a bare limit. Conversely the `g=0,b!=0` chart is retained;
for example `a=b=1,e=g=0` has `C7_b=-6`, `I7=-3/2`. No single chosen
coordinate chart is being substituted for the whole class.

## 4. Connected good strata and actual meridian transport

Define a good locus by requiring birational normalization, the prescribed
single finite cusp and infinity type, and otherwise only ordinary nodes.
There are three irreducible starting parameter spaces:

* The finite-seven/infinity-seven open subset of `M` in (10)–(12), which
  includes every exceptional curve (9).
* The finite-seven/infinity-nine family (3) with `d=-3b/2`, over `b!=0`.
* The finite-nine/infinity-seven family (3),(5), over `b!=0,-3/4`.

In each case its good locus is Zariski open. The inherited corrected
incidence proof applies with its hypotheses paid as follows. First remove
parameters with a second preimage of the finite cusp. Near `t=0`, either
the quadratic coefficient `b` of `U` or the quadratic coefficient `g` of
`V` is a unit, so additional preimages cannot approach that section.
The infinity preimage is unique. The extra-preimage incidence is closed
in the proper `P1` family, so its parameter image is closed.

Now resolve the prescribed cusp and marked infinity branch along their
fixed sequences. All final odd coefficients (4), (6), (8), or (12) are
units on their respective strata. The pullback centre ideals are fixed
powers of the normalization parameter times units, so the normalization
actually lifts after the first exclusion. The centres are the unique
successive branch sections; their descriptions on the `b` and `g`
charts agree on overlaps and therefore glue. Vanishing intermediate
even coefficients does not change a free centre into a satellite centre,
by the same marked-chart argument as in the audited predecessor.

Remove next the closed image of the nonimmersion locus. The lifted map
is then unramified and separated. Its diagonal in the fibre product is
both open and closed, making the off-diagonal double incidence proper
over the base; the pairwise-distinct triple incidence is proper as well.
Nontransverse pairs and triple incidences have closed algebraic parameter
images. Their removal leaves exactly the good locus. Nonbirational maps
are excluded by their nontransverse generic double pairs. The genus ledger
fixes the node count, rather than an extra finite census assumption.

Every nonempty resulting good locus is path connected. For the two
one-dimensional strata this is the complex line minus finitely many
points. For (10), each of its `b!=0` and `g!=0` charts is an affine
three-space with one coordinate inverted. A nonempty Zariski open subset
of such a chart is path connected: two points can be joined inside their
complex affine line after avoiding its finitely many intersections with
the excluded algebraic set. The two good chart loci intersect whenever
they are both nonempty, since they are open in one irreducible variety.
Their union is therefore path connected. This argument retains every
good exceptional member, not merely the dense normalized subfamily.

Over each good locus, the nodes form a finite etale multisection. Blow up
that entire multisection, without choosing global labels for nodes or
branches. Together with the prescribed resolutions this yields a smooth
proper family with the reduced total curve plus infinity divisor relative
normal crossing. Along any compact path in the good base, lift a real
tangent field in divisor charts and patch it tangent to every stratum.
Properness gives its flow. This identifies the actual affine complements:
all centres lie in the removed divisor, and the complex normal orientation
preserves positive curve meridians up to their access conjugacies.

Consequently, **one actual two-positive-meridian certificate in a stratum
supplies the same generation statement for every good member of that
stratum**. The map is a diffeomorphism of actual open complements, with
the strict curve and infinity divisor marked. It does not transport an
abstract infinity representation or a projective quotient in place of
that complement.

## 5. Audited representative inputs and the whole-support exclusion

The separate producer
[higher-cusp braid source](../../04-computation/planar_jc48_sep07_higher_braid.py)
certifies these representatives, all with `U=t^4+t^3+t^2`:

| Stratum | `V` | Normalized parameters |
|---|---|---|
| finite7 / infinity7 | `2t^6` | `b=1,d=e=0` |
| finite7 / infinity9 | `2t^6+3t^5-(3/2)t^3-(3/2)t^2` | `b=1,d=-3/2` |
| finite9 / infinity7 | `2t^6+(12/7)t^5-(6/7)t^3-(6/7)t^2` | `b=1,d=-6/7` |

Its independently audited exact geometry establishes derivative gcd `t`, no extra
cusp preimage, transverse double pairs without triples, and node counts
four, three, three. The six frozen rational path witnesses pass286548
always-active exact gates under both normal and optimized execution;
independent replays match the2920-byte output. The
[primary braid proof](planar_jc48_sep07_higher_braid.md) and its independent
audit pay the actual marked affine-group consequence.

With the three representative two-meridian certificates and their
geometric inputs independently accepted, Section 4 applies to all
three good strata. The actual whole-support node ledger with `N>=2`
gives generic degree `d<=2a`, where `a` is the generic retained-sheet
count. If two positive curve meridians generated the affine complement,
their transitive monodromy would satisfy

```text
2(d-a-1)>=d-1,
```

by the audited moved-cycle graph bound. But `d<=2a` makes the left side
at most `d-2`, a contradiction. Thus those certificates exclude
every curve in the declared class as a whole Keller nonproperness support.
This consumer does not need an unproved classification of arbitrary
transitive groups generated by longer cycles.

## 6. Exact evidence and scope

The [standalone source](../../04-computation/planar_jc48_sep07_higher_odd.py)
and [output](planar_jc48_sep07_higher_odd.out) verify the complete local
coefficient exhaustion, its nonzero ninth boundary, both exceptional and
regular coordinate charts, their overlap, the actual coefficient family,
and the three genus ledgers. They import no repository producer. The
connected-family proof is analytic/algebraic, not a finite parameter scan.

```sh
python3 -B 04-computation/planar_jc48_sep07_higher_odd.py
python3 -B -O 04-computation/planar_jc48_sep07_higher_odd.py
```

The preserved data are the actual normalization, finite and infinity
branches, node inventory and positive affine meridians. Normalizing by
`a` alone loses the exceptional chart; (10) repairs precisely that loss.
The local finite-eleventh search stops with the constant numerator `-44`.
The representative certificates and independent classification/transport
audit are now complete. Together with the previously audited ordinary-cusp
and finite-(2,5) classes, this excludes the entire birational (4,6)
normalization class with exactly one finite odd cusp (2,m), m>=3, and
at least two ordinary nodes as whole Keller nonproperness support.
Neither a general cusp classification beyond degree pair (4,6) nor JC(2)
is asserted. No theorem identifier is allocated.

The normal and optimized replays pass **65 always-active exact gates**
and are byte-identical to the frozen output. Source/output are frozen:

```text
source SHA256 edd6c8c77716bb695c4ed902586c19e68486119fddfbab1845c1bb8614515f73
output SHA256 c880a825cd246197e7d2a1d49f23d2cf0ad3b8065c0f9fa42c58f2a875f3d958
semantic SHA256 f6acc585ae29ec30e27fd861744c3b41776497ebbb4ead4d9d7a7a565f9f5c60
```

The independent analytic/source audit and actual representative path
certification both pass. Multiple cusps, fewer than two nodes, and other
degree pairs are outside this class theorem.
