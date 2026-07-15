# The n=12 object is cover debt on periodic Kakeya combs

*codex-2026-07-15-S14, with codex-S10 relative-ramification and exact-H6 continuation.
Companion to HYP-6820, THM-815/823/836/841/844/845/847/857/858,
and the exact Fano/chi7 flood and Hamming-six contraction audits.*

## 0. Honest verdict

The requested uniform theorem—emptiness of the sporadic branch at `n=12`—is
still open.  This session did not turn the accumulated finite banks into that
global statement.  It did, however, remove several named residuals and change
the shape of the remaining object:

1. THM-847 closes the last effective-order-at-most-twelve Hamming-five
   common-sheet language. THM-858 proves that no new language appears through
   effective order 21. All-one, all-three, and mixed one-plus-three are the
   complete list through 21 and are uniformly empty at proper lift height;
   the remaining order bank is finite but not yet classified.
2. The proper scale-one and first ramified Hamming-six faces are now exact.
   THM-857's
   `580,919,164`-node recursion over all 924 deletion roots has one covering
   terminal, the doubled AP `2[12]`, and twenty loose terminals.  It subsumes
   the earlier nonprimitive contraction and fourteen-row AP-germ closure,
   removing all 909 formerly open primitive-core rows and the exceptional
   mixed-parity branch at scale one. THM-860 makes every primitive ramified
   language finite, and THM-861 evaluates all 64 `c=2` signed-cycle roots:
   the only cover is the ordinary AP `[12]`. Transport remains open from
   `c=3` onward.
3. THM-836's local `s=5` shell is exactly four congruence classes, and the full
   class `d=11 mod 52` is uniformly impossible by one explicit divisor-grid
   witness.  The classes `15,37,41 mod 52` remain open uniformly, and the
   U-independent single-numerator endpoint-grid template is now proved
   impossible in all three.
4. THM-841's proposed toothpick breakpoint tree does not exist.  The ladder has
   a much simpler all-order totient formula.
5. The Fano/`chi_7` organization of the 21 `j=4` flood bodies is exact but not
   a symmetry quotient and not a local three-needle closure.

This is progress toward the sporadic theorem, not the theorem itself.

The H5 continuation first changes “unbounded orders” from a vague size
condition to an arithmetic incidence condition, and then removes the
unboundedness.  If `S` is any nonempty colour set and
`m_i=D_i/gcd(D_i,lcm(S^c))`, common-sheet coverage forces

```text
sum_(i in S) ceil(2m_i/13)/m_i>=1.
```

Thus no order has a private maximal prime power.  After the exact order-21
census, every live order is `{2,3,7}`-smooth and any new row has
`min D_i<=21<max D_i`.  Reapply the cut across adjacent valuation levels:
the `2,3,7` exponent ranges are at most `3,2,1`.  Comparison with the order at
most 21 gives `max D_i<=10,584`, so the presentation bank is uniformly finite.
On the 185 last order-21 residuals, both raw and
complement-conditioned order tournaments are transitive despite 574 edge
flips, while a top-prime subset cut rejects every row.  The H5 analogue of
cover debt is therefore a prime-power carrier hypergraph decorated by
complement-lcm fibres and affine owner intervals.  A pair tournament sees the
relative ranking but not the fibre on which union coverage fails.

## 1. Reframing the branch

The historical definition says that a primitive tight twelve-speed set has a
top-peel eleven-core whose maximum is strictly above the extremal `1/12`
floor.  That definition is correct but describes the packet only after the
hard event has happened.  The operational definition exposed by the Hamming
and sheet proofs is more useful.

For a current core `P`, put

```text
E(P) = {t : ||pt||>1/13 for every p in P}
D_u  = {t : ||ut||<=1/13}.
```

The state is not just `P` or its maximum.  It is

```text
(literal components of E(P), remaining labelled comb progressions,
 last speed, endpoint owners/activation bits).
```

A continuation is tight exactly when the remaining periodic danger combs
cover every component of `E(P)`.  This makes the sporadic branch a labelled
component-cover problem with arithmetic motion, rather than a search for an
exceptional list of integers.

For a component `I` and chosen combs `S`, define its overlap debt

```text
omega_I(S) = sum_(u in S) |I intersect D_u|
             - |I intersect union_(u in S) D_u|.
```

Coverage requires the union term to equal `|I|`.  The discrepancy ladder
controls the first term.  What it does not see is how much of that mass is
wasted as `omega_I`.  At radius seven the mean-density coefficient changes
sign; overlap debt, unique ownership, and component incidence are therefore
the natural next potential.  This is the exact content behind the earlier
informal phrases “owner diversity” and “Kakeya needles.”

## 2. Kakeya needles, with the analogy made precise

The needles here are one-dimensional periodic combs, not arbitrary planar
Kakeya line segments.  Each speed supplies many very short parallel translates
on the circle.  The relevant question is whether six labelled combs can place
teeth so that their union covers a disconnected residual with the required
arithmetic phases.  The analogue preserves thinness plus movable incidence;
it does not import planar Kakeya dimension estimates.

At the doubled-AP equality, division by two gives

```text
core={1,...,6},             needles={7,...,12}.
```

The strict core has twelve components, total measure `27/65`, and longest
component `1/13`.  Its component--comb graph has 34 incidences.  Four
components have zero overlap debt and a unique full owner:

```text
[7/39,12/65]   -> 11,       [14/65,3/13] -> 9,
[10/13,51/65] -> 9,         [53/65,32/39] -> 11.
```

The other component degrees are `2,3,4,6`.  Thus equality is neither a uniform
density event nor a pairwise tournament event.  It is an exact tiling with four
rigid single-owner pins and eight overlap-bearing components.

There is also a clean arithmetic contraction.  A scale-one H6 packet is

```text
([12]\R) union {r+13h_r:r in R},       |R|=6, h_r>=1.
```

If it is nonprimitive, six retained labels must all be divisible by a common
integer.  There are six multiples in `[12]` only for the divisor two, forcing

```text
R={1,3,5,7,9,11},       gcd=2,       every h_r odd.
```

Writing `h_(2i-1)=2k_i+1` and dividing by two gives

```text
{1,...,6} union {6+i+13k_i:i=1,...,6}.
```

Zero positive `k_i` is `[12]`; one through five positive coordinates are the
already closed Hamming radii one through five; six positive coordinates form
one top-half H6 chamber.  Exact longest-component recursion closes that chamber
in

```text
1, 54, 3612, 130515, 2104, 2, 0
```

states by depth, with no covering prefix.  Hence `2[12]` is the only possible
tight nonprimitive H6 packet.  Before THM-857, the remaining H6 problem was
intrinsically primitive, but “primitive” had two layers.  The odd-label row has retained
core `{2,4,6,8,10,12}`: its all-odd height fibre contracts as above, whereas
mixed height parity makes the completed packet primitive.  That exceptional
fibre must not be dropped.

For every other row, let `f` count full antipodal missing pairs.  Exactly
`2f` thirteenth-grid points lie in the retained-core safe set.  At each such
cusp, compare the owner's tooth reach with the first incoming provider tooth.
The resulting weighted handoff graph is not a tournament: a pair can carry
both directions or neither because its two cusp gauges differ.  Directed-cycle
products close six of the twenty `f=3` rows; eight more force one of ten fixed
coordinates, and exact slice trees close all ten in 3,699 states.  Six
product-`1/16` rows survived that quotient.  This is the first place where a
projective AP grid cusp, rather than a runner or whole component, is the right
local vertex—but it is a lossy local chart, not the final H6 object.

THM-857 completes the picture by refusing that quotient.  Numerically order
all six proper lifts and retain the literal open residual after each insertion.
THM-815's longest-component inequality gives a finite next-speed cap at every
prefix.  A child is certified dead either when it contains a full safe tooth,
which would force the next ordered speed below the current speed, or when an
emerging exact child component already makes the cap smaller than the least
legal future lift.  This produces the all-root census

```text
nodes by depth = 924, 83,881, 8,906,315, 559,202,706,
                 12,671,505, 53,812, 21;
covers = 1; loose terminals = 20.
```

The unique cover lifts missing labels `{1,3,5,7,9,11}` to
`{14,16,18,20,22,24}`.  An independent replay reconstructs unions of closed
danger teeth from scratch and checks each expanded child endpoint-for-endpoint.
Thus the scale-one H6 object is the labelled component--future-comb action
tree, not the antipodal count, germ-cycle product, or a runner tournament.

The first ramified fibre makes that distinction sharper.  At `c=2`, common
sheet routing forces one signed-doubling Hamiltonian cycle on the six missing
labels.  There are only 64 such roots, yet the cycle does not decide their
metric fate.  Running the literal step-26 component action produces
`41,882,982` logical nodes and one cover: replace labels `7,...,12` by
`1,3,...,11`; together with the retained even speeds this is exactly `[12]`.
Thus the scale-two sporadic branch is empty, but not because the signed cycle
is inconsistent.  The arithmetic stalk is perfectly viable on all 64 roots;
the metric fibre kills 63 roots and every non-AP height.  This is the cleanest
current example of a tournament-related carrier refining the search without
being the underlying cover object.

A separate implementation supplies a smaller crosscheck of that action tree:
all 52 `f=2` roots whose exact first-speed cap is at most 156 close in
`9,888,159` states, with zero covering prefix and 1,323 independent
closed-danger reconstructions.  It is now a subcertificate, not an open
frontier.

## 3. The deep shell is an incidence problem too

In the guarded `(13d,5d)` two-sheet branch, THM-836 originally excluded the
thin shells `s=2B-13d=1,3`.  At `s=5`, the packing width is less than four and
`B=9 mod 13`.  The two coefficient-nine owners must therefore occupy the only
signed pair among the top four residue classes:

```text
{6,7},       owners={B-3,B-2},
d mod 52 in {11,15,37,41},       omitted class={3,10}.
```

This is an all-size classification of the local owner equations.  It is not
an existence theorem for full cores.

For `d=11 mod 52`, however, the unit

```text
q=5d,       p=(45d-1)/26
```

is safe for every allowed signed-complement lift, while the folded target
value is `2/5<11/13`.  This contradicts even `E_U subset H_d`, so that whole
class is uniformly empty.  At `d=15`, divisor completeness and parity support
leave 71,644 of 1,605,632 signed lifts; 3,004 cover the `q=75` grid, four cover
the `q=195` grid, and none cover both.  That last statement is finite at
`d=15`, not uniform over its congruence class.

There is also an all-size negative result about the proof architecture.  If a
single unit numerator is required to be deep for every possible choice of the
seven free lifts, a 35-speed skeleton localizes it to four intervals around
`4/13` and `9/26`.  Neither `q=5d` nor `q=13d` meets those intervals in a unit
point for `d mod 52 in {15,37,41}`.  So the successful `d=11` column has no
uniform analogue.  The missing object cannot be a scalar numerator selected
before the lift set is known; it must be lift-dependent, multi-column, or on a
different denominator.

The first multi-column repair also fails uniformly in the relaxed model.  A
unit `q=13d` column is dangerous on at least three free raw classes; a
nonforced unit `q=5d` column is dangerous on at least two, except for the
genuine family `d=41 mod 52`, `p=+/-(45d+1)/26`, whose support is exactly
`{4}`.  A labelled matching therefore chooses one relaxed shell-admissible
lift set that kills any fixed pair of unit endpoint columns; in the singleton
case the same concrete lift `B-18` kills both signs.  The all-size proof is a
432-row affine endpoint certificate with minimum strict margin `1/6545`.
This is a sharper architecture no-go, not a shell closure: THM-772 divisor
completeness and THM-803 parity support can delete relaxed choices.

The lesson matches H6: an owner inequality locates a narrow shell, but a
cross-incidence of numerator grids, signed lifts, divisor obligations, and
parity support decides it.  Ordering the four residue classes gives a
transitive tournament in either tested gauge and loses the witness.

## 4. Toothpick self-similarity was the wrong recursion

For a Farey gap with endpoint denominators `i,j` at order `k`, adjacency gives

```text
i+j>k,       min(floor(k/i),floor(k/j))=1.
```

Only one endpoint can have a nontrivial multiple chain.  Its collision with
the other primitive endpoint lies above `lambda=1/(k+1)` unless both events
are primitive.  Therefore each THM-841 gap has at most one actual kink, at
`lambda=1/(i+j)`.  There is no divisor-refined dyadic tree.

The correct self-similarity is totient/Farey accretion.  Every rung `r>=3` is
globally linear:

```text
m_r(lambda)=1-(2lambda/r) sum_(d<=floor(k/r)) phi(d)/d.
```

The analogous all-order factorial-moment formula is also linear for
`rho>=3`; only `m_1,m_2,S_2` retain the primitive left--right overlap term.
The literal Farey-end growth agrees with the toothpick totals at
`1,2,4,8 = 1,3,11,43`, creating the mirage, but at 16 it is `159`, not `171`.
The first local recurrence failure already occurs at `n=5`.

This correction is methodologically useful: the ladder is not asking for a
recursive fractal quotient.  It is asking for a one-kink Farey stalk plus a
totient multiplicity ledger.

## 5. Fano and chi7: exact organization, destroyed symmetry

The 21 `j=4` flood bodies are exactly

```text
{8,...,14} union {a,b},       1<=a<b<=7,
```

so they are the edges of `K_7`.  Labelling the seven points by the nonzero
vectors of `F_2^3` partitions those edges into the seven Fano lines.  The
Boolean `B_3` inclusion--exclusion word is

```text
++-+--+
```

and its negative masks `{3,5,6}` form one Fano line.

That organization does not quotient the numerical work.  Among all 168
elements of `GL(3,2)`, only the identity preserves any one of the exact body
observables `r(E)`, `m(E)`, or `V1(E)`.  Three one-sided endpoint needles trace
a weak switch path in the cube, visit at most four masks, and meet at most two
points of the negative Fano line.  They cannot realize the full `chi_7`
carrier locally.  Finally every flood body has empty modular obligation set,
so the obligation/lcm filter is literally silent.  The `j=4` flood tail still
requires its exact lower-tree interval geometry; the Fano probe diagnoses why
the proposed shortcut cannot prune it.

The representation-theoretic obstruction is exact.  On the 21-dimensional
edge space of `K_7`, the orbit of `r(E)` has rank seven and is the point
potential `r_{ab}=x_a+x_b`.  The orbits of `m(E)` and `V1(E)` both have full
rank 21.  Point incidence and Fano-line incidence together have rank only 13,
leaving eight invisible edge directions; explicit four-edge curls of both
`m` and `V1` are nonzero.  Thus the still-open flood geometry lives precisely
in directions that every point/Fano marginal forgets.

Two edges can nevertheless be closed without pretending the Fano organization
is a symmetry.  For flood `(5,7)`, a fixed-`E_2` lower bound closes 181,445 of
525,362 third-level nodes before evaluating `m_3`; exact `m_3` closes another
339,348.  Only 4,569 `E_3` nodes and 28,847 bottom sweeps remain, all positive,
with exact minimum `7/858`.  The same literal computation on `(6,7)` leaves
9,560 residual `E_3` nodes and 73,323 bottom sweeps, all positive, with exact
minimum `97/4004`.  These are two of the 21 floods; the other nineteen must be
transported or recomputed because their numerical weights have no nontrivial
`GL(3,2)` stabilizer.

## 6. Tournament threads after the audit

The recurring tournament failure is now precise.  Scalar gauges often produce
a transitive tournament with score histogram `0,...,m-1`, no directed cycles,
singleton SCCs, and one Hamiltonian path even when the proof verdict is hard.
This happened in the H5 contexts, the flood bodies, the shell classes, and the
component marginals.  Transitivity says the planning observable is ordered; it
does not say the cover predicate is pairwise.

Changing the vertices remains productive.  Useful vertices in this session
were:

- residual components and remaining comb obligations;
- zero-debt owner pins;
- divisor-grid numerators and signed lift choices;
- `K_7` flood edges/Fano lines;
- odd-cycle collections in the exact H-drift formula.

The H-drift result is the tournament-side version of the same warning:

```text
E[Delta H|T]
 = (a_(n-2)-(n-1)H)/C(n,2)
 = 2(K-nH)/C(n,2),
```

where `K` is a cycle-length-weighted OCF partition function.  `H` alone loses
the endpoint coefficient and even the sign of the drift at `n=6,H=23`.
The next Krylov step is exact as well.  If `S` sums the increments over all arc
flips and `A_T(x)=sum_j a_j(T)x^j`, then

```text
S a_j=(j+1)a_(j+1)+(n-j)a_(j-1)-(n-1)a_j,
S A=(1-x)((1+x)A'-(n-1)A),
S K=a_(n-3)+a_(n-2)-C(n,2)H.
```

So the functional form is a finite tridiagonal Walsh/Krylov ladder, not an
all-`n` two-scalar law.  Through `n=7`, `(H,K)` happens to determine `S K`
because the ambient kernel direction `(1,-5,4)` is absent from the realized
tournament locus.  At `n=8` it appears in exactly one of 1,727 `(H,K)` fibres:
`(H,K)=(345,2820)` has two `K` drifts, `-240/7` and `-60/7`.  This is a precise
next-size boundary, not a heuristic failure of mean reversion.  Pointwise,
`Delta_e K` is four times the gained-minus-lost odd-cycle-forest weight after
deleting the affected cycle vertices; 6,340 incidences through `n=7` verify
that current law.

The ladder itself diagonalizes.  With `m=n-1` and
`z=(1-x)/(1+x)`, the exact all-size identity is

```text
A_T(x)=sum_r H_(2r)(T)(1-x)^(2r)(1+x)^(m-2r)
      =(1+x)^m B_T(z),
S A_T=(1+x)^m[-2z B_T'(z)].
```

Thus continuous flip-sum evolution is the radial dilation
`B(z)->B(exp(-2t)z)`.  The scalar `H`-drift looks nonlinear only after the
even radial modes are projected onto one coefficient.  This is a genuine
self-similarity, but it belongs to the tournament needle stalk—not to the
Farey violation ladder, whose toothpick recursion was refuted above.

Likewise, the proposed all-`n` black self-line law loses its first untested
case at `n=8` (`404` versus `SC(8)=176`), while only the Klein-four orbit
structure survives.  Both corrections point away from scalar class counts and
toward decorated incidence/orbit packets.

The concurrent THM-855 F6 second-moment closure sharpens the positive side of
that warning.  For the score-axis energy `x`, one arc flip satisfies the exact
pointwise law

```text
E[(Delta x)^2|T]=16(n-4)x/C(n,2)+64,
```

and the stationary variance is `2n(n-1)(n-2)`.  Thus scalar first-moment
closure can fail while a quadratic fluctuation--dissipation law survives.
That law concerns the collar event walk, however, not the cover inequality
itself; the two walls should not be conflated.

The pulled THM-856 supplies the cover-side second-order object directly.  For
remaining combs `S`, put

```text
Phi_H(E,S)=mu(E)-sum_(u in S)mu(E intersect D_u)
 +max_(trees T on S)sum_(uv in T)mu(E intersect D_u intersect D_v).
```

Hunter--Kounias gives `mu(E minus union D_u)>=Phi_H`.  At independent
densities its coefficient is `(165-22|S|)/169`, so it is positive for seven
combs and negative for eight.  This is the correct quadratic use of overlap
debt.  The same audit also proved that every single-comb schema is genuinely
empty at seven.

The initial global pair formula pulled with THM-856 was false: it omitted the
containment plateau between unequal teeth.  For `x=ga,y=gb`, `(a,b)=1`, the
correct exact law is

```text
mu(D_x intersect D_y)
 =4/169+[Q(a+b)-Q(|a-b|)]/(169ab),
Q(c)=r(13-r),  r=c mod 13.
```

The correction has both signs (`D_6 intersect D_7` has measure `2/91`, below
`4/169`) and refutes a raw `C(E)/min(x,y)` deficit bound along the legal ray
`(x,y)=(6g,7g)`, `g=1 mod 13`.  What survives is sharper: for `E` with `c_E`
components,

```text
|mu(E intersect D_(ga) intersect D_(gb))
 -mu(E)mu(D_a intersect D_b)| <= 2c_E/g.
```

Thus the Hunter graph is weighted and undirected, with edge coordinates
`(reduced ratio, common scale, mod-13 sawtooth, endpoint discrepancy)`.  It is
not a switched tournament.  The remaining seven-comb lemma is a maximum-tree
statement on those projective edge packets, followed by a mechanical-word
argument on any consecutive/AP-window residual.

The resulting unification is operation-indexed.  The H6 germ graph uses a
min-plus directed-cycle obstruction; the seven-comb wall uses a max-plus
graphic-matroid basis; the H-drift stalk uses a tridiagonal flip generator;
the common-sheet H5 bank uses complement-lcm subset cuts on a prime-power
carrier hypergraph; and the flood atlas uses the full `K_7` edge module.  In
each case a scalar score orders candidates, but the proof lives in which value
is incident to which cusp, endpoint, fibre, flip, or edge.  The closest common
object is therefore a typed operation profile over literal cover obligations,
not a tournament on runners.

## 7. The sharpened next obligations

The shortest honest route to the requested theorem is now:

1. **Ramified H6 transport:** THM-859 proves that common dilation conjugates
   THM-857's full action tree and that order `D` is exactly the number of deck
   masks created by insertion. THM-860 bounds the presentation base and
   THM-861 closes the full `c=2` fibre except for `[12]`. Evaluate THM-862's
   scale-three bank: `336` two-order-one/four-order-three, `672`
   one-order-one/five-order-three, and `496` all-order-three unit contexts.
   Then lift surviving states
   across `D>1` fibres without losing literal components, phase masks, labelled
   future progressions, last-speed order, equality ancestry, or the shortcut
   witness.  There is no remaining primitive scale-one label census.
2. **Finite ramified H5 sheets:** classify the THM-858 strip
   `min D_i<=21<max D_i<=10,584`.  Every order is `{2,3,7}`-smooth, every maximal
   prime power is shared, and all complement-lcm relative-capacity cuts hold.
   The complete bank through 21 is already reduced to the closed legacy
   languages; the next recursion should act on the decorated prime-power
   carrier hypergraph rather than enumerate the full rectangular order box.
3. **Deep `s=5`:** the uniform single-numerator endpoint-grid template is
   impossible in all three remaining classes, and every fixed pair of unit
   endpoint columns is jointly killable in the relaxed lift model.  Seek a
   genuinely U-dependent column family, three-or-more-column incidence using
   divisor/parity restrictions, a non-endpoint denominator, or the missing
   divisor/satellite obligation.  Do not extrapolate the finite `d=15` census.
4. **General two-sheet transport:** retain signed component--return sum arcs
   and endpoint ancestry outside the sixteen no-switch radius types.
5. **`j=4` flood tail:** port the exact `(5,7)`/`(6,7)` bulk screen across the
   nineteen remaining flood bodies.  Fano symmetry, `chi_7`, and local three-needle
   paths are ruled out as pruning quotients, so each transport needs a numeric
   monotonicity proof or its own exact sweep.

The object underneath these obligations is a ramified sheaf of periodic combs
over literal safe components, carrying owner, scale, residue, future-action
data, and on sheeted branches the complement-lcm fibre.  The newer theorems
meet at that description: THM-857 resolves the literal metric stalk at scale
one, THM-858/860 bound arithmetic ramification, and THM-861 resolves the first
nontrivial fibre without collapsing its sparse sheet cycle.  The older
viewpoints—Farey gaps, sheets, tournaments, Fano masks, continued fractions,
matroid circuits, tropical cycles/trees, and Kakeya needles—are useful charts
precisely to the extent that they preserve this cover-debt packet.

## 8. A recursive definition of the underlying object

The accumulated work suggests a sharper definition than “a difficult
twelve-set.”  The live object is a **ramified component-cover automaton** with
five coupled layers:

1. The projective base records common scale, effective sheet orders, unit
   classes, and prime-power carriers.  THM-858 proves that the H5 common-sheet
   part of this base is finite, but does not decide its metric fibres. THM-859
   gives effective order an intrinsic dynamic meaning: it is the cardinality
   of the deck-mask orbit of one insertion, hence the ramification degree of
   the failed scale quotient.
2. A metric fibre records the exact open components of the current strict-safe
   set and every remaining labelled arithmetic progression of comb speeds.
3. The transition labelled by `u` is the idempotent intersection
   `T_u(E)=E intersect Safe(u)`, together with numerical-order and endpoint
   ancestry.  Covering is exactly arrival at the empty fibre.
4. A proof certificate is a local obstruction to future empty fibres: a
   longest-component cap, full safe tooth, streaming cap witness, owner pin,
   complement-lcm cut, tropical cycle, or Hunter tree credit.
5. The legal operation category records whether the next move is insertion,
   peel/deletion, replacement, scale change, sheet lift, or descent.  Two states
   are continuation-equivalent only when every legal future word has the same
   terminal emptiness verdict.  This is the Myhill--Nerode right congruence of
   the partial cover action.  A proposed quotient is transportable only
   after it is proved to be a congruence for the required operations; THM-840
   shows why insertion-Markov endpoints alone fail this test for deletion.

This definition explains the historical viewpoints without identifying any
one of them with the whole problem.  Hamming shells choose a local coordinate
chart around the AP.  Sheet and dyadic descent move in the projective base.
Farey/continued-fraction addresses track endpoint genealogy.  Fourier modes
and relation-lattice or matroid circuits measure global resonance but forget
which component is hit.  AP-germ cycles are a min-plus necessary condition for
local cusp coverage; the Hunter/Kruskal tree is a max-plus rank-two lower bound
on the uncovered mass.  Exact radius-seven coverage still lives in the
higher-order cover nerve, or dually in a component-weighted fractional-cover
LP/Farkas certificate retaining edge incidence.  Tournaments impose a binary order on one projection and are
usually non-conservative: thousands of edges can switch while the terminal
cover verdict is unchanged.  Fano and `chi_7` coordinates organize seven-fold
incidence but lose the metric edge curl.

In this language THM-857 computes the entire H6 fibre over the trivial
scale-one base and finds a unique empty terminal, the equality section
`2[12]`. THM-859 proves that common dilation is an exact action conjugacy and
that `D=c/gcd(c,w)` is the ramification degree of its failure: an order-`D`
insertion acts through `D` distinct deck masks. THM-858 bounds the nontrivial
H5 ramification base, THM-860 bounds the primitive H6 base, and THM-861
computes the complete `c=2` metric fibre with unique AP terminal `[12]`. The
central missing theorem is now narrower: construct the phase-decorated
transport functor on the remaining `D>1` fibres, beginning with THM-862's
scale-three contexts, so that it commutes with the required legal
operations and preserves terminal emptiness and certificate witnesses.  Merely
matching order multisets, tournament fingerprints, component counts, or
pairwise overlaps cannot supply that transport.
