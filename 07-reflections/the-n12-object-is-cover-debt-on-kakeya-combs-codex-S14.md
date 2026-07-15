# The n=12 object is cover debt on periodic Kakeya combs

*codex-2026-07-15-S14.  Companion to HYP-6820, THM-815/823/836/841/844/845/847,
and the exact Fano/chi7 flood and Hamming-six contraction audits.*

## 0. Honest verdict

The requested uniform theorem—emptiness of the sporadic branch at `n=12`—is
still open.  This session did not turn the accumulated finite banks into that
global statement.  It did, however, remove several named residuals and change
the shape of the remaining object:

1. THM-847 closes the last effective-order-at-most-twelve Hamming-five
   common-sheet language.  All-one, all-three, and mixed one-plus-three are now
   uniformly empty at arbitrary lift height.
2. The nonprimitive scale-one Hamming-six branch contracts exactly.  Its only
   possible tight row is the doubled AP `2[12]`; every other nonprimitive row
   is loose.  Primitive H6 remains open.
3. THM-836's local `s=5` shell is exactly four congruence classes, and the full
   class `d=11 mod 52` is uniformly impossible by one explicit divisor-grid
   witness.  The classes `15,37,41 mod 52` remain open uniformly.
4. THM-841's proposed toothpick breakpoint tree does not exist.  The ladder has
   a much simpler all-order totient formula.
5. The Fano/`chi_7` organization of the 21 `j=4` flood bodies is exact but not
   a symmetry quotient and not a local three-needle closure.

This is progress toward the sporadic theorem, not the theorem itself.

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
tight nonprimitive H6 packet.  The remaining H6 problem is now intrinsically
primitive.

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
Likewise, the proposed all-`n` black self-line law loses its first untested
case at `n=8` (`404` versus `SC(8)=176`), while only the Klein-four orbit
structure survives.  Both corrections point away from scalar class counts and
toward decorated incidence/orbit packets.

The concurrent THM-810 F6 second-moment closure sharpens the positive side of
that warning.  For the score-axis energy `x`, one arc flip satisfies the exact
pointwise law

```text
E[(Delta x)^2|T]=16(n-4)x/C(n,2)+64,
```

and the stationary variance is `2n(n-1)(n-2)`.  Thus scalar first-moment
closure can fail while a quadratic fluctuation--dissipation law survives.
The seven-comb wall has exactly this shape: its linear density potential stalls
at `14/13`, but cover debt has a nonnegative quadratic part.  The concrete next
experiment is not to transplant the tournament axis blindly; it is to extract
the size-biased component/comb event kernel from the THM-815 recursion and
compute the drift of `omega_I` and `omega_I^2`.  This preserves the union-cover
predicate that the score-axis quotient would destroy.

## 7. The sharpened next obligations

The shortest honest route to the requested theorem is now:

1. **Primitive H6:** extend the component--comb recursion from the closed
   nonprimitive chamber to the 923 primitive missing-label rows, using overlap
   debt and unique-owner pins before attempting a raw full tree.
2. **Unbounded H5 sheets:** classify common-sheet presentations with an
   effective order above twelve.  The entire bounded bank is already empty.
3. **Deep `s=5`:** find uniform grid witnesses for the remaining classes
   `15,37,41 mod 52`, or exhibit which additional divisor/satellite obligation
   is missing.  Do not extrapolate the finite `d=15` census.
4. **General two-sheet transport:** retain signed component--return sum arcs
   and endpoint ancestry outside the sixteen no-switch radius types.
5. **`j=4` flood tail:** continue the exact resumable sweeps.  Fano symmetry,
   `chi_7`, and local three-needle paths are now ruled out as pruning quotients.

The object underneath these obligations is a ramified sheaf of periodic combs
over literal safe components, carrying owner, scale, residue, and future-action
data.  The older viewpoints—Farey gaps, sheets, tournaments, Fano masks,
continued fractions, and Kakeya needles—are useful charts precisely to the
extent that they preserve that cover-debt packet.
