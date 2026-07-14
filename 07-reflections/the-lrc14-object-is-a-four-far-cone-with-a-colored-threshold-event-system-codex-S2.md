# The First Open LRC14 Stratum Is A Four-Far Cone With A Colored Threshold Event System

*codex-2026-07-14-S2. A creative information-preservation session following
HYP-6780, the completed Lean THM-755 chain, HYP-6785, and the exact endpoint
sidecar audit. Updated after the exact HYP-6830 fragmentation falsifier.*

## 1. First challenge the word "four-dimensional"

There is a literal four-dimensional object at the current frontier, but it is
not the whole LRC14 moduli space.

THM-738 closes the rows with at most three speeds above 14. The first
unclosed far-count stratum therefore has a fixed nine-speed core
`C subset {1,...,14}` and four ordered far speeds. In gap coordinates those
four speeds are controlled by

```text
(a,d1,d2,d3) in Z^4,
w=(a,a+d1,a+d1+d2,a+d1+d2+d3).
```

This is the honest four-dimensional chart. There are `binom(14,9)=2002`
possible small cores. On the exactly-nine-small stratum the core is unique.
THM-741 itself quantifies over rows with at least nine small speeds, so its
2002 bodies overlap on `f<=3`; the cone atlas describes its first new rung,
not its whole quantified body as a disjoint union.

The global primitive 13-speed space modulo dilation has twelve ratio degrees
of freedom, and the `f>=5` strata have larger far cones. Any claim that the
entire problem is four-dimensional already assumes a descent theorem that the
repository does not possess. The first reframe is therefore negative and
useful:

> Do not reduce the dimension in prose. Exhibit the functor that performs the
> reduction and list the fiber information it must retain.

## 2. The outer object

For a fixed core `C`, the covering far-speed set `X_C` is a semilinear
rank-four cone chart. Covering depends only on divisibility by `2,...,14`, and
hence on residues modulo `360360`. Any nine-element subset of `{1,...,14}` has
gcd one, so primitivity is automatic on this chart. Thus `X_C` is a finite
union of residue classes intersected with the positive ordered cone. This gives
finitely many residue addresses; the chart itself remains infinite. "Cone"
does not mean an ordinary real convex cone: `X_C` need not be addition-closed.
It does carry an exact positive-integer dilation action `g -> k*g`, because a
divisor of a far speed remains a divisor after multiplication by `k`. This
action fixes `C` and scales only the four far coordinates, so it is not global
speed dilation and can change exact clearance and other metric data. The exact
covering canary with `C={1,...,9}` has `M=6/61` at far tuple `(22,26,28,60)` and `M=12/121`
after that far tuple is doubled. The action must be retained rather than
silently quotiented.

### The chart has two inequivalent directions to infinity

The far-dilation action is only the radial direction. THM-776 supplies a
transverse frequency refinement. For prime `N>110`, the ray

```text
C={1,...,9},       g_N=(15,110,N,1092N)
```

lies in `X_C`, is primitive and covering, and has no divisor packet of size
seven. Yet after peeling `1092N`, the core good set has unbounded component
count. The fixed base `{1,...,9,15,110}` has a safe interval
`[1/14,111/1540]`; the `N`-runner cuts it by at least
`N/1540-8/7` disjoint teeth.
Exact component counts are `66,104,174,310` for
`N=211,503,1009,2003`.

So the first four-dimensional object is not merely an ordered cone with a
radial boundary. Its fiber cellulation refines without bound along a
scale-free path. Coherent dilation increases the number of copies of an old
cell pattern; tooth insertion creates new walls inside an old cell. These are
different degenerations and require different compactification coordinates.
The failed complementarity `r_P<=B(c*)` collapsed them into one scalar.

The replacement suggested by the cap theorem is relational: for a named peel
`v`, retain `r_P/(v|G'_P|)=pi*v*(P)/v`, the divisor-support profile, and the
endpoint owners that created the cuts. This is not yet a finiteness theorem.
It is the smallest current proposed payload that distinguishes radial replication
from transverse refinement while remembering why THM-776 closes the whole
prime ray by the capped peel.

THM-777 identifies the exact rate law behind that closure. If the old fiber has
safe mass `mu` and `r` components, a new frequency `N` leaves at least
`6mu/7-2r/(7N)` safe mass and creates at most `N+r` components. A named interval
is only one way to certify the initial state. A peel of speed `aN` therefore sees
a bounded normalized wall load and becomes terminal once `a` crosses an explicit
state-dependent threshold. In this chart, increasing geometric complexity can
decrease proof difficulty.

The transition also iterates. Certified lower/upper states evolve by

```text
T_N(mu_lower,r_upper)
  = (6mu_lower/7-2r_upper/(7N), r_upper+N).
```

The underlying insertions commute, but the enclosure operators `T_N` do not:
their error remembers which frequency was inserted first. This apparent defect
has an exact resolution. If `x<y`, inserting `x` before `y` improves the mass
lower bound by

```text
(2/7)*(y-x)/(xy)*(x+y+r_upper/7)>0.
```

Thus adjacent exchanges make increasing frequency the unique optimal order.
In a four-far one-peel search, the 24 formal peel/order gauges collapse to four
canonical certificates, one per peel. Insertion order is therefore not new
moduli data; it is an observer choice with a canonical representative in the
certificate calculus. The nonterminal boundary must make all four certified
mass bounds or cap tests fail, slow the peel relative to wall creation, or
retain owner correlations that the two-number state intentionally erased.

That observation does not make LRC periodic. Inside one residue class,
changing height changes:

- the rational endpoint phases;
- the cyclic event order and simultaneous tie blocks;
- the safe component lengths;
- the pair-sum moduli used by peak witnesses;
- the THM-755 ratio and Bernoulli discrepancy;
- the embedded integer relations and their coefficient heights.

The point of the chart is not finite enumeration by itself. It supplies the
base over which a fiberwise family of predicate-bearing event words varies;
a uniform wall stratification and transport law have not yet been built.

For `g in X_C`, the speed vector gives a cocharacter

```text
phi_{C,g} : T -> T^13,
t |-> (s*t mod 1)_{s in S(C,g)}.
```

Let `K=[1/14,13/14]^13`. The universal safe incidence

```text
Z_C={(g,t): phi_{C,g}(t) in K}
```

is the cleanest exact formulation of the residual primitive-covering
exactly-`f=4` problem. LRC on this chart is surjectivity of `Z_C -> X_C`, after
the dilation and noncovering reductions (including THM-366). This formulation
does not subsume the separately dispatched noncovering rows.

This uses exact torus language: speeds define a cocharacter of the compact
torus and the runner motion is its image. The target `K` is a closed arc box,
not an algebraic toric subvariety, so “toric variety” would be stronger than
the construction and there is no toric proof here.

## 3. The inner object

Many hard infinite families are affine rays

```text
V(c)=cA+R.
```

Put `u={ct}`. Then exactly

```text
(c*a_i+r_i)t = a_i*u+r_i*t mod 1.
```

One two-torus field

```text
Phi(u,t)=min_i ||a_i*u+r_i*t||
```

contains every member of the family as the slope-`c` slice `u=ct`. Adding
the integer slope `c` and clearance level `lambda` produces the mixed
four-coordinate suspension `(u,t,c,lambda)`. Here `c` is discrete and
`u=ct`; at fixed `c` only `(t,lambda)` vary continuously. This is an exact
incidence description, not a second four-dimensional manifold. It is also
broader than a fixed `X_C`: a general affine family can cross far-count
strata, as the executable AP dilation and shear examples do.

This recovers several old views as literal coordinate projections:

- `R=0`: a cylinder, hence dilation invariance;
- bounded `R`: transverse shear, not negligible error;
- threshold `lambda`: the persistence direction;
- fixed `c`: the ordinary exact circle arrangement;
- varying `c`: arithmetic tomography by closed geodesic slopes.

THM-742 already proved this strand-versus-area idea for an additive cluster
chart. HYP-6815 is the arbitrary owned-offset version, with the threshold
coordinate retained.

### The first proved chart change between scales

The concurrent THM-761 result turns one part of the proposed recursion into a
theorem. For `V=cD union W`, choose a core time `t0` and inspect the sheets

```text
t_k=(t0+k)/c,  k in Z/cZ.
```

Every core runner `cd` has exactly the phase of `d` at `t0` on every sheet.
The exceptions become an inhomogeneous covering problem on the finite cycle:

```text
w t_k=(w t0)/c + (w mod c) k/c.
```

This is controlled forgetting with an exact annihilation certificate. The
core event word can be dropped for this operation because the required core
margin is sheet-invariant. The exception packet cannot be reduced to separate
runner summaries: it needs the shared sheet index, owner-labelled residues,
offsets, gcd multiplicities, and closed boundary convention.

THM-761 proves a free sheet when the summed bad-sheet budget is below `c`; in
LRC14 it closes every one-exception packet and the large-scale regime with at
most six exceptions (`c>=43` uniformly in the coprime case). Seven exceptions
can tile the sheet cycle, and the small-scale/scale-free splice remains open
in HYP-6830. This sheet chart is distinct from both global dilation and the
fixed-core far-coordinate action above.

That distinction is structural. A fixed-core far ray `C union kF` has scaled
core `F` of size four and exception set `C` of size nine. The sheet identity is
still exact, but THM-761's free-sheet bound does not fire at `r=9`. The literal
four-far chart therefore sits inside the open seven-or-more-exception residue,
not inside the newly closed high-support regime.

### Tight-core recursion reveals the missing fiber

THM-769 makes the same preservation issue exact at the 12-core boundary. At a
tight maximum `p/(13s)`, the endpoint owners lie in the on-sheet packet
`E=sU`, while every off-sheet speed is strictly interior near that maximum.
This local split is not enough to preserve tightness. The exceptions must cover
all `s` lifts over every point of the quotient core's loose set `G_U`.

The deciding object is therefore the persistent incidence field

```text
(tau,j,w) with tau in G_U, j in Z/sZ,
and w closed-dangerous on the j-th lift of tau.
```

Its first consequences are already rigid: a deep primitive tight packet needs
at least two off-sheet speeds; two force `s=2` and a persistent opposite-parity
two-colouring; the three-speed equality edge is either touched by a half-sheet
tightener or has `s=3` with persistent ownership of all three colours. A residue
multiset at one maximum, a gcd budget, and a runner tournament each forget the
quantifier over all `tau in G_U`.

THM-770 supplies the complementary shallow evidence. Across all
`13^12` full-residue packets of lift height at most twelve, endpoint-owned cell
incidence leaves only thirteen zero-defect dilates and only `{1,...,12}` after
primitive normalization. Its diagnostic tournament has all `66` pairwise
comparisons tied; the endpoint incidence hypergraph still separates those
thirteen solutions from more than `2.3e13` candidates. This is an exact warning
against demanding that a useful proof carrier be binary or pairwise. The
recursive state is a coloured hypergraph varying over a metric base, with the
sheet action and closed boundary convention retained.

### Truth compression is not transport compression

The concurrent THM-772/773/774/775 programs make the preservation hierarchy
more exact. THM-772 and THM-774 are proved; THM-773 has a stored exact audit but
remains claimed at theorem level, and THM-775's dyadic descent remains claimed.

At prime sheet size seven, THM-773's stored audit shows that the six nontrivial
finite-field moments recognize exactly the 5,040 permutation states among all
`7^7` owner-token configurations. Those moments therefore preserve the Boolean
predicate "the sheets are exactly covered." But they give the same full-grid
values on every permutation. Linear and circular tournament gauges likewise
collapse all assignments to one unlabelled node, while exact chamber movies
exhibit equal coarse states with different next owners and free-sheet futures.
Truth needs the moments; transport needs the labelled assignment, event step,
and endpoint schedule.

At binding scales two and three, THM-772 points in the opposite direction.
Eligibility at one modulus is not the predicate. The proof obligation
is simultaneous compatibility of all unit-fraction columns, including the
modulus-6/12 splice. That carrier is a hypergraph of proof obligations; its
pairwise tournament loses the joint column. THM-774 proves a second compression
for the two-sheet packet: after folding `(x,y)` to `(a,b)=((x+y)/2,(x-y)/2)`,
eligibility and opposite colour should become one `l1` distance inequality.
Its exact endpoint and measure audit shows that the folded metric preserves the
local obstruction but still not the cross-modulus transport data.

THM-775 then proposes a genuinely recursive use of the same carrier. An
imprimitive deletion should force a unique dyadic seam, lift the owner field from
`Z/2` to an exact `2+1+1` tiling on `Z/4`, and transfer divisor completeness to a
smaller quotient. If completed, the descent datum is not merely `gcd=2`; it is
the labelled seam plus the quotient's persistent unit-grid obligations.

The resulting rule is task-indexed. Preserve only moments for exact-cover truth,
the full owner stalk for future events, the obligation hypergraph for persistent
unit-grid compatibility, and the folded metric for the proposed local two-sheet
obstruction. "Keep everything" and "find one invariant" are both too coarse.

### A smaller local address than the edge table

At one time write `x_i={s_i t}`, `kappa_i=floor(14x_i)`, and
`rho_i={14x_i}`.  The 14-sector label of a directed difference factors as

```text
gamma_ij=[kappa_j-kappa_i-1_{rho_j<rho_i}] mod 14.
```

So 156 directed edge labels are not independent.  They come from thirteen
vertex potentials and the inversion matrix of one global weak order.  This is
an exact address for the sector movie and a better state than a raw binary
tournament.  It is not complete without wall flags: `x_i=13/14` is safe, but
the one-sided sector label merges it with dangerous points immediately above
the endpoint.  It also does not retain safe lengths or scale transport.

## 4. The dual object

For every coefficient vector `z in Z^13`, project the relation to

```text
(m,n)=(z dot A,z dot R) in Z^2.
```

On the two-torus its character is `exp(2*pi*i*(m*u+n*t))`. On slope `c` it
becomes `exp(2*pi*i*(c*m+n)*t)`. Therefore the relation survives exactly on
the dual line

```text
c*m+n=0.
```

This is the exact bridge among four repository languages:

```text
normalized shape        -> m
owned residue/detuning  -> n
scale slope              -> line c*m+n=0
eligible marked modes    -> full z-preimage over that line
Fourier/relation mass    -> multiplicities and weights on that preimage.
```

For finite trigonometric polynomials, integration on a slope fiber selects
exactly the marked coefficient vectors `z` above the line. Bare projected
points `(m,n)` lose multiplicity and coefficient weights. With a justified
Abel/Fejer limit the full weighted preimage can recover the
one-dimensional Haar measure of the pulled-back safe indicator. It cannot
recover closed-threshold nonemptiness or isolated equality witnesses; the
AP/Goddyn-Wong boundary cases can be lonely while safe measure is zero.

The shape alone invents relations. In the executable audit, the AP shape has
36 support-three relations `e_i+e_j-e_k`; one owned `+1` shear destroys 11 if
attached to owner 1, 9 if attached to owner 7, and 6 if attached to owner 13.
The missing information is not "a residue" in the abstract. It is the
pairing between the residue coordinate, its owner, and the relation
coefficient.

The relation lattice must also remain embedded and marked. A code weight
enumerator, lattice rank, additive energy, or successive minimum forgets how
the relation meets the thirteen coordinate facets of the safe cube. That
incidence is where the observer lives. More precisely, after owner labels and
the target anchor are fixed, the full kernel embedded in `Z^13` determines a
primitive positive speed vector: its normal is unique up to sign and
positivity selects the sign. For a nonprimitive vector it determines only the
primitive normalization, enough for `M` but not for component count.

## 5. There is no demonstrated nontrivial universal quotient

The phrase "what information must be preserved?" is incomplete until the
next question is named. The identity representation is of course universally
sufficient; the open issue is a useful compressed carrier that survives more
than one proof operation.

### To decide LRC truth at one fixed threshold

Start at `t=0`, where every runner is dangerous. Record the cyclic endpoint
events. Each simultaneous block retains entering owners and exiting owners.
If `B_-` is the dangerous-owner set just before a block, with exiting owners
`E` and entering owners `A_in`, then the exact tie and next-cell states are

```text
B_tie = B_- \ E,
B_+   = (B_- \ E) union A_in.
```

Integrating these typed updates reconstructs the open-cell cubical path and
all zero-dimensional equality states:

```text
beta(t)_i=1_{||s_i*t||<1/14} in {0,1}^13.
```

The row is lonely exactly when an open cell or exact tie state reaches the
zero vertex. For this Boolean question, metric gap lengths can be forgotten.

Owner colors and tie blocks cannot. If one runner exits when another enters,
their scalar currents cancel. The aggregate current would erase a boundary
state that may be the only equality witness.

### To preserve exact clearance

The fixed-threshold event word is insufficient. Retain the whole threshold
filtration or an exact peak witness, value, and maximality certificate. A time
witness alone supplies only a lower bound; a certified `(t*,M)` with an exact
upper bound determines `M`, but it does not determine safe measure or
component topology.

### To preserve measure and autocorrelation

Retain rational endpoint phases and gap lengths. The signed owner endpoint
divisor reconstructs the Bernoulli `B2` formula. The unweighted endpoint
tournament sees interlacing but not metric discrepancy.

### To preserve covering

Retain the divisor mask. Neither the runner tournament nor endpoint
tournament determines it.

### To preserve the capped-envelope decision

For peeled speed `v_peel` and peeled core `B`, retain the projective comparison
`v_peel*|G'_B|/r_B > 1/pi`, equivalently `pi*v_peel*|G'_B| > r_B`, or its exact
comparison bit. A single lift in the endpoint audit changed the cap status
with zero endpoint edge flips.

### To preserve deletion, peel, or scale transport

Retain owner labels and the action on the scale/residue fiber. A quotient can
be truth-safe for the current row and illegal for the next observer.
THM-765 makes the needed metric sidecar exact in a one-runner peel: safe
component length plus midpoint phase determines tooth containment, while the
gcd translation deck detects whether all lifted components can be covered.
Its hereditary-primitivity corollary gives no new pruning on exactly `f=4`:
deleting a far speed leaves the gcd-one nine-core, and deleting a small speed
leaves eight distinct integers in `[1,14]`, already of gcd one.

### To preserve a proof

Retain the next operation, certificate availability, discharge mode, and
named residual debt. This is proof state, not extra mathematical truth, but
without it local certificates cannot be composed honestly.

The hierarchy is:

```text
truth   = initial state + signed owner event blocks
metric  = truth + rational phases or threshold filtration
functor = metric + arithmetic action + marked relation embedding
proof   = functor + next observer + certificate/discharge
```

## 6. The exact canaries

The companion 552-row audit gives a compact impossibility theorem for several
tempting compressions.

Raw runner tournaments have mixed fibers for covering, `M`, cap status, and
discrepancy. Raw endpoint tournaments also mix all four. Adding the divisor
mask fixes covering. Adding the cap sign fixes the THM-755 bit. Signed
endpoint phases fix `B2`. The value of `M` remains mixed within those carrier
fibers until a peak witness, value, and maximality certificate are attached.

The cleanest pair is

```text
{1,...,13}       M=1/14
{1,...,12,26}    M=2/27.
```

They share the endpoint tournament, divisor mask, and cap sign. A quotient
which keeps order, covering, and the terminal cap route still loses the
clearance.

The affine audit supplies a different canary. For the same shape, same owned
offset, and the same `c mod 14`, scales `c=2` and `c=16` both have `M=2/15`,
but their safe measures are `5/84` and `115/1904`, with 4 and 30 components.
Even exact `M` plus residue data is not a metric-topology carrier.

The outer chart adds a complementary failure: multiplying only the four far
coordinates preserves membership in `X_C` but changes `M` on the core-nine
canary above. The relevant object is the monoid action on event fibers, not
the orbit set of gap tuples.

THM-764 turns the concurrent HYP-6820 audit into an exact task-specific
carrier. For `15<=q<=28`, a rational witness is decided by the zero-owner set
and owner-resolved signed-unit-pair deck. The `c=26` shear ray blocks every
denominator from 15 through 25 and first witnesses at 27; a gcd-incoherent
uncapped residual blocks the same window and first witnesses at 26. This deck
preserves the bounded rational-witness predicate, not full LRC or metric data.

THM-761 supplies the complementary chart: the coherent `c=26` row closes on
its native sheet cycle without a bounded-denominator theorem. The invariant is
therefore not one universal deck, but an atlas of exact carriers with declared
transitions: endpoint words, signed-pair decks, scale sheets, and blocker
complexes.

## 7. What unrelated threads contributed

### Tournament switching and tilings

The same cube with a different group folding has different invariants. The
transferable datum is the action and a gauge section, not the underlying
bijection. For LRC the analogue is scale/residue action plus the marked
observer.

### Coding theory and matroids

Identical weight enumerators can hide decomposable and indecomposable support
graphs. Relation support counts are not relation realization. LRC needs
support incidence, coefficient height, and the embedding against safe-cube
facets.

### Ising and Lee-Yang

Total odd-cycle mass and compatible cycle packings separate at higher
interaction order. A scalar moment or `p0` loses the organization of
obstructions. The analogue here is the full threshold or miss-count
filtration. Root surfaces are candidate phase classifiers in the finite banks
audited by that thread; no global classification follows, and they do not
replace an endpoint certificate.

### Resolvent folds

A symmetric inner block becomes exact only after center coupling,
antisymmetric leakage, and boundary sectors are retained. This suggests a
dictionary in which projective shape plays the inner-page role and owned
shear plus endpoint boundary data play leakage sidecars. No exact
resolvent-to-LRC map has been constructed.

### CRT and p-adic trees

Residue skeleton and valuation height are independent channels. The affine
audit gives an archimedean version of the same lesson: equal residue and owner
can retain `M` while safe length and component count move.

### Canonical metagraph addresses

HYP-6825's rooted line address can serve as a finite combinatorial base, but
it cannot by itself be the LRC quotient. The `V_N` family demonstrates the
missing stalk: bounded divisor-scale type coexists with an owner-labelled wall
word whose length grows without bound. This does not prove that the rows share
one metagraph node. It proves what any pullback from that unrelated tiling
atlas must preserve: metric gaps, threshold walls, owners, and peel-relative
load remain fiber data after canonical labeling.

### Cut/cycle topology and anti-local testing

The tested local invariant views can agree while global realizability or
higher packing differs. The tiling-cycle-moment THM-555 sharpens the warning:
aggregate cut scores and low cycle counts may survive while higher overlap and
packing data are lost. The endpoint loop and endogenous blocker complex are independently equivalent
fixed-row presentations of LRC truth, but no pairing, chain map, or gluing
theorem between them has been built.

### Observer categories and baby-Hodge holes

The observer-category thread separates anchored cyclic order, metric gap
widths, and observer placement. The LRC object is marked: an observer-blind
order class can contain both safe and unsafe placements. By analogy, THM-509's
tournament baby-Hodge holes give a different warning: a continuous relaxation
can miss integer realizability. No LRC torus-suspension hole has been exhibited;
the transfer is the discipline of carrying a realizability or carry syndrome
before a continuous certificate discharges an integer row.

### Observer-cut ledgers and proof circuits

These ledgers supply an audit discipline, not a universal quotient theorem.
Relative to a named next observer, every changed predicate should have
reconstruction, exact/coboundary cancellation, dual annihilation, family
descent, a boundary stop, or named residual debt.

## 8. A different style of mathematical thought

The conventional question asks for an invariant of a speed set. The more
productive question asks for a contract between representations:

```text
What operation comes next?
Which fibers does this quotient merge?
Which theorem-facing predicates vary inside those fibers?
What action survives on the forgotten coordinate?
Which sidecar repairs the first failure?
```

This turns failed invariants into useful output. A failure names the next
coordinate rather than merely disqualifying an analogy.

It also suggests replacing one linear proof ledger by a bicomplex:

```text
d_geom  = cross an endpoint/contact wall
d_arith = change residue, valuation, scale, or detuning depth.
```

If exact versions of the two moves commute, a quotient may descend. If they do
not, their commutator could locate an owner, holonomy, coefficient-height, or
proof-route sidecar that was erased. At present this is a schematic program:
the maps and comparison law have not been constructed, and the word
"curvature" alone proves nothing.

## 9. Recursive view of infinity

HYP-6780 showed why raw maximum speed is the wrong induction variable. The
four-cone suggests a replacement.

For packets `cD union W`, this is no longer entirely schematic. THM-761 is an
exact descent to the finite pointed cycle `Z/cZ`: the core is fiber-exact and
the exceptions become the next level's residue runners. The immediate open
boundary is explicit: the seven-exception strata outside THM-767's
commensurate event-pierce lane, excessive-gcd descent, and HYP-6830's
peel-relative replacement for the refuted claim that scale-free cores have
controlled raw good-set fragmentation. Corrected THM-767 shows that an exact
deck tiling cannot cross an event: it is chamber-locked, while the coincidence
law locates double-boundary walls without governing chamber persistence.
The projective compactification below is needed for arbitrary unbounded
sequences outside that proved sheet chart.

At infinity, keep a flag:

```text
leading projective shape A
first nonzero owned detuning R
residue/valuation address
adaptive exact-period/blocker deck
cluster tree or relation circuit
terminal certificate or smaller descendant.
```

Pure dilation is the zero-detuning cylinder. Additive clusters,
multiplicative packs, and hierarchical clocks suggest boundary faces, with
dispatches proved only on their existing theorem domains. A projective or
tropical compactification of `X_C` has not been constructed, so arbitrary
"paths at infinity" are not yet defined.

A precise target should quantify over an unbounded sequence in `X_C`, pass to
a subsequence with stabilized residue address and limiting projective face,
and then prove that the subsequence enters a known coherent face, gains a
uniform lonely margin, emits a bounded-height marked circuit, or maps to a
strictly smaller packet. This remains schematic until the compactification,
cluster boundary faces, and well-founded packet order are defined. Such a
theorem would make finite computations terminal leaves of a proof rather than
samples from an infinite box.

## 10. Concrete next lemmas

1. Extend THM-777 beyond its positive-good-state/proportional-peel face: stratify
   sequences by safe-mass decay relative to component, wall-creation, and peel rates; attach
   endpoint owners, and prove every residual state is bounded or recursively descends.
2. Classify the seven-exception `1/7` chambers left after THM-767's
   event-pierce/chamber-locking lane and formalize the terminating excessive-gcd
   descent `c -> c/g`.
3. Formalize the semilinear four-cone and colored endpoint-loop criterion.
4. Compute the event cocycle of `g -> k*g`: classify the owner blocks that
   split, merge, or reorder and find the smallest sidecar transporting truth,
   `M`, measure, and component topology along the action.
5. Construct a comparison theorem between cubical zero-state reachability and
   an empty protected edge in the endogenous pair-sum blocker complex,
   including the data transported between presentations.
6. Test a finite-jet theorem on affine rays outside THM-761's sheet regime;
   THM-764 refutes a uniform `q<=25` terminal window.
7. Prove a four-circuit localization lemma: outside known coherent faces, a
   blocker-complete point forces a bounded-height marked relation involving all
   four far coordinates.
8. Build an observability matrix over the actual `f=4` cone, not a curated row
   bank, and find the smallest sidecar portfolio separating every
   truth/metric/peel-changing fiber pair.
9. Use the existing 2002-core runner as a base case only after its output and
   coverage protocol are present and resumable.

## 11. What this session did not prove

The affine suspension is an exact reparameterization, not a compactness
theorem. The semilinear residue quotient does not make LRC residue-periodic.
The colored threshold event system has not been transported uniformly over
the cone. The projective compactification, comparison map, four-circuit
localization, and cone descent remain open. Torus, persistence, and curvature
language is retained only where an exact object or falsifiable proposal has
been stated.

The durable conclusion is narrower:

> The first unclosed exactly-`f=4` predicate-bearing object is a semilinear
> rank-four chart carrying a fiberwise family of marked torus loops. It has a
> geometry-side owner-colored threshold event presentation and an
> obstruction-side embedded relation/blocker presentation. Any useful quotient
> must preserve the action required on its forgotten fiber, and the amount of
> retained information must be indexed by the next mathematical question.
