# Source orientation, a complete boundary bound, and six-coordinate contacts

2026-10-18, clock_duality_oct15 / frontier_petersen_oct18.

**PROVED, scoped, self-audited:** the sibling covariance tensor, the
clock-preserving packet identity, the finite reduction for the fixed
frontier's boundary imbalance, and the support-set contact compiler.
**PROVED with a reproducible finite certificate:** the sharp fixed-boundary
range `-72 <= D0(N) <= 1428` for every integer `N >= 1`.
**FINITE-EXACT:** two concentration tests at every `4096 <= N < 2^20` and
every first-arrival prefix; all 570 marked six-colour patterns extend to
the mixed-quadrilateral and six-coordinate graphs tested here.
**OPEN:** the concentration inequalities at arbitrary scale/depth, universal
Collatz entry, and six-colour forcing. No new resolution of the planar
Jacobian problem is claimed in this session.

## 1. Inheritance and concept board

- Closest mechanisms: [THM-4690, first-arrival port and integer-clock
  decompositions](../../01-canon/theorems/THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md),
  the [Eisenstein source lift, E1--E4](root_eisenstein_first_arrival_20261018j.md),
  and audited [THM-4700, unit-path ring and mixed
  quadrilaterals](../../01-canon/theorems/THM-4700-the-unit-path-ring-of-the-six-colour-field-short-reciprocal-walks-and-the-heron-census-of-coefficient-triangles.md).
- Canonical hostiles: the expanding Gram direction at `N=524288` in E2;
  the failure of trace contraction at `N=1204`; and port-full pieces whose
  port-only gluings cannot restrict the extension relation (Proposition G
  in [the previous ternary note](ternary_ports_golden_defect_20261018j.md)).
- Corrected near miss: the incoming ideal-ladder prediction of persistent
  twin dissipation failed at depths 8--10. Its sign survived those tests,
  but its magnitude decayed. See [the incoming unit-path synthesis,
  section 6](unit_path_ring_kuratowski_pythagorean_yang_mills_20261018.md).
  We do not use the old persistence heuristic.
- Least-used sidecars: the oriented square of an Eisenstein charge, the
  finite stopped-source boundary term, and the exact set of attainable
  colour differences on each coefficient factor.

The live board is: **actual source orientation** (anchor), **stopped-source
scale events**, **integer-clock sibling packets**, **mixed interior
contacts** (niche), and **marked extension relations versus minors**
(wildcard). The new source calculation changes the first three together:
rotation-invariant energy misses a leading cancellation, while the lost
orientation has a bounded, completely calculable initial boundary term.
The geometric calculation changes the last two: even a closed six-way
contact is insufficient in the tested contexts; exact boundary predicates,
not connectedness or minor structure, remain the test.

## 2. The full sibling tensor

Use the original source, without changing its residues. For primitive odd
seeds `s != 5 mod 8`, let `B(s)=4s+1` and

```
c_r(s;N) = sum_{j: B^j(s)<4N} K_r(B^j(s)),  r=0,1,2,
K_r(m) = (0,-3,3)[m+r mod 3]       if m<2N,
         (-1,-1,2)[m+r mod 3]      if 2N<=m<4N.
```

Only `r=0` is the actual ROOT source. The other rows use its **same** seeds,
routes, ports, and integer clocks. Write `c=(c_0,c_1,c_2)^T` and define

```
S_N = sum_{x=1 mod 8, x<2N} c(x;N)c(2x+1;N)^T.
```

**S1 (PROVED).** Entrywise, for all integer `N -> infinity`,

```
S_N = (N/7) [ 5 -1 -4 ] + O(log(N+1)).                 (1)
              [-1 -4  5 ]
              [-4  5 -1 ]
```

Thus the leading trace is zero, although the actual row's sibling product
has leading term `5N/7`. The tensor's nonzero eigenvalues are
`+3N/sqrt(7)` and `-3N/sqrt(7)`. This is not a positive-semidefinite
contraction in every source direction.

**Proof.** A chain ending in the upper band has row vector
`U(r)=((-1,-1,2)[r+j])_(j=0,1,2)`, independently of its length. A chain
with `L` lower-band points and none in the upper band has vector
`D_L(r)=(sum_{i<L}(0,-3,3)[r+j+i])_j`; it depends only on `L mod 3`.
The two sibling shell types of THM-4690 have, for each residue `r mod 3`,
counts `N/(12*4^L)+O(1)` and `N/(6*4^L)+O(1)` respectively. Their
products are `U(r)D_L(2r+1)^T` and `D_L(r)U(2r+1)^T`.

Summing over `r`, both product matrices equal

```
L=1 mod 3: [ 9  0 -9 ]       L=2 mod 3: [ 9 -9  0 ]
           [ 0 -9  9 ]                  [-9  0  9 ]
           [-9  9  0 ]                  [ 0  9 -9 ].
```

For `L=0 mod 3` they vanish. The geometric weights are `16/63,4/63,1/63`
for positive `L` in these three classes. Multiplying by `1/12+1/6=1/4`
gives (1). There are `O(log N)` shells. The endpoint exceptions already
identified in THM-4690 change only `O(1)` bounded summands per shell.
The unused geometric tail contributes `O(1)`. This proves the error claim.
The eigenvalues follow from the displayed matrix. QED.

The source Gram itself is
`(17N/7)(3I-J_3)+O(log N)`. This follows either by the same residue shells
or by E3 together with the residue-cyclic symmetry of their leading terms.

**Orientation form.** Let `omega^2+omega+1=0`, with positive imaginary
part, and `A_s=a_s+b_s omega` be E1's exact lift, so
`c_r=2 Re(omega^r A_s)`. Equation (1) is equivalently

```
sum_x A_x conjugate(A_(2x+1)) = O(log N),
sum_x A_x A_(2x+1) = ((2-omega)/7)N + O(log N).          (2)
```

The Eisenstein norm of `2-omega` is exactly 7. Its appearance is explained
by the shell sum, not by fitting a nearby constant or by an analogy with
another occurrence of 7. The Hermitian correlation loses the main term;
the un-conjugated correlation retains it.

## 3. Folding siblings without losing their clocks

Let `V` be the fixed, forward-closed 627-seed frontier in
`first_arrival_20261018/V_closed_Y1000.txt`. Remove from the sibling sum
pairs with either endpoint in `V`, obtaining `S_N^V`. Each seed belongs
to at most one sibling pair, so each matrix entry changes by at most
`9|V|`.

**S2 (PROVED).** For each remaining pair `(x,y=2x+1)`, the two sources
either both first enter `V` after the same native depth or neither enters.
When they enter, their port is the same and their integer clocks are
`t,t+1`. Their exact polynomial contribution is

```
z^t(A_x + z A_y).                                      (3)
```

Consequently every finite first-arrival prefix can use one two-term packet
instead of two source records. At `z=-1` its charge is
`(-1)^t(A_x-A_y)`. This substitution preserves the source, stopping rule,
port, native cutoff, and both clock exponents. It does **not** identify the
polynomial (3) with a constant packet at arbitrary `z`.

Proof: the two first native edges have the same target and clocks differing
by one (THM-4690). Neither source is already stopped. Their remaining path
is therefore identical, including any failure to enter. QED.

The independent packet-source Gram is exactly

```
P_N = sum_s c(s)c(s)^T - S_N^V - (S_N^V)^T
    = (N/7) [ 24 -15  -9 ] + O(log(N+1)+|V|).          (4)
              [-15  42 -27 ]
              [ -9 -27  36 ]
```

This is also the Gram of the literal folded source: replace the two rows
by `c(x)-c(y)` and zero. In particular the actual source baseline is
`P_N[0,0]=24N/7+O(log N+|V|)`, while the trace has **no leading reduction**.

## 4. A complete bound on the initial orientation

Let `W_u(N,d)` be the exact Eisenstein sum of sources whose first arrival
at `u in V` has native depth at most `d`, with the sign of their integer
clock. Keep non-arriving sources outside this stopped sum. Set

```
E(N,d) = sum_u (2 Re W_u)^2,
T(N,d) = 6 sum_u Norm(W_u),
D(N,d) = 3E(N,d)-T(N,d) = 6 Re sum_u W_u^2.
```

At depth zero the only sources are the literal seeds already in `V`:
`D0(N)=D(N,0)=sum_(s in V)[3 c_0(s;N)^2-sum_r c_r(s;N)^2]`.
No orbit computation is needed for this number.

**B1 (PROVED finite reduction).** Put `R=65536>max(V)=63665`. All values
of `D0(N)` for positive integer `N` occur already at `1<=N<64R`.

**Proof.** A point `m=B^j(s)` changes its source contribution at the entry
event `E_j=floor(m/4)+1` and the lower-band event `L_j=(m+1)/2`.
At entry add `lambda=(−1,−1,2)` rotated by `m mod 3`; at the second event
add the rotated `kappa-lambda`, with `kappa=(0,−3,3)`.

Above `R`, entry events have `j>=2` and lower events have `j>=1`.
Their next copies satisfy exactly

```
E_(j+1)=4E_j-2,       L_(j+1)=4L_j-1.                  (5)
```

There are no ties there. Two events of the same type would identify two
primitive `B`-chains. For opposite types, `3E_j-2` and `3L_k-1` are both
divisible by 4: the former is `4^(j-1)(3s+1)`, and the latter is
`4^k(3t+1)/2`, with `t` odd. Equality of the event times would require
these multiples of 4 to differ by 1. The affine maps in (5) preserve the
order of any two distinct integer event times, including mixed types.
Neither type can occur at a block boundary `4^kR`, by the same congruences.
There is exactly one event of each type per seed in each such block.

At the block start, and then at corresponding positions in the event
order, the source rows transform by

```
c_r -> kappa(s+r) + c_(r+1).                           (6)
```

At the start this follows directly by inserting `s` and mapping each old
chain point `m` to `B(m)`; the integer threshold inequalities are exact.
Each subsequent event update is the same cyclic rotation, since
`B(m)=m+1 mod 3`. Three applications of (6) give the identity, because
the three kappa entries sum to zero. Thus the complete vector of 627
source rows repeats after three event blocks, not merely its norm.
The transient `N<R` and the first three blocks suffice. QED.

The same proof works for any finite set of primitive odd seeds, taking a
power of four `R>max(seed)+1`. Forward closure is unnecessary for this
boundary calculation. The repeated state is the full source array, so
the reduction also certifies any statistic of that array, not only D0.

**B2 (PROVED with a finite certificate).** For this exact frontier,

```
-72 <= D0(N) <= 1428  for every integer N>=1.           (7)
```

Both endpoints occur: `N=171` and `N=66475` respectively. The event sweep
below `64R=4194304`, licensed by B1, is a complete certificate. The script
also checks a fourth stable block: 1,254 complete source-vector states
repeat, and 667,128 scalar source comparisons agree. A separate C engine,
using the three integer cutoff updates rather than this event-order
reduction, checks every `1<=N<2^24` and reproduces both extrema and their
first occurrences. The range is special to the declared frontier.

## 5. Two precise remaining concentration questions

Consider, separately, the following claims for `N>=4096`, every finite
native first-arrival cutoff `d`, and this fixed `V`:

```
(OC)  D(N,d) <= D0(N).
(PB)  E(N,d) <= P_N[0,0].                              (8)
```

They are **OPEN**, not consequences of B1 or of the packet identity.
They concern the actual marked source. Neither is an operator norm bound
for arbitrary input rows, nor does it assume that all sources enter.

**CONDITIONAL consequence.** OC and (7) give
`E(N,d)<=T(N,d)/3+476`. If E2's OPEN trace bound
`T(N,d)<=(3/2)T_source(N)` also holds, E3 gives the explicit bound

```
E(N,d) <= (51/7)N + 36J + 479,
J = 1 + ceil(log_4(12N+1)).                            (9)
```

With the stronger trace coefficient 1, the right side becomes
`(34/7)N+24J+478`. PB instead directly bounds energy by the independently
computed source packet baseline (4). These are possible routes to the
existing ROOT-current/density-exponent budget. They are not universal-entry
or universal-Collatz conclusions.

**FINITE-EXACT.** An incremental C scan tests OC, PB, and the inherited
trace controls for every `4096<=N<2^20`, at all first-arrival prefix depths
(maximum 190; maximum integer clock 313), with no violation. Python
constructs the literal source and packets independently at all dyadics
through `2^20`, plus `N=1204,5997`: 523,130 packet identities, 263,168
scalar source checks, 2,228 independent route checks, and 2,219 overlapping
C/Python energy comparisons. Route caps fail explicitly; no source is
silently omitted in these finite experiments.

**Hostile controls.** Dropping `D0` fails already at depth zero:
`D(4096,0)=951>0`. Requiring monotonicity at every depth also fails:
at `N=4096`, `D(31)-D(30)=786`; at `N=5997`, the increase at depth 46
is 1047; at `N=2^20`, the increase at depth 29 is 152910. The survivor
is the cumulative boundary-referenced inequality OC, not stepwise
dissipation. A proof must pay for these genuine increases.

## 6. Mixed quadrilaterals and the completed coefficient space

We incorporated THM-4700's exact quadrilaterals with two sides of length
1 and two of length `1/|k|`, for each `k=b,c`, into the actual coefficient
copy labels. Their unit sides activate cross-coordinate circuits; their
other sides become unit edges after scaling by the corresponding
coefficient. These are closed interior contacts, not just port gluings.

**FINITE-EXACT:** the five factor sizes are `(688,130,10,10,130)`; all
570 canonical patterns on the original eight marked ports still extend.
The literal-edge solver verifies 200,628,600 factor constraints across
its saved colourings. This graph has all its axis unit edges and the
five named contact families; it is not the complete ambient unit graph.

We then executed the proposed six-dimensional completion in
[the earlier circuit note, section 3](golden_contact_circuit_20261017i.md),
using the port-preserving basis

```
(1,u,v,uv,Y,Yv),       Y=-Z/conjugate(Z), Z=d+u,
v^2=Y-(d+u)v,         u^2=a+bu+cY.                     (10)
```

The inherited degrees `[K(u):K]=3`, `[K(u,v):K(u)]=2`, and independence
of `1,u,Y` prove that these six addresses are independent over
`K=Q(zeta15,sqrt(-11))`. The missing `Yv` is retained, not quotiented
into another port or treated as an independent field extension.

Put `r2=m2-b`, `r0=(r2*a+1)/c`, `r1=(r2*b+m1-a)/c`, using the inherited
cubic coefficients `u^3=m2*u^2+m1*u+1`. Multiplication gives

```
uY = r0+r1*u+r2*Y,
u^2 v = a*v+b*uv+c*Yv,
uv^2 = r0+r1*u-a*v-(d+b)*uv+r2*Y-c*Yv.                (11)
```

The last line is an authenticated unit direction involving **all six**
coordinates. We realize its coefficient translations, the `u^2v`
circuit, and the `v`-multiple of the old Y contact, alongside all earlier
families. Every family norm is checked exactly in the cubic field and
quadratic v extension. All eight ports remain `(p,0,0,0,0,0)`.

**FINITE-EXACT:** factor sizes `(1220,250,246,250,250,250)`; eight named
families; all 570 patterns again extend. The large implicit vertex count
is a Cartesian product, not evidence of forcing. The positive witnesses
apply to precisely the declared axis and family edges. No complete
six-coordinate norm-one census, and no general colouring of the full
field, is asserted.

## 7. Exact support compression and a sharper forcing obligation

**G1 (PROVED).** Fix an additive colouring with values in `Z/q`. For a
Cartesian contact family with independent pair sets `P_i` in its active
coordinates, put

```
D_i = { colour_i(y)-colour_i(x) : (x,y) in P_i }.
```

All edges in that family are proper if and only if

```
0 not in D_1 + ... + D_r.                              (12)
```

Proof: each selection of one factor pair yields one contact edge, and its
colour difference is their sum. Conversely every element of that sumset
comes from one of these independent choices. QED. This equivalence is
specific to the Cartesian family and the additive ansatz; it does not
assert that every graph colouring is additive.

The compiler stores one support bit per attainable difference, propagates
sumsets by forward Boolean implications, and forbids the zero bit. This
is exact: present differences force every actual sum; when zero is absent,
assign precisely the true support and sum bits. It preserves every
possible boundary assignment within the additive ansatz. The six-way
family alone has 25,633,184,000 Cartesian constraints, represented by
eleven support groups. The full problem uses 16,422 Boolean variables
and 364,517 clauses. Verification recomputes the sumsets directly from
each saved colouring, independently of the SAT auxiliaries.

**G2 (PROVED).** Call a nonempty subset of `Z/6` rich when it contains
two consecutive residues. If five factors in a contact have rich
difference sets, their sum with any other nonempty factors is all of
`Z/6`; hence no additive colouring can satisfy that contact.

Proof: adding `{a,a+1}` to a nonempty proper subset increases its size
by at least one. Otherwise the subset is invariant under the generating
translation 1, so it is already all six residues. Five increases suffice.
Extra summands cannot undo coverage. QED.

This makes a sharper next obligation than simply connecting copies:
force five rich directional difference sets **simultaneously**, or find
a stronger obstruction to arbitrary colourings outside the additive
ansatz. Existing witnesses must have at least two non-rich factors in
each six-way contact. A non-rich difference set is contained in one
parity coset, or is an antipodal pair, or a singleton: this follows by
listing independent sets of the six-cycle. The mod-2 and mod-3 freedoms
are therefore an exact possible escape, not merely a numerical analogy.

Controls include positive and hostile SAT support examples, exhaustive
literal sumset comparisons over the ternary alphabet through arity six,
and an exact 63-state sumset automaton covering all `63^6` sequences of
nonempty six-colour supports. The new compiler also reproduces the old
five-coordinate family counts and all 570 extensions. It rechecks the
old solver's witnesses against the new support groups, so matching solver
answers alone are not the comparison criterion.

The independent replay checks 137,249 ternary support sequences and
17,670 old-witness/new-support comparisons. Across the 6,270 saved
six-way witness/group combinations, the numbers of rich factors are
0, 1, 2, 3 with counts 420, 5326, 448, 76 respectively. None has four;
five would be impossible by G2. This is a description of these witnesses,
not a maximum over all colourings of the graph.

### Incoming shade controls and a complete local branch cover

The meaningful later pull is `82c9e40041`: [marked snark shade
lifts](snark_boundary_transfer_20261018a.md), [generic cubic retractions
and exceptional phases](cubic_boundary_search_20261018a.md), and the
[thirteen-new-projected-pair filter](exceptional_factor_support_20261018a.md).
The last filter requires a loopless projection to its exact eight anchors;
we do not apply it just by counting physical edges in our larger product.
The cubic packet's complete-graph retractions also reinforce why our
selected-family tests must retain their narrower scope.

The shade-lift mechanism gives a concrete next refinement of G2. Its
pair labels lie in the nonzero elements of `F_2^2`, with a separate binary
shade; these pair labels are **not** `Z/3` flow values. In contrast, the
following uses the actual reductions of our additive colour difference
in `Z/6` modulo 2 and 3. The two constructions share a cycle-consistency
test, not an identification of these label spaces.

**G4 (PROVED).** Every additive six-colouring satisfying one six-way
Cartesian contact belongs to a union of 375 explicitly specified local
congruence branches. A branch chooses two distinct factors and, for each,
one of the five conditions

```
all directional colour differences = a mod p,
(p,a) in {(2,0),(2,1),(3,0),(3,1),(3,2)}.               (13)
```

Proof: G2 leaves at least two non-rich factors. The classification of
non-rich subsets above says exactly that each is contained in a coset
modulo 2 or modulo 3. There are `binomial(6,2)*5^2=375` choices. They
overlap (singletons satisfy two types), but cover every permitted additive
colouring. They are necessary branches, not sufficient colourings. QED.

Each choice gives equations `f(y)-f(x)=a mod p` on **every original
factor pair** in that directional pair set. Several contacts may share
these vertices. For a fixed global branch selection these linear
equations are consistent iff every oriented cycle has total increment
zero modulo p. A spanning-tree potential proves sufficiency; summing a
cycle proves necessity. An inconsistent system returns the literal
cycle of source edges and its nonzero sum. Isolated factor vertices
remain in the state. The remaining unit-edge and contact constraints
still have to be checked after this necessary test passes.

This is an implementable reroute: eliminate periodic-support branches
with authenticated cross-interior contacts, rather than assuming a
single generic coefficient shift. It does not yet eliminate the additive
ansatz, much less all graph colourings. Across many contacts the cover is
not a claim of only 375 global cases; choices must be made consistently
at every contact.

**FINITE-EXACT:** all 63 nonempty difference masks satisfy the classification.
The decoder handles all 570 saved six-coordinate witnesses and eleven
six-way groups, checking 6,840 original-vertex congruence systems with
415,867 equations. It retains one explicit compatible branch selection
per witness. No one local branch is common to all 570 saved witnesses in
any of the eleven groups; the most frequent branches cover between 234
and 353 witnesses. This is witness adaptivity, not a proof that no
alternative common-branch family of witnesses exists. Hostile controls
are a nonzero mod-2 triangle and a nonzero mod-3 triangle; a zero-sum
mod-3 triangle with an isolated vertex is the positive control.

There is also an exact averaging connection with the incoming snark
calculation. Averaging its fixed-pairing defect count over all fifteen
palette pairings gives a constant. Averaging our three residue phases
kills the leading sibling tensor's directional signal. Both operations
project onto an invariant statistic and discard the diagnostic carried
by the marked representation. This explains why restoring the sidecar
is useful; it supplies no transfer of a snark obstruction into a Collatz
energy estimate.

## 8. Petersen minors do not retain the boundary relation

**G3 (PROVED; abstract graph control).** Mark all ten vertices of the
Petersen graph, and subdivide every edge once. The resulting 25-vertex
graph still has a Petersen minor with those ten branch vertices retained.
But every assignment of six colours to the ten marked vertices extends:
each of the fifteen new degree-two vertices can independently choose a
colour different from its two endpoints. The original Petersen graph
does not allow equal colours at adjacent marked vertices.

In particular, the all-equal boundary assignment is a minimal witness to
the loss of information. Contracting the chosen half of each subdivided
edge recovers exactly the original edge set, so the minor survives while
the extension predicate changes. This is an abstract control; an omitted
edge between actual points at distance one cannot be omitted when claiming
a complete geometric unit graph.

The independent audit enumerates all 109,299 canonical assignments of
at most six colours to the ten ports. Exactly 5,649 extend to the
original graph; all 109,299 extend after subdivision.

The relevant substitution invariant is the Boolean boundary relation
`R_G: assignments on ports -> {extends, does not extend}`. Under gluing,
relations combine by conjunction and existential quantification over the
identified internal colours. Equal boundary relations are interchangeable
in such contexts; equal minors are not.

**Typed bridge, not a proof transfer.** Source: the exact signed source or
the marked geometric graph. Target: a smaller representation. Maps: phase
trace, minor contraction, or the exact packet/support compilers. Preserved
by the exact compilers: the intended current or selected-edge boundary
predicate, with stopping clock or factor-pair independence retained.
Destroyed by the coarse maps: source orientation or extendability data.
The cheapest decisive tests are (1), the depth-increase examples, and G3's
all-equal ports. There is no claim that Collatz transport and graph
colouring are isomorphic problems.

## 9. Reproduction and next proof targets

Scripts live in `04-computation/experiments/`; receipts and explicit
colourings in `charge_contexts_20261018k/`. Run from that experiments
directory, with numpy, python-sat, sympy and the existing numba dependency
available for the field engine:

```bash
python3 source_pair_tensor_20261018k.py --max-power 20 --small 1024 --out charge_contexts_20261018k/source_tensor.json
cc -O3 -Wall -Wextra source_concentration_20261018k.c -o /tmp/source_concentration_20261018k
/tmp/source_concentration_20261018k first_arrival_20261018/V_closed_Y1000.txt 20 prefix > charge_contexts_20261018k/source_concentration_everyN.json
python3 frontier_boundary_automaton_20261018k.py --out charge_contexts_20261018k/boundary_automaton.json
cc -O3 -Wall -Wextra frontier_boundary_audit_20261018k.c -o /tmp/frontier_boundary_audit_20261018k
/tmp/frontier_boundary_audit_20261018k first_arrival_20261018/V_closed_Y1000.txt 24 > charge_contexts_20261018k/boundary_everyN.json
python3 golden_reflection_context_20261018j.py --cubic --links 0 --mixed-quads 2 --out charge_contexts_20261018k/mixed_quads.json
python3 six_coordinate_context_20261018k.py --mixed-quads 2 --out charge_contexts_20261018k/support_compiler_control.json
python3 six_coordinate_context_20261018k.py --six --dense-six --mixed-quads 2 --out charge_contexts_20261018k/six_dense.json
python3 interface_audit_20261018k.py --directory charge_contexts_20261018k
python3 periodic_support_branches_20261018k.py --input charge_contexts_20261018k/six_dense.json --out charge_contexts_20261018k/periodic_branches.json
```

The next source target is a cumulative orientation-debt certificate for
OC that pays the observed positive increments; a stepwise sign rule is
already excluded. PB is a parallel, scalar target with a source-only
baseline. The boundary term is now completely controlled, so neither
target should spend effort bounding that term again.

The next geometric target is a small, exact context that forces the
required directional support richness or eliminates all consistent
periodic branches, checked by its marked extension relation. More copies
of an unchanged port-full context are insufficient.
The support compiler makes additional authenticated unit families cheap
to test without pretending that family membership is a complete census.
