# Indexed shared square clocks have one signed integer coordinate per component

**PROVED:** an exact compiler and counting formula for a finite network of
paths with prescribed ordered adjacent sums, interval hosts, and actual
indexed shared ports. Optional collision guards are retained explicitly.
**FINITE-EXACT:**15,276 specified clock-network cases compared against direct
enumeration, plus adversarial controls. This is an all-size rule for this
signed-affine carrier, not for unrestricted tournament wiring or arbitrary
unlabelled path attachment.

## 1. Inheritance, input, and what is new

Incoming commit `cf7c4fc327` supplies
[SQUARE_WORD_RECONSTRUCTION.md](../../session42_shared_lattices_square_colors_20261023/SQUARE_WORD_RECONSTRUCTION.md)
and its
[independent audit](../../session42_shared_lattices_square_colors_20261023/audit-square-words/REVIEW.md).
They prove the affine reconstruction of a single clock word, its interval-minus-
holes start domain, odd/even Hamiltonian fiber bounds, and the exact alternating-
block classification of two realizations. The first square-clock twin is at36.
These single-path results are inherited, not new discoveries here.

The new step is to compile **all indexed shared ports together**, propagate a
single signed integer root per connected component, and count the compatible
starts after pulling back every bound, collision hole, cycle, and absolute mark.
Marginally feasible starts cannot be chosen independently at shared ports.

Input for each named path a is:

- a finite ordered list `s_(a,1),...,s_(a,k_a)` of positive adjacent sums;
- an integer host `[1,N_a]`, with `N_a>=1`;
- the requirement that its own `k_a+1` vertices are distinct;
- shared-port equalities of actual pairs `(path a, position i)`;
- optionally, absolute indexed marks `x_(a,i)=v` and explicit inequalities
  between pairs of indexed ports.

For square clocks, each sum is the square of its supplied positive root.
The theorem also holds for arbitrary prescribed sums. Different paths may
intersect unless an inequality explicitly forbids it. No unindexed endpoint
or unspecified vertex equality is silently converted into an indexed one.
If `k_a+1=N_a`, the local distinctness and host conditions make the path
Hamiltonian; otherwise it is a partial path.

## 2. The inherited local domain

Set `c_(a,0)=0` and `c_(a,i)=s_(a,i)-c_(a,i-1)`. Every realization is exactly

```
x_(a,i) = epsilon_i t_a + c_(a,i),     epsilon_i=(-1)^i.
```

The host conditions give an integer interval `[L_a,U_a]` for `t_a`.
Two equal-parity positions collide identically iff their constants agree;
that makes the whole local domain empty. Otherwise an opposite-parity pair
excludes at most one integral start. Let `E_a` be the set of these excluded
starts inside the interval, with repetitions removed. Thus the exact local
domain is either empty or

```
D_a = [L_a,U_a] intersect Z, minus E_a.
```

This includes adjacent equal vertices: a square half-diagonal cannot be used
as a loop in the simple graph. Keeping the square label and integer interval
alone does not enforce this condition.

## 3. Signed propagation and the exact component theorem

A shared indexed port `x_(a,i)=x_(b,j)` becomes

```
t_b = eta t_a + q,
eta = epsilon_i epsilon_j,
q = epsilon_j (c_(a,i)-c_(b,j)).                       (1)
```

Form the undirected multigraph on named paths with one edge for each such
equality. Parallel edges and loops retain distinct incidence obligations.
Choose a root and spanning tree in each connected component. Propagating(1)
along the tree yields

```
t_a = sigma_a r + b_a,     sigma_a in {+1,-1},
sigma_root=1, b_root=0.                               (2)
```

Every original edge, including a non-tree edge, now requires

```
(sigma_b-eta sigma_a) r = eta b_a+q-b_b.               (3)
```

Its left coefficient is0 or+/-2. A zero coefficient requires zero right
side; otherwise the component is inconsistent. A nonzero coefficient fixes
one rational root value, which must be an integer. All such pinned values
must agree. An absolute indexed vertex mark similarly pins r, with coefficient
`epsilon_i sigma_a` equal to+/-1.

Next pull back every local domain through(2). When `sigma_a=1`, its root
interval is `[L_a-b_a,U_a-b_a]`; when `sigma_a=-1`, it is
`[b_a-U_a,b_a-L_a]`. A local hole e pulls back to
`sigma_a(e-b_a)`. Intersect all root intervals to obtain `[L,U]`, and take
the **union** of all pulled-back holes, clipped to that interval, to obtain H.

**Theorem.** If no inconsistency or empty local domain has occurred, the
compatible configurations of this connected component are in bijection with

```
([L,U] intersect Z) minus H,
```

further restricted to the unique pinned root value if one exists. Hence a
free component has exactly

```
max(0,U-L+1) - |H|
```

configurations, with H understood to be empty when `L>U`. A pinned component
has exactly one or zero configurations according as its integral pin belongs
to the displayed domain. Distinct components multiply **only when no
additional constraint couples them**.

**Proof.** Every compatible configuration satisfies the tree equations, hence
has the form(2), with the root's own start equal to r. All non-tree equations
are exactly(3); no relation is omitted. Substituting(2) into the local domain
conditions gives precisely the intersected intervals and union of holes.
Conversely any root passing these tests gives all starts by(2), then every
vertex value by the clock formula. It satisfies every original shared port,
every local bound and every local distinctness condition. The root start
recovers r, so this correspondence is injective as well as surjective. QED.

The graph/trees are computation devices, not quotients of the actual path
owners. Choosing another root or tree changes the coordinate by a signed
integer affine substitution and leaves the realized configurations unchanged.
There are polynomially many stored local holes (at most one per position
pair), and propagation uses the explicit incidence graph. No Cartesian
enumeration of all local start domains is needed for this component count.

## 4. Explicit cross-path distinctness guards

If two indexed ports lie in the same shared-port component, substituting(2)
into the forbidden equality `x_(a,i)=x_(b,j)` gives another equation
`A r=B`, with `A=0` or+/-2. To **forbid** this equality:

- when `A=0,B=0`, reject the component;
- when `A=0,B!=0`, there is no restriction;
- when `A!=0`, exclude `B/A` only if it is integral.

Thus explicit same-component distinctness adds holes or detects an impossible
identification; the same count theorem applies.

For ports in different components, the forbidden equality has the form

```
alpha r_c - beta r_d = gamma,    alpha,beta in {+1,-1}.
```

It is a joint forbidden diagonal in two component domains. The compiler
retains it as such; the simple product count is no longer asserted. For
example two singleton paths in `[1,3]` have3 choices each, but requiring their
values to differ leaves6 configurations, not9. Additional finite factor
contraction or enumeration may count this coupled model. It cannot be
replaced by deleting a fixed unary hole independently in each component.

## 5. Positive controls and the decisive failures

The receipt retains every input, root interval, hole, pin, and complete start
list for the following examples.

**Marginal feasibility does not synchronize two ports.** Take one edge of
sum4 and another of sum16 in host `[1,15]`. Requiring their first vertices to
agree is feasible. Requiring their second vertices to agree is separately
feasible. Requiring both is impossible, because the same ordered endpoints
would have to sum to both4 and16. The parallel shared edges expose a
nonzero translation inconsistency in(3).

**A fixed congruence register loses integer clocks.** For any modulus M, the
valid square edges `(1,3)` and `(1,(M+2)^2-1)` agree at both indexed ports
modulo M, but cannot agree at both ports as integers. This is an obstruction
to dropping actual clocks in favor of a fixed modulus. It does not invalidate
an exact finite cyclic-coordinate carrier with its full coordinate retained.

**An odd cycle can pin a noninteger.** Glue the end of each of three one-edge
paths to the start of the next, with sums4,9,16. Every path separately has
valid realizations in `[1,15]`, but closure pins the first value to11/2.

**An integral pin can hit a collision hole.** The same cyclic gluing with
sums4,9,9 and host `[1,8]` pins the vertices to `2,2,7`. The pin is integral
and lies in every interval, but the first edge would be a forbidden loop.
The pulled-back local hole correctly rejects it.

**A real square triangle tests the positive direction and the bounds.** Sums
25,49,36 under the same cyclic gluing give the unique indexed vertices
`6,19,30`. Host `[1,30]` accepts; `[1,29]` rejects, even though each edge still
has valid local realizations. Retaining only the affine closure equation
would miss the host failure.

**A free component really can have many choices.** Glue the end of a sum9
edge in `[1,8]` to the start of a sum16 edge in `[1,12]`. The root interval is
`[1,5]`, with hole1, leaving the four start pairs
`(2,7),(3,6),(4,5),(5,4)`. Hole1 is the sum16 half-diagonal `(8,8)`.

**Pairwise compatibility of filtered domains is insufficient.** The lattice
lane's independent audit supplied three actual local domains: a singleton
path in `[1,2]` gives `{1,2}`; a sum4 edge in `[1,3]` gives `{1,3}`; and a
sum9 edge in `[1,7]` gives `{2,3,4,5,6,7}`. Identify the start port of every
path. Each two-path assembly is feasible, with common starts1,2,3 for the
three respective pairs, but their joint intersection is empty. Thus a
pairwise-Helly theorem for pure cyclic cosets cannot be imported after these
interval and collision-hole filters. The full component intersection is
essential even when the shared-port support is a tree.

**The incoming full-word twin is repaired by any indexed mark.** The order36
clock word from the incoming certificate has exactly two starts,3 and9.
The compiler recovers both. Marking the actual value of the first realization
at any one of its36 specified positions selects start3 uniquely. An unindexed
claim that some vertex equals17 does not help: both Hamiltonian realizations
already contain17. The specific position is the needed address.

## 6. Relation to the other lanes

The incoming
[square-color packet](../../session42_shared_lattices_square_colors_20261023/square-colors/RESULTS.md)
at `cce2774f3c` independently proves the same sharp71 founder-saturation
threshold as [our report](REPORT.md), and strengthens the ancestry carrier
with exact parent-forest multiplicities, a weighted joint partition function,
and harmonic marginal probabilities. Our earliest22 joint-feasibility hostile
and pinned three-vertex proper-coloring move remain distinct additions. The
incoming feature-color scan concerns **path observation**, whereas our95
candidate screen concerns **proper graph coloring**; their success/failure
counts are not contradictory.

Here the root r plays the role of the single retained component phase, but
its domain is an integer interval with collision holes, not an arbitrary
cyclic group. A loop in the shared-port network either checks translation or
pins twice r. Actual path positions determine the signs and offsets. This is
a concrete answer to the Fibonacci address warning: local lengths or local
support sets do not identify the indexed positions at which paths are glued.
It is also a concrete low-dimensional coupled-move carrier: changing one
root may change all vertices of arbitrarily many attached paths while
preserving every shared equality. Legality still requires staying in the
root domain and retaining all optional cross-component factors.

## 7. Reproduction and scope

```
python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/square_sum/check_shared_clocks.py
python -B -O 05-knowledge/results/attachment_lattices_square_sum_20261023/square_sum/check_shared_clocks.py
```

The checker compares the compiler against an independent direct reconstruction
for all15,276 cases in these explicitly finite universes:

- every clock word with0..3 edges, roots2..5, and host1..12;
- every subset of the four possible indexed equalities between two one-edge
  paths, with sums and common host each chosen from `{4,9,16}`;
- three cyclically glued one-edge paths with each sum in `{4,9,16}` and
  independently chosen hosts1..8.

Additional controls test all positions of the inherited36 twin, absolute
marks, local permanent collisions, same-component and cross-component
distinctness, and the displayed positive/hostile cases. Normal and optimized
runs each pass30,743 checks, including comparison with
[shared_clocks_receipt.json](shared_clocks_receipt.json). These finite checks
support the implementation; the theorem is proved by the signed propagation
and exact pullback argument above.
