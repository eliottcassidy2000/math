# LRC(14): owner reuse, located LCM seams, and the affine holonomy interface

**Session:** codex-2026-07-19-S78 owner-reuse probe

**Status:** research memo.  It proves several new conditional lemmas about a
hypothetical six-comb cover and isolates one operation-level continuation
lemma.  It does **not** prove six-comb noncoverage, LRC(14), or uniform
emptiness of the `n=12` sporadic branch.

## 0. Verdict

The strongest immediate bridge is more rigid than THM-1244's two-edge debt.
THM-1198 forces every one of the six fast labels to own private mass on the
full slow gap.  An irredundant chain of individual strict teeth must therefore
visit all six labels.  Its consecutive handoffs contain a spanning tree of
five **actual interior overlaps**, and every one has the label-aware quantum

```text
gcd(u,v)/(14uv)=1/[14 lcm(u,v)].
```

Consequently a putative cover forces the fully located clock debt

```text
H_fast >= 1/c+(7/12) sum_(uv in T) 1/lcm(u,v).        (A)
```

There are no unlocated pair-period errors in (A).  This is a genuine
scale-covariant strengthening of THM-1178 and is the global rank-five
backbone to which THM-1244's protected private stalk should be attached.

Private mass also forces literal owner reuse.  THM-1244's selected owner `b`
occurs in more than `308b/(3705c)` private tooth cells.  Thus

```text
b/c >=3705/154
```

forces at least three private `b`-teeth and hence an interior private stalk.
For the full-gap private needle of any fast label `d`, the sharper condition
is simply `d>14c`.  An interior stalk has one of two exact forms:

1. distinct left/right wall owners give a rank-two `u-b-v` seam rooted on
   the same private interval; or
2. the same owner occurs at both ends, in which case the stalk is exactly a
   complete safe gap of that owner contained in one `b`-tooth, forcing the
   toothpick separation `u>6b`.

The remaining obstruction is now algebraically precise.  Chronological
handoffs are relations on the fast speeds `d_i`; THM-1240's positive
holonomy is naturally a relation on the shifted clocks `Q_i=c+d_i`.  The
affine shift can reverse the sign.  The exact mixed-circuit identity below
shows how the lost sign is paid by reciprocal tooth-centre drift.  What is
still unproved is a cover-preserving lift of an oriented private germ around
the blocker cycle.  If that lift closes, the address factors telescope and
THM-1240's product inequality gives a genuine integer descent.  If it first
fails, the desired replacement theorem should turn the interrupter into a
new disjoint located handoff, increasing the LCM debt in (A).

## 1. Setup and conventions

Let

```text
D_s={t in R: ||st||<1/14},
G=[(14k+1)/(14c),(14k+13)/(14c)],
|G|=6/(7c),
c<d_1<...<d_6,                                         (1)
```

and suppose the six strict danger combs `D_(d_i)` cover `G`.  Translate the
circle lift so that `0<=k<c`; then `G` lies strictly inside `(0,1)`.

The strict convention matters topologically.  Closed teeth may cover by a
zero-length meeting.  A finite open cover of the compact interval `G`, after
deleting redundant teeth, has strict positive consecutive overlaps.  Measure
identities are unchanged at the finitely many walls.

THM-1198 gives, for every fast label `d_i`, a private-provider set

```text
Q_i subset G intersect D_(d_i),
Q_i intersect D_(d_j)=empty for j!=i,
|Q_i|>=1/(49c).                                        (2)
```

On THM-1244's protected slowest-spoke component `K`, some label
`b in {d_2,...,d_6}` has private mass

```text
|U_b|>44/(3705c).                                      (3)
```

Finite wall sets can be removed from `Q_i,U_b` without changing their
measure.  Their positive components may therefore be treated as ordinary
open intervals with one-sided labelled wall germs.

## 2. Proved now: all six private needles force a fully located LCM tree

### 2.1 The chronological chain

Choose a finite subcover of `G` by individual open teeth and delete teeth
until it is irredundant.  Order the selected teeth by their left endpoints:

```text
J_a=(alpha_a,beta_a),              a=1,...,N.          (4)
```

Irredundancy excludes containment.  Hence, after removing impossible ties,

```text
alpha_1<...<alpha_N,
beta_1 <...<beta_N,
alpha_(a+1)<beta_a.                                    (5)
```

For `a<N`, the right endpoint `beta_a` is before the right endpoint of `G`,
and for `a>1`, `alpha_a` is after the left endpoint of `G`; otherwise later
or earlier teeth would be redundant.  Thus every consecutive overlap

```text
W_a=(alpha_(a+1),beta_a)                               (6)
```

is the full raw tooth-to-tooth overlap and lies in `int(G)`.  It is not an
arbitrary pair intersection clipped by a slow-wall endpoint.

Every private set `Q_i` contains points covered by no other label.  Therefore
the selected word of speed labels visits all six labels.  Consecutive labels
are distinct because distinct teeth of one speed are disjoint.  The graph of
unordered consecutive label pairs is consequently connected.  Choose a
five-edge spanning tree `T` in that graph and retain one witnessing handoff
for each edge.

### 2.2 Every tree edge has its own LCM clock quantum

Suppose a retained handoff goes from the tooth of speed `u`, address `n`, to
the tooth of speed `v`, address `m`.  Formula (6) gives

```text
omega
 = (14n+1)/(14u)-(14m-1)/(14v)
 = [v(14n+1)-u(14m-1)]/(14uv)>0.                      (7)
```

The positive integer numerator is divisible by `gcd(u,v)`.  Hence

```text
omega>=gcd(u,v)/(14uv)=1/[14 lcm(u,v)].               (8)
```

This qualifier is essential.  An arbitrary component of
`G intersect D_u intersect D_v` can be clipped by an endpoint owned by `c`
and need not have the label-gcd quantum.  The irredundant consecutive chain
is what removes that flaw.

### 2.3 Hunter consumes the five actual intervals

Put `H_fast=sum_i 1/d_i`.  Forest Hunter on the fixed tree `T` gives

```text
0 >= |G|-sum_i mu(G intersect D_(d_i))
       +sum_(uv in T)mu(G intersect D_u intersect D_v). (9)
```

The sharp singleton discrepancy used in THM-1237 is

```text
mu(G intersect D_d)<=|G|/7+6/(49d).                  (10)
```

Using `|G|=6/(7c)` in (9)--(10) proves

```text
H_fast
 >=1/c+(49/6)sum_(uv in T)mu(G intersect D_u intersect D_v)
 >=1/c+(7/12)sum_(uv in T)1/lcm(u,v).                (11)
```

Equivalently, with `delta=cH_fast-1`,

```text
delta>=(7c/12)sum_(uv in T)1/lcm(u,v).                (12)
```

This proves (A).  THM-1178 had the same Hunter architecture with the
phase-free quantum `1/(14uv)`.  Equation (12) promotes it to the actual
common-clock quantity `1/lcm(u,v)` on one specially located tree.  Under a
common dilation, both sides of (12) are invariant.

### 2.4 Multiplicity version

Let `m_uv` count consecutive handoffs with unordered label pair `{u,v}`.
Two distinct such handoff intervals are disjoint: otherwise the common point
would lie in two distinct teeth of at least one of the two fixed speeds.
Therefore

```text
mu(G intersect D_u intersect D_v)
 >=m_uv/[14 lcm(u,v)].                                (13)
```

The positive support of `m_uv` is connected.  Give edge `uv` weight
`m_uv/lcm(u,v)`.  Every edge of `K_6` lies in one third of its labelled
spanning trees, so a maximum-weight tree `T_*` satisfies

```text
H_fast
 >=1/c+(7/12)sum_(uv in T_*)m_uv/lcm(u,v)
 >=1/c+(7/36)sum_(u<v)m_uv/lcm(u,v).                 (14)
```

Thus owner reuse has an immediate scalar consumer: a repeated handoff is
not merely another word symbol; it adds disjoint pair mass on the same
Hunter edge.

### 2.5 Relation to THM-1244

The full-gap tree does not make THM-1244 obsolete.

* On `G`, all six fast labels appear by their global private needles, so
  rank five and five located seams are available.  But `d_1` is not safe on
  all of `G`, and the tree by itself does not place the S82 seven-wall safe
  reservoir or its deletion circuit.
* On `K`, both `c` and `d_1` are safe.  Only five blockers remain, at least
  three are forced, and THM-1244 locates a macroscopic private stalk relevant
  to the protected seven-wall problem.

The useful combined carrier is therefore a full rank-five LCM tree on `G`,
rooted at the THM-1244 private label and carrying its rank-two `K`-subforest
when the germ dichotomy below is in the distinct-owner branch.

## 3. Proved now: private mass forces address-cell reuse

A danger tooth of speed `b` has length `1/(7b)`.  If a private set of mass
`M_b` meets `r_b` distinct `b`-teeth, then

```text
r_b/(7b)>=M_b.                                        (15)
```

For THM-1244's selected owner, (3) gives

```text
r_b>308b/(3705c).                                     (16)
```

Consequently

```text
b/c>=3705/154  ==>  r_b>=3.                           (17)
```

At most two private teeth can touch the two endpoints of `K`, so (17)
forces an interior private tooth and an interior private component.

For every full-gap private needle (2), the stronger uniform count is

```text
r_(d_i)>=ceil(d_i/(7c)).                              (18)
```

In particular,

```text
d_i>14c  ==>  r_(d_i)>=3.                             (19)
```

This is literal owner reuse across address cells.  Any irredundant subcover
must contain a tooth at every such private address, since no other label can
cover the private mass there.  The low alternatives `b/c<3705/154` and
`d_i<=14c` are substantially smaller ratio charts, though not finite
absolute-address charts by themselves.

The already-proved projective ceiling `d_i/c<2345` gives a complementary
uniform extraction without a high/low split.  Each label has at most `2011`
teeth meeting `G`; endpoint sweeping therefore extracts an oriented
private interval of length `const/c`.  The useful new point in (17)--(19) is
not merely a tiny uniform component: in the high branch the same owner is
forced into at least three different tooth addresses.

## 4. Proved now: an interior private stalk is either a rooted two-seam V or
an exact toothpick

Let `J=(x,y)` be a component of the interior of a private `b`-set, contained
in an interior `b`-tooth.  At `x`, some other danger tooth must cease to be
active; otherwise private ownership would continue to the left.  At `y`,
some other tooth must begin to be active.  Because the six strict combs cover
the ambient interval, neither endpoint can be a bare boundary of the
`b`-tooth.  We may therefore choose labels and addresses such that

```text
x=(14m+1)/(14u),             left tooth u exits,
y=(14n-1)/(14v),             right tooth v enters.    (20)
```

The signs `+` at the left and `-` at the right are the oriented off-grid
germ data lost by a residue mask.

### Distinct endpoint owners

If one can choose `u!=v`, then

```text
|J|=[u(14n-1)-v(14m+1)]/(14uv)
    >=gcd(u,v)/(14uv).                                (21)
```

Immediately to the left of `x`, the `u`- and `b`-teeth overlap; immediately
to the right of `y`, the `b`- and `v`-teeth overlap.  Since these events are
interior, their pair components have the exact lower bounds

```text
omega_(ub)>=gcd(u,b)/(14ub),
omega_(bv)>=gcd(b,v)/(14bv).                          (22)
```

The two distinct edges `{u,b},{b,v}` form a rank-two forest rooted at the
macroscopic private owner.  Extend it to the full tree in Section 2.  This
is the requested direct coupling of THM-1244's private mass and its two
located seams; no anonymous tree selection is needed in this branch.

### One endpoint owner on both sides

If no distinct choice is available, the same label `u` supplies both walls.
It is dangerous immediately left of `x` and immediately right of `y`, and
safe throughout `J`.  The two teeth must be consecutive, so

```text
y-x=6/(7u).                                           (23)
```

Moreover the closed `u`-safe gap `[x,y]` lies strictly inside one
`b`-danger tooth.  Comparing widths gives

```text
6/(7u)<1/(7b),                 hence u>6b.             (24)
```

This is an exact toothpick/self-similar branch, not a failure to identify
endpoint owners.  It already pays a macroscopic ratio cut.  Any closing
argument should recurse or reroute this complete safe gap rather than
discarding it as a same-label degeneracy.

## 5. Proved now: every repeated owner carries positive address holonomy

Write a consecutive part of the irredundant chain as teeth

```text
(s_0,n_0),(s_1,n_1),...,(s_r,n_r).                   (25)
```

In the canonical lift every address satisfies `0<n_i<s_i`.  Indeed, a
fast tooth with address zero cannot meet a `c`-safe point strictly beyond
`1/(14c)`, and the reflected argument excludes address `s_i` near one.

For one handoff put

```text
h_i=n_(i+1)s_i-n_i s_(i+1).                           (26)
```

The strict ordering of both endpoints and their positive overlap gives the
exact band

```text
|s_i-s_(i+1)|<14h_i<s_i+s_(i+1),       h_i>0.         (27)
```

Its overlap and reciprocal-centre drift are complementary:

```text
omega_i
 =[s_i+s_(i+1)-14h_i]/[14s_i s_(i+1)],               (28)

h_i/[n_i n_(i+1)]
 =s_i/n_i-s_(i+1)/n_(i+1)>0.                         (29)
```

Equation (28) is the located Hunter credit.  Equation (29) is the functional
form of the chronological `H`-drift.  A handoff cannot make both quantities
small without paying through its speed/address scales.

If `s_r=s_0` is the next occurrence of a repeated label, chronology forces
`n_r>n_0`.  The path relations have the exact positive holonomy identity

```text
(n_r-n_0)s_0 product_(j=1)^(r-1)n_j
 =sum_(i=0)^(r-1) h_i
    (product_(j<i)n_j)(product_(j>i)n_(j+1))>0.       (30)
```

Equivalently the product holonomy of the addressed return is

```text
product p_i/product q_i=n_r/n_0>1.                    (31)
```

Thus the high-owner branches (17) and (19) do not merely contain several
private cells.  They contain at least two consecutive owner-return segments,
each with a canonical positive address holonomy.  Reversing such a segment
would be an integer descent, but the legality of that reverse transport is
exactly the missing metric continuation theorem.

## 6. The no-reuse branch meets the blocker cycle in a typed circuit

If the irredundant word has exactly six teeth, every label occurs exactly
once.  Its chronology is a Hamiltonian path

```text
lambda_1 -> lambda_2 -> ... -> lambda_6.              (32)
```

Choose a THM-1240 blocker functional graph.  Every selected blocker cycle
has an edge directed backward in the order (32).  A backward edge from the
later spoke label `s_r` to the earlier blocker label `s_0`, together with
the forward chronological subpath from `s_0` to `s_r`, is a typed mixed
circuit.  If the backward edge is adjacent, this is already a directed
two-cycle; a tournament completion would erase the fact that its two
directions have different meanings.

This gives a clean combinatorial dichotomy:

```text
N>6:  a repeated-owner return with positive address holonomy;
N=6:  a mixed chronological/blocker circuit.          (33)
```

The second branch is where the affine carrier shift becomes decisive.

## 7. Proved now: the exact affine mixed-circuit identity

Take the chronological subpath in Section 6, with endpoint teeth
`(s_0,n_0)` and `(s_r,n_r)`.  At the centered spoke of `s_r`, write

```text
Q=c+s_r,                 t=P/Q,
r_*=P s_0-NQ,            |r_*|<Q/14,                 (34)
```

where `N` is the blocker tooth address of `s_0`.  THM-1240/1226 give the
positive shifted-clock quantity

```text
K=P(c+s_0)-N(c+s_r)=Pc+r_*>0.                         (35)
```

On the unshifted fast speeds, however, the closing residual is

```text
R=P s_0-Ns_r=Nc+r_*=K-(P-N)c.                         (36)
```

Its sign is not controlled by (35).

Let

```text
Delta=sum_(i=0)^(r-1) h_i/[n_i n_(i+1)]
     =s_0/n_0-s_r/n_r>0.                              (37)
```

Direct expansion, equivalently the relation-cycle identity with the
chronological path inserted, gives

```text
R=N n_r Delta+(s_0/n_0)(P n_0-N n_r).                (38)
```

Hence

```text
P n_0-N n_r
 =(n_0/s_0)(R-N n_r Delta).                           (39)
```

This is the exact interface between the positioned handoff stalk and the
centered positive holonomy.  It yields a rigorous sign/invoice dichotomy:

```text
P n_0<N n_r,
or
R>=N n_r Delta
  =N n_r sum_i h_i/[n_i n_(i+1)]>0.                  (40)
```

In the second branch the affine residual pays the entire reciprocal-centre
drift.  Equations (28) and (40) expose the prospective closing mechanism:
every path edge pays either overlap mass through `omega_i` or drift through
`h_i`; the remaining task is to aggregate those complementary payments on a
cover-preserving circuit.

There is also a useful scale consequence.  If `R<=0`, then
`-r_*>=Nc`, while strict danger gives `Q>14|r_*|`.  Thus

```text
R<=0  ==>  Q>14Nc  ==>  s_r>(14N-1)c.                (41)
```

For `N>=2`, this is already `s_r>27c` and forces the global repeated-owner
branch (19).  For `N=1`, only the thin ratio window `13c<s_r<=14c` can avoid
that same conclusion.

The sign reversal is real, not hypothetical.  Take

```text
c=5, k=0, s_r=82, s_0=9,
Q=87, P=9, N=1, r_*=-6.                              (42)
```

Then `P/Q=3/29` is the centered spoke in the canonical `c`-gap and speed
`9` is dangerous there.  Nevertheless

```text
K=Pc+r_*=39>0,                 R=Nc+r_*=-1<0.         (43)
```

Therefore no proof may copy the positive `Q=c+d` holonomy sign onto a
chronological relation on the `d`-vertices.  The affine drift (36)--(40) is
not an error term to discard; it is the functional coordinate that must be
transported.

## 8. Candidate closing statements

The following two formulations are precise enough to falsify and narrow the
remaining work.

### Candidate A: closed germ lift or a new located handoff

Choose a directed centered blocker cycle `C`.  Every cycle label has a
global private needle by THM-1198; root the construction at THM-1244's
protected private stalk when that label is available.  For each edge
`i -> j`, attempt to continue the outgoing oriented germ of owner `j` from
the spoke cell into the next private address cell.

> Either the whole cycle lifts to a closed sequence of literal oriented
> private germs with positive address factors `a_i`, or the first failed
> continuation exposes an additional interior handoff interval, disjoint
> from the handoff occurrences already charged for its unordered label pair.

In the closed case the desired potential comparison is automatic:

```text
product_i (N_i a_(i+1))
 < product_i (P_i a_i),                               (44)
```

because the `a_i` telescope and THM-1226 gives
`product P_i>product N_i`.  In the failed case, (13)--(14) consume the new
handoff as additional located LCM mass.  What is open is the topological
assertion that a first failed germ really produces a *new disjoint* handoff
rather than returning on a previously charged wall with changed lift height.
That is why the full exact kernel/gcd sheet and canonical gap address must be
part of the state.

### Candidate B: complementary drift/overlap closure on a mixed circuit

For a backward mixed circuit, use the exact split (28), (38):

```text
handoff pair sum = 14 h_i +14s_i s_(i+1)omega_i,
mixed sign       = R-Nn_r sum_i h_i/(n_i n_(i+1)).    (45)
```

The candidate quantitative lemma is:

> One can select a backward circuit and a spanning tree extension so that
> either its mixed sign is compatible with the legal direction of an
> oriented private germ, or the complementary overlap terms `omega_i`, with
> repeated-edge multiplicities retained, make the located tree inequality
> (14) exceed the actual harmonic slack `H_fast-1/c`.

This is not yet proved.  It is substantially narrower than an unspecified
address-compression theorem: all quantities in (45) are exact integers or
rational interval lengths, and the only qualitative input still missing is
which germ direction is legally cover-preserving.

## 9. Obstructions and counterexamples that the lemma must survive

1. **The `q=15` blocker-complete guardrail (THM-1247).**  Full masks, a legal
   Kakeya cut, every carrier spoke blocked, every positive fast-fast sum beat
   blocked, and a nontransitive tournament coexist with a lonely phase.
   Thus Candidate A must use off-grid germ orientation and owner reuse, not
   sampled blocker incidence.
2. **The resonant external star (THM-1239 and its multi-cell extensions).**
   One `14a`-type label can block many cracks or all six sampled spokes.
   There is no scale-free blocker-degree contradiction.  Continuous private
   teeth and handoff multiplicity survive this collapse.
3. **Arbitrary-gap survival is false.**  A selected toothpick gap may be
   completely erased while another address cell contains the lonely phase.
   The transport is allowed to change gaps; it must retain the global cell
   address when it does.
4. **Primitive normalization is mandatory (MISTAKE-184).**  Multiplying a
   relation cannot manufacture new samples or `1/y` discrepancy.
5. **Low-height phase data are not an inverse theorem (MISTAKE-185).**  The
   exact high kernel and torsion sheet can freeze a finite interval even when
   every low-height detuning crosses it.
6. **The full tree is existential.**  It is the tree extracted from the
   actual chronological word, not an independently optimized numerical tree.
   The multiplicity maximum in (14) is legitimate only after the disjoint
   handoff intervals have been located.
7. **Positive shifted holonomy is not positive fast-speed holonomy.**  The
   exact row (42)--(43) rules out this tempting sign shortcut.

## 10. Exact computational experiment

The next computation should enumerate the proof object, not only speed
packets.

### 10.1 State and constraints

For each canonical gap and compact ratio row retain

```text
(c,k; d_1,...,d_6;
 irredundant tooth word (label,address);
 private components with signed endpoint owners;
 handoff lengths omega_i, drifts h_i, multiplicities m_uv;
 located maximum LCM tree;
 centered spokes (P_i,Q_i), full blocker sets and residues r_ij;
 exact primitive kernel/gcd/torsion sheets).           (46)
```

The word constraints can be imposed directly by exact inequalities

```text
alpha_i<alpha_(i+1)<beta_i<beta_(i+1),                (47)
```

together with endpoint coverage and the requirement that all six labels
occur.  This searches exact arrangement cells and can reject impossible word
types before enumerating large absolute speeds.  THM-1233 bounds each label
to at most `2011` teeth on `G`, so the occurrence carrier is uniformly
finite even though its addresses and common scale are not.

### 10.2 Primary measurements

For every feasible row report:

* `delta=cH_fast-1` and the base/multiplicity LCM debts (12), (14);
* owner occurrence counts, interior-stalk type (`u!=v` or exact `u>6b`
  toothpick), and the number of repeated-owner returns;
* all selected blocker cycles and every backward mixed circuit;
* the sign triple `(K,R,Pn_0-Nn_r)` and the exact ratio
  `R/[Nn_r Delta]` when defined;
* which path edges pay mostly through `omega_i` and which through `h_i`;
* whether a literal oriented germ lift exists across each blocker edge;
* if a lift fails, whether the interrupting handoff is genuinely new and
  disjoint for purposes of `m_uv`.

The first falsification targets should be:

1. “the THM-1244 private owner lies on every blocker cycle” — there is no
   present reason to believe this;
2. “`K>0` implies `R>0`” — already refuted by (42)--(43);
3. “some backward edge always has a prescribed mixed sign” — test both
   signs before building a proof around one;
4. Candidate A's new-disjoint-handoff conclusion;
5. whether every same-owner endpoint row really enters an already closed
   macroscopic-cut or transverse-corridor branch.

Normal and optimized exact-arithmetic runs should be byte-compared.  Any
finite census must state the scale/address bound it covers; it is a
falsifier and pattern extractor until a descent makes that bound uniform.

## 11. Tournament and alternate-carrier audit

The carrier choice changes what “cycle” means.

* **Runners with speed order.**  The tournament is transitive, with score
  histogram `(0,1,2,3,4,5)`, no directed triangles, six singleton SCCs, and
  one Hamiltonian path.  It loses every address, overlap, and private region.
* **One-tooth-per-owner chronology.**  In the `N=6` branch it gives the same
  transitive fingerprint, now with the Hamiltonian path (32).  Adding one
  backward blocker obligation creates a typed SCC of size equal to the
  intervening path segment; an adjacent backward edge gives a directed
  two-cycle.  This is not a tournament because the forward handoff and
  backward blocker edge are different observables and must both remain.
* **Repeated tooth occurrences.**  These are the faithful vertices for
  owner reuse.  Their chronological graph has return circuits with numerical
  holonomy `n_r/n_0`; quotienting occurrences to runners preserves connected
  rank and `m_uv` but destroys the address return.
* **Centered spoke obligations.**  Their functional digraph forces a cycle
  and retains the positive `P/N` product only after `(P,N,Q,r)` are attached.
  Completing unspecified pairs by speed order is a gauge and adds no proof.
* **Gaps or private components.**  They preserve literal cover truth and
  off-grid continuation, but need owner/address labels to know which wall can
  be crossed.
* **LCM clocks.**  They are the correct scalar shadow of located handoffs,
  but by themselves forget chronology and which oriented germ generated the
  overlap.
* **Fano lines, residues, and Fourier modes.**  They organize seven-way
  incidence or global mass, but THM-1247 shows that they do not reconstruct
  the off-grid owner word.

The smallest faithful carrier exposed by this probe is

```text
(full irredundant tooth-occurrence word;
 five-edge located LCM tree and multiplicities;
 one protected private stalk with signed endpoint germs;
 repeated-owner address returns;
 centered blocker cycle on shifted clocks;
 affine drift R and full exact kernel sheet;
 current and adjacent gap addresses).                (48)
```

The tournament is useful telemetry on two projections of (48), but the
proof-bearing object is their typed fibre product.

## 12. Highest-leverage next move

The new mathematics changes the priority.  There is no need to spend another
session buying generic pair-positioning errors: the six-cover already
contains a rank-five tree of actual LCM-quantized overlaps.  Nor is another
bare Fano, `chi_7`, blocker-degree, or selected-gap survival statement likely
to help.

The highest-leverage task is to prove or refute Candidate A on exact
arrangement cells, using (38)--(45) as the arithmetic ledger.  A successful
closed germ lift turns THM-1240's positive product holonomy into a
well-founded address descent because the intermediate owner addresses
cancel.  A failed lift must be charged as a new handoff occurrence and hence
as additional LCM mass.  The private-cell thresholds (17)--(19) and the exact
toothpick branch (23)--(24) ensure that high-ratio owners cannot hide in one
unoriented component.

That is the sharpened frontier: **close the oriented germ lift, or prove that
every failure spends a new located LCM seam.**
