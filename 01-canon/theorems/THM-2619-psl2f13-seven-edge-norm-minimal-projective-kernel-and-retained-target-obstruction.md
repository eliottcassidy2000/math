---
id: THM-2619
title: "PSL2(F13) seven-edge norm, minimal projective kernel, and retained-target obstruction"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  In PSL2(F13), U=[[1,1],[0,1]] and g_t=[[0,1],[-1,t]] for
  t=3,5,6 generate the full group and g_t has projective order seven.
  The forward/reverse products of the seven conjugates g_t^k U^a g_t^-k
  are central exactly when t-a/t+a is in +/-{3,5,6}; the six resulting
  five-element atlases cover every a in F13*.  This is an ordered exact
  nonabelian norm, but not a target transition.  The normalizer of the
  target deck <U> has order 78 and no element of order seven, so the norm
  necessarily cycles seven distinct C13 decks.  Its minimal faithful
  nonnegative permutation realization is the 14-state projective line.
  There the successful norm is identity and the THM-2602 twisted return
  against U^(7a) consists only of infinity; its F13 target part is empty.
  Restoring lawful adjacent 13-state target charts via the natural g_t
  intertwiner restores the original holonomy U^(7a).  All equivariant
  intertwiners form a C13 torsor, and cancellation occurs exactly after an
  externally supplied connection sum c_0+...+c_6=-7a, with 13^6 choices.
  The rank-one Bruhat cell gives an exact 13 by 13 independent-action
  square, but after the future coordinate is given the lawful opposite sign,
  its order-seven trace wall has Fourier support only on lambda=nu and hence
  only target difference zero.
  No nontrivial 2x2 real nonnegative projective representation can carry
  the order-seven clock or order-thirteen target deck.  Thus the construction
  is an exact model of the missing principal-C13 bibundle sidecar, not a
  physical LRC edge, row exclusion, or proof of LRC(14).
source: psl2f13-seven-edge-transition-audit-2026-07-28
depends_on:
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
related:
  - THM-2605-inverse-root-dipole-connection-and-mixed-square-invoice
  - THM-2613-canonical-root-diagonal-opposite-shift-section
  - THM-2615-physical-diagonal-toric-kernel-and-dipole-radon-invoice
  - THM-2607-constant-six-rail-boundary-holonomy-invoice
  - THM-2608-alternative-rail-clock-collapse-and-missing-transition-index
  - THM-2609-external-target-section-itinerary-saturation-and-root-state-no-go
  - THM-2610-chronological-paired-slice-marked-triangle-graft-and-action-axis-boundary
script: 04-computation/lrc14_psl2f13_nonabelian_norm_thm2619.py
output: 05-knowledge/results/lrc14_psl2f13_nonabelian_norm_thm2619.out
script_sha256: c69a6bc5cdd376029fc70b8efaaee6e7f588ac779a29fab529bad4d2cce38214
output_sha256: 0d4c8dc038412e368c0991cacf44636f378a64d7e5bccc6e57f4c242c87d2aaa
hash_basis: LF-normalized bytes
---

# THM-2619 -- the nonabelian norm closes only after the target deck moves

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2602 isolates the live LRC transition test.  Seven ordered nonnegative
kernels must contain a positive target-residue path

```text
q -> q-7a,                                                (1)
```

not merely have nonzero vertex marginals or trivial monodromy after a state
space has been enlarged.  The finite group `PSL_2(F_13)` contains a striking
seven-factor identity which at first appears to meet that demand.  The exact
audit below identifies both its genuine content and its precise type failure.

The gain is real: order matters, the construction has a faithful
nonnegative permutation realization, and its trace atlas covers all twelve
nonzero roots.  The failure is equally sharp: the order-seven motion cannot
normalize one retained target `C_13` deck.  Once the missing adjacent target
identification is supplied, the nonabelian cancellation disappears and the
original `7a` holonomy returns.

## 1. The exact nonabelian norm

Work in `SL_2(F_13)` modulo its centre `{+I,-I}`.  Put

```text
U^a = [[1,a],[0,1]],

g_t = [[0,1],[-1,t]],                 t in {3,5,6}.        (2)
```

For all three traces,

```text
g_t^7=-I,                                                 (3)
```

so `g_t` has projective order seven.  Each pair `U,g_t`
generates the full group of order

```text
|PSL_2(F_13)|=1092.                                      (4)
```

Define the seven conjugate parabolics

```text
A_k(a)=g_t^k U^a g_t^(-k),             k=0,...,6.         (5)
```

Their two ordered norms telescope before any entrywise expansion:

```text
N_t^+(a)=A_0 A_1 ... A_6
          =(U^a g_t)^7 g_t^(-7),

N_t^-(a)=A_6 A_5 ... A_0
          =g_t^7 (g_t^(-1)U^a)^7.                        (6)
```

The moving traces are

```text
tr(U^a g_t)=t-a,
tr(g_t^(-1)U^a)=t+a.                                    (7)
```

An element of `SL_2(F_13)` has noncentral projective order seven exactly
when its trace belongs to

```text
{+/-3,+/-5,+/-6}={3,5,6,7,8,10}.                         (8)
```

The moving matrices in (7) are never central.  Equations (3), (6), and (8)
therefore give the equivalences

```text
N_t^+(a) is central  iff  t-a in +/-{3,5,6},
N_t^-(a) is central  iff  t+a in +/-{3,5,6}.              (9)
```

For `a!=0`, the exact success atlas is

| `t` | forward `N_t^+` | reverse `N_t^-` |
|---:|:---|:---|
| `3` | `6,8,9,10,11` | `2,3,4,5,7` |
| `5` | `2,8,10,11,12` | `1,2,3,5,11` |
| `6` | `1,3,9,11,12` | `1,2,4,10,12` |

The union is every element of `F_13^*`.  Reversing the same seven factors
changes `t-a` to `t+a`; this is genuine ordered noncommutativity, not a
marginal support effect.

## 2. Why the norm must move the target deck

Let

```text
P=<U> isomorphic to C_13.                                 (10)
```

This is the abstract target-translation deck: in the standard projective
chart it acts by `q -> q+1`.  Its normalizer in `PSL_2(F_13)` is the upper
triangular Borel subgroup,

```text
|N(P)|=13*6=78.                                          (11)
```

In particular `N(P)` has no element of order seven.  Equivalently, an
order-seven clock cannot act nontrivially on a retained `C_13` deck: such an
action would give an order-seven element of

```text
Aut(C_13) isomorphic to C_12,                             (12)
```

which is impossible.

Consequently the seven subgroups

```text
P_k=g_t^k P g_t^(-k)                                     (13)
```

are distinct.  Any two intersect only in the identity.  The cancellation
in (9) is possible precisely because the factors do **not** remain in one
target translation group.  This is the group-theoretic version of the
root-deck loss in THM-2607 and THM-2610.

It also gives an exact boundary for the newly promoted THM-2605 affine
path.  That theorem's invariant phase `v=kr+q` and update
`q -> q-alpha` use one fixed affine `C_13` translation coordinate.  An
order-seven `g_t` cannot normalize that deck by (11); adjoining the PSL clock
would replace the coherent THM-2605 phase by seven different Sylow decks.
Thus the nonabelian norm cannot type THM-2605's physical path as the missing
chart/next-target path while preserving its proved affine invariant.

There is no nonabelian group extension on one `13 x 7` clock/root fibre that
repairs this.  Since `7` does not divide `|Aut(C_13)|=12`, every semidirect
product `C_13 semidirect C_7` is the direct product.  On a retained principal
target deck the seven conjugates are therefore all the same translation and
their product is `U^(7a)`, not the identity.

## 3. The faithful nonnegative realization adds a cemetery state

Let `PSL_2(F_13)` act on

```text
P^1(F_13)=F_13 union {infinity}                           (14)
```

by fractional linear transformations.  This gives faithful `14 x 14`
permutation matrices, hence deterministic nonnegative kernels.  In this
representation:

```text
U^a: one 13-cycle on F_13, fixing infinity;
g_t: two disjoint 7-cycles;
P_k: fixes f_k=g_t^k infinity and cycles the other 13 states. (15)
```

Thus (9) becomes a completely lawful ordered permutation identity on the
fourteen-state set.  It does not become a target transition.  For every one
of the thirty successful `(t,orientation,a)` triples, projective centrality
gives

```text
rho(N_t^+/- (a))=I_14.                                   (16)
```

THM-2602 asks for the twisted return after the nonzero base holonomy
`U^(7a)`.  From (16),

```text
Fix(rho(N_t^+/- (a)) rho(U^(7a)))={infinity},

Fix(...) intersect F_13=empty.                           (17)
```

So the full permutation trace is positive for exactly the wrong reason:
its sole return is the extra cemetery point.  The target-residue twisted
trace is zero.

Restricting each factor to the affine thirteen-state target alphabet does
not help.  Every successful seven-factor path has three, four, or five
starting residues that visit `infinity`; respectively ten, nine, or eight
avoid it.  The restricted product is the partial identity on those avoiding
residues.  Since `7a!=0`, its product with `U^(7a)` still has empty diagonal.
Thus the cheapest decisive test already fails on every success in (9).

### 3.1 Fourteen states are not an accidental choice

The natural action is the minimal faithful permutation action that can still
carry `U` nontrivially.  Here is a short proof specialized to the only
possible competing size.

Suppose `PSL_2(F_13)` acted on thirteen states with `U` as a 13-cycle.  The
action would be transitive and, by simplicity, faithful.  A point stabilizer
`H` would have order `1092/13=84` and would be self-normalizing.  Its Sylow
seven subgroup is unique, since its number divides `12` and is `1 mod 7`.
Choose a nonidentity seven-element `x in H`.  Its centralizer in
`PSL_2(F_13)` is the order-seven nonsplit torus, so its conjugacy class `C`
has size `1092/7=156`.

Put `m=|C intersect H|`.  Then `1<=m<=6`.  Counting incidences between `C`
and the thirteen conjugates of `H` gives

```text
13m=156f
```

for an integer `f`, hence `m=12f`, a contradiction.  Therefore no faithful
thirteen-state permutation kernel exists.  The fourteenth point in (14) is
forced, and (17) shows exactly how it corrupts the LRC return test.

## 4. Lawful adjacent target charts restore the old holonomy

The projective action contains a more faithful local description of the
type failure.  Set

```text
V_k=g_t^k U g_t^(-k),
f_k=g_t^k infinity,
Q_k=P^1(F_13)\{f_k}.                                     (18)
```

Then `Q_k` is a regular `C_13` torsor for `V_k`, and `A_k(a)=V_k^a`
preserves `Q_k`.  But an LRC edge at clock `k` must end in the next target
alphabet `Q_(k+1)`, not back in `Q_k`.

The natural adjacent-chart map is `g_t:Q_k -> Q_(k+1)`.  More generally,
all deck-equivariant bijections are

```text
T_k(c)=V_(k+1)^c g_t=g_t V_k^c,          c in F_13.       (19)
```

They form the principal `C_13` bibundle classified by proved THM-2611.  Use
the local coordinate

```text
phi_k:Q_k -> F_13,             phi_k(x)=g_t^(-k)x.        (20)
```

The lawfully typed edge

```text
E_k(c_k)=T_k(c_k)V_k^a:Q_k -> Q_(k+1)                    (21)
```

then becomes simply

```text
phi_(k+1) E_k(c_k) phi_k^(-1)=U^(a+c_k).                 (22)
```

Since `g_t^7` is central, the sevenfold target holonomy is

```text
E_6(c_6)...E_0(c_0)=U^(7a+c_0+...+c_6)    projectively.  (23)
```

Two consequences are exact.

First, the natural connection `c_k=0` restores

```text
U^(7a),                                                     (24)
```

the original THM-2602 obstruction.  This holds for every `t,a`, whether or
not the untyped raw norm in (9) happens to be central.

Second, a lawful equivariant cancellation exists exactly when

```text
c_0+...+c_6=-7a mod 13.                                  (25)
```

There are exactly `13^6=4,826,809` such connection choices.  Their
correction comes entirely from the supplied intertwiner torsor (19); the
trace condition (9) neither constructs nor selects one.  Taking
`c_0=-7a,c_1=...=c_6=0`, for example, merely inserts the full missing target
correction on one edge.

This is the precise connection to THM-2608.  The raw projective-line product
is a dense, composable object only after the target alphabets have been
coalesced into a fourteen-state space.  Restoring the next-target index and
its equivariance either gives the natural uncancelled holonomy (24), or asks
for exactly the external sidecar (25) that LRC still lacks.

## 5. The Bruhat cell is the right square and the wrong target mode

There is one positive structural connection which survives the no-go.
Let

```text
W=[[0,-1],[1,0]],

B(r,s_B)=U^r W U^s_B
        =[[r,r s_B-1],[1,s_B]].                         (25a)
```

The `169` matrices (25a) are distinct in `PSL_2(F_13)`, and

```text
U^x B(r,s_B) U^y=B(r+x,s_B+y).                          (25b)
```

Thus the rank-one Bruhat cell is an exact algebraic `13 x 13` carrier with
independent left and right `C_13` actions.  It is the group-theoretic
counterpart of the independent dipole square typed by candidate THM-2615.

The sign of the future axis is load-bearing.  THM-2615 uses an opposite
future shift, so set

```text
s=-s_B.
```

Then left multiplication translates `r` positively, right multiplication
translates `s` negatively, as required.  In these typed coordinates the
order-seven wall is

```text
R(r,s)=1_(B(r,-s) has projective order seven)
      =1_(r-s in +/-{3,5,6}).                            (25c)
```

The two moving matrices are precisely slices of this square:

```text
U^a g_t       = B(a,-t)       projectively,  so (r,s)=(a,t),
g_t^(-1)U^a  = B(t,a),                         (r,s)=(t,-a).
                                                               (25d)
```

Their trace tests `a-t` and `t+a`, respectively.  Hence the six sets in
Section 1 are the intersections of two families of coordinate lines with
one `78`-point Bruhat trace wall.

This exact square still carries no primitive target residue.  With
THM-2615's transform convention, put

```text
Rhat(lambda,nu)
 =1/169 sum_(r,s) R(r,s) zeta^(-lambda r+nu s).           (25e)
```

Writing `u=r-s`, the sum over `r` vanishes unless `lambda=nu`.  On that
diagonal,

```text
Rhat(c,c)=1/13 sum_(u in +/-{3,5,6}) zeta^(-cu) !=0.      (25f)
```

The final inequality holds for every `c`: a nonempty proper `0/1`
polynomial of degree below thirteen cannot vanish at a primitive thirteenth
root.  Therefore

```text
support(Rhat)={(c,c):c in F_13},
{lambda-nu:(lambda,nu) in support(Rhat)}={0}.             (25g)
```

THM-2615's target Radon coordinate is `q=lambda-nu`.  The full PSL
order-seven wall is consequently concentrated in the target-zero channel.
The Bruhat cell supplies the correct two-axis *ambient carrier*, but the
trace norm populates exactly the wrong Fourier line.  A useful future bridge
would need a physical response on this square which is not merely a function
of `r-s`.

## 6. A `2 x 2` nonnegative/projective kernel cannot evade the boundary

The matrices in (2) are matrices over `F_13`; reading their residues as
nonnegative real entries does not preserve multiplication or projective
centrality.  There is also an intrinsic obstruction to any honest real
`2 x 2` replacement.

**Lemma.**  If an invertible real entrywise-nonnegative `2 x 2` matrix has
odd finite projective order, then it is projectively the identity.

Indeed, for `B=[[r,s],[u,v]]>=0`, the characteristic discriminant is

```text
(r-v)^2+4su>=0,                                           (26)
```

so both eigenvalues are real.  If `B^m=lambda I` for odd `m`, their ratio is
a real `m`-th root of one and hence equals one.  A nontrivial Jordan block
cannot have a scalar positive power, so `B` itself is scalar.  In particular
neither the order-seven `g_t` nor the order-thirteen `U` has a nontrivial
two-dimensional nonnegative projective realization.

There is an even stronger factorwise form.  If seven invertible
nonnegative `2 x 2` matrices have scalar product, every cyclic product is
the same scalar and each factor has a nonnegative inverse.  A matrix and its
inverse are both nonnegative only when the matrix is monomial.  On the
positive projective interval, a two-dimensional monomial map has finite
orbits of length at most two.  It cannot permute thirteen target residues or
carry the seven-cycle clock.

Thus the literal `2 x 2` norm and a nonnegative transition kernel live in
different categories.  The first faithful nonnegative functor is the
fourteen-dimensional permutation action, and Section 3 gives its exact
hostile consequence.

## 7. Binary-relation and tournament ledger

The faithful vertices here are the fourteen projective states, with seven
clock-coloured deterministic arrows

```text
(k,x) -> (k+1,A_k(a)x).                                  (27)
```

The target predicate is the twisted return (1).  Forgetting the colour/order
destroys the norm; deleting `infinity` destroys totality; forgetting which
`P_k` acts destroys the target deck.

A tournament is not intrinsic.  The projective action is 2-transitive, so
every invariant binary relation on distinct projective states is either all
ordered pairs or none; it cannot orient one member of every pair.  The seven
fixed points `f_k` carry a directed clock cycle, not a tournament, and any
orientation of its nonadjacent pairs is external.  The faithful relational
object is therefore the coloured permutation automaton together with the
principal-`C_13` fibre sidecar (19).

The exact connection contract is:

```text
source:       ordered PSL2(F13) conjugation norm;
map:          action on P^1(F13) by permutation kernels;
preserved:    order, multiplication, positivity, norm centrality;
added/lost:   one cemetery state is added and one fixed target deck is lost;
sidecar:      adjacent equivariant C13 intertwiners plus physical ancestry,
              clock, target gauge, owner, and repair semantics;
decisive test: target-restricted diag(K U^(7a)), which is empty.             (28)
```

This also separates the latest canon cleanly.  THM-2605 supplies a coherent
positive physical affine path but no chart/next-target identification.
THM-2607 supplies an abstract rail boundary but no physical deck
intertwiner.  THM-2608 shows clock marginalization can create a false
transition prism.  THM-2609 supplies external section itineraries but no
event-state equality.  THM-2610 joins a later action character while positive
time erases the old deck.  Proved THM-2611 classifies the missing bibundle
abstractly.  The present model shows that PSL noncommutativity does not build
that physical bibundle: it changes the `C_13` deck instead.

## 8. Exact companion and scope

Run

```text
python 04-computation/lrc14_psl2f13_nonabelian_norm_thm2619.py
python -O 04-computation/lrc14_psl2f13_nonabelian_norm_thm2619.py
```

The companion uses only exact arithmetic.  It:

1. generates all `1,092` projective matrices from each pair `U,g_t`;
2. verifies the complete projective order and trace census;
3. reconstructs all six central-norm success sets;
4. computes `N(P)`, the seven distinct conjugate decks, and their fixed points;
5. realizes every factor as a fourteen-state permutation kernel;
6. checks all thirty successful norm paths, affine leakage, and the exact
   cemetery-only versus target-empty twisted return;
7. checks every one of the `13*7*3` equivariant maps (19), the natural
   holonomy (24), and representative full connection laws (23); and
8. verifies that every total connection holonomy has exactly `13^6` lifts;
   and
9. reconstructs the `169` Bruhat coordinates, the `78`-point order-seven
   wall, and its exact target-zero Fourier support (25g).

Normal and optimized executions byte-match the stored transcript after LF
normalization.  No floating point, random sampling, or inferred LRC carrier
is used.

What is **not** supplied is a physical LRC realization of `g_t`, a
same-ancestry adjacent target bibundle, a positive Boolean mass for (21), a
terminal word/owner/repair transport, or a reason to select one of the
`13^6` corrections in (25).  No scalar row is excluded and LRC(14) remains
open.

**QED pending independent hostile audit.**
