---
id: THM-1248
title: THE SLOWEST-ROOTED CENTERED BLOCKER LASSO HAS A FINITE ORIENTED ADDRESS WORD AND LIFTS TO A BOUNDED BOUNDARY-EVENT WALK
status: PROVED (exact lower-gap cocycle and translation gauge; least-positive determinant remainder and gcd sheet; uniform 1,174-symbol relative-address bank; binary chronological digits on every target-at-most-source-clock edge; exact oriented-germ and all-candidate sampled-owner reconstruction; additive and contracting affine lasso transport; proper tail length at most four and slowest-rooted multiplier below 13/19; every blocker tooth has two interior located gcd seams, with fourth-support/protected placement on source-safe/ascent walls; two-wall owner pigeonhole reducing the incidence-free/reuse-free lasso to the slowest two-cycle; canonical tooth-instance expiry tournament and piecewise-monotone tail-to-cycle lift; protected-tree tariff corollary; seven exact guardrails including realizable two-wall residuals and identical-seam retracing; dependency-free exact referee; sorry-free Lean arithmetic core)
source: codex-2026-07-19-S82 continuation with address-compression, cycle-model, and package-audit agents
depends_on: [THM-1233, THM-1237, THM-1240, THM-1244, THM-1250]
related: [THM-1198, THM-1226, THM-1242, THM-1243, THM-1246, THM-1247, THM-1252, THM-1253, THM-1254, THM-1256, HYP-7870]
script: 04-computation/lrc14_centered_blocker_address_compression_thm1248.py
output: 05-knowledge/results/lrc14_centered_blocker_address_compression_thm1248.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCCenteredBlockerAddressCompression.lean
script_sha256: 0595ced97eaf4469752c9334b93119f357f7645875fa4d8944012d6ef18c08eb
output_sha256: 97fbe205630b5cae842ebb87a25ac841cad45252ce54ec83ee0c8b4e8b7ff88f
formalization_sha256: 4b135cb3eda4e26db0db07403becc3680cc64d2c4b97e7f8a8d4d8c639486191
---

# THM-1248 — centered blocker address compression and wall export

## 1. The lower-gap blocker cocycle

Let

```text
G=G_k(c)=[(14k+1)/(14c),(14k+13)/(14c)],   0<=k<c,   (1)
```

be a complete closed safe gap of the positive integer speed `c`.  Suppose
the six strict danger combs of

```text
c<d1<d2<d3<d4<d5<d6                                  (2)
```

cover `G`.  Put

```text
t0=(k+1/2)/c,             Q_i=c+d_i.                  (3)
```

For every `i`, choose an integer `P_i` nearest `Q_i t0`, set

```text
t_i=P_i/Q_i,
epsilon_i=P_i-Q_i t0 in [-1/2,1/2],                  (4)
```

and choose one THM-1240 blocker `j=b(i)`.  Thus `j!=i` and `d_j` is
strictly dangerous at `t_i`.  There is a unique integer `N_i` such that

```text
beta_i=d_j t_i-N_i,               |beta_i|<1/14.      (5)
```

The raw tooth address `N_i` is not the correct invariant.  Define instead

```text
X_i=cP_i-kQ_i,
r_i=P_i d_j-N_iQ_i=Q_i beta_i,
M_i=k+N_i,
ell_i=P_iQ_j-M_iQ_i,
delta_i=P_j-M_i.                                      (6)
```

Direct expansion using `Q_j=c+d_j` gives the exact cocycle

```text
ell_i=X_i+r_i.                                        (7)
```

THM-1240 says

```text
Q_i/4<X_i<3Q_i/4.                                    (8)
```

Combining (5), (7), and (8) produces the sharp positive central band

```text
5Q_i/28<ell_i<23Q_i/28.                              (9)
```

In particular `0<ell_i<Q_i`; no modular choice is hidden in (9).

## 2. Gauge removal, determinant remainder, and gcd sheet

An integral translation of time by `n` changes the addresses by

```text
P_i ->P_i+nQ_i,       P_j->P_j+nQ_j,
k   ->k+nc,           N_i->N_i+nd_j,
M_i ->M_i+nQ_j.                                      (10)
```

The four quantities `X_i,r_i,ell_i,delta_i` are unchanged.  Thus the
absolute lift `N_i` is gauge, whereas `delta_i` is a genuine relative tooth
address.

Let

```text
D_ij=P_iQ_j-P_jQ_i.                                  (11)
```

Then (6) gives

```text
D_ij=ell_i-delta_iQ_i.                               (12)
```

By (9), `ell_i` is exactly the least positive residue of `D_ij` modulo
`Q_i`, and

```text
M_i=floor(P_iQ_j/Q_i),
gcd(Q_i,Q_j) divides ell_i.                           (13)
```

The second statement follows before any reduction: every common divisor of
the two clocks divides both terms in `P_iQ_j-M_iQ_i`.  Consequently the
faithful sampled edge is not just a directed pair of runners but

```text
(Q_i,Q_j; least-positive ell_i; relative delta_i;
 gcd(Q_i,Q_j)-sheet).                                 (14)
```

It preserves blocker danger, phase chronology, and exact torsion while
discarding only the integral tooth lift.

## 3. The finite relative-address bank

Substitute (3)--(6) and `ct0=k+1/2`.  All absolute addresses cancel:

```text
beta_i=delta_i-1/2+(d_j/Q_i)epsilon_i-epsilon_j,      (15)
epsilon_j=(d_j/Q_i)epsilon_i+delta_i-1/2-beta_i.      (16)
```

Writing `a_i=d_j/Q_i`, the rounding and danger bounds imply

```text
-a_i/2-1/14<delta_i<1+a_i/2+1/14.                    (17)
```

THM-1233 gives `d_j<2345c`, while `Q_i>2c`.  Hence
`a_i<2345/2`; since `delta_i` is integral, (17) yields the uniform bank

```text
-586<=delta_i<=587.                                  (18)
```

There are therefore only `1,174` possible relative address symbols on an
arbitrary blocker edge.  More sharply, if

```text
d_j<=Q_i=c+d_i,                                      (19)
```

then `a_i<=1`, so (17) lies strictly inside `(-1,2)` and

```text
delta_i in {0,1}.                                    (20)
```

Every speed-descent edge has this binary property, but (19) also includes
many ascent edges.

Put `theta_i=ell_i/Q_i`.  Equations (6) and (12) give

```text
Q_j(t_i-t_j)=theta_i-delta_i.                         (21)
```

When (20) holds, the digit records which chronological side of `t_j` the
source phase occupies, and (9) becomes the quantitative separation

```text
5/(28Q_j)<|t_i-t_j|<23/(28Q_j).                      (22)
```

Thus a binary edge is an adjacent-centered-tooth transport, not a bare
orientation.

### 3.1 The oriented tooth germ is already in the quotient

The left and right distances from `t_i` to the two walls of the blocking
`d_j` tooth are

```text
lambda_i^-=1/(14d_j)+r_i/(d_jQ_i),
lambda_i^+=1/(14d_j)-r_i/(d_jQ_i).                  (22a)
```

Their signed imbalance is

```text
lambda_i^--lambda_i^+
 =2r_i/(d_jQ_i)=2(ell_i-X_i)/(d_jQ_i).              (22b)
```

The needed vertex coordinate is itself in the common orbit (40):

```text
X_i=Q_i/2+c epsilon_i=(Q_i+E_i)/2.                  (22b')
```

Thus the off-grid orientation identified by the q=15 guardrail is not an
additional raw-address coordinate: it is reconstructed exactly from the
full vertex-edge quotient `(E_i,Q_i;ell_i,delta_i)`.  The smaller edge tuple
(14) alone intentionally omits the vertex coordinate `X_i`.

Moreover `theta_i` lies strictly between zero and one and `delta_i` is an
integer.  Equation (21) therefore gives the exhaustive side rule

```text
delta_i<=0  => t_j<t_i and the path selects the left wall,
delta_i>=1  => t_i<t_j and the path selects the right wall.               (22c)
```

For a binary edge, `0` and `1` are literally its two oriented germ symbols.
For a large digit, its sign class still selects the wall and its magnitude
retains the inter-spoke address displacement.

## 4. Exact cycle invoices

Let

```text
i_1->i_2->...->i_m->i_1,              2<=m<=6,       (23)
```

be a directed cycle of the selected blocker map.  Summing (21) around the
cycle gives the additive invoice

```text
sum_s delta_s/Q_(s+1)=sum_s theta_s/Q_(s+1).          (24)
```

After division by `sum_s 1/Q_(s+1)`, the right side is a weighted average
with every `theta_s in (5/28,23/28)`.  Therefore every cycle contains both

```text
some delta_s<=0,                    some delta_r>=1.  (25)
```

The positive remainders also have a multiplicative form.  In the canonical
gap gauge, `P_i>0`.  Moreover `N_i>0`: since `t_i>1/(14c)` and `d_j>c`,
one has `d_jt_i>1/14`, so (5) excludes `N_i=0`; hence `M_i=k+N_i>0`.
Each identity

```text
P_iQ_j=M_iQ_i+ell_i>M_iQ_i                           (26)
```

can now be multiplied around (23), giving the positive address holonomy

```text
prod_s P_s>prod_s M_s.                               (27)
```

This product statement is canonical-gauge data; (12), (15), and (24) are
the gauge-invariant content.

## 5. The blocker cycle is a uniformly contracting affine local system

Equation (16) transports the centered rounding errors along every edge:

```text
epsilon_(s+1)=a_s epsilon_s+b_s,
a_s=d_(s+1)/(c+d_s),
b_s=delta_s-1/2-beta_s.                              (28)
```

After one circuit,

```text
epsilon_1=A epsilon_1+B,
A=prod_s a_s=prod_s d_s/(c+d_s),
B=sum_s b_s w_s,
w_s=prod_(r>s) a_r.                                  (29)
```

Every factor is positive and less than one.  Since `m>=2` and THM-1233
gives `d_s/(c+d_s)<2345/2346`,

```text
0<A<(2345/2346)^2,
1-A>1-(2345/2346)^2=4691/5503716.                    (30)
```

Thus every fixed choice of slopes and beta/delta data has a unique affine
fixed point.  For the actual spoke errors, `|epsilon_1|<=1/2`, so

```text
|B|<=(1-A)/2.                                        (31)
```

Eliminating the strict blocker errors from (29)--(31) yields the finite-word
exclusion target

```text
|sum_s (delta_s-1/2)w_s|
 <(1/14)sum_s w_s+(1-A)/2.                           (32)
```

The sampled cycle is therefore neither an unbounded address walk nor a free
combinatorial cycle.  It is a length-at-most-six word over (18), carried by
a uniformly contracting affine error system and the exact gcd residues
(13).

## 6. Every blocker tooth exports two phase-located gcd seams

Every directed cycle (23) has a speed-ascent edge `i->j` with `d_i<d_j`;
otherwise speed would decrease strictly around a cycle.  At `t_i`, runner
`j` lies strictly inside one danger tooth `J` by construction.  At `t_j`,
runner `j` has depth greater than `1/4`, so the segment from `t_i` to `t_j`
leaves `J`.  Let `w` be its first boundary point.

The source remains quantitatively safe there.  Since the width of `J` is
`1/(7d_j)`,

```text
|w-t_i|<1/(7d_j),
||d_iw||>1/4-d_i/(7d_j)>3/28>1/14.                   (33)
```

This wall event is in fact protected.  Let `S_i` be the closed `d_i`-safe
component through `t_i` and put

```text
K_i=G intersect S_i.                                 (34)
```

If `Delta_i=||d_it_i||`, the distance from `t_i` to either endpoint of
`S_i` is at least

```text
rho_i=(Delta_i-1/14)/d_i>5/(28d_i).                  (35)
```

The centered-spoke calculation used in THM-1244 shows that the distance to
either endpoint of `G` is still larger than `rho_i`.  Every point of the
whole target tooth `J` is less than its width from `t_i`, and

```text
1/(7d_j)<1/(7d_i)<5/(28d_i)<rho_i.                   (36)
```

Hence the closure of `J` lies in `int(K_i)`.

At `w`, the carrier `c`, source `d_i`, and target `d_j` are all safe (the
target is equality-safe).  The assumed six-cover therefore forces a third
fast label

```text
h notin {i,j}                                         (37)
```

to be strictly dangerous at `w`; equivalently, `h` is a fourth support when
the carrier is counted.  Its open tooth contains `w`, so it overlaps `J` on
a nonempty open interval `Omega_ijh` adjacent to `w`.  By (36), this overlap
lies wholly in `int(K_i)`.  Its endpoints are rational tooth endpoints.  A
positive endpoint difference has an integer numerator divisible by
`gcd(d_j,d_h)`, and therefore

```text
mu(Omega_ijh)>=gcd(d_j,d_h)/(14d_jd_h).               (38)
```

The source-safe ascent is the strongest version, but seam existence does not
need source safety.  Section 8.2 proves that the closure of **every** selected
blocker tooth lies in `int(G)`.  At either of its two walls, the carrier and
target are safe, so the strict six-cover chooses some fast owner `h!=j`
(possibly `h=i`).  Its open tooth contains the wall and overlaps the target
tooth on a positive inward interval.  Both endpoints are fast-tooth
endpoints, hence the same divisibility argument gives

```text
mu(Omega_jh)>=gcd(d_j,d_h)/(14d_jd_h).               (38a)
```

Thus every actual blocker edge exports two opposite, located gcd seams.
On a source-safe wall the owner is also distinct from the source; on an
ascent, (36) places the seam in the protected source component.  This is the
coverage-sensitive information absent from the six sampled phases: the
functional graph has become a two-sided wall-event hypergraph with genuine
metric credits.

## 7. The slowest orbit anchors one positioned Hunter credit

Apply the same argument to the selected outgoing edge `d1->d_j`.  It is
automatically a speed ascent.  Here `K_1` in (34) is exactly THM-1244's
protected slowest-spoke component `K`, so (38) is a located overlap between
two of the five active blocker labels inside `K`.

THM-1244 supplies a rank-two forest of two distinct projected handoff edges
inside `K`.  At least one of those edges differs from the wall edge
`{d_j,d_h}`.  Adjoin that edge to `{d_j,d_h}`.  Two distinct edges form a
forest, so THM-1244's forest can be reselected with one prescribed credit:

```text
one of its two positioned Hunter edges is the wall seam (38).             (39)
```

The quantitative scalar and two-unlocated-edge debts of THM-1244 remain
valid, since each selected overlap still contributes at least
`gcd(d2,...,d6)/(14d6^2)`.

Starting the functional orbit at `d1` connects this anchored seam to its
eventual cycle.  The initial edge need not itself lie on that cycle, and the
cycle seam from Section 6 need not lie in the slowest component `K`.
THM-1250 sharpens (39): the graph of actual overlaps on the complete gap is
connected, so any prescribed forest of wall-seam label pairs extends to a
five-edge located spanning tree.  Repeated chronological handoffs are charged
by its multiplicity-averaged Hunter debt.

Section 8 closes the previously existential *component-path* gap: the whole
slowest-rooted lasso lifts canonically to actual danger-tooth instances.  What
remains is narrower.  Directions can reverse at centered spokes, an owner
label can recur in different teeth, and the shifted `c+d_i` holonomy does not
automatically orient the `d_i` handoff tree.

## 8. The whole slowest-rooted lasso and its boundary-event lift

### 8.1 A four-edge proper tail with exact affine transport

Write the distinct part of the functional orbit as

```text
v_0=d1, v_(s+1)=b(v_s),
v_0,...,v_(L+m-1) distinct,          v_(L+m)=v_L,     (39a)
```

where the eventual cycle has length `m>=2`.  Since there are six labels,

```text
L+m<=6,                              0<=L<=4.         (39b)
```

Thus the proper tail has at most four edges.  Put

```text
E_s=2c epsilon_s,        theta_s=ell_s/Q_s.
```

Clearing (16) with (7) and `Q_(s+1)=c+d_(s+1)` gives the exact residue step

```text
Q_s E_(s+1)-Q_(s+1)E_s
     =2c(Q_s delta_s-ell_s).                         (39c)
```

On the proper tail it telescopes to

```text
E_L/Q_L-E_0/Q_0
 =2c sum_(s<L)(delta_s-theta_s)/Q_(s+1).             (39d)
```

There is a second normalization which exposes genuine contraction.  Set

```text
z_s=epsilon_s/d_s,       rho_s=d_s/(c+d_s).
```

Then every tail edge obeys

```text
z_(s+1)=rho_s z_s
 +(delta_s-1/2-beta_s)/d_(s+1).                      (39e)
```

If the tail is nonempty, its first source is `d1`.  THM-1198 gives
`d1/c<13/6`; every later `rho_s<1`.  Hence its homogeneous multiplier is

```text
R_tail=prod_(s<L)rho_s<13/19.                        (39f)
```

With `W_s=prod_(s<r<L)rho_r`, iteration of (39e) yields the finite invoice

```text
|sum_(s<L) W_s(delta_s-1/2)/d_(s+1)|
 <1/(2d_L)+R_tail/(2d1)
   +(1/14)sum_(s<L)W_s/d_(s+1).                      (39g)
```

One can also choose the blocker map to minimize the number of distinct
vertices in this `d1`-orbit.  Then no proper tail edge `i->j` is reciprocated
by danger of `d_i` at `t_j`; otherwise resetting `b(j)=i` closes a two-cycle
and strictly shortens the orbit.  In exact quotient coordinates,

```text
||d_i t_j||
 =||(D_ij+X_j)/Q_j||
 =||(ell_i-delta_i Q_i+X_j)/Q_j|| >=1/14.            (39h)
```

So a minimal lasso is a word of length at most six with at most four explicit
one-way tests, not an unbounded functional tail.

### 8.2 Selected-wall source, target, and owner germs are quotient data

Let

```text
sigma_i=+1 if delta_i>=1,       sigma_i=-1 if delta_i<=0.
```

The selected target-tooth wall and its gauge-free displacement are

```text
w_i=(14N_i+sigma_i)/(14d_j),
w_i-t0=E_i/(2cQ_i)-r_i/(d_jQ_i)+sigma_i/(14d_j).     (39i)
```

The source residue at that wall is exactly

```text
||d_iw_i||=||A_i/(14d_jQ_i)||,
A_i=7d_j(Q_i-E_i)+d_i(sigma_iQ_i-14r_i).             (39j)
```

Thus `||d_iw_i||>=1/14` is a finite rational predicate of the quotient.  At
every selected wall, the carrier and equality-safe target already force an
owner `h!=j` and hence a located target--owner seam in `G`.  If (39j) says the
source is safe, then `h!=i` as well and the event is a fourth support.  The
ascent argument of Section 6 is the strongest special case, placing the
whole seam inside the protected source component.  Replacing `sigma_i` by
`-sigma_i` reconstructs the opposite wall and its owner predicate, so both
wall seams are quotient-visible.

There is no carrier-boundary clipping hidden here.  The centered spoke has

```text
dist(t_i,partial G)
 =[3/7-|E_i|/(2Q_i)]/c>5/(28c)>1/(7d_j),
```

so the closure of the entire blocked `d_j` tooth lies in `int(G)`.  Its
target--owner overlap is therefore a genuine interior gcd-quantized seam.

The opposite target leg is also sharper than the original central band.  Put

```text
C_i=14d_jQ_j sigma_i(t_j-w_i)
   =-7sigma_i E_j+(14|delta_i-1/2|-1)Q_j.            (39k)
```

The integer coefficient of `Q_j` is at least six.  Since `|E_j|<=c`,
`Q_j=c+d_j`, and `d_j>c`,

```text
C_i>=-7c+6Q_j=6d_j-c>5d_j.                          (39l)
```

The source lies strictly inside the target tooth, so the source-to-wall leg
is positive.  Adding it to (39l) gives

```text
|theta_i-delta_i|=Q_j|t_i-t_j|>5/14.                (39m)
```

For binary edges this strengthens (9) to

```text
delta_i=0 => 5/14<theta_i<23/28,
delta_i=1 => 5/28<theta_i<9/14.                      (39n)
```

Now retain the centered state `(Q_h,E_h)` of every candidate owner `h`.  Its
wall residue is also reconstructed without the absolute tooth address:

```text
||h w_i||=||B_(h,i)/L_(h,i)||,
B_(h,i)=7d_jQ_i(c-E_h)+7hE_i d_j-14ch r_i
        +ch sigma_i Q_i,
L_(h,i)=14c d_jQ_i.                                  (39o)
```

Consequently strict sampled ownership is not information missing from the
*full* six-vertex quotient; it is the exact predicate
`dist(B_(h,i),L_(h,i)Z)<L_(h,i)/14`.  If it holds, choose the nearest owner
tooth integer `H_(h,i)` and set

```text
alpha_(h,i)=h w_i-H_(h,i) in (-1/14,1/14),
f_(h,i)=(1/14-sigma_i alpha_(h,i))/h,
b_(h,i)=(1/14+sigma_i alpha_(h,i))/h.                (39p)
```

These are the exact forward expiry and backward overlap margins in the
travel direction; `f+b=1/(7h)`, and the adjacent target--owner overlap has
length

```text
min(1/(7d_j),b_(h,i)).                               (39q)
```

The natural sampled vertex is therefore the signed seam germ
`(j,h,w_i,sigma_i,alpha_(h,i);b_(h,i),f_(h,i))`, not a runner alone.

### 8.3 The expiry tournament canonically transports every wall to the next

Fix one segment from `w_i` to `t_j` and use the monotone coordinate
`u=sigma_i t`.  At a current point, take as vertices the active *tooth
instances*, not their runner labels.  Orient a pair toward the tooth with the
later forward endpoint; break an exact endpoint tie by `(label,tooth index)`.
The pairwise observable is the endpoint difference.  This expiry tournament
is transitive.  With `p` active instances its score multiset is
`(0,1,...,p-1)`, it has no directed cycle, `p` singleton SCCs, and one
Hamiltonian path after the lexicographic tie-break: the sorted endpoint list.

Let `U=sigma_i t_j` and let `e` be the latest active forward endpoint.  If
`e>=U`, stop at `U`; the active tooth contains `t_j`.  If `e<U`, advance to
`e`.  Strict coverage supplies a successor containing `e`, so its left
endpoint is strictly below `e` and its forward endpoint is strictly above
`e`; the two open teeth overlap positively.  The endpoints increase
strictly, so the greedy interval-cover walk terminates.  The THM-1233 ratio
bound used in THM-1250 gives at most `2011` teeth of each fast label meeting
`G`.  One fixed travel direction uses at most `6*2011=12,066` forward
endpoints.  Keeping both oriented endpoints in one direction-independent
bank gives the conservative bound

```text
6*2*2011=24,132                                      (39r)
```

boundary instances.  The lasso has at most six segments, so the conservative
two-orientation event-occurrence cap is `144,792` even when directions
alternate.

At `t_j`, append the prescribed next blocker tooth if the greedy chain ended
in another active tooth; the two open teeth overlap at the spoke.  That
prescribed tooth contains the outgoing segment up to the next selected wall
and reaches the wall in its closure; the strict next-wall owner and overlap
from (39q) cross the boundary.  Iterating stitches the complete `d1` tail and
eventual cycle into a continuous, piecewise-monotone walk through actual
danger-tooth instances and their wall closures.  This proves
undirected/component-level tail-to-cycle transport.  Its turns occur exactly
when the signs `sigma_i` reverse.  A repeated owner label may occur in a
different tooth, so the lift does not yet supply a decreasing gain or a
closed oriented holonomy.

### 8.4 Two seams, owner reuse, and the exact finite residual

If `L=0`, the slowest edge is already a cycle ascent.  If `L>0`, the first
edge `d1->b(d1)` and some cycle ascent are distinct functional edges, and
both export seams.  Write their target--owner label pairs as `e_0,e_C` and
their lengths as `omega_0,omega_C`.  If the pairs are distinct, they form a
two-edge forest.  Hunter on the full gap gives

```text
H>=1/c+(49/6)(omega_0+omega_C)
 >=1/c+(7/12)[1/lcm(e_0)+1/lcm(e_C)].                (39s)
```

Adjoin these interior wall edges to THM-1250's connected located
chronological handoff graph.  Any prescribed forest extends to a spanning
tree in their union, so both edges lie in one fully located tree.  If the
pairs coincide, the lasso has instead exposed a repeated projected owner
edge.  To decide whether it is the same physical overlap component one must
compare the unordered tooth-instance pairs

```text
{(d_j,N_i),(h,H_(h,i))},                             (39t)
```

allowing target and owner roles to swap; the owner address `H` alone is not
enough.

There is a sharper two-wall pigeonhole.  Put `v=L+m`, the number of distinct
lasso labels and functional source edges.  For each of the two wall
occurrences of each blocker tooth, choose one strict owner.  If no chosen
owner lies on the lasso and no outside owner label is reused, these choices
inject `2v` occurrences into only `6-v` outside labels:

```text
2v<=6-v.                                              (39t')
```

Since `m>=2`, this is possible only for `v=2`, hence only for the
slowest-rooted shape `(L,m)=(0,2)`.  Every lasso with at least three vertices
already forces lasso-owner incidence or repeated outside ownership.  This is
an incidence/reuse split, not yet extra measure: different wall occurrences
can still retrace one physical seam.  At an owner wall there is an integer

```text
R_(h,i)=h(14N_i+sigma_i)-14d_jH_(h,i),
|R_(h,i)|<d_j,               gcd(h,d_j) divides R_(h,i). (39u)
```

If the same `h` owns two walls, then

```text
||h(w_i-w_r)||<1/7,                                  (39v)
H_(h,i)-H_(h,r)
 =h[(w_i-t0)-(w_r-t0)]-[alpha_(h,i)-alpha_(h,r)] in Z. (39w)
```

The last integer is the exact sampled tooth displacement.  Finally the wall
and shifted-clock torsion sheets intersect only through the carrier:

```text
gcd(gcd(d_j,h),gcd(c+d_j,c+h))=gcd(c,d_j,h).         (39x)
```

This isolates the live arithmetic: a turn/reuse invoice must couple signed
tooth displacement to the `c+d_i` holonomy, not merely count owner labels.

### 8.5 Protected-tree tariff and its exact scalar frontier

There is a useful post-THM-1250 scalar corollary.  Put

```text
S=c sum_i 1/d_i,
x=d1/c,
B(x)=(9x+2)/[3x(1+x)],
w_uv=c/lcm(u,v).                                     (39y)
```

For the actual two-edge protected blocker forest `F`, THM-1244 gives

```text
S>=B(x)+(7/12)w(F),                                  (39z)
```

while extending `F` inside the union of `F` with THM-1250's connected
chronological handoff graph gives a tree `T` containing `F`, and hence

```text
S>=1+(7/12)[w(F)+w(T\F)].                            (39aa)
```

Here `F` ranges over two-edge forests on the five blocker labels
`{d2,...,d6}`, while

```text
R(F)=min_{T spanning all six fast labels, F subset T} w(T\F).
```

Eliminating the unknown protected forest gives the exact tariff

```text
S>=Phi:=B(x)+min_F[(7/12)w(F)
       +((7/12)R(F)-(B(x)-1))_+].                    (39ab)
```

The base is uniformly strict because

```text
B(x)-258/247
 =[(13-6x)(129x+38)]/[741x(x+1)]>0.                 (39ac)
```

This genuinely eliminates new arithmetic cells.  For

```text
c=52,311,
(d_i)=(98,700,111,300,747,300,1,046,220,1,743,700,2,615,550), (39ad)
```

every fast-pair lcm is `100c`, while the common fast gcd is only `10`.
Here

```text
S=117/100,
B=26659/22950,
S-B=77/9180,
Phi-S=301/91800>0.                                   (39ae)
```

Thus the base estimate passes but the positioned two-edge tariff closes the
cell outside the common-gcd branch.

The tariff also marks its own limit.  On the ordered speed path,
`w_(i,i+1)<=c/d_i-c/d_(i+1)`, so one unforced low-weight tree has total weight
less than one; this minimum-tree scalar cannot force a large invoice.  More
sharply, for

```text
(d_i)=(c+1,2c-1,4c-1,20c-1,50c-1,75c-1)            (39af)
```

for the explicit value of `c` below, all fast pairs are coprime and

```text
1+1/2+1/4+1/20+1/50+1/75=11/6=B(1).                (39ag)
```

The exact referee evaluates `Phi` at
`c=10,673,449,340,400` and the packet still passes.  Hence the current
protected-tree tariff does not settle the tangent frontier; phase, offset,
or typed-word information is still essential.

## 9. What address data remains

The quotient has removed the unbounded tooth lift, not the actual speed
packet or its torsion.  Within the sampled centered subsystem, the remaining
unbounded phase-address coordinate can be written

```text
E_i=2c epsilon_i
   =2cP_i-(2k+1)Q_i in [-c,c] intersect Z,
E_i==-(2k+1)Q_i                         (mod 2c).      (40)
```

Thus all spokes share one odd multiplier `2k+1` on the modulus `2c`.  The
finite `delta` word, the contracting errors, the exact `gcd(Q_i,Q_j)` sheet,
and this common residue orbit reconstruct every selected wall, candidate
owner residue, and signed owner expiry through (39i)--(39q).  Absolute scale
and gcd/torsion remain deliberately visible sidecars; no claim that all
arithmetic is finite is intended.

What is not present in an orbit-edge tuple is the all-candidate incidence
state.  What is not present even in the full sampled hypergraph is the
chronological tooth-instance word *between* walls.  Coverage restores a word
by the expiry walk of Section 8.3.  Downstream, THM-1252--1256 make one
deletion-minimal chronological word coherent for every blocker, so its raw
handoffs are disjoint and arbitrary greedy retracing no longer controls the
proof.  THM-1256 further shows that actual blocker marks are phase/word
aligned.  The remaining loss is therefore not a mysterious address or a
free turn choice: it is the well-founded continuation/metric conversion of
one aligned typed cell.

## 10. Seven exact guardrails

The finite relative bank cannot be replaced by a bound on absolute tooth
addresses.  For every `c>=27`, set `k=c-1` and take the fast packet

```text
{c+1,c+2,c+3,c+4,2c,3c}.                             (41)
```

At the six centered spokes, the first four labels may point to `2c` and
`2c<->3c` is an exact blocker two-cycle with relative word `(0,1)`, while
the absolute tooth addresses grow linearly with `c`.  The full seven-speed
packet is primitive.  Nevertheless

```text
t*=1-2/(5c),             min_v ||vt*||=1/5.           (42)
```

These packets are globally lonely, not six-comb covers.  They prove that
primitivity, compact ratios, sampled blocker obligations, and cycle
contraction alone cannot control the absolute lift.

Nor can every ascent digit be forced binary.  At `c=1,k=0`, the packet

```text
(1;2,16,17,34,35,2343)                               (43)
```

satisfies all THM-1233 adjacent-ratio bounds and splits into the centered
two-cycles

```text
2<->2343,              16<->17,              34<->35. (44)
```

On `2->2343`,

```text
(Q_i,P_i,N_i,M_i,delta_i,ell_i)=(3,2,1562,1562,-390,2), (45)
```

while the three cycle words are `(-390,1),(1,0),(1,0)`.  Yet `t=1/6` gives
the seven-speed packet depth `1/6`.  This is again a globally lonely
guardrail, not a cover counterexample.  The full alphabet (18), or a
coverage-sensitive replacement for it, is genuinely necessary.

A third guardrail shows that a short tail need not carry a monotone phase or
wall potential.  At `c=1,k=0`, the packet

```text
(1;2,3,4,5,6,7)                                      (45a)
```

admits the functional word

```text
2->3->4->5->6->7->6,
tail delta=(1,0,1,0),       tail E=(-1,0,-1,0,-1),
tail walls=(5/14,27/56,29/70,41/84).                 (45b)
```

The speeds increase throughout the proper tail, but phase sides and wall
coordinates alternate.  The packet is lonely at `t=1/8`.

Fourth, orbit-local data do not determine sampled ownership.  The local
word `2->3->4->5->4` is identical in

```text
(1;2,3,4,5,6,7),
(1;2,3,4,5,6,14),
(1;2,3,4,5,6,28).                                   (45c)
```

At its first wall `5/14`, the first packet has no extra owner, while the
other two have owners `14` and `28`.  All three are globally lonely.  The
full quotient distinguishes them through `(Q_h,E_h)` and (39o); the smaller
orbit quotient does not.

Fifth, take `c=2,k=0`, orbit labels `{3,4,5,6,14}`, and one external label
`h in {13,19,39}`.  Every centered vertex has a unique blocker and the forced
slowest lasso is

```text
3->5->14->4->6->4,                 h->4.             (45d)
```

The five orbit rows have

```text
(delta,ell,w)=
(1,2,3/14),(0,4,55/196),(0,8,13/56),
(1,2,5/28),(0,4,13/56),                              (45e)
```

and source-wall depths

```text
(5/14,79/196,1/4,2/7,11/28).                         (45f)
```

For `h=39`, label `14` owns `3/14` and `39` owns every other distinct lasso
wall.  Nevertheless those three `39` events lie in tooth addresses

```text
(w,H,alpha)=(55/196,11,-11/196),
            (13/56,9,3/56),
            (5/28,7,-1/28).                          (45g)
```

Repeated owner label is not repeated tooth component.  The full packets are
lonely at `t=1/8` (depth `1/8` in the `h=39` case).  Thus even unique
blockers, a nonempty tail, exact source-safe walls, and owner payment at every
sampled wall do not constitute a continuum proof.  The active tooth-instance
word between those samples is essential.

Sixth, the exceptional two-wall pigeonhole branch is arithmetically
realizable without coverage.  At `c=2,k=0`, the slowest-rooted two-cycle

```text
4<->6,       outside labels (15,19,28,43)
```

has centered phases `(1/6,1/4)`, binary edge data

```text
(delta,ell)_(4->6)=(1,2),
(delta,ell)_(6->4)=(0,4),
```

and its four target-tooth walls are owned uniquely, in wall order, by
`19,28,43,15`.  No owner lies on the lasso and no outside label repeats.
The packet is lonely at `t=1/8` with depth `1/8`.

Nor does distinct two-wall ownership force the ascent digit to be binary.
The packet

```text
(2;4,18,11,13,25,29)
```

has the same centered phases, unique wall owners `25,29,13,11`, and

```text
(delta,ell)_(4->18)=(2,2),
(delta,ell)_(18->4)=(0,10).
```

It too is lonely at `t=1/8`.  Thus (39t') is the exact combinatorial split,
not a scalar contradiction; the coherent cover word is essential.

Seventh, even the greedy event lift need not create a new Hunter occurrence
at a turn.  At `c=1,k=0`, choose

```text
(1;2,3,4,5,12,13),
2->3->12->13->4->5->12.                              (45h)
```

The adjacent edges `3->12` and `12->13` have opposite travel signs and walls

```text
w_-=83/168,                    w_+=85/182.            (45i)
```

On the incoming sweep, the greedy handoff from the `2`-tooth to the
`13`-tooth is

```text
(13/28,85/182),                 length 1/364.         (45j)
```

At the outgoing wall, target `13` and owner `2` have *exactly the same*
overlap: `H_2=1`, `alpha_2=-6/91`, and (39q) again gives `1/364`.
The stitched walk immediately retraces one located seam in the opposite
abstract direction.  The packet is lonely at `t=1/8`.  Thus an arbitrary
greedy continuation closes component transport but does not, by itself,
prove a new or disjoint multiplicity credit.  THM-1252--1256 subsequently
remove this free-choice loss by selecting every blocker from one coherent
deletion-minimal word; the example remains the exact guardrail explaining
why that coherence is necessary.

## 11. Tournament and alternate-vertex audit

For the runner tournament, use the pairwise observable

```text
D_ij=P_iQ_j-P_jQ_i=Q_iQ_j(t_i-t_j),                  (46)
```

orient by increasing lexicographic `(t_i,i)`, and use the runner index to
break `D_ij=0` ties.  The tie Hamiltonian path is the same sorted order.  It
is transitive, with score histogram `(0,1,2,3,4,5)`, no directed cycles,
six singleton SCCs, no edge flips after the gauge is fixed, and one
Hamiltonian path.  It loses the central remainder, relative digit, blocker
truth, gcd sheet, and wall ownership.

The proof-bearing directed relation is instead

```text
i -> j  iff d_j is dangerous at the selected centered phase t_i,          (47)
```

decorated by (14).  It is a loopless functional digraph rather than a
tournament: every vertex has outdegree one.  Within the `d1` component, the
eventual cycle is its unique nontrivial SCC and the proper-tail vertices are
singleton SCCs; other components may have their own cycles.

At a wall-to-spoke segment, the useful tournament instead has active *tooth
instances* as vertices.  Its observable is the difference of forward expiry
endpoints, its switch is the travel sign `sigma_i`, and its tie Hamiltonian
path is the endpoint order.  It is locally transitive; choosing the
Hamiltonian maximum is exactly the interval-cover greedy step.  Concatenated
over the lasso, these local tournaments retain tooth-component identity and
produce the piecewise-monotone boundary-event walk.  The macro turns prevent
their local orders from collapsing to one global transitive order.

We challenged runners, gaps, fixed sections, tooth addresses, phase
determinants, residues, gcd sheets, wall crossings, safe components, overlap
seams, private owners, tooth instances, boundary events, and proof
obligations as vertices.  The smallest faithful carrier found here is

```text
(c; all (Q_i,E_i); finite blocker word; positive ell residues; gcd sheets;
 derived signed owner-germ hypergraph; greedy boundary-event walk).        (48)
```

It preserves blocker truth, phase chronology, sampled ownership, oriented
germs, and tooth-component continuation while discarding the global absolute
lift.  THM-1252--1256 refine (48) to one marked deletion-minimal tooth word:
the chronological and centered-phase orders are two transitive tournaments,
and actual blocker edges lie in their agreement relation.  The still-missing
operation is an aligned typed-cell descent or metric conversion.  Another
runner tournament or unforced MST scalar destroys precisely that datum.

## 12. Verification and scope

The dependency-free exact referee replays (7), (9), (12)--(18), (21)--(32),
the tail identities and constants (39c)--(39g), sharp target clearance
(39k)--(39n), source and all-candidate owner reconstruction (39j),(39o),
signed owner margins, the tariff (39ab), the exact excluded lcm cell, the
pairwise-coprime tangent guardrail, all seven non-cover guardrails, and the
tournament methodology using only integers and exact rational arithmetic.
Normal and optimized outputs are byte-identical.  It does not computationally
certify the continuum wall or greedy interval-cover arguments; those remain
explicit paper topology providers.

The Lean module kernel-checks the lower-gap cocycle and translation gauge,
determinant and gcd identities, recentered blocker identity, central band,
global and binary digit bounds, binary phase separation, the digit-selected
germ side, exact germ imbalance, reconstruction of `X_i`, centered-residue
edge transport, normalized tail transport, wall-coordinate reconstruction,
the `>5d_j` target-clearance arithmetic and its `>5/14` phase consequence,
the universal endpoint-gcd quantum, the two-wall owner pigeonhole, reciprocal blocker residue, tail
length at most four, the slowest-rooted `13/19` cap, two-seam Hunter
rearrangement, owner pigeonhole, selected-wall source residue, all-candidate
owner-residue reconstruction, oriented owner-germ width, protected-tree
tariff combination, the strict `258/247` base, affine triangle transport,
positive triangle holonomy, and the wall margin.  Arbitrary-cycle
telescoping, nearest-integer selection, circle-norm interpretation, greedy
interval continuation, and gcd-quantized overlap topology remain paper
providers.  There are no proof placeholders or `native_decide` calls.

Frozen hashes are

```text
source         0595ced97eaf4469752c9334b93119f357f7645875fa4d8944012d6ef18c08eb
output         97fbe205630b5cae842ebb87a25ac841cad45252ce54ec83ee0c8b4e8b7ff88f
formalization  4b135cb3eda4e26db0db07403becc3680cc64d2c4b97e7f8a8d4d8c639486191
```

THM-1248 now proves uniform relative-address compression, a finite
slowest-rooted lasso, two located wall seams per blocker tooth, exact sampled
owner germs, and a canonical tooth-instance lift from the slowest wall
through the tail and eventual cycle.  The two-wall pigeonhole leaves only a
slowest-rooted two-cycle if neither lasso incidence nor outside-owner reuse
occurs.  It does not prove six-comb noncoverage or LRC(14).

Downstream THM-1252--1256 choose every blocker on one deletion-minimal tooth
word, pay every handoff, remove the arbitrary-retracing loss, and force actual
blocker phase/word alignment.  The principal live target is therefore no
longer address compression, owner discovery, tree placement, existential
tail transport, or a generic turn invoice.  It is a well-founded
complete-gap/address descent for the aligned typed cell, or a metric
conversion beyond its already-counted handoffs.  THM-1247 and the guardrails
show why sampled incidence or an unforced scalar tree alone is insufficient.
