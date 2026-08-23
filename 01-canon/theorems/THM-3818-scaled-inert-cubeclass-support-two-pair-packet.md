---
id: THM-3818
title: "Scaled inert cube classes recover the support-two pair packet"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the 5,855
  admissible support-two ratios from THM-3793, the scaled cube value recovers
  the common scale and primitive ratio for every positive common scale: the
  5,855 primitive cube sums have distinct rational-cube classes.  The placed
  primitive covector recovers
  the oriented labelled runner pair and opposite exposed facets, and the
  thirteen ambient residues modulo the pair sum recover exactly the 1/14-
  safety mask on the pair-sum time lattice.  There are 456,690 unoriented
  supports and 913,380 oriented assignments.  The support-two covector lies
  inside the rank-eleven bounded-relation span and cannot be its twelfth rank
  increment.  The full decoder-row matroid is exactly graphic: rank eleven
  forces a disconnected pair graph, while connected incidence gives the sharp
  finite bound max(n_i)<=355^12.  Common exposed facets obey a separate
  two-level orientation law.  In the exact two-component rank-eleven branch,
  the orthogonal integer lattice has two scale coordinates; if the full
  height-91^6 support-three code still has rank eleven, every bounded crossing
  row is forbidden and every internal pair has height at most 91^6.  Full
  packets recover both scales unless one component is a singleton.  This is
  an all-scale address theorem, not an LRC exclusion: the grid can be empty
  while an off-grid loneliness time exists.
source: root + lrc_reversible_address / incoming-signal extension, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (thm3791-hostile-audit, 2026-08-23).  The
  all-scale decoder, 94/5,855 census, support/orientation quotient, literal
  signed covector, eleven-dimensional opposite facets, inclusive residue
  boundary, same-packet off-grid information loss, AP13 hostile, and split/
  exponent-three collision controls were independently rederived.  The audit
  caught the provisional factor-two half-atlas defect recorded as MISTAKE-457;
  the repaired primary and independent assertion-free companions enumerate
  all 913,380 oriented assignments in 5,504,236 and 2,752,607 active gates.
  A separate factor/scan audit proves all 5,855 rational-cube classes distinct,
  checks 17,137,585 pairs, and retains taxicab and arrival hostiles.  A
  17,596-gate sidecar finds one square, eight triangular and no
  square-triangular primitive addresses.  Normal/-O/frozen streams agree for
  both packet companions and the arithmetic sidecars.
  A graphic-matroid extension then checks the diagonal incidence equivalence,
  all 5,855 ratios, 245,220 projective ratio triangles, 46,136 full-decoder
  triangle circuits, every subset of a mixed five-speed hostile, the sharp
  thirteen-vertex path, and independent rank/facet controls.  Its 36,461,514
  active requirements agree under normal, optimized, and frozen replay.
  The two-component extension adds 1,742 assertion-free exact requirements:
  it checks the quotient slope atlas, a safe 2+11 row with both scales
  recovered, and a 12+1 singleton packet fibre of size 14,077,914,720,208.
  Normal, optimized, and frozen streams agree.
depends_on:
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
  - THM-3793-inert-prime-sum-all-scale-two-cube-singleton
  - THM-2052-finite-height-forces-high-rank-bounded-relation-code
related:
  - THM-778-centered-christoffel-endpoint-skew-product
script: 04-computation/lrc14_scaled_inert_cubeclass_pair_packet_thm3818.py
output: 05-knowledge/results/lrc14_scaled_inert_cubeclass_pair_packet_thm3818.out
script_sha256: 282e79f865d943d6b7435d3375d21e28de1f7fe2ed438979cc913c1c1ce338ca
output_sha256: 9ca40f84ce8b9ce02cb22197832700699c5674a39578acc50898e276eaa29cdd
semantic_sha256: 03be043c6db1e7679e0bb858ee4e419f482133b0d459c2ee0bbfbf7237b0237f
independent_script: 04-computation/lrc14_scaled_inert_cubeclass_pair_packet_independent_audit_thm3818.py
independent_output: 05-knowledge/results/lrc14_scaled_inert_cubeclass_pair_packet_independent_audit_thm3818.out
independent_script_sha256: aaab79b1a08cad7eed01d23412c1c286662b63679d14b7e4a22dbc96592cc2c2
independent_output_sha256: a859f2863660389caf706c7892cf8f564712519539dea744d11437cb4391aeb5
independent_semantic_sha256: 998891530735604ea2ede50d0296cf003dabdfb777d585431192759c48f6989c
cubeclass_script: 04-computation/lrc14_scaled_inert_cubeclass_rational_class_extension_thm3818.py
cubeclass_output: 05-knowledge/results/lrc14_scaled_inert_cubeclass_rational_class_extension_thm3818.out
cubeclass_script_sha256: e37ae5e8adf3414081e76d56eb01607172a6364703dca76b41ac096c9f6d77c3
cubeclass_output_sha256: 6199147668a19c0a2ed403b14d889caf65078c675068f000856202cb68a861a8
cubeclass_independent_script: 04-computation/lrc14_scaled_inert_cubeclass_rational_class_extension_independent_audit_thm3818.py
cubeclass_independent_output: 05-knowledge/results/lrc14_scaled_inert_cubeclass_rational_class_extension_independent_audit_thm3818.out
cubeclass_independent_script_sha256: 0bb6b60238335695addce0010182f5222da4c4a86c6f59912bd0f36123d063ca
cubeclass_independent_output_sha256: 37aa867d653d7123435f9841cdc772c3cd71c289162006f80aeb6079a17ad7bc
cubeclass_independent_semantic_sha256: 86c481b03158ba3cb7024ef8739fde640be331d9856ce5d4fc0c1b3b4fcc06cb
square_triangular_script: 04-computation/lrc14_cubeclass_square_triangular_sidecar_thm3818.py
square_triangular_output: 05-knowledge/results/lrc14_cubeclass_square_triangular_sidecar_thm3818.out
square_triangular_script_sha256: 0050033a2782e0adf4ab86aaca74db43b7e879caf48efdc0f2fcb7da31a9e963
square_triangular_output_sha256: 673a764dd14e1b06c549a80cf72e747376d0fba6d28f17da3ac778335f3f9b15
square_triangular_semantic_sha256: e80c76bcbbf9a83c6cb87983ababe96906c1160581aaee8079b6325051eb1e53
graphic_script: 04-computation/lrc14_cube_covector_graphic_matroid_extension_thm3818.py
graphic_output: 05-knowledge/results/lrc14_cube_covector_graphic_matroid_extension_thm3818.out
graphic_script_sha256: 1abba30659a81b8524b5417ff1c7c2a7915a042b365efead087a838c4c8d9497
graphic_output_sha256: 6c985414db35c3ba05f20d6e800eba5efa2f394a4826ccb772e30d7c119be64e
graphic_semantic_sha256: c98a37e4f7bf777bf6c3027bdcf4f6bc20daa83226f382132bf68f9b3ac9165d
two_component_script: 04-computation/lrc14_two_component_decoder_quotient_extension_thm3818.py
two_component_output: 05-knowledge/results/lrc14_two_component_decoder_quotient_extension_thm3818.out
two_component_script_sha256: 2ce52b1a45c3c4416e367f8177b288f7f9dc7ac651bf58ee036137e17f2208b2
two_component_output_sha256: df4b1c86df9f9871997c41bc9d236429d0ca61fae4125bfe0ce51bfff78004e6
two_component_semantic_sha256: 3c937069964a6e2134e7da6c6009d1b3b71eecdbc2d5260807ea34bffe69ab3b
hash_basis: raw LF bytes
---

# THM-3818 -- the scaled cube address retains one exact finite time grid

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Packet and statement

Let `n=(n_1,...,n_13)` be a row of positive integer speeds.  Choose distinct
labels `i,j` and suppose, after orienting the pair, that

```text
n_i=gp,                 n_j=gq,                 1<=p<q,
gcd(p,q)=1,             p+q<=356.                         (1)
```

Put

```text
D=g(p+q)=n_i+n_j,       M=n_i^3+n_j^3,
a=q e_i-p e_j,          rho_l=n_l mod D.                  (2)
```

Assume every prime divisor of the primitive sum `p+q` is `2 mod 3` and has
exponent at most two.  The common scale `g` is arbitrary.  Then the packet

```text
(M,a,D,rho_1,...,rho_13)                                  (3)
```

recovers exactly:

1. the common scale `g` and primitive ratio `(p,q)`;
2. the oriented labelled pair `(i,j)`, its primitive relation `a`, and the
   opposite exposed facets selected by the signs of `a`; and
3. the complete pair-sum-lattice safety set

```text
K_D(n)={k in Z/DZ : ||k n_l/D||>=1/14 for every l}.        (4)
```

There are exactly `5,855` primitive ratios satisfying these hypotheses.  They
give `5,855*C(13,2)=456,690` unoriented label supports and
`5,855*13*12=913,380` oriented assignments at every allowed scale.  This is an
exact decoder for the displayed packet, not for the other eleven speeds,
their integer quotients by `D`, or any time outside `(1/D)Z/Z`.

## 2. The cube value recovers scale and ratio

First quotient by rational cubes.  On the finite set of all `5,855`
admissible primitive ratios, prime-exponent reduction modulo three and an
independent maximal-cube-divisor scan prove

```text
(p^3+q^3)/(r^3+s^3) in (Q_(>0)^x)^3
  iff (p,q)=(r,s).                                         (5a)
```

Since `M=g^3(p^3+q^3)`, its rational-cube class recovers `(p,q)` and then
`g=(M/(p^3+q^3))^(1/3)`.  Thus finite-atlas decoding works for arbitrary
positive `g`.

There is also a table-free divisor decoder on the stronger physical
all-inert subbranch.  If every prime divisor of `D=g(p+q)` is `2 mod 3`,
apply THM-3793 to

```text
x=gp,                  y=gq,
gcd(x,y)=g,            (x+y)/g=p+q.                       (5)
```

The hypotheses above are precisely its all-inert pair-sum condition and its
cube-free primitive-quotient condition.  Hence `M` has exactly one positive
distinct two-cube representation:

```text
M=(gp)^3+(gq)^3.                                           (6)
```

The representation is effectively recoverable.  For every possible sum
`d=x+y`, necessarily `d|M`, and

```text
xy=(d^3-M)/(3d).                                           (7)
```

Enumerating the positive divisors of `M`, retaining integral positive (7)
with square discriminant, therefore finds the unique pair.  Taking its gcd
recovers `g`; division by `g` recovers `(p,q)`.  This proves the all-scale
table-free decoder without bounding the inert prime powers in `g`.

## 3. The covector recovers placement and the facet pair

The vector `a` in (2) is primitive and satisfies

```text
a dot n=q(gp)-p(gq)=0,                ||a||_1=p+q<=356.    (8)
```

It has one positive coordinate `q` at `i` and one negative coordinate `-p`
at `j`, so its signed support recovers the oriented placement.  On the
projected lonely-runner zonotope of THM-3743, pairing with `a` is the same as
pairing before projection.  Thus its maximum fixes the `i` cube coordinate
at the upper endpoint and the `j` coordinate at the lower endpoint; its
minimum fixes the reverse pair.  The other eleven coordinates are free.
Because the positive speed vector cannot lie in that free-coordinate
subspace, projection preserves dimension eleven, so these are opposite
exposed facets.  No owner or facet information is inferred from `M` alone;
it is retained by `a`.

## 4. Ambient residues recover the pair-sum grid

For every label `l` and residue class `k mod D`, let

```text
r_(l,k)=k rho_l mod D,          0<=r_(l,k)<D.              (9)
```

Then

```text
||k n_l/D||=min(r_(l,k),D-r_(l,k))/D.                    (10)
```

Consequently (4) is given by the integer-only rule

```text
k in K_D(n)
iff 14 min(r_(l,k),D-r_(l,k))>=D for every l.            (11)
```

This proves that the thirteen residues, together with `D`, recover the whole
finite safety schedule.  Adding arbitrary multiples of `D` to any ambient
speed leaves (11) unchanged, which is exactly the information destroyed by
the residue quotient.

## 5. Failure boundary

The schedule is not a loneliness certificate.  In the safe arithmetic
progression row

```text
n=(1,2,...,13),              (p,q,g,D)=(1,4,1,5),        (12)
```

runner `5` lies at the origin at every time `k/5`, so `K_5(n)` is empty.
Nevertheless `t=1/14` has

```text
||l/14||>=1/14                 for 1<=l<=13.              (13)
```

Thus even the complete pair-sum grid can miss a valid off-grid time.  The
cube packet also destroys the other eleven exact speeds, other relations,
global wall chronology, arrival, owner, and every off-grid value.

The arithmetic hypotheses are sufficient, not necessary.  The controls

```text
1729=1^3+12^3=9^3+10^3                       (split prime),
515375=15^3+80^3=54^3+71^3                   (5^3 pair sum) (14)
```

show the two first failed mechanisms when the inert and exponent conditions
are dropped.  No converse or classification outside the stated packet is
claimed.

## 6. The full decoder system has a weighted graphic matroid

Form the full decoder graph `G_dec(n)` on the thirteen labels: `{i,j}` is an
edge exactly when the reduced ratio and actual gcd scale satisfy Section 1,
so the complete packet exists for that pair.  Orient every edge from its
smaller to its larger speed.  For an edge put

```text
g_ij=gcd(n_i,n_j),       ell_ij=lcm(n_i,n_j),
a_ij=(n_j/g_ij)e_i-(n_i/g_ij)e_j.                    (15a)
```

Let `Delta_n=diag(n_1,...,n_13)`.  Directly,

```text
a_ij Delta_n=ell_ij(e_i-e_j).                         (15b)
```

All diagonal entries and row scales are nonzero.  Right multiplication by
`Delta_n` and row rescaling are invertible over `Q`, so for every selected
edge set `E`, including isolated labels in the component count,

```text
rank_Q span{a_e:e in E}=13-c(G),
{a_e:e in E} is independent iff E is a forest.       (15c)
```

The minimal dependencies are exactly simple cycles.  Traverse such a cycle,
let `epsilon_e=+/-1` compare its traversal with the lower-to-higher gauge,
and put `L=lcm_e ell_e`.  Its explicit circuit is

```text
sum_(e in cycle) epsilon_e (L/ell_e)a_e=0.            (15d)
```

Dividing the coefficients by their gcd gives the primitive circuit.
Reversing a packet orientation changes only a row sign and cannot add rank.

### 6.1 Component shapes and the sharp connected terminal

For a connected decoder component `C` with `t` vertices, put
`g_C=gcd(n_i:i in C)`.  Every primitive edge coefficient is at most `355`.
Choose a spanning tree.  Deleting one column from its weighted incidence
matrix gives by leaf elimination a product of `t-1` edge coefficients; its
signed maximal minors generate the primitive one-dimensional kernel.  Hence

```text
max_(i in C) n_i/g_C <=355^(t-1).                     (15e)
```

For a primitive connected graph on all thirteen labels this is

```text
max_i n_i <=355^12
 =4006270200351760530820556640625.                    (15f)
```

The bound is sharp inside the decoder class.  The unique ratio with a
coefficient `355` is `(1,355)`, and its sum `356=2^2*89` is admissible.
Thus

```text
(1,355,355^2,...,355^12)                              (15g)
```

has decoder graph exactly the consecutive path and attains `(15f)`.
Conversely, equality forces all selected cofactor entries to be `355`, hence
every tree edge to have ratio `1:355`.  Distinct speeds force one vertex at
each depth, so the tree is the path `(15g)`, up to label permutation and
common allowed inert scale.

Every decoder row belongs to THM-2052's `W`.  Therefore

```text
dim W=11  ==>  G_dec(n) is disconnected.              (15h)
```

If the decoder rows themselves span dimension eleven, `(15c)` gives exactly
two components.  A further edge inside a component is a cycle dependency;
an edge joining them would give rank twelve.  Conversely, a connected graph
already supplies twelve independent rows and enters the finite rank-twelve
terminal, now with the direct bound `(15f)`.  No row in that box is excluded
here.

### 6.2 Facet orientation is independent of row rank

For every selected edge, its maximizing face fixes the smaller-speed cube
coordinate at `13/14` and the larger-speed coordinate at `1/14`.  Hence

```text
intersection_(e in E) F_e^+ is nonempty
iff no label occurs with both signs
iff the lower-to-higher orientation has no directed path of length two.
                                                               (15i)
```

For nonempty compatible `E`, all incident coordinates are fixed and all
others are free.  Projection is injective on this direction space because
the positive vector `n` has a nonzero fixed coordinate.  Therefore

```text
dim intersection_(e in E) F_e^+=13-|V(E)|.           (15j)
```

The same holds for minimizing faces.  A cycle is compatible precisely when
it is even and its speeds alternate local minima and maxima; every triangle
circuit is incompatible.  Rank and facets are genuinely orthogonal:

```text
speeds (2,3,8), edges (2,3),(3,8):
  independent forest, but no common maximizing face;

speeds (1,2,3,9), edges (1,3),(2,3),(2,9),(1,9):
  rank three with a common dimension-nine face,
  6a_13-3a_23+a_29-2a_19=0.                          (15k)
```

### 6.3 Exact triangle census

The `5,855` primitive ratios contain exactly `245,220` sorted projective
triples whose three pair ratios remain in the atlas.  Requiring all three
actual gcd scales to have only inert prime divisors leaves exactly `46,136`
full-decoder triangle circuits.  The unique smallest by maximum speed and
then total speed is

```text
(2,3,8),                 4a_23+a_38-3a_28=0.          (15l)
```

This is the first literal mechanism by which repeated cube-addressed rows
fail to add rank.  The graph forgets cube values, pair-sum grids, owner,
arrival, cross-component scale coupling, off-grid time, and loneliness.

### 6.4 The exact two-component quotient

Put

```text
Q=91^6=567869252041,             B=Q^2.               (15m)
```

Assume that the decoder rows have rank eleven, so their graph has exactly
two connected components `I,J`.  Write the primitive component shapes as
`u=(u_i)_(i in I)` and `v=(v_j)_(j in J)`.  Then the full integral
orthogonal lattice is

```text
V_dec^perp intersect Z^13
  ={s u direct-sum t v : s,t in Z}.                   (15n)
```

Indeed, the internal edge rows force each component restriction to be a
rational multiple of its primitive positive kernel vector, and Gauss's
lemma makes the two multiples integral.  Thus the positive primitive speed
rows in this lattice are exactly those with `s,t>0` and `gcd(s,t)=1`.

Let `W=W_(Q,3)(n)`.  In the counterexample box `sum_i n_i<=B`, THM-2052
gives `dim W>=11`, while `V_dec` is an eleven-dimensional subspace of `W`.
There are therefore only two cases:

```text
dim W=12:  the rank-twelve finite terminal;
dim W=11:  W=V_dec.                                  (15o)
```

In the second case, let `c` be any support-at-most-three row of height at
most `Q` that meets both components, and put

```text
P_c=sum_(i in I)c_i u_i,       R_c=sum_(j in J)c_j v_j.
```

Neither partial sum can vanish: with at most three occupied coordinates,
vanishing on both sides would force a nonzero one-coordinate relation on
one component.  Hence

```text
c dot (s u direct-sum t v)=0
  iff s P_c+t R_c=0
  iff s/t=-R_c/P_c.                                  (15p)
```

Every such crossing row lies outside `V_dec`, so none occurs in the
unresolved `dim W=11` branch.  Equivalently, the projective scale `s:t`
avoids the finite atlas of slopes in (15p).  Applying THM-2052's guaranteed
triple relation to two labels in one component and one in the other now
forces the outside coefficient to vanish.  Consequently every same-component
primitive pair has height at most `Q`:

```text
max(u_i/gcd(u_i,u_i'),u_i'/gcd(u_i,u_i'))<=Q,
max(v_j/gcd(v_j,v_j'),v_j'/gcd(v_j,v_j'))<=Q.         (15q)
```

Thus every mixed pigeonhole triple degenerates to an internal pair relation;
this is why bounded support-three rows cannot add rank without entering the
terminal.

The complete cube packets sharpen the remaining quotient.  Connected edge
covectors recover each primitive component shape.  If both components have
at least two vertices, one packet modulus `D_ii'` and one `D_jj'` recover
the two scales exactly:

```text
s=D_ii'/(u_i+u_i'),             t=D_jj'/(v_j+v_j').   (15r)
```

This reconstructs the speed row, but not owner/phase arrival or loneliness.
If one component is a singleton, no packet has an internal modulus there;
the other packets see its speed only through congruences modulo

```text
L=lcm{D_e:e is an edge of the nonsingleton component}. (15s)
```

This loss is real.  The two-nontrivial-component row with

```text
u=(1,3),
v=(1,3,4,9,10,16,19,22,24,33,40),
(s,t)=(1,2^42)                                      (15t)
```

has decoder component sizes `2+11`, `W=V_dec`, sum
`796046418509828`, and both scales are packet-recoverable, yet it is lonely
at `t=1/14`.  For the singleton hostile, take the twelve-vertex shape

```text
(1,4,...,4^11),          L=22906142720.              (15u)
```

The singleton speeds

```text
4763632550625148929,     4763632573531291649          (15v)
```

give identical complete labelled packet collections.  Below `B`, the same
dominance-preserving packet fibre contains exactly `14,077,914,720,208`
positive candidates.  The component quotient is therefore arithmetically
rigid precisely when neither component is a singleton; even then it is not
an LRC certificate.

## 7. Exact verification and scope

Replay the new quotient extension with

```text
python -B 04-computation/lrc14_two_component_decoder_quotient_extension_thm3818.py
python -B -O 04-computation/lrc14_two_component_decoder_quotient_extension_thm3818.py
```

and compare after LF normalization with
`05-knowledge/results/lrc14_two_component_decoder_quotient_extension_thm3818.out`.

The companion independently enumerates every coprime `p<q`, `p+q<=356`,
checks complete coordinate fibres and the divisor decoder, and reconstructs
all 156 oriented covectors for every admissible ratio.  It also tests arbitrary
inert powers in 32 scale controls and compares (11) with direct rational
evaluation.  Normal and `python -O` streams match the frozen output.
The graphic companion checks all `5,855` ratios, every ratio pair, all
`46,136` full-decoder triangle circuits, all edge subsets of a mixed
five-speed hostile, the diagonal incidence law, facet compatibility, and the
sharp thirteen-vertex path in `36,461,514` active requirements.  Its normal,
optimized, and frozen streams also agree.
The two-component companion checks the integral quotient, crossing-slope
atlas, bounded pair-clique consequence, a safe `2+11` row with both scales
recovered, and the exact `12+1` singleton fibre in `1,742` active
requirements.  Its normal, optimized, and frozen streams agree as well.

The rational-class extension independently factors and scans every primitive
cube sum, then exhausts all `17,137,585` unordered ratio pairs and scaled
taxicab/arrival hostiles.  Its separate audit reaches the same zero-collision
answer.  The square/triangular sidecar exhausts its stated Pell-capped shells;
these classical labels add no LRC predicate.

The precise connection ledger is

```text
source:      support-two Graver branch of THM-3743
target:      cube fibre, pair-sum grid, weighted graph and oriented faces
map:         pair packet -> (M,a,D,residues) -> a diag(n)
preserved:   scale, ratio, placement, grid, row matroid, circuits, facets
destroyed:   ambient quotients, cross-component scale, off-grid loneliness
sidecar:     full speed row and owner/phase chronology for any LRC use
next test:   classify the finite forbidden crossing-slope atlas together
             with facet/owner/phase arrival; packets alone are exhausted. (16)
```

In particular, (3) supplies no exclusion of a hypothetical LRC(14)
counterexample and proves no case of the Lonely Runner Conjecture.  There is
also a sharp rank boundary.  The covector `a` has support two and
`||a||_1<=356`; THM-3743, Section 4, therefore puts it automatically inside
THM-2052's bounded-relation span `W`.  A cube-addressed row can refine the
internal rank-eleven code, but it can never be the outside-`W` twelfth
relation that triggers the rank-twelve finite terminal.  Section 6 now
classifies all support-two incidence: an unresolved rank-eleven decoder graph
has at least two components, exactly two when its rows span eleven.  In that
two-component branch, bounded crossing rows are exactly a finite forbidden
slope atlas, all internal pairs are already bounded, and the full packet
either recovers both scales or retains the exact singleton congruence fibre.
The live problem is the crossing-slope complement together with owner/phase
arrival—not another support-two count.  **QED.**
