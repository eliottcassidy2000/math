---
id: THM-801
title: The seven-term tiling recursion lifts to exact three-face Cech descent, pure cubic colour interaction, and tournament/Smith Mobius curvature
status: PROVED (three-face tiling and line descent, colour laws, interval-face Mobius and Smith identities) + FINITE-EXACT (Omega/S2 codec through n=8; node/fibre classifications through n=7)
source: codex-2026-07-15-S12/S13/S11
depends_on: [THM-442, THM-553, THM-785, THM-790, THM-796]
related: [THM-549, THM-550, THM-781, THM-805, THM-809, HYP-2685, HYP-3234, HYP-6825, HYP-6870, HYP-6880]
verification:
  - 04-computation/mobius_cech_metagraph_codec_codex_S12.py
  - 05-knowledge/results/mobius_cech_metagraph_codec_codex_S12.out
  - 05-knowledge/results/mobius_cech_metagraph_codec_codex_S12.json
  - 04-computation/mobius_cech_n8_frontier_codex_S13.cpp
  - 05-knowledge/results/mobius_cech_n8_frontier_codex_S13.out
  - 04-computation/three_sorted_metagraph_continuation_minimization_codex_S11.py
  - 05-knowledge/results/three_sorted_metagraph_continuation_minimization_codex_S11.out
---

# THM-801 — Möbius/Čech descent for the merged metagraph

The historical recursions

```text
full:       A+B+C-D-E-F+G,
even half:  A+B-C,
odd half:   A+B-C+D-E-F+G
```

are scalar shadows of restriction diagrams.  At the level of tilings, lines,
and tournaments, their overlap slots retain phases and rooted-path data that
the scalar signs cannot see.  This theorem lifts the full seven-term chart to
an exact three-face descent, identifies the missing middle face, and isolates
two irreducible interactions:

1. a pure cubic blue/black colour interaction;
2. the complete `B2` Möbius curvature of legal marked-path faces, whose
   tournament and electrical readings are respectively boundary cyclicity and
   the two coordinates of the transitivity-flow cloud.

## 1. The full `B3` cover is an exact tiling descent

Write a tile as `(a,b)`, with `1<=b<a<=n` and `a-b>=2`, and put

```text
(r,c)=(a-b-1,b).
```

Thus the staircase carrier is `r,c>=1`, `r+c<=n-1`.  THM-553's three
size-`n-1` subtriangles give restriction maps

```text
d_A: a<n,       (a,b) -> (a,b),
d_B: a-b>=3,    (a,b) -> (a-1,b),
d_C: b>=2,      (a,b) -> (a-1,b-1).                           (1.1)
```

Here `A` and `C` are respectively THM-796's high- and low-end deletion faces.
The new face `B` contracts the interval-root gap coordinate.  It is a valid
lower staircase tiling, but it is not induced deletion of a fixed tournament
vertex.

The three face sets cover every tile.  Their pairwise intersections are three
copies of the size-`n-2` staircase and their triple intersection is the
size-`n-3` staircase.  Therefore

```text
|S_n|=3|S_(n-1)|-3|S_(n-2)|+|S_(n-3)|,                       (1.2)
```

which is THM-442's `A+B+C-D-E-F+G` identity.  More strongly, a triple of face
tilings glues to a unique upper tiling exactly when it agrees on every pair
overlap.  This is the ordinary sheaf property of binary assignments on a
union, but here it gives an executable inverse, not only a count.

All three restrictions commute with the all-tile complement `kappa`.

## 2. The triple overlap is the complement-line phase witness

Let `E_n=T_n/<kappa>` be the literal complement-line set.  Apply (1.1) to a
line and retain the three bare lower lines.  For `n>=6`, this is a bijection

```text
E_n  ~=  {
  (e_A,e_B,e_C) in E_(n-1)^3 :
  the three corresponding pair-overlap lines agree
}.                                                               (2.1)
```

Proof: choose arbitrary endpoints of the three lower lines.  Agreement of the
bare `ij` overlap lines gives a relative phase `delta_ij in F_2`.  Since the
triple overlap `S_(n-3)` is nonempty for `n>=6`, restricting all three phase
comparisons to it gives

```text
delta_AB+delta_AC+delta_BC=0.                                  (2.2)
```

Hence there are endpoint phases `s_A,s_B,s_C` with
`s_i+s_j=delta_ij`, unique up to simultaneous complement.  The coherently
oriented faces glue to one tiling; the remaining simultaneous choice is
exactly the upper `kappa` quotient.  This proves (2.1).

The small cases expose why `G` matters.  At `n=5` the triple overlap is empty:
there are 64 pairwise-compatible bare triples, split evenly into 32 zero- and
32 nonzero-holonomy triples, while only the 32 zero-holonomy triples come from
upper lines.  At `n=4`, four upper lines collapse onto one bare face triple.

This contrasts with THM-796's two-end cover

```text
S_n=S_(n-1)^A union S_(n-1)^C union {(n,1)}.
```

The faces miss the apex, so their bare-line map retains a two-sheeted apex
phase.  The full `B3` cover includes the apex in `B`, and its nonempty `G`
overlap kills the phase holonomy.

## 3. Expanded line/node tensor and exact finite codecs

For the apex-zero endpoint `t` of an upper line, let `(u,u')` be the ordered
upper merged-node pair and `(a,a'),(b,b'),(c,c')` the corresponding ordered
node pairs on the three faces.  Define

```text
Omega_n(u,u';a,a';b,b';c,c';UABC)                             (3.1)
```

to count lines with those data and their upper/three-face colour word.
Dropping the `B` pair and its colour, then ordering `(C,A)` as `(low,high)`,
recovers THM-796's `Xi_n` exactly.

The complementary folded `B2` sidecar uses THM-553's mirror-crossing clock
`tau=a+b-1`.  For every canonical half representative with `tau<n`, retain
the ordered bit pair `(t_s,t_(sigma s))` and count its four values separately
on each `tau` layer.  On the fixed layer `tau=n`, retain the zero/one counts.
Call the resulting raw integer vector `S2(e)`.  Then

```text
e is blue  <=>  every off-diagonal 01/10 count in S2(e) is zero. (3.2)
```

A coarse `B3` sidecar `S3(e)` counts one-bits in each of the seven local
membership strata `A,B,C,AB,AC,BC,ABC`.  The finite exact classifications are:

| n | lines | `Xi` cells | `Omega` cells | `Omega+S2` cells |
|---:|---:|---:|---:|---:|
| 4 | 4 | 4 | 4 | 4 |
| 5 | 32 | 32 | 32 | 32 |
| 6 | 512 | 509 | 510 | 512 |
| 7 | 16,384 | 16,031 | 16,308 | 16,384 |
| 8 | 1,048,576 | -- | -- | 1,048,576 |

At `n=7`, the missing gap face removes 277 of `Xi`'s 353 collision excess.
The remaining 76 double collisions are all separated by `S2`.  This proves
finite injectivity of `Omega+S2` through `n=8`, not an all-size or Markov-state
claim.

THM-809 proves the `n=8` row without constructing upper `n=8` node classes.
The strictly weaker lower-only key `Lambda_8` uses the three ordered lower
node pairs, `UABC`, and the `S2` mirror-layer populations.  Its collision
excess falls

```text
418 -> 252 -> 148 -> 74 -> 52 -> 0 -> 0
```

as `tau=3,...,7,fixed` are adjoined.  The 418 initial doubletons are a new
genealogy: none has both endpoint faces `A,C` equal and none uses any of the
four recursively stable `n=7` phase-square lines as a face.  Thus static
injectivity now holds through `n=8`, while continuation completeness remains
open.

The independent bounded continuation minimizer localizes the join.  Recursive
labelled refinement of `Xi_7` leaves only the phase/reflection square

```text
{0x12ca,0x12cb,0x146c,0x146d}.
```

Within each phase pair the `A,C` endpoint-face lines agree, while the `B` gap
lines are respectively `0x115/0x114` and `0x0c3/0x0c2`.  Phase one makes the
gap line a black loop at node 33; phase zero makes it a black cross-line
`21--33`.  Hence `Omega_7` itself separates this entire recursively stable
two-face residual; `S2` is needed for other, static `Omega` collisions.  This
is the finite-exact instance of Section 2's Cech statement: the two-face phase
bit has become ordinary incidence on the missing face.

Disintegrating literal lines inside their 6,126 coloured upper node-pair
fibres gives:

| n=7 sidecar | subcells | collision excess | fully separated fibres |
|---|---:|---:|---:|
| `S2` | 16,212 | 172 | 5,962 |
| `S3` | 15,016 | 1,368 | 5,071 |
| `(S2,S3)` | 16,368 | 16 | 6,110 |

Loops are counted once in this literal-line table; they are doubled only when
passing to endpoint incidence.

## 4. The missing gap-face node correspondence

For `s in {A,B,C}` define

```text
D_n^s(u,v)=#{t in F_n(u): pi_(n-1)(d_s t)=v}.                 (4.1)
```

Each row sums to `|F_n(u)|`; each column sums to
`2^(n-2)|F_(n-1)(v)|`, because every face omits `n-2` tile bits.  Reflection
swaps `A` and `C`, so `D^A=D^C`.  It fixes the gap direction and does not force
`D^B=D^A`.

| n | nodes | support cells `A/B/A+B` | primitive cells `A/B/A+B` |
|---:|---:|---:|---:|
| 4 | 3 | `2/2/3` | `3/2/3` |
| 5 | 10 | `7/5/10` | `10/7/10` |
| 6 | 34 | `34/23/34` | `34/30/34` |
| 7 | 272 | `264/244/269` | `272/262/272` |

Thus the gap support resolves five of the eight support-twin pairs left by
endpoint deletion at `n=7`, although the weighted/primitive endpoint row is
still sharper.  This is a direct tournament classification: every row entry
counts fixed-Hamiltonian-path presentations whose high endpoint, gap, or low
endpoint chart lands in a specified lower converse class.

## 5. Four-role colour law

Put

```text
M_n=C(n-1,2),                 r_n=(M_n+floor((n-1)/2))/2,
T=2^(M_n-1),                  U=2^(r_n-1),
L=2^(r_(n-1)+n-3),           Q=2^(n-3+floor((n-2)/2)),
J=2^(n-3),                    R=2^(2 floor((n-2)/2)).          (5.1)
```

Let `X,A,B,C` be the events that the upper, high, gap, and low lines are blue.
Then

```text
|X|=U,                 |A|=|B|=|C|=L,
|X cap B|=U,           |X cap A|=|X cap C|=J,
|A cap B|=|A cap C|=|B cap C|=Q,
|A cap B cap C|=R,     |X cap A cap B cap C|=J.               (5.2)
```

Every intersection containing `X` and at least one endpoint face also has
size `J`.  Equivalently, `X` implies `B`, and under `X` the endpoint-face
colours agree.  Boolean Möbius inversion gives the ten possible atoms:

| word `XABC` | count |
|---|---:|
| `BBBB` | `J` |
| `BKBK` | `U-J` |
| `KBBB` | `R-J` |
| `KBBK`, `KBKB`, `KKBB` | `Q-R` each |
| `KBKK`, `KKKB` | `L-2Q+R` each |
| `KKBK` | `L-2Q+R-U+J` |
| `KKKK` | `T-3L+3Q-R` |

All omitted words have count zero.  Constraint-graph components prove (5.2):
one face has `L` solutions; any two faces give parity-constant gap diagonals
and `Q` solutions; all three faces have `2 floor((n-2)/2)` free components;
upper symmetry forces the gap face and collapses every diagonal to the `J`
locus.  The audit agrees exactly for `n=4..7`.

## 6. The three-role colour defect is one pure cubic

Now forget the gap face and let `X,Y,Z` be upper-, low-, and high-face blue.
Set

```text
a=beta_n=U/T,             b=beta_(n-1)=L/T.
```

THM-796 gives the intersection moments

```text
p_X=a, p_Y=p_Z=b,
p_XY=p_XZ=p_XYZ=ab,       p_YZ=b^2.                           (6.1)
```

Therefore the exact probability generating polynomial is

```text
P_n(x,y,z)
 = ((1-a)+ax)((1-b)+by)((1-b)+bz)
   + ab(1-b)(x-1)(y-1)(z-1).                                 (6.2)
```

Indeed the true and product laws have the same one- and two-role marginals,
so their difference vanishes when any variable is one and must be a multiple
of the displayed cubic.  Comparing the `xyz` coefficient gives `ab(1-b)`.

Consequences:

```text
cum(X,Y,Z)=E[(X-a)(Y-b)(Z-b)]=ab(1-b),                        (6.3)
T^3 cum(X,Y,Z)=U L (T-L).
```

If `R` is an exact blue-role set, Boolean Möbius inversion separately gives

```text
P(exactly R)-P_ind(exactly R)=(-1)^(3-|R|)ab(1-b).            (6.4)
```

Encoding nonempty subsets by binary masks `1,...,7`, the signs are
`++-+--+`: the Legendre word `chi_7(1),...,chi_7(6)` with a positive top slot.
This is an actual constant-amplitude realization of HYP-3234's odd chart.
All connected Eisenstein cyclic imbalances vanish because the singleton
corrections agree and the pair corrections agree.

Conditioning converts the invisible cubic into compensating `B2` charges:

```text
Cov(Y,Z | X=1)=b(1-b),
Cov(Y,Z | X=0)=-a b(1-b)/(1-a).                               (6.5)
```

Their weighted sum is zero, explaining global pairwise independence.  Boolean
subset inversion in (6.4) and the partition-lattice cumulant in (6.3) are
related here but are not the same operation; the seven-term union formula is
not itself a cumulant.

## 7. The legal interval-face Möbius operator

Deleting an internal vertex can destroy the inherited marked Hamiltonian
path.  The legal marked-path faces are the contiguous intervals `[i,j]` of
the path, ordered by prefix/suffix deletion.  Below the full interval this is
a product of two chains.  The Möbius function of a chain vanishes across a
gap of at least two, so its entire top primitive is supported on four faces:

```text
Omega f(t)
 = f_n(t)-f_(n-1)(d_L t)-f_(n-1)(d_H t)
   +f_(n-2)(d_H d_L t).                                      (7.1)
```

Thus `Omega` is not an ad hoc finite difference.  It is the full Möbius
primitive of the legal recursive face poset.  It must be distinguished from
the full Boolean vertex-deletion lattice of unmarked tournaments and from the
Boolean tile cube used for Walsh/root packets.

For any additive substructure count, `Omega` selects exactly the copies that
meet both path endpoints.  In particular, put

```text
q(t)=Omega C3(t).
```

Then

```text
q(t)=#{v in {2,...,n-1}: T[{1,v,n}] is cyclic} >= 0,          (7.2)
C3_n(t)=C3_(n-1)(d_Lt)+C3_(n-1)(d_Ht)-C3_(n-2)(core)+q(t).   (7.3)
```

Equation (7.3) is the positive recursive form of the transitive-to-distributed
coordinate.  Since `E4=n(n^2-1)/3-8C3`,

```text
Omega E4=2(n-1)-8q.                                          (7.4)
```

For the apex-zero endpoint `s_0(e)` of each complement line,

```text
sum_e z^q = 2^(M_n-2n+7)(3+z)^(n-4).                         (7.5)
```

If `q_0=q(s_0(e))` and `q_1=q(kappa s_0(e))`, the stronger edge law is

```text
sum_e z^(q_0) w^(q_1)
 = 2^(M_n-2n+5)(1+w)^2(3+zw)^(n-4).                          (7.6)
```

Every internal endpoint triple contributes `3+zw`; complement preserves its
cyclicity.  The two boundary triples give `(1+w)^2`, so
`q_1-q_0 in {0,1,2}`.  With `r=r_n`, the blue sub-polynomial is

```text
2^(r-n+2)(1+w^2)(3+z^2w^2)^((n-4)/2)                         (n even),
2^(r-n+2)(1+w^2)(1+zw)(3+z^2w^2)^((n-5)/2)                  (n odd). (7.7)
```

Black is total minus blue.  Reflection acts freely on black lines, preserves
`(q_0,q_1)` and the projected merged-node pair, and therefore makes every
coefficient of every fixed-node-pair black refinement even.

The node polynomial

```text
K_u(z)=sum_(t in F_n(u)) z^(q(t))                             (7.8)
```

counts Hamiltonian-path presentations by cyclic triangles meeting their two
path endpoints.  It separates `3/3, 9/10, 30/34, 238/272` nodes at
`n=4,5,6,7`; adjoining total `C3` raises the last count to `249/272`.  It is a
compact horizontal fingerprint, weaker than THM-796's full primitive face row.

The common-core endpoint-node pair is an objectively related missing `Xi`
coordinate because

```text
q=C3(u)-C3(l)-C3(h)+C3(core).
```

Adding it changes `Xi_7` from 16,031 to 16,110 cells and reduces collision
excess from 353 to 274.

## 8. Smith curvature recovers the two-dimensional flow cloud

Raw-bit conventions matter.  THM-796's atlas bit `t_(xy)=0` means the arc
`x->y`; THM-790 uses the orientation bit

```text
b_(xy)=1-t_(xy).
```

Use the latter in the canonical horizontal Smith current of the staircase:

```text
J_n^h(t)=sum_((x,y) in S_n) b_(xy)/(n-1-y),
J_n^v(t)=J_n^h(sigma t).                                     (8.1)
```

Put

```text
lambda=e_1-e_n,              epsilon=e_1+e_n-(n-2),
a=b_(n,1),                   D=(n-2)(n-3).
```

Then for every tiling and every `n>=4`,

```text
D Omega J^h = a(n-2)+e_1-(n-2),
D Omega J^v = a(n-2)-e_n,                                    (8.2)
```

and hence

```text
D(Omega J^h-Omega J^v)=epsilon,                              (8.3)
D(Omega J^h+Omega J^v)=lambda+(2a-1)(n-2).                   (8.4)
```

Proof: common-core tiles cancel in (7.1).  A top-exclusive tile cancels
against its shifted low face.  A bottom-exclusive tile contributes
`1/(n-2)-1/(n-3)=-1/D`; the apex contributes `1/(n-2)`.  Reflection fixes the
apex and swaps the legs, giving (8.2).

Thus the antisymmetric primal/dual Smith curvature is exactly THM-785's
transverse black defect.  The symmetric curvature, after subtracting the apex
EMF, is THM-790's longitudinal leg current, with
`Delta E4=8 lambda` and `Delta C3=-lambda`.  The full two-dimensional
`(lambda,epsilon)` line cloud is an electrical Möbius-curvature cloud.

Blue implies zero antisymmetric curvature, but the converse is false: balanced
black lines exist.  Reflection changes the sign of `epsilon` while preserving
the merged node pair.  Hence every black node-pair fibre has a sign-symmetric
`epsilon` histogram and even zero multiplicity.  Signed `epsilon` cannot by
itself cause THM-785's oriented quotient drift; only its absolute strata and
their unequal fibre allocation can correlate with that drift.

THM-811 subsequently gives the complete joint polynomial of
`(Delta C3,q0,q1,epsilon)`.  It proves that B3 determines the two linear Smith
currents while `q` is their positional quadratic overlap, with
`Cov(q_i,epsilon^2)=-(n-4)/8`.  Its curvature node polynomial is complete
through `n=7`, but its black edge carrier is not: reflection-orbit B2/B3 data
remain necessary for positional identity.

## 9. What an exact recursive fibre must preserve

THM-796 writes a tiling as

```text
t=(core, low-leg word, high-leg word, apex).
```

For a node `u` and fixed literal core `c`, define the Boolean function

```text
I_(u,c)(p_L,p_H,a)=1[reconstruct(c,p_L,p_H,a) in F_n(u)].     (9.1)
```

Its anchored Boolean Möbius coefficients on
`L disjoint-union H disjoint-union {a}` are invertible.  Therefore the family
of core-keyed coefficient stalks is an exact node-to-tilings and
tiling-to-node carrier.  Its bidegrees distinguish same-leg, cross-leg, and
apex-containing interactions; `D_n`, `K_u`, `Omega_n`, `S2`, and `S3` are
different marginals or compact projections of this complete stalk.

This definition is exact but exponentially large.  HYP-6880 asks whether a
smaller transported subfamily remains injective and continuation-useful.  A
static finite codec is not automatically a recursive state: THM-796's strong-
lumpability failures still require the literal lower tiling/path orbit or an
equivalent action on the stalk.

## 10. Historical chart reconciliation and preservation boundary

- Full `A+B+C-D-E-F+G` is literally the `B3` face cover in Sections 1–2.
- `A+B-C` is the scalar `B2` union chart.  THM-796 categorifies its two endpoint
  faces only after retaining the uncovered apex bit.  THM-550's even half-
  tiling recurrence is a different mirror-folded carrier, not this same set.
- In the corrected odd half chart the generators are `(A,D,B)`, with `C` an
  overlap and `D` a corner of equal scalar size.  The colour-event assignment
  `A=low`, `B=high`, `D=upper` realizes its `++-+--+` sign word through (6.4),
  but is a chart isomorphism, not an identification of the geometric slots.

The assumption challenged here is that recursive tournament vertices must be
literal vertices or arcs.  Useful vertices include interval faces, gap-
contracted roots, complement-line phases, colour events, boundary cyclic
triples, Smith currents, and proof obligations.  The gap face preserves a
lower tournament tiling and exact tile-overlap ownership; it destroys induced-
subtournament ancestry.  None of these quotients preserves the LRC loneliness
predicate, whose metric scale, owner, threshold, carry, and continued-fraction
transport still require a separate stalk.

Concurrent THM-805 computes the unmarked staircase Smith network's Tutte
polynomial as `prod_(k=1)^(n-2)(x-1+[k]_y)`.  That is the natural scalar base
for a marked deletion-contraction refinement by the `A/B/C` ownership slots
and `q` curvature.  It does not identify the cell-level series-parallel Smith
network with the class-level wiggly network: the latter's exact `n=7`
potential has adjacent-level concordance failures, and the incoming dual-side
audit refutes simple resistance reciprocity.  Equations (8.3)–(8.4) are local
face identities and do not require either global concordance claim.

Tournament Analysis uses the information carriers

```text
node_pair_colour -> B3_address -> B2_address -> B2_B3_join
 -> C3_boundary_pair -> Xi -> Omega -> Omega_B2_join -> exact_line
```

as vertices.  The pairwise observable is the number of unordered literal-line
pairs separated by a carrier; retention and retention per log-cell are the
two gauges.  At `n=7` both tournaments are transitive (score histogram
`{0:1,...,8:1}`, zero directed triangles, nine singleton SCCs, one Hamiltonian
path) and 35 edges flip between gauges.  This is telemetry about information
economy, not an importance ordering.
