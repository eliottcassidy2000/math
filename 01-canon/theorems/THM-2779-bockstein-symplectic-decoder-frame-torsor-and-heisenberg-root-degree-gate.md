---
id: THM-2779
title: "Bockstein--symplectic decoder/frame torsor and Heisenberg root-degree gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The THM-2771
  augmentation decoders and the THM-2772 normalized transverse frames are
  exact affine F13-torsors with the same gauge-invariant -aN7 transgression.
  Their determinant face is the central phase of H_13: it survives on the
  thirteen-root coefficient space, dies after forgetting phases to Boolean
  root permutations, and first admits a faithful permutation carrier in
  degree 169, exactly the endpoint-origin fibre size of THM-2772.  All
  62,377,224 normalized endpoint squares pass the exact coefficient gate.
  After a frame and two-digit section are chosen, this action is exactly the
  carry-suppressed H_13 action on Z/169, whose central step is +13 and whose
  +1 odometer differs at the carry wall.  A flat dilation-natural
  target-to-root map must vanish; a graded digit lift is the sharp survivor.
  No physical
  same-ancestry lift, semantic edge, row exclusion, or LRC(14) conclusion is
  proved.
source: root/bockstein-symplectic-twin-transgression-2026-07-28
audit: >
  root/audit-2779-2026-07-28 independently checked the decoder/frame torsors,
  Heisenberg and two-digit odometer actions, permutation-degree gate,
  dilation boundary, THM-2625 import signs and validity gates, both 28,392
  endpoint-edge banks per field, the 62,377,224-square reduction, witnesses,
  Pluecker identities, characteristic-zero inference, all four LF hashes,
  and both companions under normal/-O/stored replay: ACCEPT.
depends_on:
  - THM-2771-joint-c7-c13-right-wing-mixed-spectrum-and-commuting-square-no-go
  - THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction
  - THM-2697-filtered-affine-handoff-germ-category-and-base-signature-holotopy-boundary
  - THM-2620-endpoint-pair-parabolic-transvection-and-translation-gauge-boundary
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
related:
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2695-secondary-kummer-bockstein-picard-divisibility-spectrum-and-prime-alignment-boundary
  - THM-2763-carrier-equivariant-endpoint-address-extension-and-gauge-obstruction
script: 04-computation/lrc14_bockstein_symplectic_heisenberg_gate_thm2779.py
output: 05-knowledge/results/lrc14_bockstein_symplectic_heisenberg_gate_thm2779.out
script_sha256: 4c6a58c80ddd4be0fd9bdd297b310df054bbc08996eb223d519d3cce6b8ed13a
output_sha256: f7c96259777a3ab4a5e46cac8666181ae77a3be2e440cee8785997507706791a
secondary_script: 04-computation/lrc14_bockstein_symplectic_endpoint_gate_thm2779.py
secondary_output: 05-knowledge/results/lrc14_bockstein_symplectic_endpoint_gate_thm2779.out
secondary_script_sha256: 004e06c617f9305e2f0bc30871926e3faa7843f47dcf63af1fd8a892e63101e4
secondary_output_sha256: a89c00a3830ee9ff282cc5e4557d41293af5d6f0e7feabd5d3c7e7808591e754
hash_basis: LF-normalized bytes
---

# THM-2779 -- the twin transgressions have the same gauge space, but not yet the same carrier

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2771 and THM-2772 independently reached the same seven-chart correction:

```text
right-wing coefficient Bockstein     -> -a on every chart;
endpoint determinant Mobius face     -> -a on every chart.             (1)
```

This theorem determines why that agreement is rigid, and why it is still
not a physical closure.  The two constructions have isomorphic
thirteen-point spaces of null choices.  Their common invariant is the
central commutator of the finite Heisenberg group.  A thirteen-root
coefficient vector can retain that commutator as a phase, but a Boolean
permutation of thirteen roots cannot.  The smallest faithful permutation
carrier has `13^2=169` points, precisely the endpoint-origin fibre already
present in THM-2772's universal pullback.

The result converts two vague missing maps into two finite gates:

```text
coefficient gate: choose one point of an F_13 decoder/frame torsor;
physical gate:    lift the central phase to a 169-point same-ancestry
                  permutation carrier, or retain it as a lawful phase;
scale gate:       keep the graded root digit when a dilation edge occurs. (2)
```

The first gate is solved algebraically below.  The last two remain physical.

## 1. The complete augmentation-decoder torsor

Let

```text
k=F_13,
A=k[C_13]=k[u]/(u^13-1)=k[epsilon]/(epsilon^13),
epsilon=u-1,
N_13=1+u+...+u^12=epsilon^12.                              (3)
```

Use THM-2771's intrinsic Bockstein target collapse

```text
S_beta=(0,0,8,0,0,0,0,0,9,9,9,9,8) in the u basis.       (4)
```

It has a simple augmentation zero:

```text
S_beta=epsilon V_beta,                    V_beta a unit.   (5)
```

The printed decoder is

```text
K_beta=(3,5,5,5,7,12,2,8,2,9,8,11,0),                    (6)

S_beta K_beta=epsilon.                                    (7)
```

The complete solution set of `(7)` is

```text
Dec(S_beta)
 ={K_beta+cN_13:c in F_13}.                               (8)
```

Indeed, if `K` and `K'` solve `(7)`, then

```text
epsilon V_beta(K-K')=0.
```

Since `V_beta` is a unit,

```text
K-K' in Ann(epsilon)=k epsilon^12=kN_13.                  (9)
```

Conversely `epsilon N_13=0`, so every point in `(8)` is a solution.
Thus `Dec(S_beta)` is an affine `F_13`-torsor, not a unique inverse.
The ambiguity is exactly the target norm mode.

This gauge does not change the final chart invoice.  If `B_e` denotes one
physical chart row before target decoding, replacing `K_beta` by
`K_beta+cN_13` adds

```text
c aug(B_e)N_13                                             (10)
```

to that row.  Summing its target-zero column over charts adds

```text
c aug(S_beta)=0.                                           (11)
```

Hence THM-2771's decoded column may move locally, but its uniform chart
convolution remains

```text
C*N_7=-N_7,                 N_7=1+x+...+x^6.              (12)
```

The decoder gauge loses provenance and cannot change the holotopy invoice.

## 2. The complete normalized-frame torsor

Let

```text
V=F_13^2,                         omega(v,w)=det(v,w).      (13)
```

Fix a nonzero target direction `t`.  The normalized transverse frames are

```text
Fr(t)={s in V:omega(s,t)=1}.                               (14)
```

Choose one `s_0` in `Fr(t)`.  Nondegeneracy of `omega` gives

```text
Fr(t)={s_0+ct:c in F_13}.                                  (15)
```

Thus `Fr(t)` is another affine `F_13`-torsor.  There are thirteen frames
over each of the `168` nonzero `t`, hence

```text
168*13=2,184                                               (16)
```

normalized ordered frames, agreeing with THM-2772.

For any endpoint origin `(L,R)` and any `s in Fr(t)`, put

```text
F_(epsilon,eta)=omega(L+epsilon s,R+eta t),
epsilon,eta in {0,1}.                                      (17)
```

Bilinearity gives

```text
F_00-F_10-F_01+F_11=omega(s,t)=1.                          (18)
```

The result is independent of `(L,R)` and of the frame gauge `s->s+ct`.
For a fixed marker `a!=0`, the seven identical faces therefore give

```text
c_j=-a,                    sum_(j in F_7)c_j=-7a.           (19)
```

This is exactly the THM-2591 invoice.

## 3. The coefficient twin-transgression theorem

After choosing basepoints `K_beta` and `s_0`, equations `(8)` and `(15)`
admit the equivariant bijection

```text
K_beta+cN_13  <->  s_0+ct.                                 (20)
```

Both sides are translated freely by `c in F_13`.  Both gauges disappear
under their boundary maps:

```text
S_beta(K_beta+cN_13)=epsilon,
omega(s_0+ct,t)=1.                                         (21)
```

Multiplication by `-a` and uniform extension over the seven chart edges
therefore sends both sides to the same point

```text
-aN_7.                                                     (22)
```

So the THM-2771 and THM-2772 transgressions are not merely numerically
equal.  They have isomorphic homotopy fibres over their common normalized
boundary.  Their ambiguity is one norm/longitudinal gauge in each case.

The basepoint choice in `(20)` is load-bearing.  THM-2771 supplies the
coefficient basepoint `(6)`.  An oriented target basis supplies `s_0`, but
canon has not transported that basis to the physical root deck on the same
Boolean ancestry.  THM-2695 gives the general warning in another coefficient
sequence: equal Bockstein values without a naturality square do not identify
the underlying lifts.  Equation `(20)` is an algebraic torsor
identification, not that missing physical square.

## 4. The determinant face is a Heisenberg central phase

For any odd prime `p`, define the finite Heisenberg group

```text
H_p=F_p^3,

(x,y,z)(x',y',z')
 =(x+x',y+y',z+z'-yx').                                  (23)
```

With

```text
X=(1,0,0),             Y=(0,1,0),             Z=(0,0,1),
```

one has

```text
[X,Y]=Z,                 |H_p|=p^3.                       (24)
```

More generally, the central exponent in a commutator is

```text
xy'-x'y=det((x,y),(x',y')).                               (25)
```

Thus THM-2772's determinant Mobius face is the exponent of the Heisenberg
commutator.

This commutator already acts on a `p`-dimensional **coefficient** space.
On basis vectors indexed by `r in F_p`, let

```text
T e_r=e_(r+1),                 M e_r=zeta_p^(-r)e_r.       (26)
```

Then

```text
TM=zeta_p MT.                                              (27)
```

The central element is visible as a phase.  After forgetting phases to the
underlying Boolean root permutation, however, `M` fixes every root and the
commutator in `(27)` disappears.  This is the exact distinction between the
coefficient success in `(22)` and the missing physical root-deck lift.

## 5. Thirteen roots are impossible; 169 endpoint origins are minimal

No faithful permutation action of `H_p` has degree less than `p^2`.

To prove the lower bound, let `H_p` act on a set of size `<p^2`.  Every orbit
has `p`-power size, hence size `1` or `p`.  On a `p`-point orbit the image is
a transitive `p`-subgroup of `S_p`.  Since

```text
v_p(p!)=1,
```

that image has order `p`, is cyclic, and kills `[H_p,H_p]=<Z>`.  Fixed
orbits kill it as well.  Therefore `Z` acts trivially on the whole set, so
the action is not faithful.

The bound is sharp.  On `F_p^2` define

```text
rho(x,y,z)(r,w)=(r+x,w+z-yr).                             (28)
```

Equation `(23)` is exactly the composition law for `(28)`.  The action is
transitive, and an element fixing every `(r,w)` has successively
`x=z=y=0`; hence it is faithful.  Therefore

```text
mu_perm(H_p)=p^2,                  mu_perm(H_13)=169.      (29)
```

THM-2772's pullback has exactly `169` endpoint origins over every faithful
carrier address.  After choosing oriented coordinates, `(28)` gives an
abstract faithful Heisenberg action on that origin fibre.  This cardinal
match is now structural rather than accidental: `169` is the first
permutation degree capable of retaining the determinant commutator.

There is a sharper intrinsic description using THM-2620.  Fix a normalized
frame `(s,t)` and the target sector `q=s`.  Write an endpoint origin as

```text
R=ws+vt,                         Delta=det(s,R)=v.          (30)
```

On these `169` origins define

```text
X=T_(s,1): R |-> R+det(s,R)s,
Y=tau_t:    R |-> R+t,
Z=tau_s:    R |-> R+s.                                    (31)
```

In `(v,w)` coordinates,

```text
X(v,w)=(v,w+v),       Y(v,w)=(v+1,w),       Z(v,w)=(v,w+1),

[X,Y]=Z.                                                   (32)
```

Thus the THM-2620 transvection and target translation generate exactly a
faithful `H_13` action on the endpoint-origin fibre.  Each fixed determinant
sector is one central `Z`-cycle.  At `Delta=1`,

```text
L=R+q=R+s=Z(R),                                           (33)
```

so the central edge is literally the ordered endpoint edge `R->L`.
THM-2625 gives nonzero canonical current at all `169` vertices for every
nonzero `s`.

This is an algebraic target-dipole-to-**endpoint-deck** map, not yet the
THM-2542 physical root deck.  THM-2625's weights are nonzero but are not
proved invariant, covariant, or even phase-compatible under `(31)`.  The
frame gauge `s->s+ct` also changes the actual central direction `Z` while
leaving the scalar determinant invoice unchanged.  An equivariant
decoder/frame torsor identification on one Boolean ancestry remains the
missing sidecar.

### 5.1 Every normalized square is occupied coefficientwise on the canonical control

The THM-2625 exact control is stronger than vertexwise support.  Let
`P_L,Q_R` be its separate left/right endpoint factors.  In each of the two
certified finite-field specializations, exhaustive reconstruction gives,
for every nonzero `d in V` and every `L,R in V`,

```text
P_L-P_(L+d) !=0,             P_L+P_(L+d) !=0,
Q_R-Q_(R+d) !=0,             Q_R+Q_(R+d) !=0.             (33a)
```

There are `168*169=28,392` source edges and the same number of target edges
per field.  Consequently, for every one of the `2,184` normalized frames
and every pair of endpoint origins, all four entries of

```text
v=(P_0Q_0,P_1Q_0,P_0Q_1,P_1Q_1)                           (33b)
```

and all four Hadamard coordinates

```text
D_0=(P_0+P_1)(Q_0+Q_1),
D_1=(P_0+P_1)(Q_0-Q_1),
D_2=(P_0-P_1)(Q_0+Q_1),
D_3=(P_0-P_1)(Q_0-Q_1)                                   (33c)
```

are nonzero, while

```text
D_0D_3=D_1D_2.                                            (33d)
```

Thus all

```text
2,184*169*169=62,377,224
```

normalized endpoint squares pass the full coefficient-side
Segre--Pluecker gate in both exact fields.  A nonzero image proves the
corresponding characteristic-zero cyclotomic coefficient is nonzero.

The gain is exact common-array structure: the four corners belong to one
marked THM-2625 triangle, not four unrelated sector witnesses.  They are
still four **endpoint coefficients**, not the four absent/present Boolean
carrier-allocation states on one physical support.  Nor does `(33a)` prove
that the current is equivariant under the Heisenberg action `(31)`.

### 5.2 The 169-point carrier is the two-digit odometer model

The degree `169` coincidence has an exact arithmetic explanation.  Put

```text
Omega_2=Z/169Z
```

and choose the standard digit section

```text
iota(v,w)=v+13w,                0<=v,w<13.                (33e)
```

Conjugating the endpoint operations `(31)` by `iota` gives

```text
X(n)=14n mod169,
Y(v+13w)=(v+1 mod13)+13w,
Z(n)=n+13 mod169.                                      (33f)
```

Indeed,

```text
14(v+13w)=v+13(w+v) mod169.                            (33g)
```

Thus `(33f)` is a faithful permutation action of `H_13` on the two-digit
cyclic address set, and

```text
[X,Y]=Z.                                                (33h)
```

The point is not that the odometer extension splits.  If

```text
A(n)=n+1 mod169
```

is the actual cyclic successor, then

```text
A=Y                 when v!=12,
A=Z o Y             when v=12.                         (33i)
```

There are `156` nonwrapping addresses and `13` carry addresses.  The central
step `Z=+13` is precisely the correction paid by the chosen low-digit
section at the carry wall.

This is, up to the inherited unit normalization of the central/kernel
generator, the two-digit pushout of THM-2657's length-six nonsplit odometer
extension.  THM-2657 records class `7` in its fixed root normalization,
whereas `(33i)` records class `1`; multiplication by `7` identifies them.
It identifies the **abstract form** of the endpoint Heisenberg
carrier with a root-depth address model after choosing both a normalized
frame and the digit section `(33e)`.  It does not identify THM-2625 endpoint
origins with physical LRC cylinders on one ancestry.  In particular, one
may not call `+13` a physical arm/root move in the marked current until that
intertwiner is supplied.  Equation `(33i)` also shows exactly why replacing
the carry-suppressed `Y` by physical `+1` would be false.

## 6. Dilation forces a graded root map

THM-2697 identifies one further obstruction.  Suppose a proposed
target-role-to-root map treats a target label as unchanged by a dilation
edge, while collapsing the physical root filtration to one `F_13` root
torsor.  Dilation induces

```text
D_target=identity,                    D_root=0.             (34)
```

Any linear natural intertwiner `pi` would satisfy

```text
D_root pi=pi D_target,
0=pi.                                                        (35)
```

Thus a nonzero flat, dilation-natural `pi` is impossible under these
hypotheses.

The sharp survivor is THM-2697's graded digit chain

```text
F^0/F^1 -> F^1/F^2 -> ... -> F^5/F^6.                     (36)
```

Between consecutive graded copies, dilation transports the digit rather
than acting by zero.  A scalar choice at one depth propagates through all
six grades, giving thirteen graded intertwiners, twelve nonzero.  Therefore
a future physical lift must do at least one of the following:

1. retain the root-depth object `(36)`;
2. give the target role its own nontrivial dilation action;
3. work in the inverse-branch affine germ groupoid rather than the collapsed
   forward root torsor.

Equation `(35)` is conditional on the fixed-target action in `(34)`; it is a
typing gate, not a no-go against every possible enriched carrier.

## 7. Exact source, target, loss, and cheapest next tests

The proved candidate mechanism is

```text
source:
  THM-2771 intrinsic coefficient Bockstein;
  THM-2772 determinant square and endpoint-origin fibre;

map:
  decoder/frame torsor equivalence (20);
  Weyl/Heisenberg central-phase identification (27);

target:
  the common coefficient correction -aN_7;

preserved:
  target orientation, normalized determinant, seven-chart invoice,
  all decoder/frame gauge classes, the endpoint Heisenberg action;

lost:
  positive Boolean ancestry, semantic owner edge, local chart provenance,
  a physical root-deck identification, semantic ancestry, H-equivariance of
  the current weights, root-depth provenance.                              (37)
```

Every proposed completion now has three cheap decisive tests.

1. **Phase or 169-point carrier.**  Does it lawfully retain the central
   phase in `(27)`, or exhibit a same-ancestry 169-point action conjugate to
   `(31)`?  A bare thirteen-root permutation cannot pass.
2. **Naturality.**  Does its target/root map commute with every dilation it
   crosses?  If it uses `(34)`, `(35)` forces it to be zero.
3. **Semantic attachment.**  Does one nonzero right-wing current and one
   THM-2542 chart edge occupy the same lifted carrier?  Full coefficient
   support or an abstract action on the fibre is not enough.

Passing all three would turn `(22)` into the physical mixed square requested
by THM-2591.  This theorem passes only the coefficient test.

## 8. Exact companion and scope

Run

```bash
python 04-computation/lrc14_bockstein_symplectic_heisenberg_gate_thm2779.py
python -O 04-computation/lrc14_bockstein_symplectic_heisenberg_gate_thm2779.py
python 04-computation/lrc14_bockstein_symplectic_endpoint_gate_thm2779.py
python -O 04-computation/lrc14_bockstein_symplectic_endpoint_gate_thm2779.py
```

Each normal/optimized pair byte-matches its stored transcript:

```text
05-knowledge/results/lrc14_bockstein_symplectic_heisenberg_gate_thm2779.out
05-knowledge/results/lrc14_bockstein_symplectic_endpoint_gate_thm2779.out
```

The primary no-`assert`, integer-only companion verifies:

- convolution rank twelve and kernel `<N_13>`;
- all thirteen decoders in `(8)`;
- all `2,184` normalized frames and `371,293` origin/frame Mobius faces;
- the order-`2,197` Heisenberg law, commutator, faithful 169-point action,
  and all `4,826,809` action-pair compositions;
- the `169` endpoint transvection/translation commutators, all thirteen
  determinant-sector central cycles, and the thirteen `Delta=1` endpoint
  edges;
- the `169` two-digit odometer conjugacies, including the `156/13`
  nonwrap/carry split and the central `+13` correction;
- the thirteen Weyl phase checks and Boolean centre loss;
- the flat-zero versus graded-thirteen dilation intertwiner boundary.

The secondary exact companion imports and reruns THM-2625's dual
Lucas/order validity gates, reconstructs both endpoint-factor transforms,
checks all `28,392` nonzero source edges and all `28,392` nonzero target
edges in each field, and thereby certifies the `62,377,224` normalized
endpoint squares in Section 5.1.

An independent hostile audit reconstructed both affine torsors, the
Heisenberg commutator and minimum permutation degree, the two-digit
odometer conjugacy and carry wall, and the dilation-natural zero gate.  It
also rederived the secondary endpoint DFT signs and factorization, checked
the two exact edge banks and universal-square reduction, replayed both
companions under normal and optimized Python against their stored outputs,
and verified all four LF-normalized hashes.

No physical same-ancestry endpoint lift, root-deck map, semantic transition,
positive current, scalar-row exclusion, or LRC(14) conclusion is proved.

QED.
