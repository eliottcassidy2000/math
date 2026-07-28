---
id: THM-2772
title: "Carrier-allocation pullback, central K4 fibre, and mixed-face obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  THM-2763's
  extended carrier address and THM-2625's endpoint plane
  have the same cardinality but forget different coordinates; their
  canonical pullback over the old target quotient has size 13^6.  The
  Boolean bare/source/target/both allocation lift has fibres 4,2,1 according
  to zero harmonics, with a central abstract K4.  After an ordered endpoint
  frame, the 2,185 admissible determinant sectors are exactly zero plus all
  nonconstant affine F13 fields on that K4, and all 2,184 nonzero-target
  sectors have formal THM-2757 charge but zero Boolean mixed curvature.  On
  a coherent endpoint-translation parallelogram, however, the determinant
  mixed face is exactly det(s,t); normalized oriented frames give c_j=-a
  and pay the seven-clock Cech invoice algebraically.  Existing canon does
  not construct one common current over those four varying target sectors or
  identify carrier allocation with the endpoint translations.  On
  one genuinely common endpoint atom the third Segre-Hadamard coordinate is
  the required mixed face.  No such physical common atom, root-deck map,
  physical Cech correction, row exclusion, or LRC(14) conclusion is proved.
source: root/carrier-allocation-holotopy-pullback-2026-07-28
audit: >
  arm-wing-contingency-2026-07-28 (independent pullback/type audit, complete
  affine/Segre/Pluecker and virtual-face census replay, symplectic endpoint
  parallelogram audit, dependency/scope audit, normal/-O/stored/hash replay:
  ACCEPT)
depends_on:
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2620-endpoint-pair-parabolic-transvection-and-translation-gauge-boundary
  - THM-2625-canonical-endpoint-current-full-transvection-sector-survival
  - THM-2757-marked-reference-opposite-edge-clutch-transgression
  - THM-2763-carrier-equivariant-endpoint-address-extension-and-gauge-obstruction
related:
  - THM-2374-binary-allocation-complete-subcube-dirichlet-spectrum
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2634-endpoint-pair-two-carry-cospan-and-single-carry-no-go
  - THM-2756-opposite-edge-projectors-parity-cancellation-and-integral-clutch
  - THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks
  - THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile
  - THM-2770-tree-incidence-a-d-weyl-clutch-and-four-vertex-fan-dichotomy
  - THM-2771-joint-c7-c13-right-wing-mixed-spectrum-and-commuting-square-no-go
script: 04-computation/lrc14_carrier_allocation_pullback_k4_segre_thm2772.py
output: 05-knowledge/results/lrc14_carrier_allocation_pullback_k4_segre_thm2772.out
script_sha256: 10b681f575fb51eb4b1af4bc909fba89846b85bd5da36fc069dac97ae2ebe409
output_sha256: 9cd1e82634822e4997d9966a8252b566d099050e42593a47c49b8a0f387190e5
hash_basis: LF-normalized bytes
---

# THM-2772 -- allocation pullback, central `K4`, and the mixed face

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2763 replaces the false collapsed-carrier quotient by the exact address
`(r,k,l)`, retaining the total source and target carrier harmonics.  THM-2625
retains instead an allocated endpoint pair `(L,R)` and fills every admissible
determinant sector.  Both finite objects have `13^4=28,561` elements.  That
cardinality match is not a map: the first forgets endpoint origin, while the
second forgets carrier harmonics and allocation state.

This theorem computes their canonical common pullback and then adds the
Boolean carrier-allocation bits.  The zero-harmonic fibre is the first literal
abstract four-state object naturally adjacent to THM-2757:

```text
bare, source-only, target-only, both.                         (0)
```

It also finds a surprising exact dictionary.  Once an ordered endpoint frame
is supplied, THM-2625's `2,185` admissible determinant sectors are precisely
the zero affine field and every nonconstant affine `F_13` field on the four
allocation states.  Every nonzero-target sector is formally charged by
THM-2757's marked opposite-edge transform.  But every one of those affine
fields has zero mixed rectangle.  Thus determinant charge is not the
transition-dependent mixed face required by THM-2591.

There is a sharp positive exception.  If the two allocation toggles are
realized on one common endpoint atom as independent source and target
translations, the determinant parallelogram has mixed curvature `det(s,t)`.
Seven normalized oriented frames would give the exact correction `c_j=-a`
and pay THM-2591's Cech invoice.  The theorem proves this symplectic identity
and isolates the missing physical lift: the four corners vary the target
sector, while current canon supplies neither their common current nor the
translation meaning of the Boolean allocation bits.

The genuine mixed face appears only after four coefficient values are formed
on one common endpoint atom.  Its exact formula is the third coordinate of a
Segre--Hadamard factorization.  Existing canon does not prove that coefficient
survives on the carrier/endpoint pullback, so this remains a proof-complete
finite theorem under audit rather than an LRC closure.

## 1. The universal carrier/endpoint pullback

Put

```text
V=K_13/L_13 isomorphic F_13^2.                              (1)
```

Use THM-2763's distinguished `e_0` section.  Its full carrier quotient has
the typed coordinate form

```text
G_full isomorphic V direct-sum F_13,k direct-sum F_13,l,

[(r,k,l)] |-> (q,k,l),
q=[r-(l-k)e_0] in V.                                      (2)
```

Here `k,l` are total source/target carrier harmonics.  THM-2620's allocated
endpoint plane is

```text
E=V x V,
pi_E(L,R)=L-R=q.                                          (3)
```

The canonical common object is the fibre product

```text
P=G_full x_V E
 ={(q,k,l,L,R):L-R=q}
 isomorphic {(q,k,l,R)}.                                  (4)
```

Therefore

```text
|G_full|=13^4,
|E|=13^4,
|P|=13^6=4,826,809.                                      (5)
```

Every carrier address has exactly `13^2=169` endpoint origins in `P`.  The
determinant sector is now a lawful function

```text
Delta=det(L,R)=det(q,R).                                  (6)
```

The construction is universal, not merely a count.  If a set `X` has maps

```text
f:X->G_full,                g:X->E
```

whose `V` projections agree, then

```text
x |-> (f(x),g(x))
```

is the unique map `X->P` commuting with both projections.  Thus `(4)` is the
minimal common object in the pullback sense.  A chosen section may use fewer
coordinates, but that section is exactly new endpoint-origin/allocation data.

The equality of the two `13^4` cardinalities does not supply it.  In
particular, THM-2763's common-harmonic kernel has size thirteen, whereas the
common endpoint-translation fibre in `(3)` has size `169`.  They are not the
same torsor.

## 2. The Boolean allocation lift

Let

```text
epsilon_S,epsilon_T in {0,1}                              (7)
```

record whether the carrier is present at the source and target endpoints.
Absence forces zero harmonic, while presence permits zero or nonzero
harmonic.  Over the harmonic plane define

```text
A={(epsilon_S,epsilon_T,k,l):
     epsilon_S=0 implies k=0,
     epsilon_T=0 implies l=0}.                            (8)
```

Projection to `(k,l)` has the exact fibre census

```text
(k,l)=(0,0):                    size 4 on one point;
exactly one of k,l is zero:     size 2 on 24 points;
k,l are both nonzero:           size 1 on 144 points.     (9)
```

Hence `A` has `196` points.  Define the target-extended allocation lift

```text
A_tilde=V x A.                                             (9a)
```

It has

```text
169*196=33,124                                             (10)
```

points.  At `(k,l)=(0,0)`, the four states are exactly `(0)`.  With the bare
state marked, the three perfect matchings are the source toggle, target
toggle, and diagonal pairing.  This is the precise abstract `K4` required by
THM-2757 in every fixed `q` fibre.

The projection

```text
A_tilde -> G_full,       (q,epsilon_S,epsilon_T,k,l) |-> (q,k,l), (10a)
```

is not injective.  At zero harmonic it forgets whether a carrier was absent
or present with its zero Fourier mode.
Consequently even THM-2763's full address is not a factor-allocation label.
No physical nonvanishing of all four central states is inferred here; `(9)`
is the exact structural reason a separate allocation sidecar is necessary.

## 3. The determinant-sector/affine-square dictionary

Choose an ordered basis `(s,t)` of `V`.  The choice is an endpoint-allocation
frame; it is not supplied canonically by THM-2763.  For every admissible
THM-2625 sector put

```text
a=det(q,s),                       b=det(q,t),

v_00=Delta,
v_10=Delta+a,
v_01=Delta+b,
v_11=Delta+a+b.                                      (11)
```

The map

```text
q |-> (a,b)                                               (12)
```

is an isomorphism: if both determinants vanish, `q` lies in the spans of
both basis vectors and hence is zero.  Thus `(11)` is a bijection

```text
admissible (q,Delta) sectors
  <-> {zero affine field} union
      {nonconstant affine F_13 fields on {0,1}^2}.         (13)
```

Indeed the count is

```text
1+13(13^2-1)=1+2,184=2,185.
```

The twelve structurally impossible labels `q=0,Delta!=0` would be the
constant nonzero fields and are absent from both sides.

Mark `00`, order the arms as `10,01,11`, and apply THM-2757's opposite-edge
difference.  Direct substitution gives

```text
D(v)=(-2b,-2a,0).                                        (14)
```

Because two and three are units in `F_13`, the primitive `C3` projection of
`D(v)` is nonzero exactly when `q!=0`.  Hence all `2,184` nonzero-target
sectors are formally charged, independently of `Delta`; the sole neutral
profile is `(q,Delta)=(0,0)`.

However

```text
v_00-v_10-v_01+v_11=0.                                  (15)
```

The determinant label is affine on the allocation square and has zero
Boolean curvature.  Formula `(14)` is a first-order marked-reference charge,
not THM-2591's mixed-square transition correction.

THM-2625 proves that every sector in `(13)` carries nonzero canonical current.
It does not transport one current coefficient through the four allocation
states `(11)`.  Multiplying a nonzero sector current by four scalar labels is
not a literal carrier allocation.

### 3.1 Endpoint-parallelogram transgression: an exact conditional bridge

The affine fixed-`q` dictionary is not the only way the determinant can meet
an allocation square.  Suppose a future physical construction retains one
common endpoint atom `(L,R)` and identifies the source and target carrier
toggles with endpoint translations by `s,t in V`:

```text
L_epsilon=L+epsilon s,       R_eta=R+eta t,
epsilon,eta in {0,1}.                                      (15a)
```

This hypothesis is strictly stronger than anything proved by THM-2625 or
THM-2763.  On such a square put

```text
F_(epsilon,eta)=det(L+epsilon s,R+eta t).                  (15b)
```

Bilinearity gives the division-free identity

```text
F_00-F_10-F_01+F_11=det(s,t).                             (15c)
```

Thus the determinant itself supplies a canonical finite symplectic, or
Heisenberg-commutator, mixed face.  It is independent of the endpoint origin
`(L,R)`.  This is the missing term in `(15)`: that formula kept `q` fixed,
whereas the four corners of `(15a)` have

```text
q_(epsilon,eta)=q+epsilon s-eta t.                         (15d)
```

Consequently `(15c)` crosses four generally different target sectors rather
than living inside one fixed `(q,Delta)` cell.

For an ordered endpoint frame, `d=det(s,t)` is nonzero.  There are exactly

```text
2,184
```

ordered frames of each prescribed `d in F_13^*`, hence `26,208` ordered
bases in total.  In particular there are `2,184` normalized oriented frames
with `d=1`.

This produces an exact candidate for THM-2591's root-deck map.  Fix its marker
`a in F_13^*`.  If the seven clock faces admit coherently oriented normalized
frames, define

```text
c_j=-a(F_(j,00)-F_(j,10)-F_(j,01)+F_(j,11))=-a.           (15e)
```

Then, more strongly than the cyclic condition alone,

```text
a+c_j=0 for every j,
sum_(j in F_7)c_j=-7a.                                    (15f)
```

So the finite cochain produced by `(15e)` pays the exact Cech invoice `(27)`.
More generally, frames of areas `d_j` give

```text
sum_j c_j=-a sum_j d_j,                                   (15g)
```

and flatten the seven-cycle exactly when `sum_j d_j=7`; a constant area
`d` gives only the scaled invoice `-7ad` and requires `d=1` or an additional
lawful normalization.

The source, map, and target of this conditional mechanism are now precise:

```text
source:   one common current on the four endpoint corners (15a),
          including their four q-values and determinant sectors;
map:      the determinant Mobius boundary (15c), followed by -a;
target:   the transition cochain c in (15e), with invoice (15f);
preserve: clock, semantic word, common Abel boundary, endpoint order,
          carrier gauge, and the identification of allocation toggles
          with the two endpoint translations.                            (15h)
```

Every clause in `(15h)` is load-bearing.  If `s,t` are dependent, `(15c)` is
zero.  If they are independent but have unnormalized area `d`, the correction
is scaled as in `(15g)`.  Most importantly, separate nonzero currents in the
four THM-2625 sectors do not constitute one common current on `(15a)`, and
THM-2763's scalar carrier harmonics `(k,l)` are not yet identified with the
endpoint translation vectors `(s,t)`.  Hence `(15c)--(15f)` give an exact
positive algebraic bridge and its cheapest decisive physical test, not a
constructed LRC carrier or a row exclusion.

## 4. The coefficient-level Segre--Hadamard identity

Suppose a stronger future construction retains one common endpoint atom.
Write its source amplitudes as `P_0,P_1` for bare/carried and its target
amplitudes as `Q_0,Q_1`.  The four allocated current coefficients are

```text
v=(P_0Q_0,P_1Q_0,P_0Q_1,P_1Q_1).                       (16)
```

They satisfy the rank-one Segre equation

```text
v_00v_11-v_10v_01=0.                                    (17)
```

Over every commutative ring, THM-2757's marked difference factors exactly as

```text
D_1=(P_0+P_1)(Q_0-Q_1),
D_2=(P_0-P_1)(Q_0+Q_1),
D_3=(P_0-P_1)(Q_0-Q_1).                                 (18)
```

The third coordinate is precisely the Boolean Mobius/ANOVA mixed face:

```text
mu=D_3=v_00-v_10-v_01+v_11.                             (19)
```

Retain also the invariant Hadamard coordinate

```text
D_0=v_00+v_10+v_01+v_11=(P_0+P_1)(Q_0+Q_1).            (19a)
```

The four coordinates satisfy the sharp Pluecker/Segre relation

```text
D_0D_3=D_1D_2.                                          (19b)
```

Thus one first-order charge is insufficient, but the three-coordinate gate

```text
D_0!=0,                   D_1!=0,                  D_2!=0 (19c)
```

forces the mixed face `D_3!=0` over `F_13`.  This is the cleanest positive
conditional bridge in the theorem: a common atom with nonzero invariant,
source contrast, and target contrast automatically pays the Boolean mixed
invoice before any Cech/root-deck map is applied.  Exact exhaustion gives

```text
20,736 / 28,561 quadruples satisfying (19c),
14,400 / 20,736 when all P_0,P_1,Q_0,Q_1 are nonzero.    (19d)
```

This identifies the exact coefficient that a physical four-state carrier
must preserve.  Merely observing primitive marked charge is insufficient.
The sharp all-endpoint-factors-nonzero hostile over `F_13` is

```text
(P_0,P_1,Q_0,Q_1)=(1,1,1,2),

v=(1,1,2,2),             D=(11,0,0),             mu=0.  (20)
```

A two-sided mixed positive control is

```text
(P_0,P_1,Q_0,Q_1)=(1,2,1,3),

v=(1,2,3,6),             D=(7,9,2),              mu=2.  (21)
```

Exhaustion of all `13^4=28,561` endpoint quadruples gives

| marked profile | mixed face | count |
|---|---:|---:|
| charged | nonzero | `24,192` |
| charged | zero | `3,744` |
| flat | nonzero | `144` |
| flat | zero | `481` |

The `144` flat-but-mixed cases are the sharp boundary

```text
P_1=Q_1=0,                  P_0Q_0!=0,
```

for which `D` is a nonzero constant vector.  Restricting to all four endpoint
amplitudes nonzero gives

| marked profile | mixed face | count |
|---|---:|---:|
| charged | nonzero | `17,424` |
| charged | zero | `3,168` |
| flat | zero | `144` |

Thus on the all-nonzero locus `mu!=0` implies marked charge, but marked charge
still does not imply `mu!=0` by `(20)`.

The sharp flatness criterion behind the census is

```text
D_1=D_2=D_3
 iff
 P_1(Q_0-Q_1)=0 and Q_1(P_0-P_1)=0.                     (22)
```

It follows from `D_1-D_3=2P_1(Q_0-Q_1)` and
`D_2-D_3=2Q_1(P_0-P_1)`.

### 4.1 Virtual mixed spectrum versus four-corner co-support

There is a second exact boundary after Fourier allocation.  Let `a,b` be
independent source/target carrier twists in `F_13`.  Write `p(a),q(b)` for
the carried source/target endpoint arrays and `p_bare,q_bare` for the two
no-carrier endpoint values.  Define the four coefficient-side allocation
arrays

```text
B(a,b)=p_bare q_bare,
P(a,b)=p(a) q_bare,
Q(a,b)=p_bare q(b),
H(a,b)=p(a)q(b).                                         (22a)
```

These are raw endpoint-factor products.  They do not assert four physical
events on one common support.  With unnormalized independent DFT and hats
for the one-dimensional transforms,

```text
Bhat(k,l)=169 p_bare q_bare 1_((k,l)=(0,0)),
Phat(k,l)=13 q_bare phat(k) 1_(l=0),
Qhat(k,l)=13 p_bare qhat(l) 1_(k=0),
Hhat(k,l)=phat(k)qhat(l).                                (22b)
```

Thus the fourfold primal support is contained in the neutral corner
`(0,0)`.  It is exactly that singleton only when the four central factors
`p_bare,q_bare,phat(0),qhat(0)` are all nonzero.

The additive allocation face is nevertheless

```text
Omega=B-P-Q+H
     =(p_bare-p(a))(q_bare-q(b)).                         (22c)
```

For every mixed character `k,l!=0`, equations `(22b)` give

```text
Omegahat(k,l)=Hhat(k,l)=phat(k)qhat(l),

Bhat(k,l)=Phat(k,l)=Qhat(k,l)=0.                          (22d)
```

THM-2763's exact diagonal-fibre certificate gives full thirteen-element DFT
support for each of its source and target carrier arrays.  Hence the
coefficient-side tensor in `(22a)` has

```text
Omegahat(k,l)!=0             for all 12^2=144 mixed pairs. (22e)
```

At every one of those addresses, however, only the `both` corner is present
among the four transforms.  The virtual face equals that one surviving
corner; it is not a common-co-support `K4`.  This is stronger than saying the
mixed coefficient might cancel: full mixed spectral survival can occur while
three of the four allocated vertices are absent at each mixed address.

There is a sharp warning for the Pluecker gate `(19c)`.  At a mixed address
the transformed four-vector is

```text
(Bhat,Phat,Qhat,Hhat)=(0,0,0,h),                  h!=0.
```

Its Hadamard coordinates are

```text
(D_0,D_1,D_2,D_3)=(h,-h,-h,h),                           (22f)
```

so it passes `(19c)` and `D_0D_3=D_1D_2` perfectly.  Therefore the Pluecker
gate is a positive bridge only when imposed on a **common atom before Fourier
allocation**.  Applied after transforming the four corners separately, it
also accepts the one-corner virtual-face hostile.

The raw THM-2763 twist support `{0,2,10}` for each endpoint array is support
in the **dual evaluation variable**.  Its Cartesian square is not a primal
four-corner atom.  Also, THM-2763 certifies `q_bare!=0` and both carried DFTs
at zero, but does not track the corresponding no-source-carrier `p_bare` for
this four-state packet.  Therefore this theorem does not assert that the
neutral corner is physically occupied by all four states.

Equation `(22d)` is the sharp THM-2771-to-THM-2772 stopping mechanism:
mixed Fourier support is a lawful virtual Mobius face, while a literal marked
carrier requires common vertex co-support before Fourier projection.  The
current THM-2771 file is only a reserved candidate and is not a dependency;
the elementary transform identities `(22a)--(22d)` and THM-2763's proved
one-dimensional support are sufficient here.

## 5. The exact source/target/map invoice

This theorem makes the next positive bridge precise.

### Source

A seven-clock coefficient array on the allocation pullback

```text
A_tilde x_(G_full) P                                       (23)
```

must retain on one common physical endpoint atom:

```text
clock and semantic word/owner;
the two allocation bits;
the exact carrier address (r,k,l);
the endpoint pair (L,R), q=L-R, and Delta=det(q,R);
the marked triangle and one common Abel boundary.          (24)
```

### Map

Before determinant/Radon summation, take the allocation Mobius boundary

```text
mu_j=(P_(j,1)-P_(j,0))(Q_(j,1)-Q_(j,0)).                 (25)
```

Then a separately supplied carrier-equivariant functional must give

```text
pi_j:mu_j |-> c_j in F_13.                               (26)
```

The endpoint-parallelogram functional `(15c)--(15e)` is now an exact algebraic
candidate for `(26)`: on a normalized oriented translation square it sends
the determinant mixed face to `-a`.  Neither THM-2625 nor THM-2763 constructs
the common coefficient square `(25)`, identifies its allocation toggles with
those endpoint translations, or proves that this candidate acts on the
physical root deck.

### Target and exact Cech test

Let `g_j=a` be THM-2542's seven chart transitions.  The corrected cochain
`g+c` is a coboundary exactly when its cyclic sum vanishes.  Therefore a
root-trivializing mixed face must pay

```text
sum_(j in F_7)c_j=-7a.                                   (27)
```

Necessity is THM-2591's invoice.  Sufficiency at the finite cochain level is
also immediate: on a cycle, a one-cochain has zero sum iff successive partial
sums define a vertex potential.  This does not make the `c_j` physical.

### Preserved and lost data

The proposed map must preserve the common endpoint atom, factor allocation,
endpoint order, carrier representative gauge, target difference, determinant
sector, clock, and semantic word.  It may reselect an exact triangle only if
those identifications survive.  It may not first sum over endpoint origins,
forget the zero-harmonic allocation bits, replace the common atom by separate
marginals, or substitute the first two coordinates of `(18)` for `mu`.

## 6. Cheap decisive tests

Every incoming physical construction can now be audited in this order.

1. **Allocation fibre.**  At `k=l=0`, does it distinguish all four states in
   `(0)`?  If not, it is not a marked carrier `K4`.
2. **Endpoint origin.**  Is one `R` retained before determinant/Radon
   pushforward?  If not, it does not live on `(23)`.
3. **Common atom.**  Do the four coefficients arise from the same two source
   and two target amplitudes as in `(16)`?  Separate marginals are blocked by
   THM-2634.
4. **Mixed face.**  Does `(19)` survive in one sector and clock?  THM-2757
   charge alone fails by `(20)`.
5. **Co-support.**  Is the mixed face represented by all four allocation
   vertices on one address, or only by the virtual one-corner mechanism
   `(22d)`?
6. **Cech invoice.**  After the root-deck map, does `(27)` hold?  Zero cyclic
   sum is another vertex gauge.
7. **Carrier covariance.**  Are the carrier and all endpoint factors
   transported together under THM-2763's representative gauge?  Its
   fixed-carrier hostile rules out a shortcut.
8. **Symplectic normalization.**  If the allocation square is realized by
   endpoint translations, are its vectors independent and normalized so that
   `det(s,t)=1`, or do their seven areas at least sum to `7`?  Otherwise
   `(15g)` pays the wrong Cech class.

## 7. Cross-frontier comparisons, not dependencies

Three nearby results clarify the scope but are not used in the proof.

- THM-2768 proves that the free product `C2*C3` becomes the finite `A4/S4`
  frame only after adding a load-bearing face relation.  Here the Boolean
  allocation vertices and marked `C3` charge likewise do not create the
  mixed face `mu` for free.
- THM-2769 proves that full global `S4` monodromy and a square product do not
  remove local even-parity ramification.  Here full address and full formal
  sector charge likewise do not force local mixed curvature.
- THM-2770 gives an
  invertible `A3=D3` incidence clutch but not a signed coefficient selector.
  Here `(14)` is an invertible frame-dependent label transform, not a selector
  for one physical coefficient on `(23)`.

These are mechanism comparisons only.  No JC, graceful-tree, or additional
monodromy statement enters this theorem.

## 8. Exact evidence and stopping boundary

Run

```text
python 04-computation/lrc14_carrier_allocation_pullback_k4_segre_thm2772.py
python -O 04-computation/lrc14_carrier_allocation_pullback_k4_segre_thm2772.py
```

Both modes byte-match the stored transcript.  The dependency-free companion
uses explicit exception gates and checks:

1. all `169` harmonic pairs and the complete allocation-fibre census;
2. all `2,185` admissible sector profiles, their inverse dictionary, marked
   differences, and zero mixed rectangles;
3. every one of the `28,561` endpoint-origin determinant cells;
4. every increment pair at one endpoint origin and every endpoint origin for
   one normalized frame in `(15c)`, together with the exact ordered-frame
   area census `2,184` for each nonzero area and the seven-clock invoice;
5. every one of the `13^4` Segre endpoint quadruples, the complete and
   all-nonzero censuses, Pluecker identity and gate census, flatness
   criterion, hostile `(20)`, and control `(21)`; and
6. an exact `F_53` two-dimensional Fourier control with corner support sizes
   `1,13,13,169`, singleton fourfold intersection, and all `144` virtual
   mixed addresses, including the `144/144` one-corner Pluecker hostile; and
7. the `13^4`, `13^4`, `13^6`, and `169` pullback cardinalities.

The universal property in Section 1 and the Cech iff in Section 5 are proved
directly and do not rely on enumeration.  Normal, optimized, and stored
outputs match after LF normalization, and the frontmatter pins their exact
hashes.

The exact proved boundary is

```text
one supplied endpoint-allocation frame
  -> a bijection from all admissible THM-2625 sectors
     to affine marked-K4 labels
  -> formal THM-2757 charge on every q!=0 sector;

but every determinant label has zero mixed curvature,
and G_full forgets both zero-harmonic allocation state and endpoint origin.
                                                                    (28)
```

The qualified positive exception is `(15c)`: a determinant on a coherent
endpoint-translation parallelogram has mixed curvature `det(s,t)` and a
normalized seven-clock family pays the Cech invoice exactly.  That square
varies `q` and is not supplied by the existing current constructions.

Not proved are a physical occupation of the four central allocation states,
a common endpoint atom for all four currents, survival of `mu` in one
determinant sector, the root-deck functional `(26)`, the cyclic invoice
`(27)` on a physical current, realization of `(15a)` by carrier allocation,
promotion of a virtual face `(22d)` to common co-support, semantic
arrival, all-`91`-unit current, row exclusion, or LRC(14).

The independent hostile audit rederived the pullback and its universal
property, repaired the allocation lift to the typed object
`A_tilde=V x A`, checked every affine, Segre--Hadamard, Pluecker, virtual-face,
and symplectic formula and census, and replayed normal and optimized modes
against the stored output and declared hashes.  No remaining defect was
found.

QED.
