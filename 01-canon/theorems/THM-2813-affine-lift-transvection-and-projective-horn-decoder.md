---
id: THM-2813
title: "Affine-lift transvection, fixed-sheet isotropy, and projective horn decoder"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The thirteen full-depth affine lifts of THM-2807 form one explicit
  C13 transvection torsor.  Their common depth-two shadow is the
  fixed-line Heisenberg shear of slope ten, while their relative top-graded
  action is (r,h)->(r,h+tr).  The low-digit-seven sheet is the exact fixed
  locus, every off-sheet point is a free thirteen-orbit, and one adjacent
  normal jet sharply recovers t.  THM-2803's coordinates 6,7 separately
  decode the determinant fibre under an explicit common-scalar/origin
  hypothesis, giving a conditional 2+1 decoder for 169 abstract flags.
  THM-2806's fixed Rees D3 has no typed lift action and does not by itself
  define the normal jet.
  No physical allocation-to-endpoint map, root/Cech map, row exclusion,
  or LRC(14) conclusion is claimed.
source: root/lrc-affine-lift-horn-decoder-2026-07-28
audit: >
  root/lrc-holotopy-allocation/audit-transverse-flag-2026-07-28
  (independent affine-composition, fixed/free orbit, depth-two and top-jet,
  projective-decoder, THM-2806 typing, normal/-O/stored/hash and docs audit:
  ACCEPT)
depends_on:
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
  - THM-2803-endpoint-current-determinant-fibre-projective-nonflatness
  - THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go
  - THM-2807-positive-graded-address-two-simplex-and-allocation-lift-boundary
related:
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
script: 04-computation/lrc14_affine_lift_transvection_horn_decoder_thm2813.py
output: 05-knowledge/results/lrc14_affine_lift_transvection_horn_decoder_thm2813.out
script_sha256: 255ce911a18f33d6b700213d6362886970b12170c324d39bed82a69a51b63e83
output_sha256: 32f61740a0e7e73384c3d1ff1a83ba30d4cd3ebeca5d228b57e0bf2510928d58
hash_basis: LF-normalized bytes
---

# THM-2813 -- the address simplex is a fixed horn for a transvection

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2807 proves that one positive address triangle has thirteen full-depth
affine lifts, all agreeing on its low-digit-seven sheet.  This could have
been an unstructured thirteenfold ambiguity.  It is not.  Relative to one
lift, the other twelve are the powers of one fixed-line transvection:

```text
A_t(y)=y+t*13^5(y-7 mod13)                       mod13^6.       (1)
```

Its first normal action is the standard Heisenberg shear.  The selected
simplex lies in the fixed line, so no amount of on-sheet address data can
choose `t`; one adjacent sheet does so sharply.

There is a complementary coefficient coordinate.  THM-2803's all-minors
theorem makes the two endpoint values at central coordinates `6,7` a
projective decoder for the determinant fibre.  Under one explicit
common-scalar/origin hypothesis, those two values recover `delta`, while
the third THM-2807 vertex repeats coordinate `7` and remains blind to `t`.
The normal jet supplies the missing coordinate and yields `13^2=169`
abstract flags.

This is an exact address/coefficient flag theorem, not the missing physical
map.  THM-2806's central Rees residue is also one, but its marked selector is
not identified with the address lift torsor and no `t` action on it is
defined.  Equality of those two scalars does not construct the boundary map
between them.

## 1. The thirteen affine lifts

Put

```text
p=13,               M=p^6=4826809,              S=p^5=371293,

n_0=3454614,        n_+=3454627,                 n_a=4143978.   (2)
```

THM-2807 gives the unique depth-five exponent class

```text
k_0=23098 mod p^4
```

and its thirteen depth-six lifts

```text
k_t=k_0+t p^4,
m_t=14^(k_t) mod M,
v_t=(1-m_t)n_0 mod M,
g_t(y)=m_t y+v_t,                         t in F_13.       (3)
```

Their exact values satisfy

```text
m_0=2652079,                 v_0=352469,

m_t=m_0+tS,                  v_t=v_0+6tS                 mod M,

g_t(n_0)=n_0,                g_t(n_+)=n_a.                (4)
```

The first equality in the second line follows either from THM-2807's exact
lift calculation or directly from

```text
(1+p)^(p^4)=1+p^5                              mod p^6.   (5)
```

Since `m_0=1 mod p`, multiplication by `m_0^(-1)` does not change a
coefficient already divisible by `S`.  Therefore

```text
m_t m_0^(-1)=1+tS                              mod M.     (6)
```

Now define the genuinely relative map

```text
A_t=g_t g_0^(-1).                                         (7)
```

It is not enough to compare `g_t(y)-g_0(y)` at a common input; equation
`(7)` compares the two lifts as transformations.  Both `g_t` and `g_0` fix
`n_0`.  Equations `(6)--(7)` therefore force

```text
A_t(y)
 =n_0+(1+tS)(y-n_0)
 =y+tS(y-7 mod p)                              mod M.     (8)
```

Equivalently its affine pair is

```text
A_t=(1+tS,-7tS).                                          (9)
```

Because `S^2=0 mod M` and `pS=0 mod M`,

```text
A_t A_u=A_(t+u),                  A_t^p=1.                (10)
```

Thus the lift set is a principal `C_13` torsor under the explicit
transvection group `(9)`.

## 2. The common depth-two shear is the transgression ten

Every `g_t` has the same reduction modulo `p^2`, because the lift correction
in `(4)` is divisible by `p^5`.  Since

```text
14^(k_0)=1+p(k_0 mod p)=1+10p                    mod p^2,
```

and the map fixes the low-digit-seven line, the common quotient affine pair
is

```text
g_t mod p^2=(1+10p,-7*10p)=(131,104)              mod169. (11)
```

Write a two-digit address as

```text
y=v+p w,                    v,w in F_13.                  (12)
```

Then `(11)` is

```text
(v,w) |-> (v,w+10(v-7)).                                  (13)
```

THM-2779's endpoint Heisenberg action has generators

```text
X(v,w)=(v,w+v),                 Z(v,w)=(v,w+1).           (14)
```

Hence the common depth-two address action is exactly

```text
(Z^(-7)X)^10.                                             (15)
```

The two quotient vertices of THM-2807 lie on `v=7`:

```text
n_0=(7,6),                  n_+=n_a=(7,7)          mod169. (16)
```

Thus the quotient horn is fixed pointwise, while its normal line is acted
on with slope ten.  This recovers THM-2807's independent edge calculation

```text
4079=10 mod13.                                            (17)
```

The number ten is therefore the first normal transvection of the collapsed
horn, not merely a matching residue.  This remains an address statement:
THM-2779 explicitly does not identify its endpoint-origin plane with the
positive physical rail sheet.

## 3. Fixed-sheet isotropy and the sharp normal jet

For a nonzero `t`, equation `(8)` gives

```text
A_t(y)=y
 iff p^6 divides t p^5(y-7)
 iff y=7 mod p.                                           (18)
```

Therefore the exact fixed locus and free complement are

```text
L={y:y=7 mod13},                     |L|=13^5=371293,

|M\L|=12*13^5,
number of free C_13 orbits off L=12*13^4=342732.          (19)
```

Define the transverse and top-digit coordinates

```text
r=y-7 mod13,
h=floor(y/13^5) mod13.                                    (20)
```

The four middle address digits disappear under this quotient, and `(8)`
becomes

```text
(r,h) |-> (r,h+t r).                                      (21)
```

This is the standard parabolic shear `X^t` from `(14)`, after shifting the
fixed line from `v=0` to `v=7`.  On the full depth-six address group the
unit relative lift is

```text
A_1=O^(-7*13^5) X_full^(13^4).                            (22)
```

Let

```text
H(y)=floor(y/13^5) mod13.                                 (23)
```

Its exact cocycle is

```text
H(A_t y)-H(y)=t r.                                        (24)
```

Every on-sheet address observable is therefore blind to `t`, because every
`A_t` fixes every point of `L`.  On the adjacent sheet `r=1`,

```text
H(A_t y)-H(y)=t,                                          (25)
```

so one transverse address recovers the lift uniquely.  This is sharp in
both senses:

1. zero transverse addresses cannot distinguish any of the thirteen lifts;
2. every off-sheet point has a free thirteen-orbit, attaining the sharp
   thirteen-state equivariant cost classified abstractly by THM-2611.

The minimal missing address sidecar is therefore one **normal jet**, not
another invariant of the fixed simplex.

## 4. The projective two-value determinant decoder

Fix a nonzero central direction `s`.  In THM-2803's notation write

```text
j_delta(w)=J_s(w s+delta t_s),           delta,w in F_13, (26)
```

where `det(s,t_s)=1`.  Its all-minors theorem applies in particular to the
coordinate pair `w=6,7`:

```text
j_delta(6)j_epsilon(7)-j_delta(7)j_epsilon(6) !=0
                                            for delta!=epsilon. (27)
```

Every displayed coefficient is nonzero.  Hence

```text
delta |-> [j_delta(6):j_delta(7)] in P^1                 (28)
```

is injective.  The exact inherited universe is

```text
168 directions * C(13,2)=13104
```

nonzero specialized minors in each certified field for this one coordinate
pair.  Either field proves the corresponding characteristic-zero
cyclotomic minor nonzero; the second is an arithmetic control.

For a physical interpretation impose the following hypothesis explicitly.

> **Common-scalar/origin hypothesis.**  One endpoint current is attached to
> the positive horn with one endpoint origin and one nonzero scalar
> `lambda`, so its first two vertex values are
>
> ```text
> (C_0,C_+)=lambda(j_delta(6),j_delta(7)).                 (29)
> ```

Under this hypothesis `(28)` recovers `delta`.  No theorem cited here
constructs `(29)`.

The third THM-2807 vertex has the same endpoint quotient coordinate as the
second:

```text
n_+=n_a mod169.
```

A quotient-only common-origin scalar model therefore has

```text
(C_0,C_+,C_a)
 =lambda(j_delta(6),j_delta(7),j_delta(7)).               (30)
```

It is independent of the full-depth lift `t`.  Consequently the horn map

```text
(delta,t)
 |-> [j_delta(6):j_delta(7):j_delta(7)]                  (31)
```

has exactly thirteen points in every fibre.  Adjoining the normal jet `(25)`
gives

```text
(delta,t)
 |-> ([j_delta(6):j_delta(7)],t),                         (32)
```

which is injective on all `13^2=169` abstract flags.

Equations `(29)--(32)` are the sharp **2+1 decoder**:

```text
two common-scalar on-sheet coefficients decode delta;
one adjacent-sheet normal jet decodes t.                 (33)
```

The flag cardinality agrees with THM-2779's minimal faithful Heisenberg
permutation degree.  Cardinality and the shear formula do not canonically
identify this flag with its endpoint-origin fibre; a physical equivariant
map is still required.

## 5. THM-2806's Rees face does not itself define the normal jet

THM-2806 proves on its literal fixed selector

```text
(B,P,Q,H)=w(1,delta_0,delta_0,delta_(0,0)),               (34)

(v00,v10,v01,v11)=w(169,13,13,1),

(D0,D1,D2,D3)=w(196,168,168,144),

(v_13(v00/w),v_13(v10/w),v_13(v01/w),v_13(v11/w))
 =(2,1,1,0),                                              (35)

D3/v11=144=1 mod13.                                      (36)
```

The tempting comparison is

```text
normalized Rees face=1,
normal jet of A_1=1.                                     (37)
```

But THM-2806 also proves:

1. the only pre-Fourier point supporting all four allocation states is
   `(0,0)`, where its raw mixed face is zero;
2. `D3` is the sum of `144` bare-only complement points;
3. literal allocation has zero endpoint step and is projectively scalar;
4. the Rees profile lives on the fixed marked selector.

Every `A_t` fixes the THM-2807 low-digit-seven address sheet pointwise.
THM-2806's Rees face lives on a separately fixed marked selector.  No proved
map identifies that selector with this address sheet or defines an `A_t`
action on it, so physical `t`-invariance is not yet typed.

For the sharp scalar-information hostile, repeat the single residue in
`(36)` over a merely formal `t` index:

```text
(1,1,...,1)                         over t in F_13,       (38)
```

whereas the normal jets in `(25)` are

```text
(0,1,...,12).                                             (39)
```

The companion checks `(38)!=(39)` directly, but `(38)` is not asserted to
be a physical torsor action.  It says only that one scalar residue cannot
encode thirteen jets.  Under any future covariant identification placing
the Rees object on the fixed address sheet it would become invariant;
without such an identification it has no `t` coordinate at all.  The Rees
`D3` is an internal allocation filtration, not by itself the adjacent-sheet
cocycle `(24)`.  Reducing both objects to the displayed scalar one forgets
the boundary map needed to compare them.

A valid positive bridge would have to construct, on one common physical
atom with `r!=0`, a gauge-covariant map

```text
beta_normal:
  allocation Rees mixed face
       -> [H(A_t y)-H(y)]/r,                              (40)
```

while retaining endpoint origin, clock, word/owner, allocation flags, and
the physical address generator.  If `(40)` sent `(36)` to one, it would
select `A_1` uniquely.  THM-2806 supplies neither a nonzero raw common-atom
face nor off-sheet transport, so it does not supply `(40)`.

## 6. First failed physical implication

The first invalid implication is

```text
equal positive whole cylinders + nonzero central Rees D3
  does not imply
one common endpoint current with an adjacent-sheet normal derivative.     (41)
```

The failures occur at distinct typed maps:

- THM-2807 transports whole weighted cylinders, not every hidden factor or
  allocated endpoint atom.
- THM-2803 supplies coefficient profiles, not their attachment to that
  physical rail sheet.
- THM-2806 supplies a fixed-selector allocation filtration, not a transverse
  address motion.

The cheapest decisive physical test is correspondingly small:

1. retain one nonzero raw common allocation atom at an address `y` with
   `y!=7 mod13`;
2. retain its endpoint origin and the same scalar in `(29)`;
3. compare the selected affine lift with the base lift on that atom;
4. verify

   ```text
   A_1(y)-y=13^5(y-7 mod13);                              (42)
   ```

5. verify that the raw mixed face, not only its central sum, maps to the
   normal jet through `(40)`.

At `y=8 mod13`, equation `(42)` is the minimal `+13^5` test.  Failure at
step one is a co-support/allocation obstruction; failure at step four is
carry holonomy; failure at step five leaves a virtual internal Rees face.

## 7. Source, target, map, and loss ledger

| item | exact content |
|---|---|
| source | THM-2807's positive address horn and thirteen affine lifts; THM-2803's determinant profiles; THM-2806's fixed Rees square |
| address map | `g_t g_0^(-1)` followed by the low/top-digit quotient `(20)` |
| coefficient map | the projective pair at coordinates `6,7` |
| target | the abstract `169`-element determinant/lift flag `(32)` |
| preserved | address generator, fixed low-digit sheet, quotient coordinates `6,7`, relative origin gauge, and conditionally one common scalar |
| destroyed | four middle digits, hidden factors, positive off-sheet support, semantic owner/source, root/Cech data, and physical endpoint allocation |
| missing sidecar | one common off-sheet allocation atom and the normal boundary map `(40)` |
| cheapest test | one residue-eight atom with displacement `+13^5` |

The common depth-two slope ten is a real normal action, and the top lift
index is a real graded action.  Neither supplies the semantic or root map.
Even success in `(40)` would select only the affine lift and determinant
flag.  It would not pay the seven-chart Cech invoice, preserve a canonical
source-labelled endpoint, exclude a relation row, or prove LRC(14).

## 8. Exact companion

Run

```text
python 04-computation/lrc14_affine_lift_transvection_horn_decoder_thm2813.py
python -O 04-computation/lrc14_affine_lift_transvection_horn_decoder_thm2813.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_affine_lift_transvection_horn_decoder_thm2813.out.
```

The dependency-free companion uses explicit exception gates and no Python
`assert`.  It verifies:

1. all thirteen exact powers, multipliers, translations, and horn images;
2. `g_t g_0^(-1)` as affine pairs, not only pointwise differences;
3. all `169` group products;
4. all `2,197` common depth-two action entries and all `2,197` relative
   low/top-digit action entries;
5. the exact fixed/free counts and all adjacent-sheet jets;
6. the formal `13`-to-one and `169`-point flag arithmetic inherited from
   THM-2803's proved all-minors theorem; and
7. THM-2806's `169/13/13/1` raw support, raw-flat sole common point,
   central Hadamard vector, and the formal one-scalar-versus-normal-jet
   hostile without asserting a physical lift action on the Rees selector.

It deliberately does not rerun THM-2803's expensive immutable cyclotomic
endpoint reconstruction.  The only use of that theorem is the proved
all-minors implication `(27)`.

**QED.**
