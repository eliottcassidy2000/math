---
id: THM-3448
title: "Weighted Keller maps with cyclic Jelonek inertia of every length"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  In the THM-3438
  weighted family, the C=0 Newton polygon and the
  ordinary discriminant component give the exact generic Jelonek anatomy.
  The bounded inverse primitive has index one along C, and an explicit
  degree-five determinant-one member has three escaping sheets in one C3
  inertia orbit.  A subfamily realizes component-generic C_ell inertia for
  every ell>=2, while no weighted quartic can realize C3.  Clean affine-
  primitive clearing formulas require explicit residue separation.
source: codex2 weighted-Keller boundary synthesis, 2026-08-15
audit: independent Newton, reconstruction, primitive-residue, degree-five, infinite-family, and hostile-quartic derivation audit
depends_on:
  - THM-3438-weighted-lift-keller-degree-spectrum
  - THM-3441-weighted-quartic-jelonek-components-and-inertia-parity
related:
  - THM-3440-weighted-lift-cyclic-infinity-torsor-and-7x13-character-grid
  - HYP-9027-twojet-disc-jelonek-odd-exponent-law
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
script: 04-computation/jc_weighted_cyclic_jelonek_inertia_thm3448.py
output: 05-knowledge/results/jc_weighted_cyclic_jelonek_inertia_thm3448.out
script_sha256: a0410cf3d921fce6ddaa801c332c2cad10fae2a279efbada2aadb35bc498ebff
output_sha256: 768b2b1a75b937daa161bcad9f39310829bccbc843799425309ea2d0fe087fbe
semantic_sha256: dbb99005e17dc8e2d69d514244640530476b74f01a723d5e5bd418bf42b27a80
hash_basis: LF-normalized bytes
---

# THM-3448 -- weighted Keller maps with cyclic Jelonek inertia of every length

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem closes the genuine-C3 lane left open by THM-3441, but only for
the three-dimensional weighted Keller family.  It does not classify arbitrary
Keller boundaries or touch `JC(2)`.

## 1. Setup and exact component statement

Use the notation of THM-3438.  Thus, over a characteristic-zero field,

```text
R(w)=integral_0^w p(s) ds,
T(w)=R(w)-Pw+cQ,
P=BC,                  Q=AC^2,                         (1)
gamma=(P-p(w))/c=-partial_w(T)/c.                      (2)
```

The source reconstructs from a root of `T` by

```text
x=C/gamma,
y=(w-gamma)/C,
z=gamma[gamma(gamma-1+a)-a w]/(b C^2).                (3)
```

Put

```text
n=deg R,              mu=ord_0 R=ord_0 p+1.           (4)
```

Because `p(0)=0`, `p(1)=-c!=0`, and `R(1)=0`,

```text
2<=mu<=n-1,           w=1 is a simple root of R.       (5)
```

Assume throughout Sections 1--5 that every nonzero root of `R` is simple.
Let `Delta(P,Q)=Disc_w T`.  Its reduced branch curve is irreducible, and
after pullback there is an irreducible reduced polynomial `L_p` such that,
up to a nonzero scalar,

```text
Delta(BC,AC^2)=C^mu L_p(A,B,C).                        (6)
```

The exact nonproperness set is

```text
S_(F_p)= V(L_p),                  if (n,mu)=(3,2),
S_(F_p)= V(C) union V(L_p),       otherwise.           (7)
```

At the generic point of `C=0`, the component ledger is

| case | escaping sheets | inertia | raw `w` order | `w`-index |
|---|---:|---|---:|---:|
| `mu=2` | `n-3` | identity `1^n` | `2` | `1` |
| `mu>2` | `n-2` | `(mu-1) 1^(n-mu+1)` | `mu` | `1` |

Here a cycle of length one means no permutation.  At the generic point of
`L_p=0`, exactly two sheets escape, inertia is one transposition, the raw
order is one, and the `w`-index is zero.

No component in `(7)` is finite affine branching: THM-3438 gives
`det J(F_p)=bc!=0`.  Both are source-infinity phenomena.

## 2. The `C=0` Newton polygon and escaping roots

Fix a generic target point `(A,B,0)` and use `C` as transverse parameter.
The controlling terms of

```text
T_C(w)=R(w)-BCw+cAC^2                              (8)
```

give the Newton vertices

```text
(0,2),             (1,1),             (mu,0).          (9)
```

If `mu=2`, the three points are collinear.  The two roots issuing from zero
both have valuation one.  Writing `w=CW`, their residual quadratic is

```text
r_mu W^2-BW+cA,                                       (10)
```

where `r_mu` is the initial coefficient of `R`; away from its discriminant,
both roots reconstruct finitely by `(3)`.

If `mu>2`, `(9)` has two sides.  It gives

```text
one root with v_C(w)=1,
mu-1 roots with v_C(w)=1/(mu-1).                       (11)
```

The first root reconstructs finitely.  On the other roots,

```text
w^(mu-1)~B C/r_mu,
gamma~-(mu-1)BC/c.                                    (12)
```

Hence `y` and, generically, `z` have pole order

```text
(mu-2)/(mu-1),                                        (13)
```

and the `mu-1` roots form one tame inertia cycle.

Now let `r` be a nonzero root of `R`.  Simplicity gives `p(r)!=0`, and the
corresponding root is an ordinary power series in `C`.  Put

```text
g_r=-p(r)/c.                                           (14)
```

At that branch, `gamma->g_r`.  The designed root `r=1` has `g_1=1` and
reconstructs finitely.  Every other nonzero root escapes: from `(3)`, its
leading `y` and `z` residues are

```text
y~(r-g_r)/C,
z~Z_r/C^2,
Z_r=[g_r^2(g_r-1+a)-a r g_r]/b.                       (15)
```

These two residues cannot vanish simultaneously.  Indeed, `r=g_r` and
`Z_r=0` would force `r=g_r=1`, the excluded designed root.

There are `n-mu-1` roots of the type `(15)`.  Combining them with `(11)`
proves the escaping-sheet counts in the table.  It also proves that `C=0`
is generically nonproper except when `(n,mu)=(3,2)`.

In that exceptional cubic, the root near one and both roots from `(10)` are
finite.  Moreover the residual discriminant in `(10)` is, up to a unit,
`L_p(A,B,0)`.  Therefore every point of `C=0` outside `L_p=0` has a finite
local inverse algebra, so `C` is not a separate Jelonek component.

## 3. Raw discriminant, inertia, and the universal index one

The root differences in the zero cluster compute the exact order of `(6)`.
For `mu=2`, the one pair contributes order two.  For `mu>2`, differences
among the `mu-1` cyclic roots contribute

```text
2 binom(mu-1,2)/(mu-1)=mu-2,                           (16)
```

and differences between the valuation-one root and those cyclic roots
contribute two more.  Thus

```text
v_C(Disc_w T)=mu.                                      (17)
```

The normalized local cover has one tame cycle of length `mu-1`, so its
maximal-order discriminant defect is `mu-2`.  Discriminants of an order and
its normalization differ by the square of the order index.  Consequently

```text
mu=(mu-2)+2 i_w,              i_w=1.                   (18)
```

Thus the bounded primitive `w` always carries one unit of monogenic index
along `C`.  For `mu=2`, its entire `C^2` factor is an index square.  For
`mu>2`, the same index square sits beside genuine cyclic inertia.

This sharpens THM-3441's quartic identity `2=0+2*1`; it does not identify a
finite critical branch.

## 4. The ordinary discriminant component and exact Jelonek set

The incidence ramification equations are

```text
T(r)=0,               p(r)=P.                          (19)
```

They parametrize the branch curve by

```text
(P,Q)=(p(r),q(r)),
q(r)=[r p(r)-R(r)]/c.                                  (20)
```

At a generic point,

```text
dQ/dP=r/c.                                             (21)
```

The tangent slope recovers `r`; hence the curve is irreducible, generically
reduced, and has one generic double root.  Since `(A,B,C)->(P,Q,C)` is
birational over `C!=0`, its closure is the irreducible component `L_p` in
`(6)`.

Transversely to `L_p`, the two roots satisfy `w-r~+-sqrt(t)` and
`gamma~sqrt(t)`.  Because `C` is a unit there, `(3)` gives
`x~t^(-1/2)`.  Exactly those two sheets escape and are exchanged.  Thus

```text
sigma_L=(12)1^(n-2),
v_L(Disc_w T)=1,             i_(w,L)=0.                (22)
```

Off `C L_p=0`, the monic algebra `(T)` is finite and separable and both `C`
and `gamma=-T_w/c` are units.  Formula `(3)` then reconstructs every source
coordinate inside that finite algebra.  Together with the exceptional-cubic
argument in Section 2 and the escaping branches above, this proves `(7)`.

## 5. Generic affine primitives and the residue gate

Let

```text
Theta=lambda_x x+lambda_y y+lambda_z z                (23)
```

be an affine source primitive, and let `q_Theta` be its monic degree-`n`
polynomial over the target function field.  The clean intrinsic clearing
formulas require the following Zariski-open conditions at the chosen generic
point of `C=0`:

1. `lambda_z!=0`, and the `Z_r` in `(15)` for
   `r notin {0,1}` are nonzero and pairwise distinct;
2. when `mu>2`,

   ```text
   lambda_y+lambda_z a(mu-1)B/(bc) !=0;                (24)
   ```

3. the residual values on the finite sheets are distinct.

The pole profile is then

```text
H=n-mu-1 roots of order 2,
mu-1 roots of order (mu-2)/(mu-1) if mu>2,
two finite roots if mu>2, and three if mu=2.            (25)
```

The minimal scalar clearing order and the cleared-discriminant exponent are

```text
rho_C=2n-mu-4,                                         (26)

E_C=-v_C(Disc q_Theta)
   =2n(n-1)-mu^2-2mu-4.                                (27)
```

Indeed, `(26)` is the sum of the pole orders in `(25)`.  For `(27)`, the
pole-two pairs contribute `2H(n+mu)`, and the remaining cyclic pairs
contribute `mu^2-4`; their sum is exactly the right side.  In particular,

```text
E_C == mu == mu-2  (mod 2),                            (28)
```

so primitive parity equals the sign of the `C`-inertia.

Along `L_p`, require `lambda_x!=0`.  Two roots have pole order `1/2` and the
other `n-2` are finite.  There are `2n-3` pairs involving an escaping root,
so

```text
rho_L=1,                 E_L=2n-3.                     (29)
```

This again has the parity of the transposition.

The residue gate is essential, not cosmetic.  Let `r` be any root of

```text
2r^4-4r^3+5r^2-6r-5=0                                 (30)
```

and take

```text
R(w)=w^2(w-1)(w-r)/(r-1),              c=1.            (31)
```

This is a valid quartic weighted seed with simple nonzero roots and
`kappa!= -2`, but exact reduction modulo `(30)` gives `Z_r=0` while
`r-g_r=r+r^2!=0`.  Its remaining nonzero branch therefore has pole order one,
not two, for a generic affine primitive.  The correct ledger is

```text
(rho_C,E_C)=(1,6),                                     (32)
```

not `(2,12)`.  The inertia, escaping count, raw order, index, and Jelonek set
remain unchanged.  This is the equality/failure boundary for `(26)--(27)`.

## 6. The first Keller `C3` component

Choose

```text
p(w)=4w^3-5w^4,
R(w)=w^4-w^5,
q(w)=3w^4-4w^5,
b=c=1,                 a=-7/6.                         (33)
```

With

```text
u=1+xy,                 gamma=1-(7/6)xy+x^2z,          (34)
```

the weighted map is

```text
F_5=( [u+3u^4 gamma^2-4u^5 gamma^3]/x^2,
      [1+4u^3 gamma^2-5u^4 gamma^3]/x,
      x gamma ).                                       (35)
```

Both quotients cancel polynomially.  Exact expansion gives

```text
det J(F_5)=1,
total coordinate degrees=(17,16,4),
z-degrees=(3,3,1).                                     (36)
```

Its inverse polynomial is

```text
f(w)=w^5-w^4+Pw-Q,             P=BC, Q=AC^2.           (37)
```

The exact discriminant is

```text
Disc_w f=C^4 L_5,                                      (38)

L_5=3125A^4C^4-2500A^3BC^3+256A^3C^2
    -50A^2B^2C^2-36AB^3C+256B^5C-27B^4.               (39)
```

The polynomial `L_5` is irreducible and `L_5(A,B,0)=-27B^4`, so a generic
`C`-meridian with `B!=0` avoids the other component.

Set `C=s^3` and `w=sW`.  Then

```text
s^(-4)f(sW)=sW^5-W^4+BW-As^2,
residual=W(B-W^3).                                     (40)
```

The three nonzero roots satisfy `W^3=B` and form one cycle under `s`-monodromy.
On them,

```text
gamma~-3Bs^3,
x->-1/(3B),
y~W s^(-2),
z~-(7/2)BW s^(-2).                                    (41)
```

Thus exactly three sheets escape and

```text
sigma_C=(123)1^2.                                      (42)
```

For `lambda_y-(7/2)B lambda_z!=0`, all nine root pairs involving an
escaping sheet have `C`-valuation `-2/3`.  Hence the complete degree-five
ledger is

| component | escaping | inertia | raw order | `w`-index | `rho` | `E` |
|---|---:|---|---:|---:|---:|---:|
| `C=0` | 3 | `C3` | 4 | 1 | 2 | 12 |
| `L_5=0` | 2 | transposition | 1 | 0 | 1 | 7 |

The separate `Q=infinity` place has THM-3440's five-cycle, while THM-3438
gives global monodromy `S_5`.  These are three distinct boundary/cover
phenomena, not one cyclic global cover.

## 7. Every cyclic length and the sharp quartic boundary

For every integer `ell>=1`, take

```text
p_ell(w)=(ell+1)w^ell-(ell+2)w^(ell+1),
R_ell(w)=w^(ell+1)(1-w),
q_ell(w)=ell w^(ell+1)-(ell+1)w^(ell+2).               (43)
```

Then

```text
p_ell(1)=-1,
p_ell'(1)=-2(ell+1),
a_ell=-(2ell+1)/(2ell),
n=ell+2,                 mu=ell+1.                     (44)
```

All THM-3438 seed gates hold.  For every `ell>=2`, `C=0` is a genuine
Jelonek component with

```text
ell escaping sheets in one C_ell inertia orbit,
v_C(Disc_w)=ell+1,              i_w=1,
rho_C=ell-1,                    E_C=(ell-1)(ell+3).     (45)
```

There are no additional nonzero roots, so the pole-two residue condition is
vacuous; the cyclic amplitude condition is a nonempty Zariski-open condition.
This realizes component-generic cyclic inertia of every length in explicit
three-dimensional determinant-one Keller maps.

For a weighted quartic, `n=4` and `(5)` forces `mu<=3`.  The `C`-cycle has
length at most two, and the only other generic component has transposition
inertia by Section 4.  Hence

```text
no weighted quartic has component-generic C3 inertia;  (46)
degree five is the first weighted realization.         (47)
```

The degree-five witness is `z`-cubic.  It therefore does not settle the
narrower quartic/two-jet Keller-`C3` question outside this weighted family.

## 8. Connection and loss ledger

| field | exact content |
|---|---|
| source | THM-3438 weighted seed `p`, inverse primitive `w`, and target pullback `P=BC,Q=AC^2` |
| target | the reduced Jelonek components and their tame inertia cycles |
| preserved | generic field degree, source reconstruction, discriminant order, inertia, and pole valuation |
| raw primitive loss | `w` can identify finite sheets and contributes the universal index square |
| affine-primitive sidecar | all pole residues, cyclic amplitude, and finite residual separation |
| cheapest positive `C3` | degree-five seed `(33)` |
| cheapest hostile to automatic `E` | quartic algebraic seed `(30)--(32)` |
| quartic survivor | no `C3` inside the weighted family; arbitrary quartic/two-jet boundary remains open |

The theorem refutes any remaining claim that genuine cyclic Keller boundary
inertia must be a transposition, while preserving the corrected law

```text
cleared-discriminant parity = infinity-inertia sign.   (48)
```

It gives no `JC(2)`, D5, LRC, or classification of all Keller maps.

## 9. Exact companion

The deterministic companion checks the general Newton/index/pair invoices
for all `2<=mu<n<=12`, six symbolic discriminant pullbacks, the exceptional
cubic residual algebra, full expansion and Jacobian of `(35)`, `(38)--(42)`,
the first sixteen cyclic seeds, and the exact quartic residue hostile.  Normal
and optimized replay commands are

```text
python -B 04-computation/jc_weighted_cyclic_jelonek_inertia_thm3448.py
python -B -O 04-computation/jc_weighted_cyclic_jelonek_inertia_thm3448.py
```

Normal and optimized runs reproduce the pinned transcript byte-for-byte after
LF normalization.  The structural proof, residue hostile, and exact replay
have passed independent audit.  This completes the proof.  QED.
