# LRC address/chamber session: determinant survival and the missing common atom

**STATUS: READY — FINITE-EXACT candidate + exact structural no-go.**

No tracked file was edited.  The two companions in this directory replay
byte-identically under ordinary and optimized Python and use explicit
`require` gates with no `assert` nodes.

## 1. Inheritance and live concept board

The closest proved mechanisms are:

- THM-2763, which replaces the false old-`169` fixed-carrier bank by
  `G_full` and its representative gauge;
- THM-2772, which proves the determinant-ready pullback
  `G_full x_V E` and isolates the missing common-current four-corner lift;
- THM-2754 and the all-`81` wing audit, which provide the physical
  clock/label chamber;
- the source-clocked natural sheet, which is a fixed-physical-clock slice of
  the one-sided carrier family and has gain `11` on all `81` endpoint labels.

The canonical hostile is THM-2763's fixed-carrier gauge failure.  The closest
corrected near miss is MISTAKE-313: omitting the source-one clock manufactured
a false wing scalar.  The least-used relevant sidecar is reserved THM-2773's
pure-`B` equal-mass/nontranslate face.

The board used here was:

```text
G_full / G_delta;
allocated endpoint plane E;
physical all-81 chamber;
source-clocked one-sided bank;
factor first-failure origin;
common physical ancestry;
determinant sector.
```

## 2. The four finite objects must remain distinct

Put `V=K_13/L_13 isomorphic F_13^2`.

THM-2763 gives

```text
G_full = {(r,k,l): r.W+k-l=0}/(L_13,0,0)
       isomorphic V x F_13 x F_13,
|G_full|=13^4.
```

With the canonical `e_0` section, write its coordinates as `(q,k,l)`, where

```text
q=[r-(l-k)e_0] in V.
```

Forgetting the common carrier harmonic gives

```text
G_delta isomorphic V x F_13_(k-l),
|G_delta|=13^3.
```

THM-2620's allocated endpoint plane is a different object,

```text
E=V x V,                    (L,R) |-> q=L-R,
|E|=13^4.
```

The determinant-ready object is their pullback

```text
P=G_full x_V E
 ={(L,k,R,l)}
 isomorphic {(q,k,l,R)},
|P|=13^6.                                                   (1)
```

The projection `P -> G_full` forgets the common endpoint origin `R` and has
fibre `V` of size `169`.  For fixed `q != 0`, the function

```text
Delta=det(L,R)=det(q,R)
```

takes every value exactly thirteen times on this fibre.

The dual typing is equally important.  The left bank lives on

```text
(ell,a) modulo (W,1),
```

and the conjugate right bank uses

```text
(ell,b) modulo the analytic gauge (ell,b)~(ell+W,b-1).
```

The physical carrier translations are `-a/13` at the source and `+b/13` at
the target.  Their product descends under

```text
(ell,a,b)~(ell+W,a+1,b-1).                              (2)
```

Independent inverse transforms of the two `13^3` banks produce the six
coordinates `(L,k,R,l)` in (1).  This is an actual transform construction of
the pullback, not a cardinality identification.

Factor allocation is finer again.  A monolithic selected carrier needs only
its total harmonic `k` or `l`; an all-factor mask must also retain the
individual Fourier origins whose sum is that total.  The wing difference is
especially delicate because

```text
one-sided minus common
 = one-sided intersect complement(translated full section),
```

and the complement is a union.  The lawful atomic refinement disjointizes it
by the first failed translated factor.  The order used in the exact audit is

```text
shifted E3, shifted clock, shifted q1, shifted q2,
shifted c2, shifted c3.                                  (3)
```

## 3. A genuine selected-wing determinant sector survives

The broad `A` chamber is not the cheapest factor-labelled probe.  At
`(s,t,e)=(8,3,3)`, the source coefficient sees only first-failure `E3`, while
the target coefficient sees both `E3` and `c2`.  The gain-`2` `B` cell
`(2,3,3)` has the same coefficient-level mixed-origin obstruction.

The gain-`2` `D` cell

```text
(s,t,e)=(0,10,3)
```

is factor-pure at the first-failure resolution: both source and target
coefficients come only from shifted-`E3` failure.  Its exact raw values are

```text
source = 2853968755527296447040,
target = 5707937511054592894080 = 2 source.              (4)
```

It is not a physical translate:

```text
piece counts:  source 3, target 2;
weighted mass: source 2188067024426754240,
               target 1458711349617836160.
```

Insert those two weighted carrier factors into THM-2763's fixed triangle

```text
(X,m,Y)=(12*13^4,1,38*13^4)
```

before target difference.  The complete separately allocated dual banks have
support

```text
left  641/2197,
right 510/2197,
```

and pass the representative gauge (2).

After the independent three-dimensional inverse transforms, the actual
THM-2620 determinant sector

```text
q=(1,0),       Delta=1,       k=l=0                       (5)
```

has all thirteen endpoint edges nonzero.  Its sector sum in THM-2625's first
certified exact field is

```text
64970359711488738 mod 352341050142921841 != 0.            (6)
```

Because this specialization is a proved cyclotomic ring homomorphism, (6)
certifies characteristic-zero nonvanishing.  This is the requested actual
determinant-sector test after a lawful gain-`2` factor-labelled selector.

It does not make the source and target carriers translates, positive
endpoint-current weights, or a physical transition.

## 4. The exact missing physical lift: the `11` corner is absent

The all-bank companion reconstructs every lawful `(s,t,e)` cell:

```text
81 labels x 7 physical clocks = 567 cells.
```

In every cell it checks the literal object decomposition

```text
one-sided = common disjoint-union wing.
```

Exactly `193/567` source wings and `193/567` target wings are physically
nonempty.  Nevertheless, after pulling the target wing back to the source
chart,

```text
source wing intersect pullback(target wing) = empty
```

in all

```text
567/567 cells.                                           (7)
```

This is the object-level anatomy.  If `A` and `B` denote the two aligned
one-sided carrier supports, then the common carrier is `A intersect B`, while
the two wings are the opposite exclusive pieces

```text
A minus B,                B minus A.
```

Their intersection is identically empty.  Hence the wing pair cannot supply
one physical atom carrying the Boolean allocation state `(1,1)`.  An
independently formed endpoint Fourier product `P_1 Q_1` may be nonzero—as
(6) proves—but it is the virtual mixed face of THM-2772, not a common
physical K4 corner.

This also explains the source-clocked gain coincidence.  The natural
source-clocked carrier is the fixed-`e` slice of the one-sided family, and

```text
one-sided = common + wing.
```

The A chamber and the source-clocked bank can both show gain `11` after
clock quotient while the added source/target wing pieces remain mutually
exclusive.  Equality of quotient gains does not restore common ancestry.

Selected first-failure boundaries sharpen (7):

- `A=(8,3,3)`: coefficient source is `E3`-only, coefficient target is
  `E3+c2`;
- gain-`2` `B=(2,3,3)`: the same coefficient-level mixed origin, despite
  equal total source/target physical mass;
- gain-`2` `D=(0,10,3)`: both are `E3`-only, but they are still disjoint,
  nontranslate, and mass-unequal.

Thus changing from A to D removes the first-factor ambiguity but not the
common-atom obstruction.

## 5. Smallest remaining sidecar

A future positive lift must not use the two wing cofibres as its two carried
arms.  The minimal still-lawful source object must retain one ancestry atom
`xi` together with

```text
clock, semantic word/owner, and external (s,t,e);
four allocation amplitudes P_0(xi),P_1(xi),Q_0(xi),Q_1(xi);
the exact address (r,k,l);
the endpoint pair (L,R);
two endpoint displacements s_xi,t_xi in V;
all four varying target sectors
  q, q+s_xi, q-t_xi, q+s_xi-t_xi;
the first-failure/factor origin at both endpoints.         (8)
```

The determinant mixed face is then lawfully

```text
det(L,R)-det(L+s_xi,R)-det(L,R+t_xi)
 +det(L+s_xi,R+t_xi)=det(s_xi,t_xi).
```

The present computation shows:

```text
coefficient pullback + determinant sector:  EXISTS;
common wing atom:                           IMPOSSIBLE;
atom-dependent one-sided/common lift:       OPEN.
```

Reserved THM-2779 proposes a distinct next-layer Heisenberg obstruction:
a literal permutation carrier preserving the symplectic commutator needs
`169`, not `13`, states.  That result should remain separate.  The current
no-go occurs earlier, before any faithful root-deck action is available.

## 6. Reproduction and evidence

From the repository root:

```bash
python3 .scratch/lrc_address_chamber_subagent/all81_common_atom_no_go.py
python3 -O .scratch/lrc_address_chamber_subagent/all81_common_atom_no_go.py

python3 .scratch/lrc_address_chamber_subagent/extended_d_sector_probe.py
python3 -O .scratch/lrc_address_chamber_subagent/extended_d_sector_probe.py
```

Normal and optimized outputs match byte-for-byte.

LF-normalized SHA-256:

```text
all81 script   a043248e907f402f283f370f194f4bc3be7d0f72afa531eb321563b537b5ad1e
all81 output   150b05b0887ab38044243ee3175786579d1404c09cef26e23a5841ca636ecba8
sector script  0aed15f738255316d59b546f1659d322197b30f355c15e4ac9cadc3658b0303d
sector output  cfdc4d3b4e466f4fdffa8deaa621afa55ddabbdfd42c3119427a333367411d9d
```

The hash block above must be refreshed if either final script is edited.

## 7. Scope

This session proves one exact selected-wing determinant-sector
nonvanishing and one all-bank physical common-atom no-go.  It excludes no
scalar row.  It does not construct the one-sided/common four-corner lift,
the endpoint-dipole-to-root-deck map, a physical Heisenberg action, or
LRC(14).

