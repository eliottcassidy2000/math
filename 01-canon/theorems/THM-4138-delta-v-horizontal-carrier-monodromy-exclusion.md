---
id: THM-4138
title: "Delta-V horizontal-carrier monodromy exclusion"
status: >
  PROVED RELATIVE TO THM-4120/4122/4130/4134 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. The degree-16 and degree-15 finite horizontal BC
  carriers left by THM-4134 are impossible. The BC divisor forces a unique
  nodal horizontal image through a rank-one Mordell--Weil classification;
  puncture-avoiding vanishing loops, the two BC transpositions, and an
  orbit-merger capacity inequality then contradict transitivity. Hence the
  theta-only exact-M=8 Delta_V=0 wall is empty. Other collision walls,
  M>=9, other cells, JC(2), and DC(2) remain open.
source: codex-frontier-synthesis-creative-20260825w
depends_on:
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
related:
  - THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy
  - THM-4053-jc2-live-max-eight-trichotomy-and-eisenstein-survivor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
script: 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138.py
output: 05-knowledge/results/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138.out
independent_audit_script: 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138_independent_audit.out
script_sha256: 61744bea33439e2a1be80ed86e9fa5a0ca8bb82b90be1575af05b0840afed9e2
output_sha256: f8f72753eb04eb6e1a00706f6ff1177f287e6c783084de5d7816cd6003e0f97c
independent_audit_script_sha256: 107c9d0d7ce86f68e9ea64ccf0d882917b288d8cc9c342da28a27be68adb72fe
independent_audit_output_sha256: 88a3a9d6368017800a046e2705ad069e29db9483334b0e395c2ace8202bd77fb
semantic_sha256: 187c7e26f7b8ec16fa3a286502a3d87cc7edf6e4a8b6804047f0e5b46fa88449
independent_semantic_sha256: 187c7e26f7b8ec16fa3a286502a3d87cc7edf6e4a8b6804047f0e5b46fa88449
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone SymPy audit checks the BC degree tower, the normalized
  elliptic surface, exact multiples through 12P, the all-n height ledger,
  forced horizontal curve and target contacts, marked branch points, and
  sharp support-capacity hostiles. Normal, optimized, and two hash-seeded
  streams byte-match the frozen output.
independent_audit: >
  ACCEPT. A clean-room standard-library referee imports no primary code. It
  implements rational functions and the elliptic group law from scratch,
  independently recovers the Mordell--Weil and contact data, locates the BC
  points on the standard second vanishing cycle, audits the pushed-loop
  surface relation and fixed-sheet carrier, exhausts all 331,776 ordered
  four-generator systems on four letters for the orbit-merger inequality,
  and checks both identity cases and sharp one-unit hostiles. All replay
  modes agree and the semantic digests match.
---

# THM-4138 -- Delta-V horizontal-carrier exclusion

**PROVED RELATIVE TO THM-4120/4122/4130/4134 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.**

THM-4134 leaves two lower-degree branches because the quadratic BC punctures
may map to finite points on a horizontal component of the nonproperness set.
That carrier is much more rigid than an arbitrary moving branch value. Its
normalization forces one explicit nodal cubic, and its two moving punctures
cost only two transpositions in the monodromy generation budget. The affine
critical points force more fixed sheets than that budget can connect.

## 1. Inherited residuals

On the theta-only exact-`M=8` wall `Delta_V=0`, THM-4134 reduces every
hypothetical nonautomorphic Keller pair to one of

```text
generic D!=0:  degree n=16, affine critical length L=19,
               packet (7,5,3,2,2,1);

secondary D=0: degree n=15, affine critical length L=18,
               packet (7,3,2,2,2,2,1).                 (1)
```

In each row the two BC punctures have index two and finite image. Every other
boundary puncture maps to the target origin `O`. We prove that neither row in
`(1)` exists.

## 2. The BC divisor is a quadratic target section

After absorbing a nonzero scalar, the normalized BC divisor is

```text
B=A1_r,                         q=a^3/2+r^2.             (2)
```

Let `S` be its irreducible horizontal image in the nonproperness set and let
`S~` be the normalization. THM-4122 gives `S~=A1`. The boundary map factors

```text
B --h--> S~ --q_S--> A1_q,
deg(h) deg(q_S)=2.                                         (3)
```

Both maps are polynomial: each normalization has a unique place at infinity,
and `(2)` has no finite pole. If `deg(q_S)=1`, the affine coordinates of
`S~` would give a nonzero `k(q)`-point on the target pencil, contradicting
THM-4120's exact calculation

```text
E_q(k(q))={O}.                                             (4)
```

Therefore

```text
deg(h)=1,                         deg(q_S)=2.              (5)
```

The degree-one map identifies the affine parameter on `S~` with an affine
linear function of `r`; its target coordinates are consequently polynomials
in `r`, not merely rational functions.

## 3. Mordell--Weil classification after quadratic base change

Scale the target coordinates and `r` to obtain

```text
E_t: y^2=x^3-3x+2+t^2.                                    (6)
```

Its discriminant and singular fibres are

```text
Delta=-432t^2(t^2+4),
I2 at t=0,          I1 at t=+-2i,          IV* at infinity. (7)
```

Indeed, in the `u=1/t` integral model,

```text
a4=-3u^4,             a6=u^4+2u^6,
(v(c4),v(c6),v(Delta))=(4,4,8).                           (8)
```

This is a rational elliptic surface. Its Neron--Severi rank is ten and the
reducible-fibre root lattice is `A1+E6`, so Shioda--Tate gives Mordell--Weil
rank one. Put

```text
P=(1,t).                                                   (9)
```

The section meets the nonidentity `I2` component and a nontrivial `IV*`
component, with local height contributions `1/2` and `4/3`. Since `P.O=0`,

```text
<P,P>=2-1/2-4/3=1/6.                                     (10)
```

The trivial lattice has absolute discriminant `2*3=6`, while the rational
surface's Neron--Severi lattice is unimodular. If `P=mG+torsion` for a
primitive free generator `G` and torsion order `tau`, the determinant formula
gives

```text
1/6=<P,P>=m^2 tau^2/6.                                   (11)
```

Thus `m=tau=1` and

```text
MW(E_t)=ZP,                         MW_tors=0.             (12)
```

## 4. Polynomial sections and the forced horizontal curve

For nonzero `j`, Shioda's height formula gives

```text
(jP.O)=(j^2-12+3 epsilon_2(j)+8 epsilon_3(j))/12,         (13)
```

where `epsilon_2(j)=1` for odd `j` and `epsilon_3(j)=1` when `3` does not
divide `j`. This intersection is entirely at finite `t`. A nonmultiple of
three meets a nonidentity `IV*` component, while

```text
3P=(1+4t^2/9,-t(8t^2+27)/27)                             (14)
```

reduces on the additive identity component to `(4/9,-8/27)`, with nonzero
additive parameter `-3/2`. No nonzero multiple of that parameter reaches the
zero section in characteristic zero.

Equation `(13)` vanishes exactly for `|j|=1,2,3`. Hence the complete list of
polynomial sections is

```text
+-P,            +-2P,            +-3P.                   (15)
```

The first hostile is

```text
4P=((16t^2+81)/(4t^2),
    (8t^4+216t^2+729)/(8t^3)),                            (16)
```

which already has a finite pole. The polynomial coordinate-degree pairs for
`P,2P,3P` are `(0,1),(0,1),(2,3)`. THM-4122 requires the positive intrinsic
pair `(2rho,3rho)`, so only `+-3P` survives, with `rho=1`.

Undoing the scaling gives

```text
U=a/2+16r^2/(9a^2),
V=+-(r+64r^3/(27a^3)).                                   (17)
```

Both signs have the same irreducible image

```text
S: V^2=(U-a/2)(U+a/4)^2.                                 (18)
```

With normalization `U=a/2+z^2`, `V=z(z^2+3a/4)`,

```text
q|S=a^2(8a+9z^2)/16.                                    (19)
```

The curve `S` has an ordinary node at

```text
(U,V,q)=(-a/4,0,5a^3/64),                                (20)
```

whose Hessian determinant is `3a`. It passes smoothly through the second
target node `(a/2,0)` at `q=a^3/2` and is transverse to both local nodal
branches there. It misses the first target node `(-a/2,0)` at `q=0`, meeting
that fibre instead at two smooth points with `U=-7a/18`.

## 5. Puncture-avoiding vanishing loops

Scale `a=1` and use THM-4130's reference fibre `q_*=1/4`. The two finite BC
branch values are

```text
Q_+-=(1/18,+-11i/54).                                    (21)
```

They lie exactly on the standard second real vanishing cycle. Avoidance is
therefore not automatic.

For the first vanishing cycle, transport from `q=0` to `q_*` by a small
complex detour around the horizontal collision value `5/64`; a fibrewise
loop detour returns to the standard first cycle, which misses `Q_+-`. For the
second cycle, use one boundary-parallel push-off inside its Milnor annulus.
Push only near `Q_+-`, preserving the unique transverse intersection with the
first cycle. Denote the resulting pair by `delta_0,delta_1`.

A thin regular neighborhood of `delta_0 union delta_1` is a once-punctured
torus, and its complementary disk contains `O,Q_+,Q_-`. Thus, up to
orientation and meridian order,

```text
pi_1(E_(q_*)-{O,Q_+,Q_-})
 = <delta_0,delta_1,mu_O,mu_+,mu_- |
      [delta_0,delta_1] mu_O mu_+ mu_-=1>.                (22)
```

Eliminating `mu_O` shows that

```text
delta_0, delta_1, mu_+, mu_-                             (23)
```

are a complete free generating set.

## 6. Fixed sheets and complete monodromy

Let `r_i` be the number of affine preimages of the target node `o_i`. The
Keller condition makes the map a local biholomorphism at every such point.
In target Morse coordinates the base equation is `uv=q-q_i`, and its pullback
has the same equation. Each puncture-avoiding loop in Section 5 consequently
has one closed degree-one lift in every affine inverse neighborhood. Distinct
affine preimages give distinct sheets.

At `o_1`, the two intersections of `S` with a nearby Milnor annulus are the
two square-root points on the standard circle; the chosen inner or outer
parallel circle avoids both and retains those local closed lifts. Away from
`S_F`, the Keller map is proper, quasifinite, and etale, hence finite etale;
global transport preserves and separates the injected sheets. Therefore, for

```text
X=rho(delta_0),                   Y=rho(delta_1),
#Fix(X)>=r_0,                     #Fix(Y)>=r_1.            (24)
```

The only finite ramification values are `Q_+,Q_-`. Their single index-two BC
preimages give transpositions

```text
T_+=rho(mu_+),                    T_-=rho(mu_-).           (25)
```

Removing finitely many points from the connected projective source curve
does not disconnect it, so the unbranched cover of the three-punctured target
is connected. Its monodromy is transitive, and `(22)` shows that it is
generated by

```text
X,Y,T_+,T_-.                                              (26)
```

## 7. Orbit-merger obstruction

> **Orbit-merger lemma.** If permutations `g_1,...,g_s` generate a transitive
> action on `n` letters, then
>
> ```text
> sum_j max(|supp(g_j)|-1,0) >= n-1.                      (27)
> ```

Start with `n` singleton orbit blocks. A nonidentity permutation of support
size `d` can meet and merge at most `d` current blocks, decreasing the block
count by at most `d-1`; the identity decreases it by zero. Reaching one block
requires at least `n-1` mergers. This proves `(27)`.

Put `A=supp(X)`, `D=supp(Y)`. From `(24)` and the inherited identity
`r_0+r_1=L`,

```text
|A|+|D| <= 2n-L.                                         (28)
```

Each transposition in `(25)` contributes one merger unit.

For the generic residual `(n,L)=(16,19)`, equation `(28)` gives support sum
at most `13`. If both `X,Y` are nonidentity, the total merger capacity is at
most `13`; if exactly one is the identity, it is at most `14`; if both are
identities, it is two. Every value is below the required `n-1=15`.

For the secondary residual `(n,L)=(15,18)`, the corresponding maxima are
`12,13,2`, all below `n-1=14`. Both rows contradict transitivity.

The identity case is sharp and load-bearing: an `(n-2)`-cycle together with
two adjacent transpositions connects all `n` letters and attains capacity
`n-1`. The inherited fixed-sheet bound wins by exactly one additional unit.
This proves that both residuals in `(1)` are impossible. **QED.**

## 8. Consequence, scope, and replay

THM-4134 already excludes the degree `20,19` full-boundary branches. Sections
1--7 exclude its only remaining degree `16,15` horizontal branches. Hence

```text
the theta-only exact-M=8 Delta_V=0 wall is empty.          (29)
```

This does not cross either other collision wall, treat `M>=9`, prove entry
into this reduced cell, exclude another cell, or prove `JC(2)` or `DC(2)`.

Replay the primary and clean-room audits with

```text
python3 -B 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138.py
python3 -B -O 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138.py
PYTHONHASHSEED=271828 python3 -B 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138.py
python3 -B 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138_independent_audit.py
python3 -B -O 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138_independent_audit.py
PYTHONHASHSEED=8675309 python3 -B 04-computation/jc23_delta_v_horizontal_carrier_monodromy_exclusion_thm4138_independent_audit.py
```

All six streams byte-match their frozen outputs, and the independent semantic
digests agree. **QED.**
