---
id: THM-4139
title: "Rational three-cycle, order-six lift, and horizontal-carrier separation"
status: >
  PROVED RELATIVE TO THM-4134 + FINITE-EXACT + VERIFIED + INDEPENDENTLY
  AUDITED, WITH RESERVED CANDIDATE NARROWED. The complete rational
  preperiodic graph of
  x |-> x^2-29/16 consists of eight quarter-integers and one exact 3-cycle,
  so there is no rational 6-cycle. This is the unique centered monic
  quadratic over Q with a nondegenerate AP-supported 3-cycle. The unique
  projective interpolant has a determinant-one lift B with B^3=-I and
  B^6=I. The normalized cycle parameters do give three rational points on
  the a=-48,q=0 nodal target fibre, but they are disjoint from THM-4134's 3P
  horizontal section; the pre-reserved horizontal-carrier identification is
  REFUTED. Separate exact censuses type period six modulo 43^2 and for the
  squaring map at 63. No new planar-Jacobian consequence follows.
source: codex-planar-jc-dynamics-20260825
depends_on:
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
related:
  - THM-4138-delta-v-horizontal-carrier-monodromy-exclusion
script: 04-computation/quadratic_29_16_cycle_audit.py
output: 05-knowledge/results/quadratic_29_16_cycle_audit.out
independent_audit_script: 04-computation/quadratic_29_16_cycle_audit_independent.py
independent_audit_output: 05-knowledge/results/quadratic_29_16_cycle_audit_independent.out
script_sha256: 167748a82f6969bbce7f6fe2cb23b8fab3475a2ae8b7739e840156c2d5e4b7ee
output_sha256: 50f7826af4a25f6aac6b1d07525f8c8c0ee4278e33eff84d02f0f1b2f101e5f7
independent_audit_script_sha256: 3c4d9542a8382bbb02fbae7ea1de78c6df0fbac74abe7833b933051bfec9c302
independent_audit_output_sha256: c9687ca3b914def3094f649b489494634505f0b2864d83395aef1d0226f519d7
hash_basis: raw LF bytes
primary_audit: >
  PASS. Dependency-free exact Fraction arithmetic rebuilds the complete
  rational graph, both projective lifts, the nodal normalization and 3P
  separation, all finite functional graphs, and all exponent-doubling
  cycles. Normal, optimized, and hash-seeded replays byte-match the frozen
  output.
independent_audit: >
  ACCEPT. A separate script imports no primary code, uses direct first-return
  censuses, checks the 43-adic tube residue by residue, independently
  evaluates the zero-fibre and horizontal-section coordinates, and rebuilds
  the cyclotomic factor by elementary integer polynomial division. Normal
  and optimized replays byte-match its frozen output.
---

# THM-4139 -- rational three-cycle order-six lift

**PROVED RELATIVE TO THM-4134 + FINITE-EXACT + VERIFIED + INDEPENDENTLY
AUDITED, WITH RESERVED CANDIDATE NARROWED.**

Put

```text
f(x)=x^2-29/16,
T(u)=(u^2-29)/4,                     u=4x.               (1)
```

The theorem has five typed layers: rational dynamics, arithmetic-progression
rigidity, a projective order-six lift, finite-ring period-six controls, and a
roots-of-unity explanation of the Mersenne exception. A proposed
horizontal-carrier identification from the reservation is tested in the
actual THM-4134 target pencil and fails; only a weaker zero-fibre
normalization incidence survives.

## 1. Complete rational preperiodic graph

The rational preperiodic set is exactly

```text
PrePer(f,Q)={a/4 : a in {-7,-5,-3,-1,1,3,5,7}}.         (2)
```

Its complete directed graph is

```text
-7/4 ->  5/4 -> -1/4 -> -7/4,
 7/4 ->  5/4,
-5/4 -> -1/4,
 3/4 -> -5/4,             -3/4 -> -5/4,
 1/4 -> -7/4.                                           (3)
```

Thus the displayed orbit is the only rational periodic orbit. In
particular,

```text
f has no rational point of exact period 6.                (4)
```

### Denominator and escape proof

Let `x` be rational and preperiodic. For an odd prime `p`, if `v_p(x)<0`,
then `v_p(-29/16)>=0`; the two terms have unequal valuations and

```text
v_p(f(x))=2v_p(x).                                       (5)
```

This tends strictly to minus infinity, so no odd prime divides the
denominator of `x`.

At `p=2`, write `s=v_2(x)`. If `s<-2`, then
`v_2(f(x))=2s<s`. If `s>-2`, the constant term dominates and
`v_2(f(x))=-4`, after which the preceding escape applies. The only possible
bounded value is `s=-2`. Then `x=a/4` with `a` odd, and

```text
a^2-29 = 4 (mod 8),                                     (6)
```

so the next iterate again has 2-adic valuation `-2`. Every rational
preperiodic point therefore has the form `a/4` with `a` odd.

Let

```text
R=(2+sqrt(33))/4.                                       (7)
```

This is the positive root of `r^2-r-29/16=0`. If `|x|>R`, then

```text
f(x)=x^2-29/16>|x|,                                     (8)
```

and the positive forward orbit thereafter increases to infinity. Since `R`
is irrational, a rational preperiodic `a/4` must have
`|a|<2+sqrt(33)`, hence `|a|<=7`. Exact substitution gives `(3)`, completing
the all-`Q` proof without a height cutoff or an unproved cycle conjecture.

## 2. Unique arithmetic-progression 3-cycle

Let `f_c(x)=x^2+c` over `Q`, and suppose an exact 3-cycle has underlying set

```text
{m-d,m,m+d},                         d!=0.               (9)
```

Changing the sign of `d` if necessary, orient it as

```text
m-d -> m+d -> m -> m-d.                                 (10)
```

Subtracting consecutive dynamical equations gives

```text
-4md=d,                         d(2m+d)=d.               (11)
```

Therefore

```text
m=-1/4,                         d=3/2,
c=m-d-m^2=-29/16.                                     (12)
```

So `(3)` is the unique nondegenerate AP-supported 3-cycle for a centered
monic quadratic over `Q`.

### Marked-cycle chart

The following marked rational 3-cycle chart provides an independent view:

```text
D(t)=2t(t+1),                              t!=0,-1,

p_0=(t^3+2t^2+t+1)/D(t),
p_1=(t^3-t-1)/D(t),
p_2=-(t^3+2t^2+3t+1)/D(t),

c(t)=-(t^6+2t^5+4t^4+8t^3+9t^2+4t+1)
      /(4t^2(t+1)^2).                                    (13)
```

Direct substitution gives `p_0 -> p_1 -> p_2 -> p_0`. The three possible AP
equations factor as

```text
2p_1-p_0-p_2 =  (t-1)(t^2+t+1)/(t(t+1)),
2p_0-p_1-p_2 =  (t+2)(t^2+t+1)/(t(t+1)),
2p_2-p_0-p_1 = -(2t+1)(t^2+t+1)/(t(t+1)).                (14)
```

When `t^2+t+1=0`, all three marked points coincide, so this cyclotomic factor
is a collision rather than an exact cycle. The rational noncollision
representatives are

```text
t in {1,-1/2,-2};                                       (15)
```

they give the same point set and `c=-29/16`, with cyclically shifted
markings. Relabeling acts by

```text
rho(t)=-1/(t+1),
A=[[0,-1],[1,1]],              A^3=-I,        A^6=I.     (16)
```

### Exact `3-4-5` specialization

The AP locus selects `t=1` up to relabeling. Then

```text
(t+1,t)=(2,1),
((t+1)^2-t^2, 2t(t+1), (t+1)^2+t^2)=(3,4,5).             (17)
```

The cycle denominator `D(1)=4` is Euclid's even leg, and the denominator of
`c(1)` is its square. The exact specialized identities are

```text
16=4^2,                 29=5^2+2^2,
7=4^2-3^2,              (-7,-1,5) is an AP.              (18)
```

This is a coordinate-level explanation of the `3-4-5` appearance, not a
claim that Pythagorean triples generate general quadratic cycles.

## 3. Unique projective interpolant and its order-six lift

Three distinct ordered points determine a unique element of `PGL_2`. The
one cyclically interpolating `(3)` is

```text
M(x)=(4x-13)/(16x+12),                                  (19)

M(-7/4)=5/4,        M(5/4)=-1/4,        M(-1/4)=-7/4.
```

Its determinant-one lift is

```text
B=[[1/4,-13/16],[1,3/4]],
det B=1,              tr B=1,
char_B(lambda)=lambda^2-lambda+1,
B^3=-I,               B^6=I.                            (20)
```

Thus `M` has projective order three while `B` has linear order six. This is
the exact central-sign mechanism: passage from `SL_2` to `PGL_2` destroys
`-I`. It is not a rational period-six orbit of `f`.

## 4. The `a=-48` zero fibre and the failed carrier identification

THM-4134's target pencil is

```text
E_q: V^2=U^3-(3a^2/4)U+q-a^3/4.                         (21)
```

At `q=0` it is nodal:

```text
V^2=(U-a)(U+a/2)^2.                                     (22)
```

Its normalization is the exact weight-`(2,3)` map

```text
nu_a(r)=(r^2+a, r(r^2+3a/2)).                           (23)
```

Normalize the AP cycle by `r=4x+1`. Its ordered values become

```text
-6 -> 6 -> 0 -> -6.                                     (24)
```

For `a=-48`, `(23)` sends them to

```text
nu(-6)=(-12,216),
nu( 6)=(-12,-216),
nu( 0)=(-48,0),                                         (25)
```

three rational points on `E_0`. The conjugated projective interpolant is

```text
N(r)=2(r-6)/(r+2),
C=[[1/2,-3],[1/4,1/2]],          C^3=-I.                (26)
```

This proves the weak zero-fibre incidence anticipated by the reservation.
It does not give an automorphism of the nodal cubic. The two conductor
preimages of the node satisfy

```text
r^2=-3a/2=72,                                           (27)
```

whereas the fixed points of `N` satisfy `r^2=-12`. Since an order-three
permutation of a two-point set must be trivial, `N` does not preserve the
conductor pair and cannot descend through the normalization.

More decisively, THM-4134's polynomial horizontal section is `3P` after

```text
q=a^3/2+rho^2,
3P=(a/2+16rho^2/(9a^2),
    -rho-64rho^3/(27a^3)).                              (28)
```

At `a=-48,q=0`, one has `rho^2=-a^3/2`, hence

```text
3P=(56/3,(5/27)rho).                                    (29)
```

Its `U` coordinate differs from both `-12` and `-48` in `(25)`. Therefore

```text
zero-fibre normalization incidence:              PROVED;
identification with THM-4134's horizontal 3P:     REFUTED. (30)
```

The pre-reserved phrase "horizontal-carrier intersection" is narrowed by
`(30)`. No BC image, monodromy, or nonproperness conclusion follows from
`(25)`.

## 5. Finite-ring positive and hostile controls for `T`

All claims in this section are **FINITE-EXACT**.

### Modulo `63`: hostile control

The complete periodic census of `T` on `Z/63Z` is

```text
(5,62,56),              (14,26,20),              (35,47,41). (31)
```

These are three cycles of length three and there is no other cycle. In
particular, `T` has no exact 6-cycle modulo `63`. The local censuses are one
3-cycle modulo `7` and one modulo `9`; CRT gives the three phase matchings in
`(31)`. This is a hostile control against reading the exponent in `2^6-1`
as a period for an unrelated quadratic operation.

### Modulo `43^2`: finite depth-six lift

The multiplier of the rational 3-cycle is

```text
(-7/2)(5/2)(-1/2)=35/8 = -1 (mod 43).                   (32)
```

Let `(b_0,b_1,b_2)=(-7,5,-1)`. For `a_0 in F_43`, Taylor expansion gives

```text
T(b_i+43a_0)=b_(i+1)+43(b_i/2)a_0             (mod 43^2),
T^3(b_i+43a_0)=b_i-43a_0                      (mod 43^2). (33)
```

The three points with `a_0=0` form the old cycle. Each of the other
`3*42=126` points has exact period six: reduction forces a multiple of three,
`(33)` excludes three, and applying it twice gives six. The residue tube is

```text
one 3-cycle + twenty-one 6-cycles.                        (34)
```

The exhaustive census of all `1849` residues confirms that these are every
exact 6-cycle modulo `43^2`. One representative is

```text
36 -> 779 -> 85 -> 1799 -> 1080 -> 1762 -> 36.           (35)
```

For `F=T^3`, the `126` noncentral points form `63` two-cycles; `T` groups
them three at a time into the twenty-one cycles in `(34)`. That `63` is the
count `3(43-1)/2`, not the Mersenne mechanism below. No compatible
`Q_43`-cycle and no rational lift of `(35)` is claimed.

## 6. Mersenne exception and squaring-map period six

This section concerns the different map `S(z)=z^2`. One has

```text
2^6-1=63=3^2*7.                                         (36)
```

The prime `3` already divides `2^2-1`, and `7` already divides `2^3-1`.
Thus exponent six introduces no primitive prime divisor.

Exact period six nevertheless occurs. The nonzero sixth dynatomic factor is

```text
Psi_(S,6)(z)
 = ((z^63-1)(z-1))/((z^7-1)(z^3-1))
 = Phi_9(z) Phi_21(z) Phi_63(z),                         (37)
```

where the right-hand factors are ordinary cyclotomic polynomials. Moreover,

```text
ord_9(2)=6,              ord_21(2)=lcm(2,3)=6,
ord_63(2)=lcm(ord_9(2),ord_7(2))=6.                      (38)
```

Therefore

```text
deg Psi_(S,6)=phi(9)+phi(21)+phi(63)=6+12+36=54,         (39)
```

giving nine exact 6-cycles. Equivalently, exponent doubling
`e |-> 2e (mod 63)` has

```text
one 1-cycle, one 2-cycle, two 3-cycles, nine 6-cycles.   (40)
```

The mechanisms differ by conductor. Conductor `9` uses increased 3-adic
depth; conductor `21` synchronizes local orders two and three by CRT;
conductor `63` combines them. "New depth without new support" is therefore
one mechanism, not a universal explanation.

## 7. Scope, validity gates, and replay

- **PROVED:** the all-`Q` graph and no-period-six result; unique AP cycle;
  both projective lifts; `3-4-5` specialization; zero-fibre incidence and
  horizontal-section separation; the no-new-prime fact; and `(37)`.
- **FINITE-EXACT:** all `63` residues for `T`, all `1849` residues for `T`,
  the `129`-point lift tube, and all `63` exponent classes for doubling.
- **REFUTED:** equality of the points `(25)` with THM-4134's horizontal `3P`
  section. The stronger reserved carrier interpretation does not survive.
- **NOT CLAIMED:** a rational period-six point for another quadratic, a
  nontrivial `Q_43` period-six cycle, or a new Keller-map obstruction.

The source of `(25)` is a one-variable cycle, the target is the normalization
of one nodal member of the THM-4134 pencil, and the map is
`x |-> r=4x+1 |-> nu(r)`. It preserves three-point incidence but destroys the
quadratic dynamics and fails the conductor-pair descent test `(27)`. There is
no proved preserved predicate involving the finite BC carrier. Accordingly,
this theorem neither strengthens THM-4138 nor closes any additional planar
Jacobian case; `JC(2)` and `DC(2)` remain open.

Reproduce with

```bash
python3 04-computation/quadratic_29_16_cycle_audit.py
python3 -O 04-computation/quadratic_29_16_cycle_audit.py
python3 04-computation/quadratic_29_16_cycle_audit_independent.py
python3 -O 04-computation/quadratic_29_16_cycle_audit_independent.py
```

Both normal and optimized replays byte-match the frozen outputs. **QED.**
