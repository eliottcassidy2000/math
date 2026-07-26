---
id: THM-2366
title: "Retained-probe target covariance and sub-thirteen budget bridge"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Retaining
  every THM-2364 blocker-probe shift inside THM-2365's lawful target
  tensor gives an exact dichotomy for the signed fully mixed colour:
  either a nonzero target survives with the same deep/probe colours, or
  its target response is exactly J(s,t)=zeta^(at)J(0,0), the conditional
  inverse-character line. A single rational nonnegative circulant tensor
  saturates the pure and fork THM-2364 floors while remaining entirely
  target-null, so mixed probe colour alone cannot break the obstruction.
  Independently, any lawful co-shift mass budget below the sharp factor
  thirteen forces positive H-drift and hence a canonical nonzero target;
  a transferred two-cover budget would give explicit drift and
  unit-deep target-energy floors. The budget transfer is OPEN. A
  noninverse probe-coloured survivor is a derived u-coloured target
  current; on a pure word its probe can be retyped as the actual factor,
  but the canonical u-sum may still cancel. No
  scalar-row exclusion, ledger decrement, or LRC(14) consequence follows.
source: codex-2026-07-25-retained-probe-target-covariance
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2364-anchored-corner-forces-mixed-deep-blocker-colour
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
related:
  - THM-2361-familywise-fixed-colour-cone-and-offdiagonal-phase-boundary
  - THM-2362-thirteen-shift-successor-statistic-and-role-jet-floor
script: 04-computation/retained_probe_target_covariance_thm2366.py
output: 05-knowledge/results/retained_probe_target_covariance_thm2366.out
script_sha256: c63bf1dc25e345b19f55887fdfab913ec898804a57ea63457d6770e5773a78c8
output_sha256: 5790e5120b2adb3dcff9c15d895be5287f57384ec9b21ea8b3ed95c87c6ab26c
hash_basis: working-tree bytes (LF)
---

# THM-2366 -- retain the probe colours through the lawful target action

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2364 forces one signed coefficient which is nonzero in the deepest
probe colour and in every named blocker-probe colour. THM-2365 supplies a
lawful two-coordinate target action, but its full target sum is an exact
circulation. The correct composition is obtained by retaining all probe
shifts before transforming the target coordinates.

The result has a sharp boundary:

```text
mixed probe colour + failure of inverse covariance
  -> a nonzero derived target fibre retaining every probe colour;

exact inverse covariance
  <- compatible with every sharp THM-2364 marginal.
                                                               (1)
```

A separate positive route avoids phase recovery entirely. If the lawful
target co-shift inherits any total-mass budget strictly below its trivial
thirteen sheets, the circulant branch is impossible.

## 1. The retained-probe lawful tensor

Put

```text
p=13,

zeta=exp(2*pi*i/p).
```

Use one THM-2365 delayed word and its lawful target translates
`F_(s,t)`, with

```text
r,s,t in F_p.
```

The deepest danger shift is

```text
Delta_r(x)=D_(c_3,r)(x).
```

Let the nonempty terminal blocker word have labels

```text
sigma={alpha_1,...,alpha_d},

d in {1,2},
```

and clock `R`, where `p|R`. For

```text
v=(v_1,...,v_d) in F_p^d
```

define

```text
K(r,s,t,v)
 =integral_T
    1_(F_(s,t))(x) Delta_r(x)
    product_(j=1)^d D_(R c_(alpha_j),v_j)(x) dx.    (2)
```

This definition is made first as a finite nonnegative table of overlap
integrals. For a pure word (`d=1`), later scalar-cover redundancy can
retype its probe as the actual named-factor physical response. For a
fork (`d=2`), it is a coefficient-derived duplicate-probe overlap bank,
not a canonically available physical measurement. The deepest complement
in `F_(s,t)` gives, pointwise almost everywhere,

```text
K(t,s,t,v)=0                                      (3)
```

for every `s,t,v`, up to the null strict-open endpoints.

At `(s,t)=(0,0)`, the delayed word already contains every anchored
blocker factor in (2). Thus

```text
K(r,0,0,0)=H_(2365)(r,0,0),
```

and the full `(r,v)` slice is exactly the THM-2364 deep/blocker-probe
table.

## 2. Full transform and lawful target

Use the normalized finite transform

```text
Khat(a,b,h,u)
 =p^(-(d+3))
  sum_(r,s,t,v)
   K(r,s,t,v)
   zeta^(a r+b s+h t+u.v).                         (4)
```

Here

```text
a,b,h in F_p,

u=(u_1,...,u_d) in F_p^d.
```

The endpoint/deep part is exactly THM-2365's transform. Every added
probe multiplier contributes through a speed `R c_(alpha_j)`, and
`R=0 mod p`. Therefore the probes are target-neutral and the actual
target vector remains

```text
q=(b,a+h) in F_p^2.                                (5)
```

The diagonal zero (3) gives the retained-probe target-line cancellation

```text
sum_(a in F_p) Khat(a,b,q_2-a,u)=0                 (6)
```

for every `b,q_2,u`. Indeed, summing `a` in (4) produces
`p 1_(r=t)`, where (3) kills the table.

## 3. The exact conditional inverse-character line

Fix a deep colour `a` and probe tuple `u`, and collapse only those
coordinates:

```text
J_(a,u)(s,t)
 =p^(-(d+1))
  sum_(r,v)
   K(r,s,t,v)zeta^(a r+u.v).                       (7)
```

Then

```text
Khat(a,b,h,u)
 =p^(-2)sum_(s,t)
   J_(a,u)(s,t)zeta^(b s+h t).                     (8)
```

THM-2364 supplies at least one tuple

```text
a!=0,

u_j!=0 for every j,
```

such that, if `rho=mu(F_(0,0))`,

```text
Re J_(a,u)(0,0)
 <=-11^d rho/(13^(d+1)12^(d+1))
 <0.                                               (9)
```

Choose one such tuple once and for all. Fourier inversion in `(s,t)`
gives the exact equivalence

```text
Khat(a,b,h,u)=0 for every (b,a+h)!=(0,0)

iff

J_(a,u)(s,t)=zeta^(a t)J_(a,u)(0,0)
             for every s,t.                        (10)
```

Thus exactly one of the following occurs:

```text
noninverse branch:
  some Khat(a,b,h,u)!=0 has q=(b,a+h)!=0;

inverse branch:
  the complete response is the nonzero conditional character
  J(s,t)=zeta^(a t)J(0,0).                          (11)
```

The two cheapest sufficient tests for the first branch are

```text
J_(a,u)(1,0)!=J_(a,u)(0,0),

or

J_(a,u)(0,1)!=zeta^a J_(a,u)(0,0).                 (12)
```

They are not necessary; any failure of (10) is decisive.

### 3a. Exact nonzero-target energy

Put

```text
c_0=Khat(a,0,-a,u),

E_nz
 =sum_((b,a+h)!=(0,0))|Khat(a,b,h,u)|^2,

Delta
 =13^(-2)sum_(s,t)
   |J_(a,u)(s,t)-zeta^(a t)J_(a,u)(0,0)|^2.       (12a)
```

Parseval and orthogonal projection onto the one inverse character give

```text
E_nz
 =13^(-2)sum_(s,t)|J_(a,u)(s,t)|^2-|c_0|^2,

Delta=E_nz+|c_0-J_(a,u)(0,0)|^2.                  (12b)
```

The last difference is the value at `(0,0)` of a sum of `168`
nonzero-target characters. Cauchy--Schwarz therefore yields

```text
Delta/13^2<=E_nz<=Delta,                           (12c)

max_((b,a+h)!=(0,0))|Khat(a,b,h,u)|
 >=sqrt(Delta)/(13 sqrt(168)).                     (12d)
```

Both energy constants are sharp. The upper equality occurs when the
noninverse part vanishes at `(0,0)`; the lower equality occurs when all
`168` nonzero-target coefficients are equal as complex numbers.

### 3b. The fixed-colour envelope criterion

Because `K` is nonnegative, define

```text
M(s,t)
 =13^(-(d+1))sum_(r,v)K(r,s,t,v)
 >=|J_(a,u)(s,t)|.                                 (12e)
```

If the inverse branch holds, then for every fixed `s_0`,

```text
sum_t M(s_0,t)
 >=13|J_(a,u)(0,0)|.                               (12f)
```

Write the inherited THM-2364 floor as

```text
kappa_d
 =11^d/(13^(d+1)12^(d+1)).
```

Equation (9) gives `|J(0,0)|>=kappa_d rho`. Hence either of

```text
sum_t M(s_0,t)<13|J_(a,u)(0,0)|,

sum_t M(s_0,t)<13 kappa_d rho                     (12g)
```

forces the noninverse branch and the energy floor (12c). The second is a
convenient sufficient estimate, not an equivalence.

This fixed-colour criterion is different from Section 6's uncoloured
`C<13` budget. Merely bounding (12e) by `C rho` or `C M(0,0)` with
`C<13` does not contradict covariance: the known mixed coefficient is
only a `kappa_d` fraction of the word mass.

## 4. What a noninverse coefficient lands

For fixed `u`, first form the finite probe projector

```text
P_u(x)
 =13^(-d)sum_v
   product_j D_(R c_(alpha_j),v_j)(x)zeta^(u.v).   (12h)
```

It is a bounded BV function. After this finite projection, the
`m`-then-`X` sums are ordinary absolutely convergent sums. Only the
internal danger-probe expansion defining `P_u` remains Abel: use one
joint Poisson-Abel smoothing there, then bounded-product `L1`
convergence recovers the finite overlap table (2) at the indicator
boundary.

After this finite `u`-projection, THM-2365's absolutely convergent
deep-multiplier and endpoint-frequency extraction applies. A noninverse
coefficient in (11) yields exact integers `m,X` with

```text
gcd(m,91)=1
```

and a nonzero derived fixed-triangle target fibre

```text
A^(u)_(X,m)(q)!=0,

q!=0.                                               (13)
```

The deep multiplier is a thirteen-unit because `m=a mod 13`, and every
live danger coefficient kills nonzero multiples of seven. Each displayed
probe residue is likewise nonzero modulo thirteen, and its live Abel
multipliers are separately prime to `91`.

For a pure word, scalar-cover redundancy permits the one probe colour to
be retyped as the actual named-factor response. For a fork, `u` records
two auxiliary duplicate-blocker probes in a coefficient-derived overlap
bank, not a canonically available physical measurement. It is not a
Bockstein, relation, or canonical word colour.

At the indicator boundary, finite inversion at `v=0` sums all
`u`-coloured currents back to the canonical uncoloured current: every
extra anchored copy collapses by `D^2=D`. This is not an identity at
finite Poisson radius. Nonzero `u`-fibres may cancel in that boundary
sum. Therefore (13) is **not** a proof that THM-2334's canonical
`A_(X,m)(q)` is nonzero. A lawful `u`-selector or a common cone would be
needed. THM-2361 does not supply it: these endpoint currents are
off-diagonal and retain its exact terminal-phase hostile.

## 5. One sharp compatibility tensor

The inverse branch is not a formal possibility only. Put

```text
mu=6/91,

G(0)=0,

G(r)=1/182                 for r!=0,

b(0)=1,

b(v)=1/12                  for v!=0.                (14)
```

Define the rational nonnegative tables

```text
H(r,s,t)=G(r-t),

K_d(r,s,t,v_1,...,v_d)
 =G(r-t)product_(j=1)^d b(v_j).                    (15)
```

They obey the diagonal zero, are independent of `s`, and every
`r`-row sum of `H` is `mu`. Their H-drift is zero. Nevertheless, for
every `a!=0`,

```text
C(a)=13^(-1)sum_r G(r)zeta^(a r)
    =-1/2366.                                       (16)
```

For every fully mixed tuple, the base coefficients are

```text
d=1:
  J_(a,u)(0,0)=-11/369096,

  S_1=sum_(a,u!=0)J_(a,u)(0,0)
     =-66/15379;

d=2:
  J_(a,u_1,u_2)(0,0)=-121/57578976,

  S_2=sum_(a,u_1,u_2!=0)J_(a,u_1,u_2)(0,0)
     =-726/199927.                                  (17)
```

These are exactly the sharp THM-2364 coefficient and corner values at
mass `mu`. At the same time,

```text
J_(a,u)(s,t)=zeta^(a t)J_(a,u)(0,0),               (18)
```

so every one of them remains on the zero-target inverse line.

The table is an abstract finite-tensor hostile, not a claim that one
canonical LRC row realizes it. It proves that nonnegativity, diagonal
zero, every sharp mixed marginal, and the complete target circulation
are mutually compatible.

There is also a quantitative obstruction to a naive positive
intertwiner. THM-2365's target action has the two generators

```text
T_(delta_eta,delta_ell)H(r,s,t)
 =H(r+delta_ell,s+delta_eta,t+delta_ell).
```

The residual `G(r-t)` is constant along both. An eta-visible delayed
blocker would compare hypothetically with the `delta_eta`-generator; a
delayed `c_3` blocker would compare with the diagonal
`delta_ell`-generator. Thus either pure role has probe-profile `l1` mass
`2`, while its corresponding lawful target orbit has mass `13`. Any
positive intertwiner carrying the hostile probe action to that
target-action generator while preserving a nonzero anchor requires
amplification at least

```text
13/2.
```

For a fork, the product probe action has mass `4`, while the full
`(delta_eta,delta_ell)` target orbit has mass `169`, requiring
amplification at least

```text
169/4.                                             (19)
```

Thus no Markov, sub-Markov, or mass-nonincreasing positive map can
identify the corresponding probe and target actions on a nonzero hostile
tensor. This is a norm obstruction to a **hypothetical** intertwiner, not
the assertion that canon already supplies such a map.

## 6. A sub-thirteen mass budget forces canonical target drift

Return to THM-2365's canonical uncoloured tensor `H` and put

```text
G(r)=H(r,0,0),

M_s=sum_r H(r,s,0),

M_0=sum_r G(r)>0.                                  (20)
```

If `D_H=0`, THM-2365 gives

```text
H(r,s,t)=G(r-t),
```

and consequently

```text
M_s=M_0                  for every s.              (21)
```

Therefore each of the following is a sufficient canonical target test:

```text
some M_s!=M_0;

some nonzero Fourier mode of s -> M_s;

sum_(s in S)M_s<|S|M_0 for some nonempty S;

sum_s M_s<=C M_0 for any C<13;

sum_s M_s<13M_0.                                  (22)
```

The non-strict threshold is exactly `C<13`; the strict inequality at
`C=13` is also decisive. The hostile (15) has equality at `13`, so the
threshold is sharp. This is an abstract lawful-`H` criterion valid for
every word. A direct comparison with a particular delayed blocker role
requires matching that blocker to the correct target-action generator.

For an eta-visible delayed blocker, the especially attractive missing
transfer is the one-factor thirteen-shift cover budget

```text
sum_s H(r,s,0)<=2H(r,0,0)          for every r.    (23)
```

It would imply (22) with `C=2`, hence `D_H>0` and THM-2365's canonical
nonzero target landing. For a delayed `c_3` role the corresponding
hypothetical comparison uses the diagonal target orbit

```text
sum_(delta_ell)
 H(r+delta_ell,0,delta_ell)<=2H(r,0,0).            (23a)
```

For a fork, the full product comparison is

```text
sum_(delta_eta,delta_ell)
 H(r+delta_ell,delta_eta,delta_ell)
 <=4H(r,0,0).                                      (23b)
```

Each residual orbit has respectively `13` or `169` identical values, so
these inequalities would likewise kill it. None is presently canonical.

The eta-axis implication is quantitative already on the thirteen masses.
Put

```text
S=sum_s M_s,

Mtilde(b)=13^(-1)sum_s M_s zeta^(bs).
```

For fixed `M_0>0` and total mass `S`, nonnegativity and Cauchy--Schwarz
on the other twelve entries give the sharp variance bound

```text
sum_(b!=0)|Mtilde(b)|^2
 >=(13M_0-S)^2/2028
 >=(13-C)^2M_0^2/2028,                            (23c)

max_(b!=0)|Mtilde(b)|
 >=(13-C)M_0/156.                                 (23d)
```

Equality in the first inequality occurs exactly when the other twelve
masses are equal. If `B(a,b,h)` is THM-2365's full transform, then

```text
Mtilde(b)=13sum_h B(0,b,h).
```

Cauchy--Schwarz in `h` loses at most `13^3`, so (23c) already implies
the H-drift bound below.

Equivalently, define the lifted functional

```text
Lambda(H)
 =sum_(r,s)H(r,s,0)-13sum_rH(r,0,0).
```

The coefficient vector of `Lambda` has squared Euclidean norm

```text
13*12^2+13*12=2028=12*13^2,                        (24)
```

and is orthogonal to every circulant table `g(r-t)`. This gives the same
estimate directly in THM-2365's normalized counting measure:

```text
D_H
 >=(13-C)^2 M_0^2/(12*13^5)                        (25)
```

under the full-axis budget in (22). Its nonzero-target energy carried by
deep colours `a!=0` is at least

```text
(13-C)^2 M_0^2/(12*13^6).                          (26)
```

For the hoped-for `C=2` transfer these become

```text
D_H>=121 M_0^2/4455516,

unit-deep target energy
 >=121 M_0^2/57921708.                             (27)
```

No per-`r` estimate is required: any aggregate certificate in (22) is
enough.

## 7. Why the budget transfer is still open

For an added blocker probe, the elementary root count gives an `l1`
budget of two. Its matched lawful THM-2365 target generator is different:
the eta-visible role uses the `delta_eta/s` action, while a delayed
`c_3` role uses the diagonal
`delta_ell:(r,t)->(r+delta_ell,t+delta_ell)` action. Either generator
co-shifts the complete present endpoint packet through a quotient-dual
direction. Other present roles move simultaneously. Their intersections
can expand even when the named danger alone has a two-cover count.

Thus (23) is not supplied by THM-2362 or THM-2364. The hostile
amplification factors (19) show that no mass-preserving positive
identification of the two shift systems can prove it. One must instead
exploit canonical overlap geometry, a signed/charged intertwiner, or
directly compute a nonconstant `M_s`.

For the eta axis, the cheapest rowwise object is only thirteen rational
masses:

```text
(M_s)_(s in F_13).                                 (28)
```

A single unequal entry, nonzero difference, or nonconstant finite Fourier
mode settles the THM-2365 branch without recovering terminal component
phase. For a pure `c_3` role, use instead the thirteen diagonal-orbit
masses

```text
N_(delta_ell)
 =sum_r H(r+delta_ell,0,delta_ell).
```

For a fork, the corresponding finite test is the full `13^2`
`(delta_eta,delta_ell)`-orbit. These are alternative lawful tests, not
identifications with the probe shifts.

## 8. Scope on the 165-row frontier

On every first-depth-one row, THM-2364 supplies the tuple in (9) and
THM-2365 supplies the lawful target tensor. Hence the retained-probe
dichotomy (11) and the compatibility boundary (15)--(19) are uniform on
the proof architecture.

What is not yet uniform is the positive side of (22). Canon does not
prove a lawful co-shift budget below thirteen, a nonconstant `M_s`, or a
selector preventing cancellation among probe colours. Accordingly:

```text
noninverse K/J branch
  -> a nonzero derived u-coloured target fibre;

sub-thirteen H budget
  -> a nonzero canonical target fibre;

current canon
  -> neither disjunct is forced on every row.       (29)
```

No scalar row is excluded, the ledger remains `165`, and LRC(14) remains
open.

## 9. Exact companion

The dependency-free companion uses exact cyclotomic and `Fraction`
arithmetic to:

- enumerate the `2,197` target-label decompositions and verify the
  diagonal-zero/root-orthogonality mechanism behind the `169` target
  lines;
- verify two-dimensional inversion on one rational hostile table and
  inverse-character support for the representative deep colour `a=5`;
- check the coefficient-space equality patterns attaining both bounds
  in (12c), together with the hostile `kappa_d` scale;
- check all twelve hostile deep colours and both pure/fork coefficient
  and corner values;
- check the arithmetic amplification ratios `13/2` and `169/4`;
- enumerate the mass-functional coefficients, their circulant
  annihilation, and squared norm `2028`;
- check the `C=2` specializations of (23c)--(27); and
- sweep separate live probe/deep multipliers through their `91`-unit
  boundary.

Run

```bash
python3 04-computation/retained_probe_target_covariance_thm2366.py
python3 -O 04-computation/retained_probe_target_covariance_thm2366.py
```

Both transcripts must match

```text
05-knowledge/results/retained_probe_target_covariance_thm2366.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Three independent hostile passes rederived the transform and covariance
identities, energy bounds, envelope, compatibility tensor, sharp mass
taxes, quantitative functional, convergence typing, and physical-scope
boundary; normal, optimized, and stored transcripts agree exactly. QED.
