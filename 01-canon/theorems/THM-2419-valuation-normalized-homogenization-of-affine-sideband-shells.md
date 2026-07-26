---
id: THM-2419
title: "Valuation-normalized homogenization of affine sideband shells"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For a
  primitive integer speed row w and nonzero sideband X, the affine
  address shell Lambda_X={a:a.w=X} has the canonical finite quotient
  Lambda_X/|X|Lambda, naturally bijective with ker(w mod |X|). When
  its finite-fibre Abel boundary limits exist, a nonzero Abel sideband
  has a surviving finite kernel fibre without choosing an origin. If
  13|X and every address has all-unit
  residue q, that fibre reduces to q. Adjoining the artificial observer
  speed X homogenizes every address with coefficient -1; valuation
  normalization gives the alternative pair (13^d,-X/13^d). A Bezout
  section contracts to an all-unit mod-13 physical relation, but changes
  exact labels by XLambda and need not preserve current amplitude.
  Self-differences lose charge; a same-X residue-zero reference is the
  sharp sufficient sidecar. No physical observer, terminal transport,
  row exclusion, or LRC(14) conclusion is claimed.
source: codex-2026-07-26-affine-sideband-homogenization
depends_on:
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2416-zero-current-or-bounded-sideband-prony-dichotomy
related:
  - THM-2284-thirteen-adic-anchored-rank-three-plucker-lift
  - THM-2325-prescribed-target-gain-full-lattice-91-unit-needle-bank
  - THM-2408-endpoint-prony-resultant-clock-separation-and-shared-node-boundary
  - THM-2412-delta-exponential-and-central-newton-layer-split
  - THM-2413-prime-index-affine-drift-and-twin-center-weld
  - THM-2418-alternating-base-thirteen-septimal-carry-matrix-and-rank-one-boundary
script: 04-computation/lrc14_affine_sideband_homogenization_thm2419.py
output: 05-knowledge/results/lrc14_affine_sideband_homogenization_thm2419.out
script_sha256: 7928b82df194f8ce2508017205f431860091f63170867ffd5f705884b251f8a0
output_sha256: beba29d94110c4750fd7ae792c0f835fbf7c677ff595a04a9987f3dc68d2b5be
hash_basis: working-tree bytes (LF)
---

# THM-2419 -- homogenizing a sideband without losing its unit residue

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2416 turns a vanished full-coordinate mean into a nonzero physical
sideband.  Its atomic addresses do not lie in the relation lattice; they
lie in an affine shell

```text
Lambda_X={a in Z^n:a.w=X},             X in 13Z\{0}.          (1)
```

The right replacement for a prematurely chosen relation is the finite
base-free quotient of this shell.  It retains the all-unit residue and
the Abel mass.  Choosing a Bezout origin then produces an ordinary
relation algebraically, but no theorem says that this chosen relation
inherits the sideband coefficient.

## 1. The affine shell is a finite lattice torsor

Let

```text
w=(w_1,...,w_n) in Z^n,           gcd(w_1,...,w_n)=1,

Lambda={r in Z^n:r.w=0}.                                  (2)
```

Fix `X!=0`, put `M=|X|`, and suppose `Lambda_X` in (1) is nonempty.
The lattice `M Lambda` acts freely on `Lambda_X` by translation.  Define

```text
K_M={u in (Z/MZ)^n:u.w=0 mod M}.                           (3)
```

Reduction modulo `M` gives a canonical map of torsors

```text
pi_X:Lambda_X/(M Lambda) -> K_M.                           (4)
```

It is a bijection.  Injectivity is immediate: if `a-a'=M b`, then

```text
0=(a-a').w=M(b.w),
```

so `b in Lambda`.  For surjectivity, choose a Bezout vector

```text
z.w=1.                                                     (5)
```

Given `u in K_M`, choose an integer lift `u_0`.  Write

```text
u_0.w=M t.
```

If `X=sigma M`, `sigma in {+1,-1}`, then

```text
a=u_0+M(sigma-t)z                                         (6)
```

satisfies `a.w=X` and reduces to `u`.  Therefore

```text
Lambda_X/(M Lambda) ~= K_M,                               (7)
```

independently of the temporary choice (5).  In particular

```text
|K_M|=M^(n-1).                                            (8)
```

## 2. Abel mass survives finite pushforward

At Abel radius `0<rho<1`, let

```text
S_rho=sum_(a in Lambda_X)c_rho(a)                         (9)
```

be absolutely convergent.  Push it through (4):

```text
C_rho(u)
 =sum_(a in Lambda_X, a=u mod M)c_rho(a),       u in K_M. (10)
```

Because `K_M` is finite, (10) is equivalently recovered by finite Fourier
inversion from the character twists of `S_rho`.  Suppose the boundary
limits `C(u)=lim_(rho->1-)C_rho(u)` exist and

```text
sum_(u in K_M)C(u)
 =lim_(rho->1-)S_rho!=0.                                 (11)
```

Then some finite fibre satisfies

```text
C(u)!=0.                                                 (12)
```

No endpoint, Bezout vector, or atomic representative has been selected.

In the THM-2334/THM-2416 application, the `K_M`-character twists are the
same bounded-product Abel currents with a finite root-of-unity weight.
The boundary argument of THM-2334 applies to every one of these finitely
many twists, and finite Fourier inversion then supplies the individual
limits assumed above.  Absolute convergence at a fixed `rho<1` alone
would not imply those boundary limits in an arbitrary series.

For THM-2416, `M` is divisible by thirteen and every atomic address in
the selected sideband has one residue

```text
a=q mod 13,                 q in (F_13^*)^n.             (13)
```

Thus every surviving `u` in (12) reduces to `q` under

```text
K_M -> K_13.                                               (14)
```

The canonical survivor is a charged finite kernel fibre, not an exact
integer relation address.

## 3. Two exact homogenizations

Every `a in Lambda_X` has the tautological observer homogenization

```text
widetilde w=(w,X),             widetilde a=(a,-1),

widetilde a.widetilde w=0.                                (15)
```

When `13|X`, (13) implies

```text
widetilde a mod13=(q,-1) in (F_13^*)^(n+1).                (16)
```

Thus the extended ten-coordinate lattice has a genuine homogeneous
all-unit relation.

There is also a valuation-normalized form.  Write

```text
d=nu_13(X)>=1,              Y_d=X/13^d,        13 does not divide Y_d.
                                                                  (17)
```

Then

```text
widehat w=(w,13^d),          widehat a=(a,-Y_d),

widehat a.widehat w=0,       widehat a mod13=(q,-Y_d).      (18)
```

Among pure powers `13^e` used as the extra speed, `e=d` is the unique
choice for which the compensating coefficient `-X/13^e` is both integral
and a unit modulo thirteen: `e<d` leaves a zero residue and `e>d` is not
integral.

Both (15) and (18) are algebraic observers.  Neither extra speed is a
factor in the original physical packet.

## 4. Bezout contraction is exact but not canonical

Choose `z` as in (5).  Then

```text
r_z=a-Xz in Lambda,                                      (19)

r_z=q mod13.                                             (20)
```

Equation (20) uses the additional hypotheses `13|X` and (13).  Under
those hypotheses every affine atom algebraically determines an all-unit
mod-thirteen physical relation after a section is chosen.  If

```text
z'=z+lambda,              lambda in Lambda,
```

then

```text
r_(z')=r_z-X lambda.                                     (21)
```

Consequently only the class

```text
[r_z] in Lambda/(M Lambda)                               (22)
```

is base-free.  It is another model of the torsor point (7).

The smallest hostile already occurs at

```text
w=(1,1),            X=13,           a=(1,12).
```

The sections `z=(1,0)` and `z'=(0,1)` give

```text
r_z=(-12,12),        r_(z')=(1,-1),                     (23)
```

differing by `13(1,-1)`.  Both reduce to `(1,12)` modulo thirteen,
but they are different exact labels.  Algebraic existence of (19) does
not transport the coefficient `c_rho(a)` to either chosen label.

## 5. The sharp reference sidecar

Two addresses on the same shell have a canonical difference:

```text
a,b in Lambda_X        =>        a-b in Lambda.             (24)
```

Self-correlation of the THM-2416 packet uses two addresses with the same
residue `q`, so

```text
a-b=0 mod13.                                               (25)
```

It restores homogeneity but erases the charged residue.

By contrast, suppose a compatible same-`X` reference packet supplies

```text
b.w=X,                 b=0 mod13.                            (26)
```

Then

```text
a-b in Lambda,         a-b=q mod13.                          (27)
```

Thus a same-shell residue-zero reference is sufficient to contract the
affine current while retaining its full mod-thirteen charge.  Equations
(25)--(27) are sharp: without labelled cross-packet phase, a difference
inside one residue class is necessarily neutral.

This source--reference composition is formally analogous to the
operation-cospan witness composition in THM-2413: the hidden
affine witnesses compose only after a distinct compatible witness is
retained.  The analogy is explanatory, not a dependency or an LRC
consequence.

## 6. The remaining mod-seven coordinate

Write `X=13Y_1`, where `Y_1=X/13`; this is distinct from the
valuation-normalized unit `Y_d=X/13^d` in (17) when `d>1`.

- If `7|Y_1`, then `91|X`, so the canonical finite fibre in `K_M` reduces
  to a fibre in `K_91`.  This does not make its coordinates units
  modulo seven.
- If `7` does not divide `Y_1`, the scalar `Y_1 mod7` is the exact
  one-dimensional affine defect.  The observer coordinate in (15) or
  (18) cancels the equivalent nonzero residue of `X` algebraically, but
  no physical Boolean factor realizes that observer.

This is the same carry-scale boundary represented by the rank-one
projector plus charged leakage in candidate THM-2418.  That reflection
analogy, and the fixed-layer/signed-leakage analogy with THM-2412, are
cross-frontier explanations only; none of these related theorems is
used in the proof.

## 7. Exact scope

This theorem proves:

- the canonical finite torsor/kernel identification (7);
- finite Abel pushforward and a surviving charged kernel fibre;
- two all-unit mod-thirteen observer homogenizations;
- exact but section-dependent physical contraction;
- neutral self-differences and the sufficient residue-zero reference
  sidecar; and
- the `7|Y_1` CRT split.

It does not prove:

- that the artificial observer is a legal LRC factor;
- that a Bezout-selected physical relation inherits sideband amplitude;
- existence of a compatible residue-zero reference packet;
- mod-seven coordinate units, terminal phase transport, a row exclusion,
  or LRC(14).

The scalar ledger remains `165`.

## 8. Exact companion

The dependency-free companion:

- exhausts small primitive two-coordinate rows and verifies the bijection
  `Lambda_X/(M Lambda) ~= K_M`;
- checks both homogenizations and uniqueness of the valuation-normalized
  pure-power observer;
- verifies the hostile (23), section ambiguity (21), neutral
  self-difference, and charged reference difference;
- checks finite pushforward mass and the exact kernel size; and
- audits both branches of the mod-seven split.

Run:

```bash
python3 04-computation/lrc14_affine_sideband_homogenization_thm2419.py
python3 -O 04-computation/lrc14_affine_sideband_homogenization_thm2419.py
```

Both commands must reproduce

```text
05-knowledge/results/lrc14_affine_sideband_homogenization_thm2419.out
```

with the LF hashes in the frontmatter.
