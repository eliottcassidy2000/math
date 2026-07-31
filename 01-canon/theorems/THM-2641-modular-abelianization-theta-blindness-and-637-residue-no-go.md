---
id: THM-2641
title: "Modular abelianization theta blindness and the 637-residue no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In PSL2(Z), the
  parabolic T generates the cyclic abelianization
  C6.  For every odd M prime to three, right multiplication by T^(2M)
  preserves the mod-two Farey state and the full level-M residue while
  generating the ternary C3 fibre.  At the live LRC levels M=lcm(49,91)=637,
  T^1274 is therefore invisible to every cited mod-2/mod-49/mod-91 address
  but shifts the C6 class by two.  Conversely the three C3-conjugates of T
  have one abelianized class and exhaust the three mod-two theta directions.
  Thus neither absolute nor relative modular abelianization canonically
  glues a Farey edge to its theta channel.  The missing post-91 spectral
  refinement is modulo 546 (or 3822 with full 637 ancestry), and a Boolean
  finite-step hostile has energy in every mod-91 class but only the zero
  mod-six class.  No LRC row is excluded.
source: root-2026-07-28-modular-connector-probe
depends_on:
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
related:
  - THM-2419-valuation-normalized-homogenization-of-affine-sideband-shells
  - THM-2424-coprime-common-root-crt-and-unit-residue-spectrum
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
script: 04-computation/modular_abelianization_theta_blindness_thm2641.py
output: 05-knowledge/results/modular_abelianization_theta_blindness_thm2641.out
script_sha256: 13340571a6fd59f963b8ad048524b319394a4635662806b8bab705bfeeb2a178
output_sha256: bc3abd6f387e18a8e7363f676a810963a47d21abe988297bbac7c4a182e2f8bf
hash_basis: LF-normalized bytes
---

# THM-2641 -- the current residue skeleton cannot see the ternary lift

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  THM-2632 leaves
one tempting possibility: identify the Farey
six-edge address `{L,R} x {three theta channels}` with the cyclic modular
abelianization `C_6`.  This theorem proves that neither natural
identification survives.  The failure persists after retaining all current
mod-49 and mod-91 LRC address data.

## 1. Free factors and the parabolic generator

Use the standard presentation

```text
Gamma=PSL_2(Z)=<S,C | S^2=C^3=1>,                         (1)

S=[0 -1;1 0],        C=[0 -1;1 1],        T=[1 1;0 1].   (2)
```

In `PSL_2(Z)`, one has `T=SC`.  Therefore

```text
Gamma^ab=C_2 x C_3=C_6,                                  (3)
```

and the class of `T` is a generator.  After choosing an orientation write

```text
chi(T)=epsilon in {1,-1} subset C_6.                      (4)
```

Every conclusion below is invariant under reversing `epsilon`.

## 2. The level-M ternary-invisibility lemma

Let `M` be odd with `3` not dividing `M`.  For arbitrary `G in Gamma` put

```text
G_j=G T^(2Mj),                         j=0,1,2.            (5)
```

Since `T^n=[1 n;0 1]`, one has

```text
G_j=G mod 2,                 G_j=G mod M,                 (6)

chi(G_j)=chi(G)+2M epsilon j.                             (7)
```

The element `2M epsilon` has order three in `C_6`.  Hence (7) exhausts the
entire `C_3` fibre while (6) fixes both the Farey parity state and the
level-`M` state.  Equivalently,

```text
ker(Gamma -> PSL_2(Z/2M)) -> C_3                         (8)
```

is surjective.  Indeed, reduction modulo two followed by the sign of
`S_3=PSL_2(F_2)` is `chi mod 2`, so the kernel lands in the even `C_3`
subgroup; the element `T^(2M)` maps to its generator.  In particular, no
finite tower whose odd modulus has only
the primes `7` and `13`, even after adjoining the mod-two Farey state,
determines the ternary coordinate in (3).

For the live packet levels,

```text
M=lcm(49,91)=637,
2M=1274=26*49=14*91=2 mod 6.                              (9)
```

The exponent `1274=lcm(2,49,91)` is the least positive parabolic shear that
is invisible at all three residue levels.  The nonnegative unimodular flanks

```text
T,                  T^1275,                  T^2549       (10)
```

are identical modulo `2`, `49`, and `91`, while their abelianized classes
are

```text
epsilon,             3 epsilon,               5 epsilon. (11)
```

Thus they have the same mod-two parent and the same `(L/R,theta)` address
for either fixed Farey child, as well as every quotient label that factors
through those mod-49/mod-91 parent residues.  Equation (10) does **not** say
that the exact speeds, Gram data, endpoint measures, or terminal word are
preserved; those are precisely possible extra sidecars.

## 3. Abelianization also forgets theta in the opposite direction

Conjugation is invisible in an abelian quotient.  The three integral
conjugates

```text
T,
C T C^(-1)       =[ 1 0;-1 1],
C^2 T C^(-2)     =[ 2 1;-1 0]                             (12)
```

therefore all have class `epsilon`.  Modulo two they are the three distinct
transvections in `GL_2(F_2)`.  Their unique fixed nonzero directions are,
in one gauge,

```text
(1,0),                    (0,1),                    (1,1), (13)
```

which are exactly the three THM-2632 theta channels.

This proves two complementary no-go statements:

```text
relative edge class chi(T): constant across all theta channels;
absolute parent class chi(G): changes under T^2 although theta does not.  (14)
```

Abelianization forgets the conjugating frame that theta needs.  The two
six-sets in THM-2632 are therefore not canonically isomorphic from modular
data alone.

## 4. Exact spectral state spaces

THM-2424 supplies positive residue energy of the form

```text
E_m=sum_(n=m mod 91) |Fhat(n)|^2.                         (15)
```

To couple it to the full modular lift one would need the strictly finer bank

```text
E_(m,c)=sum_(n=m mod 91, n=c mod 6) |Fhat(n)|^2,          (16)
```

where a single global gauge assigns `c` to the Farey edge.  Since `6` and
`91` are coprime, the post-91 state is

```text
C_6 x Z/91 = Z/546.                                      (17)
```

If the mod-49 ancestry has not already been collapsed into the parent cell,
then

```text
Z/49 x_(Z/7) Z/91 = Z/637,
C_6 x Z/637        = Z/3822.                              (18)
```

If the Farey `C_2` coordinate stays external and only the missing `C_3` is
adjoined, the corresponding moduli are `273` and `1911`.

The distinction is load-bearing.  Let

```text
F(x)=1_{{6x}<1/2}.                                       (19)
```

This is a rational Boolean finite-step function.  Summing its six translated
half-intervals gives, for nonzero frequency,

```text
Fhat(n)!=0  iff  n=6k with k odd.                         (20)
```

Because `6` is invertible modulo `91`, every residue `m mod 91` has an odd
integer `k` with `6k=m mod 91`: if the least residue has even parity, add
`91`.  Hence every bank (15) is positive, while all nonzero Fourier energy
of (19) lies in the single class

```text
c=0 mod 6.                                               (21)
```

Thus all-class mod-91 survival does not imply any prescribed mod-546 class.

## 5. Strongest survivor and next test

A lawful same-parent six-root singleton profile would have nonzero Fourier
magnitude in every `C_6` character.  Its root location, however, is encoded
in phase, not energy.  Even (16) is therefore an address state rather than a
sheet identification.  The cheapest positive experiment is:

1. construct a lawful same-parent six-root profile from the physical
   ancestor;
2. tensor it with THM-2424's mod-91 common-root packet;
3. retain the complex endpoint/reference phase;
4. test every globally admissible dihedral gauge, rejecting the connector if
   one Farey edge receives two `C_6` labels.

The source is the physical integer ancestor; the target is the refined
mod-546 or mod-3822 packet; the desired preserved predicate is one common
sheet with its complex current.  Reduction preserves residue class and
destroys the parent origin and conjugating frame.  The needed sidecar is a
common endpoint/reference phase plus one globally fixed affine
theta-to-`C_3` gauge.

No such profile is constructed here.  No LRC row is excluded, and no
Jacobian, Dixmier, knot, or graceful-tree conjecture follows.

## 6. Exact reproduction

Run

```bash
python 04-computation/modular_abelianization_theta_blindness_thm2641.py
python -O 04-computation/modular_abelianization_theta_blindness_thm2641.py
```

The companion verifies the free-factor matrices, the least level `1274`,
all `83` odd `M<250` with `3` not dividing `M`, the three explicit parents,
the conjugate theta triple, all modulus identities, and one exact frequency
in each of the `91` residue classes of the Boolean hostile.  Normal and
optimized runs must agree after LF normalization and end in
`ALL CHECKS PASSED`.
