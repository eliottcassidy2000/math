---
id: THM-2579
title: "Socle-flat target torsor and integral difference filling"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED (two ways).
  On the integral 13-colour augmentation lattice, every circulant filter
  acts on THM-2571's Cayley cokernel only through the sum of its coefficients
  modulo 13.  Translations fix the class and the punctured twelve-shift orbit
  negates it (hence preserves nontriviality), while every shift difference,
  the full thirteen-shift orbit, and one extra global factor 13 kill it.  For
  each nonzero owner colour, the
  thirteen canonical target-colour carry profiles are one constant nonzero
  cokernel coset: every pairwise target difference, every integral
  target-augmentation-zero combination, and every unnormalized full target
  Fourier numerator has an integral Cayley primitive although each absolute
  profile is obstructed.  Division by 13 is outside this integral functor:
  after global primitive reduction, a normalized DFT can restore the class.
  This proves that integral relative data lacks an absolute reference while
  leaving the normalized Abel-normal route genuinely open.  The result
  is integral coefficient algebra; it supplies no positive filling,
  semantic endpoint, scalar-row exclusion, or LRC(14) conclusion.
source: root-holotopy-2026-07-28-socle-flat-torsor
depends_on:
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
  - THM-2571-deep-colour-cayley-filling-bockstein-and-norm-curvature-split
related:
  - THM-2567-deep-coloured-duty-replica-cycle-and-augmentation-cancellation
  - THM-2572-deep-augmentation-parseval-energy-and-nonlinear-holotopy-obstruction
  - THM-2573-logarithmic-abel-normal-and-common-endpoint-jump-pairing
  - THM-2574-oriented-tooth-component-holonomy-and-fixed-frequency-descent
script: 04-computation/lrc14_socle_flat_target_torsor_thm2579.py
output: 05-knowledge/results/lrc14_socle_flat_target_torsor_thm2579.out
script_sha256: 5392102223d65c0156557224f081152cb3d5205e22681731935f576a1b5fe8fb
output_sha256: 3954af663ac162f8ec194319cd05f5509d3e4d0c73ded300921813ece0a63545
hash_basis: LF-normalized bytes
---

# THM-2579 -- target differences fill, but the absolute class does not

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED (two ways).**

THM-2571 identifies a nonzero integral first-moment class on every target
colour of the canonical digit-diagonal carry tensor.  Surprisingly, this
does not give thirteen independent obstructions.  The 13-cyclotomic socle
makes them all the **same** obstruction.

The exact shape is

```text
thirteen absolute target profiles:       one nonzero cokernel coset;

every target difference/integral DFT numerator:  a Cayley boundary.
                                                                    (1)
```

Thus difference data can be complete while the absolute charged reference
is still missing.  This is the integral counterpart of the recurring
endpoint-phase boundary: a relative atlas does not choose its origin.

## 1. Circulant filters act only through augmentation

Let

```text
M=Z[zeta_7,zeta_13],

Lambda_M={a=(a_m)_(m in F_13) in M^13 : sum_m a_m=0}.       (2)
```

Let `P` be cyclic translation in the colour index.  For
`f=(f_0,...,f_12) in M^13`, define

```text
T_f=sum_(d=0)^12 f_d P^d,

epsilon(f)=sum_d f_d.                                      (3)
```

Write

```text
beta(a)=sum_(m=0)^12 m a_m                 in M/13M.       (4)
```

Then for every `a in Lambda_M`,

```text
beta(T_f a)=epsilon(f) beta(a).                             (5)
```

### Proof

For every shift `d`,

```text
beta(P^d a)=beta(a)-d sum_m a_m=beta(a)       mod 13M.     (6)
```

Linearity gives (5). QED.

THM-2571 proves the exact sequence

```text
0 -> Lambda_M --C_tau--> Lambda_M --beta--> M/13M -> 0.    (7)
```

Consequently the action induced by the entire integral circulant algebra on
the Cayley cokernel is the one-dimensional augmentation representation:

```text
[T_f a]=epsilon(f)[a].                                     (8)
```

This immediately gives the exact filter table

```text
P^d:                         preserves the class;

P^d-P^e:                     kills the class;

I+P+...+P^12:                kills the class;

P+P^2+...+P^12:              multiplies it by 12=-1;

C_tau:                       kills the class.              (9)
```

The full orbit has augmentation `13=0 mod 13`; the punctured orbit has
augmentation `12`, a unit.  The last line is the quotient form of the fact
that `C_tau` produces boundaries.

Equation (5) is valid for arbitrary integral coefficients, not only Boolean
orbit banks.  It also shows why an unweighted rotation census cannot be
silently substituted for a carrier-weighted filter: the physical operation
must first be proved to be the same circulant `T_f` on one common lattice.

## 2. The canonical target family is one nonzero coset

Use THM-2571's primitive digit-diagonal carry cycles

```text
c^(kappa,b)=(c_m^(kappa,b))_(m in F_13),

kappa in F_7^*,                  b in F_13.                 (10)
```

They satisfy

```text
beta(c^(kappa,b))
 =Omega Y(zeta_7^kappa) !=0,                               (11)
```

where

```text
Omega=sum_(j=0)^11 (j+1)zeta_13^j,

Y(z)=6z+5z^2+z^3+12z^4+8z^5+7z^6.                        (12)
```

The right side of (11) is independent of `b`; `Y` is a unit in the
septimal factor, while `Omega` is the nonzero 13-cyclotomic socle.  Therefore,
for fixed `kappa`,

```text
[c^(kappa,0)]=[c^(kappa,1)]=...=[c^(kappa,12)] !=0

                  in Lambda_M/C_tau Lambda_M.              (13)
```

In particular, for every `b!=b'`,

```text
c^(kappa,b)-c^(kappa,b') in C_tau Lambda_M,                (14)

B_tau(c^(kappa,b)-c^(kappa,b')) in Lambda_M.               (15)
```

Thus every one of the `6*C(13,2)=468` canonical pairwise target differences
has an **integral** sawtooth primitive, even though all `78` absolute
profiles have nonintegral class.

For integer coefficients `n_b`, (11) gives the sharper criterion

```text
beta(sum_b n_b c^(kappa,b))
 =(sum_b n_b) Omega Y(zeta_7^kappa).                        (16)
```

Because a nonzero scalar in `F_13` is a unit, the combination is integrally
fillable exactly when

```text
sum_b n_b=0 mod 13.                                       (17)
```

For coefficients in the full cyclotomic ring, zero target augmentation is
still sufficient; additional annihilation can occur because `Omega` is a
socle element, not a unit.

## 3. Unnormalized Fourier numerators fill; normalization can restore

For `q in F_13`, form the **unnormalized integral numerator**

```text
d^(kappa,q)=sum_(b in F_13)
  zeta_13^(qb)c^(kappa,b).                                 (18)
```

When `q!=0`, its coefficient sum is zero exactly.  When `q=0`, its
coefficient sum is `13`.  Hence (11) gives in both cases

```text
beta(d^(kappa,q))=0,

d^(kappa,q)=C_tau B_tau d^(kappa,q),

B_tau d^(kappa,q) in Lambda_M.                             (19)
```

All `6*13=78` unnormalized target Fourier numerators are therefore
integrally fillable.  This does not say they vanish as coefficient vectors;
it says their common absolute torsion has cancelled.

The adjective **unnormalized** is load-bearing.  A standard target DFT uses
`d^(kappa,q)/13`.  Division by `13` is not an endomorphism of `Lambda_M`.
If the numerator has global coefficient content `13`, primitive reduction
can divide that content and restore a nonzero Bockstein.  If it does not,
the rational normalized coefficient retains numerator (18) over denominator
`13`, and its primitive integral class remains zero.  Thus (19) gives a
content-sensitive normalization boundary, not a no-go for normalized DFTs.

## 4. Singleton model and sharpness

The lawful mass-one singleton from THM-2571 gives a transparent model.  Put

```text
a_m^(b)=zeta_13^(b+m),                  b,m in F_13.        (20)
```

Then

```text
sum_m a_m^(b)=0,

beta(a^(b))=zeta_13^b Omega=Omega mod 13M.                 (21)
```

Each `B_tau a^(b)` has exact denominator `13`.  But

```text
a^(b)-a^(b')
 =(zeta_13^b-zeta_13^b')(zeta_13^m)_m                     (22)
```

contains the missing factor `zeta_13-1`, so its sawtooth primitive is
integral.  The companion verifies all thirteen absolute denominators and
all `78` pairwise primitives directly.

The same model sharply realizes the normalization boundary.  Its
unnormalized target Fourier numerator is

```text
sum_b zeta_13^(qb)a^(b)
 =0,                         q!=-1,
 =13a^(0),                   q=-1.                          (22a)
```

Every numerator in (22a) fills integrally.  But the normalized `q=-1`
coefficient is exactly `a^(0)` again, with Bockstein `Omega` and primitive
denominator `13`.  Hence a normalized DFT can recover precisely the absolute
class that its integral numerator killed.

Three boundaries are sharp:

1. multiplying a primitive carrier by an extra global `13` kills `beta`;
2. retaining twelve translated copies preserves the class with sign, while
   adding the thirteenth kills it;
3. one absolute profile is obstructed, while its difference from any other
   profile fills.

The first is why THM-2571 uses one global primitive clearing.  The second is
why a punctured orbit and a complete orbit have opposite integral behavior.
The third shows that no amount of relative consistency alone selects the
absolute class.

## 5. Socle interpretation and the missing sidecar

Modulo `13`, write `epsilon=zeta_13-1`.  Then

```text
Z[zeta_13]/13 isomorphic to F_13[epsilon]/(epsilon^12),

Omega=-epsilon^11.                                        (23)
```

Therefore

```text
zeta_13^u Omega=Omega,

(zeta_13^u-zeta_13^v)Omega=0.                              (24)
```

The first identity makes the obstruction universal across target colour;
the second makes every target difference forget it.  Universality and
target blindness are the same algebraic fact.

This clarifies, without closing, the relation to THM-2573.  Its logarithmic
Abel normal uses the normalized target DFT
`13^(-1)sum_s`, not merely its integral numerator.  Equation (22a) proves
that this normalization can restore a killed 13-primary class.  Therefore
THM-2579 gives no no-go for the Abel normal.  A composition would have to
show that the lawful whole-layer normal has global 13-content in its
unnormalized numerator and then perform primitive reduction on the **same**
physical carrier.  Current canon supplies neither fact.  An absolute charged
reference or equivalent oriented component trivialization remains another
possible route; THM-2574 reserves such a source but is not used here.

THM-2572's positive energy does not replace this reference.  Energy can
survive after the linear Bockstein cancels, and it is invariant under phase
information that an absolute class would need.  Thus the three layers remain
distinct:

```text
rational filling,       integral absolute class,       nonlinear energy.
                                                                    (25)
```

## 6. Exact companion and stopping boundary

Run

```bash
python3 04-computation/lrc14_socle_flat_target_torsor_thm2579.py
python3 -O 04-computation/lrc14_socle_flat_target_torsor_thm2579.py
```

Both executions must reproduce

```text
05-knowledge/results/lrc14_socle_flat_target_torsor_thm2579.out
```

byte-for-byte.  The dependency-free exact referee checks:

- all `4095` nonzero Boolean augmentation profiles and the exact
  `315/3780` integral-fillability split;
- `704340` functorial controls across all translations, all ordered shift
  differences, the full/punctured orbit sums, and the Cayley derivative;
- all thirteen singleton absolute profiles and their `78` pairwise integral
  fillings;
- THM-2571's closed `Omega*Y` factorization and its `78` absolute-class,
  `468` pair-difference, and `78` unnormalized-Fourier consequences (the
  raw canonical carrier is not reconstructed a second time);
- all thirteen singleton Fourier numerators and the exact normalized
  `q=-1` restoration hostile;
- the punctured/full-orbit boundary and the extra-factor-13 hostile.

There are `1414840` explicit checks, none implemented with `assert`.

Two independent immutable audits rederived the circulant quotient action,
the `468` pairwise primitives, the `78` unnormalized Fourier primitives, and
the singleton normalization hostile.  Both replayed normal and optimized
executions byte-for-byte against the stored transcript and independently
matched the LF hashes above.  The audits required the explicit
normalized-versus-unnormalized scope, the punctured-orbit sign, and the
statement that the closed canonical evidence replays THM-2571's factorization
rather than reconstructing its raw carrier.  Documentation checks pass.

The theorem proves an integral target-torsor law and identifies the exact
loss in passing from absolute profiles to differences, together with the
sharp normalized-DFT escape.  It does not
construct a positive Cayley primitive, identify the abstract target index
with a semantic THM-2305 endpoint, make the logarithmic normal lawful on the
same carrier, supply an owner/root reference, exclude a scalar row, or prove
LRC(14). **QED.**
