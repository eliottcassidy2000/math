---
id: THM-2593
title: "Charged target-section atlas and minimal C91 holonomy trivialization"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED (two audits).
  After any chosen affine identification of THM-2542's root-deck sheet with
  THM-2585's target-shift index, the thirteen everywhere-unit Bockstein
  factors form a faithful multiplicative local system on the same C91
  mapping torus.  The additive gauge h=-q kills the pulled-back THM-2542
  root-deck
  transition, while the multiplicative gauge Y_q^-1 kills the exact unit
  transport Y_(q+a)Y_q^-1.  The additive seven-edge deck displacement has
  order 13, and the multiplicative skew cocycle has least return period 13,
  so degree 13 is their common minimal cyclic trivializing cover.  All 546
  pulled-back owner/clock/sheet Bockstein coefficient profiles survive.  A
  correction -a on t edges of one
  fixed seven-edge lifted base-loop path leaves holonomy (7-t)a and cannot
  close unless t=7.  The affine
  identification and a physical common carrier remain unconstructed; this
  is not a semantic arrival, row exclusion, or proof of LRC(14).
source: root-holotopy-2026-07-28
depends_on:
  - THM-2542-seven-chart-cech-holonomy-and-c91-arrival-obstruction
  - THM-2585-saturated-normalized-target-projector-and-bockstein-noncommutation
related:
  - THM-2551-horizontal-transfer-transverse-projector-bicomplex-boundary
  - THM-2590-boolean-bockstein-and-theta-selector-incidence-spectrum
  - THM-2591-theta-zero-selector-cech-coboundary-and-c91-holonomy-no-go
  - THM-2592-fallback-rail-digit-diagonal-pullback-and-primitive-bockstein
script: 04-computation/lrc14_charged_target_atlas_c91_holonomy_thm2593.py
output: 05-knowledge/results/lrc14_charged_target_atlas_c91_holonomy_thm2593.out
script_sha256: 5e516aa3e3ddef475589db05fc43fd0c9d1047c2c0ca231f3172ad0508fdd11b
output_sha256: 541f0272737299cea546b3937b18b56e145ead944cccca4bae5ff8ddfca8bcc3
hash_basis: LF-normalized bytes
---

# THM-2593 -- charged target sections on the minimal C91 cover

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED (two audits).**

THM-2542 finds a nonzero additive root-deck local system over the seven
owner-clock charts.  THM-2585 finds thirteen translation-permuted target
sections whose first Bockstein factors are all units.  These two facts fit
more tightly than an external product: after choosing an affine
identification of their two `F_13` torsors, the unit factors themselves carry
a second local system on the same minimal `C_91` cover.

The two systems are trivialized by the same sheet coordinate in two different
senses:

```text
THM-2542 root-deck transition:  a  + d(-q) = 0;
Bockstein-factor transition:     Y_(q+a)Y_q^-1,
normalized by:                  Y_(q+a)^-1 [.] Y_q = 1.     (1)
```

This is a coefficient-level double trivialization.  It does not canonically
identify the root and target torsors, and it does not construct a physical
semantic edge between their underlying Boolean packets.

## 1. The two torsors and the chosen comparison

Fix a nonzero marker step

```text
a in F_13^*.                                                (2)
```

THM-2542 puts the constant transition `g_k=a` on the oriented seven-cycle of
charts.  Its mapping torus has vertices

```text
(k,q) in F_7 x F_13
```

and successor

```text
sigma(k,q)=(k+1,q+a).                                      (3)
```

The projection `(k,q)->k` is the connected degree-thirteen cover of the
seven-cycle, and (3) is one cycle of length ninety-one.

Independently, THM-2585 gives one globally primitive carrier `x` and thirteen
literal target-shift sections `D^(kappa,q)`.  Their first-Bockstein factors
are

```text
beta(D^(kappa,q))=Omega Y_q(zeta_7^kappa),                 (4)

Y_q in R_7=F_13[z]/(Phi_7(z)),
```

where every `Y_q` is a unit.  Carrier translation satisfies

```text
D^(kappa,q)(T_a x)=D^(kappa,q+a)(x).                       (5)
```

To compare (3) and (5), choose one affine torsor identification.  Changing
its origin replaces every `q` by `q+q_0`; changing its orientation replaces
`a` in (5) by a nonzero scalar multiple.  The proof is identical for every
such choice.  No theorem currently makes one of these choices canonical on a
physical LRC packet.

## 2. Additive trivialization on the cover

Pull the THM-2542 root-deck transition back along (3) and put

```text
h(k,q)=-q.                                                  (6)
```

On every lifted edge,

```text
g_k+h(sigma(k,q))-h(k,q)
 =a-(q+a)+q=0.                                             (7)
```

Thus the root-deck local system is explicitly trivial on the `C_91` cover.
After one lift of the seven-edge base loop, the sheet changes by `7a`; after
`n` such lifted traversals its additive deck displacement is

```text
7na in F_13.                                                (8)
```

It vanishes exactly when `13|n`.  This recovers the minimal-cover statement
of THM-2542 with an explicit sheet gauge that will also normalize the charged
atlas below.

## 3. A faithful unit local system on the same cover

The thirteen rows `Y_q` displayed in THM-2585 are pairwise distinct.  This
is an exact extra property of the canonical carrier, not a consequence of
unitness alone.  Define the edge unit

```text
u_(q,a)=Y_(q+a)Y_q^-1 in R_7^*.                            (9)
```

For owner colour `kappa`, let

```text
sigma_kappa:R_7->R_7,              z |-> z^kappa,

u^(kappa)_(q,a)=sigma_kappa(u_(q,a)).                       (9a)
```

Multiplication by (9a) transports the actual owner-coloured Bockstein factor
at one vertex to the factor at its successor:

```text
u^(kappa)_(q,a) beta(D^(kappa,q))
 =beta(D^(kappa,q+a)).                                    (10)
```

At the universal `R_7` level, the vertex gauge `v_q=Y_q^-1` trivializes this
transport:

```text
v_(q+a) u_(q,a) v_q^-1
 =Y_(q+a)^-1 [Y_(q+a)Y_q^-1]Y_q=1.                        (11)
```

Along any path the unit products telescope.  In particular, after `n`
traversals of the seven-edge base loop lifted from sheet `q`, the
multiplicative deck transport is

```text
U_n(q)=Y_(q+7na)Y_q^-1.                                   (12)
```

Because the thirteen `Y_q` are distinct, `U_n(q)=1` exactly when `13|n`.
For `n<13`, (12) transports between distinct vertices of the cover; it is not
a closed loop upstairs, nor is `U_n(q)` the `n`th power of one fixed unit.
Thus (8) has additive order thirteen, while the skew multiplicative cocycle
(12) has least return period thirteen.  Only the ninety-one-edge path closes
on `C_91`.  The degree-thirteen cover is not merely sufficient for both local
systems: it is their common minimal cyclic trivializing cover.

This is the useful holotopy match:

```text
additive obstruction       q -> q+a       killed by -q;
multiplicative charge      Y_q -> Y_(q+a) killed by Y_q^-1;
common return time         13 clock loops = 91 edges.       (13)
```

## 4. Charge is preserved at every lifted vertex

For every `q in F_13` and every nonzero owner colour
`kappa in F_7^*`, THM-2585 proves (4) is nonzero and `Y_q` is a unit.  Repeating
the thirteen sections over the seven clock vertices therefore gives

```text
7*13*6=546                                                  (14)
```

nonzero pulled-back Bockstein coefficient profiles.  The exact product-basis
support histogram
is the THM-2585 histogram repeated seven times:

```text
support 48: 56 profiles,
support 60: 224 profiles,
support 72: 266 profiles.                                  (15)
```

Thus passing to the minimal holonomy-killing cover loses neither a sheet nor
an owner colour.  What it loses is descent to one preferred base section:
one traversal of the seven-edge base circuit advances to a different sheet
of the faithful thirteen-element atlas.

## 5. Every proper edgewise correction on one lifted base circuit leaves a residue

The double trivialization also gives a sharp coverage invoice.  Fix one
seven-edge base-loop path and one of its lifts, hence one common marker step
`a` and sheet orbit.  Suppose a future mixed-square construction supplies the
needed root correction `-a` on exactly `t` of those seven chart edges and
supplies zero correction on the others.  The corrected cyclic sum is

```text
7a-ta=(7-t)a.                                               (16)
```

For `0<=t<=7`, `a!=0`, and characteristic thirteen,

```text
(7-t)a=0  iff  t=7.                                        (17)
```

In particular, literal corrections `-a` on only three of those base edges
leave residual holonomy `4a!=0`.  This says nothing about an arbitrary
three-cell or attaching carrier: positive cells are vertex data until a
lawful boundary map turns them into edge cochains.  Nor does it say that
every future correction must be edgewise `0` or `-a`; a genuinely mixed
2-cell may distribute other corrections whose total is `-7a`.  Equation
(17) is the sharp statement for the natural section-by-section cancellation
supplied by (7).

The same boundary explains why a positive carrier over a proper clock arc is
important but not yet a cycle-level closure.  Local gauges always exist on a
tree.  The obstruction is paid only when the seventh edge closes the loop.

## 6. Sharp scope and failure controls

Three pieces are load-bearing.

1. **Common translation atlas.**  Equation (5) comes from one fixed primitive
   carrier.  THM-2585's same-cokernel hostile shows that arbitrary independent
   representatives do not admit normalized integral sections.
2. **Unit factors.**  Without unitness, the ratios in (9) need not exist.
   Without pairwise distinctness, their skew seven-edge transport need not
   have least return period thirteen.  A constant unit atlas is the sharp
   collapsed control.
3. **A physical torsor identification.**  The equality of two alphabet sizes
   does not identify THM-2542's root deck with THM-2585's target-shift
   position.  This theorem grants an affine identification and proves what it
   would buy.  It does not construct that map inside one Boolean integral.

The construction can be realized on an external coefficient product while
all semantic vertical edges are empty; THM-2542's equal-mass hostile remains
unchanged.  THM-2551 likewise shows that product transfer preserves the
transverse kernel.  A physical use therefore still needs a common-ancestry
mixed square whose root and target labels obey (3)--(5), or a positive
semantic path living directly on the twisted `C_91` carrier.

The nearby selector candidate THM-2591 and promoted physical pullback
THM-2592 address two different parts of that debt.  They are related evidence,
not dependencies of this candidate, and no conclusion here relies on
THM-2591's promotion or on retyping THM-2592's positive cells as edges.

## 7. Exact companion

Run

```bash
python3 04-computation/lrc14_charged_target_atlas_c91_holonomy_thm2593.py
python3 -O 04-computation/lrc14_charged_target_atlas_c91_holonomy_thm2593.py
```

Both executions reproduce the stored transcript byte-for-byte after LF
normalization.  The dependency-free checker works in
`F_13[z]/(Phi_7)` and verifies:

- all thirteen slice determinants and pairwise distinct factors;
- all `78` owner/sheet charges and the lifted support histogram (15);
- all twelve nonzero mapping-torus steps and all `91` vertices per step;
- `1,092` additive and `1,092` multiplicative gauge identities;
- all `2,028` based seven-loop additive-order and cocycle-return checks; and
- all `96` proper-edge residual checks, including the three-edge `4a` wall.

There are `7,716` explicit checks; none is implemented with `assert`.

Two independent hostile audits rederived the additive gauge, the
owner-automorphism unit transport, the telescoping skew cocycle, its least
return period, the `546`-profile support census, and the literal-edge residual
invoice.  They required the cover-versus-base and vertex-versus-edge scope
repairs now stated above.  Normal, optimized, and stored transcripts agree.
On Windows the raw CRLF hashes are `1ba43358...f376` and
`1a0bcd72...c3d4`; after the declared LF normalization they are exactly the
full frontmatter hashes.

## 8. Stopping boundary

The proved implication is conditional and coefficient-level:

```text
chosen affine root/target sheet identification
 + THM-2585 common charged atlas
 -> simultaneous additive/multiplicative trivialization on C91
 -> no loss of Bockstein charge.                            (18)
```

The missing input is the first line of (18) on a lawful physical packet that
also retains semantic owner/arrival meaning.  The selected `q` is still a
target-shift position, not a nonzero target character or relation address;
`Y_q` is a Bockstein coefficient factor, not a nonnegative endpoint.  No
all-165 statement, scalar-cover contradiction, row decrement, or proof of
LRC(14) follows. **QED.**
