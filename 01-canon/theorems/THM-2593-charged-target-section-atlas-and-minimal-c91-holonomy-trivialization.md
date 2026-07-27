---
id: THM-2593
title: "Charged target-section atlas and minimal C91 holonomy trivialization"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  After any chosen affine identification of THM-2542's root-deck sheet with
  THM-2585's target-shift index, the thirteen everywhere-unit Bockstein
  factors form a faithful multiplicative local system on the same C91
  mapping torus.  The additive gauge h=-q kills the pulled-back old-head
  transition, while the multiplicative gauge Y_q^-1 kills the exact unit
  transport Y_(q+a)Y_q^-1.  Both seven-edge monodromies have exact order 13,
  so degree 13 is their common minimal cyclic trivializing cover.  All 546
  owner/clock/sheet Bocksteins survive.  A correction on t of the seven base
  edges leaves holonomy (7-t)a and cannot close unless t=7.  The affine
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
script_sha256: ae6b0f81c08663497c33727cefd5aacd05f583dc005d0d92777a8b4d5e5aecce
output_sha256: 861325a6cf4c70e3d2ea31bee902a04dfb35d91ff5ed0eb93e5bdca76288029d
hash_basis: LF-normalized bytes
---

# THM-2593 -- charged target sections on the minimal C91 cover

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-2542 finds a nonzero additive root-deck local system over the seven
owner-clock charts.  THM-2585 finds thirteen translation-permuted target
sections whose first Bockstein factors are all units.  These two facts fit
more tightly than an external product: after choosing an affine
identification of their two `F_13` torsors, the unit factors themselves carry
a second local system on the same minimal `C_91` cover.

The two systems are trivialized by the same sheet coordinate in two different
senses:

```text
old-head root transition:       a  + d(-q) = 0;
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

Pull the old-head transition back along (3) and put

```text
h(k,q)=-q.                                                  (6)
```

On every lifted edge,

```text
g_k+h(sigma(k,q))-h(k,q)
 =a-(q+a)+q=0.                                             (7)
```

Thus the root-deck local system is explicitly trivial on the `C_91` cover.
After one seven-edge base loop, the sheet changes by `7a`; after `n` base
loops its additive holonomy is

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

Multiplication by (9) transports the Bockstein factor at one vertex to the
factor at its successor:

```text
u_(q,a)Y_q=Y_(q+a).                                       (10)
```

The vertex gauge `v_q=Y_q^-1` trivializes this transport:

```text
v_(q+a) u_(q,a) v_q^-1
 =Y_(q+a)^-1 [Y_(q+a)Y_q^-1]Y_q=1.                        (11)
```

Along any path the unit products telescope.  In particular, after `n`
seven-edge base loops the multiplicative monodromy based at sheet `q` is

```text
U_n(q)=Y_(q+7na)Y_q^-1.                                   (12)
```

Because the thirteen `Y_q` are distinct, `U_n(q)=1` exactly when `13|n`.
Hence both (8) and (12) have exact order thirteen.  The degree-thirteen cover
is not merely sufficient for both local systems: it is their common minimal
cyclic trivializing cover.

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

nonzero lifted Bockstein profiles.  The exact product-basis support histogram
is the THM-2585 histogram repeated seven times:

```text
support 48: 56 profiles,
support 60: 224 profiles,
support 72: 266 profiles.                                  (15)
```

Thus passing to the minimal holonomy-killing cover loses neither a sheet nor
an owner colour.  What it loses is descent to one preferred base section:
one seven-edge circuit advances through the faithful thirteen-element atlas.

## 5. Every proper base-edge correction leaves residual holonomy

The double trivialization also gives a sharp coverage invoice.  Suppose a
future mixed-square construction supplies the needed root correction `-a`
on exactly `t` of the seven base chart edges and supplies zero correction on
the others.  The corrected cyclic sum is

```text
7a-ta=(7-t)a.                                               (16)
```

For `0<=t<=7`, `a!=0`, and characteristic thirteen,

```text
(7-t)a=0  iff  t=7.                                        (17)
```

In particular, even a three-edge attaching carrier leaves residual holonomy
`4a!=0`.  This does not say that every future correction must be edgewise
`0` or `-a`; a genuinely mixed 2-cell may distribute other corrections whose
total is `-7a`.  Equation (17) is the sharp statement for the natural
section-by-section cancellation supplied by (7).

The same boundary explains why a positive carrier over a proper clock arc is
important but not yet a cycle-level closure.  Local gauges always exist on a
tree.  The obstruction is paid only when the seventh edge closes the loop.

## 6. Sharp scope and failure controls

Three pieces are load-bearing.

1. **Common translation atlas.**  Equation (5) comes from one fixed primitive
   carrier.  THM-2585's same-cokernel hostile shows that arbitrary independent
   representatives do not admit normalized integral sections.
2. **Unit factors.**  Without unitness, the ratios in (9) need not exist.
   Without pairwise distinctness, their seven-edge monodromy need not have
   order thirteen.  A constant unit atlas is the sharp collapsed control.
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

The nearby selector and pullback candidates THM-2591 and THM-2592 address two
different parts of that debt.  They are related evidence, not dependencies
of this candidate, and no conclusion here relies on their promotion.

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
- all `2,028` based seven-loop monodromy/order checks; and
- all `96` proper-edge residual checks, including the three-edge `4a` wall.

There are `7,716` explicit checks; none is implemented with `assert`.

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
LRC(14) follows. **QED (candidate pending independent audit).**
