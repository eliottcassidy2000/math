---
id: THM-3450
title: "Marked D5 carrier isomorphism and the full-germ margin obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  After linked C91,
  Puiseux, root, and amplitude
  markings, THM-2512's rational doubly-centred interaction is encoded by one
  primitive coefficient and maps isomorphically to THM-3443's Keller
  mode-one formal line by multiplication with a unit.  This is a calibrated
  characteristic-zero carrier map, not an H1 or physical-current map.  The
  full Keller germ cannot factor through the centred interaction alone: its
  first nonprimitive terms occur exactly at orders 8 and 14 in the two ANOVA
  margin sectors.
source: codex2 D5 marked-mode synthesis, 2026-08-15
audit: independent Fourier/CRT, cyclotomic-dimension, formal-module, torsion, gauge, and order-fourteen recurrence audit
depends_on:
  - THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy
  - THM-3431-d5-secondary-h1-descent-defects-and-valuation-persistence
  - THM-3443-weighted-lift-infinity-character-eigenamplitude
related:
  - THM-3440-weighted-lift-cyclic-infinity-torsor-and-7x13-character-grid
  - THM-3437-derived-boundary-jet-euler-conservation-and-prufer-recovery
script: 04-computation/d5_marked_carrier_full_germ_obstruction_thm3450.py
output: 05-knowledge/results/d5_marked_carrier_full_germ_obstruction_thm3450.out
script_sha256: 17ef1da1413d011c788d1ff000158fa163b5c08ed284c4aa34acc2630b1642f4
output_sha256: 6398cf8c57667618de68bad0667abc1b1ebdabafda93c3da4d52fd9bf2d68b68
semantic_sha256: ac9cdaf18e641572ccac24731d0959177038a7aa95a0559872a7a500708a2f2a
hash_basis: LF-normalized bytes
---

# THM-3450 -- marked D5 carrier isomorphism and the full-germ margin obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This is the maximal lawful D5 comparison currently available.  It constructs
an explicit isomorphism only after the relevant character lines and gauges
have been marked.  It simultaneously proves why neither the mod-13 `H^1`
class nor the whole Keller Puiseux germ can pass through that one line.

## 1. The LRC primitive coefficient is the whole centred interaction

Let `A=(A_(ell,r))` be a rational `7 x 13` response table as in THM-2512,
and let `I` be its doubly-centred ANOVA interaction.  Fix the degree-91
markings of THM-3440/3443 and put

```text
eta_7=zeta_91^78,              eta_13=zeta_91^14.       (1)
```

Thus `eta_7,eta_13` are the linked normalized roots for the two CRT factors.
Define

```text
lambda_L(A)
 =1/91 sum_(ell in C7,r in C13) I_(ell,r) eta_7^ell eta_13^r
 =Ahat(1,1).                                            (2)
```

The main effects vanish at mixed frequencies, so `A` and `I` give the same
coefficient in `(2)`.  THM-2512's exact cut factorization gives, for every
`tau in F_13^*`,

```text
lambda_L(A)
 =Psi^A_(tau,1)(12,1) / [91 K_(12 tau,1)].             (3)
```

The denominator is nonzero.  With the conventional factor roots

```text
zeta_7=zeta_91^13,             zeta_13=zeta_91^7,       (4)
```

the same algebraic character is labelled `(6,2)`, and

```text
lambda_L(A)
 =Psi^A_(tau,1)(11,6) / [91 K_(11 tau,6)].             (5)
```

Equations `(3)` and `(5)` are root-gauge readdressings of one character, not
two unrelated coefficients.

There is a stronger rational statement.  Zero-margin rational arrays form

```text
Q[C7]_0 tensor Q[C13]_0
  ~=Q(zeta_7) tensor_Q Q(zeta_13)
  ~=Q(zeta_91),                                        (6)
```

and both sides have dimension `6*12=72`.  Evaluation at `(eta_7,eta_13)` is
an isomorphism.  Its 72 Galois conjugates are precisely all mixed Fourier
coefficients.  Therefore

```text
lambda_L(A), together with its marked Q(zeta_91) structure,
determines the entire rational interaction I.           (7)
```

It does not recover the grand mean or either ANOVA main effect discarded by
double centring.

For a lawful rational anchored table, THM-2512 also gives

```text
lambda_L(A)=0
 iff I=0
 iff A is on the delta-plus-six-replicas branch.        (8)
```

Thus `(2)` is a complete rational replica detector, not merely one sampled
coordinate.

## 2. The explicit marked carrier isomorphism

Put

```text
E=Q(P_0,rho,zeta_91),          rho^91=91,
t=s^91,                        R=E[[t]].                (9)
```

Let `e_(1,1)` denote the LRC character line under `(1)`, and let `e_1` denote
the marked Keller character-one line.  THM-3443 supplies

```text
u_K(t)=Xhat(1;s) in R^times,       u_K(0)=rho/91.       (10)
```

After identifying the two marked deck characters, define

```text
Phi_m:R e_(1,1) -> R e_1,
Phi_m(f(t)e_(1,1))=f(t)u_K(t)e_1.                       (11)
```

Multiplication by the unit `u_K` proves:

```text
Phi_m is a strict t-adically filtered C91-equivariant isomorphism.  (12)
```

The characteristic-zero LRC response maps explicitly as

```text
A |-> lambda_L(A)u_K(t)e_1.                            (13)
```

By `(7)--(8)`, this preserves zero/nonzero, replica/nonreplica, and—with the
Galois structure retained—the full doubly-centred rational interaction.
It does not preserve positivity, main effects, owner/deep ancestry, cut
labels, clock semantics, or physical time.

The subscript `m` is essential.  Filtered equivariant isomorphisms between
the two already-identified rank-one lines form the torsor

```text
Isom^fil_(C91)(R_chi1,R_chi1)=R^times.                 (14)
```

Equation `(11)` chooses one element of that torsor by the Keller amplitude
calibration.  It is explicit, but not unmarked or canonical.

## 3. Two exact no-go theorems before the full germ

### The coefficient-object obstruction

THM-3431's secondary LRC class lies in

```text
H^1(C13;F13),                                          (15)
```

whose additive group has exponent thirteen.  The additive group of `R` is
torsion-free and uniquely thirteen-divisible.  Every additive map in either
direction is therefore zero:

```text
Hom_Add(H^1(C13;F13),R)=0,
Hom_Add(R,H^1(C13;F13))=0.                             (16)
```

Consequently `(11)` maps the characteristic-zero response amplitude `(2)`,
not the mod-13 `H^1` class.  Calling it an `H^1` bridge would be a type error.

### The independent-gauge obstruction

Reverse the Keller inertia generator while holding the LRC roots and
generator fixed.  The target character changes from `chi_1` to `chi_(-1)`.
Since 91 is odd and the coefficient field splits `C91`,

```text
Hom_(C91)(R_chi1,R_chi(-1))=0.                         (17)
```

A diagonal reversal of both sides restores the labelled line, but a
target-only reversal kills every equivariant cross-map.  Thus a nonzero map
requires synchronized inertia, branch-zero, Puiseux-root, factor-root,
active-versus-pullback, and amplitude gauges.

These two obstructions are independent: even perfectly synchronized gauges
cannot turn the exponent-thirteen class into a characteristic-zero amplitude.

## 4. Deck semivariance and the full-germ grading

For the normalized Keller branches, THM-3443 proves

```text
X_j(s)=zeta_91^j X_0(zeta_91^(-j)s),                  (18)
Xhat(k;zeta_91 s)=zeta_91^(1-k) Xhat(k;s).             (19)
```

Therefore the coefficient of Puiseux order `m` can occur only in

```text
k=1-m mod 91.                                          (20)
```

This is a Rees-type coupling between filtration order and deck character.
Only mode one descends to `E[[t]]`; the full normalized germ moves through
other character sectors.

The rational group algebra decomposes as

```text
Q[C91]
 ~= Q
  direct_sum Q(zeta_7)
  direct_sum Q(zeta_13)
  direct_sum Q(zeta_91),                               (21)
```

with dimensions

```text
1+6+12+72=91.                                          (22)
```

The doubly-centred interaction in `(6)` is only the final primitive summand.
The two middle summands are exactly the ANOVA margin sectors removed by
double centring.

## 5. The order-eight and order-fourteen obstructions

The first missing sectors are structural, not a special target sample.
Write

```text
W_0(s)=rho V(q),              q=s/rho.                 (23)
```

For the degree-91 seed, through order `q^87` the transformed inverse equation
is exactly

```text
1-V^91+(91/90)qV^90=0.                                 (24)
```

The lower seed terms first enter at orders 88 and 89, and `P_0` first enters
at order 90.  Moreover the exact differential identity is

```text
X_0(s)=s^(-90)x_0(s)=rho[V(q)-qV'(q)]/91.              (25)
```

Solving `(24)` recursively gives the already-known order-two term

```text
[s^2]X_0=-1/(16380 rho),       mode 90=(6,12),         (26)
```

which is the cheapest obstruction to putting the **whole germ** in one fixed
character line.  Modes through order seven remain primitive mixed characters.
At order eight one obtains

```text
[s^8]X_0
 =-418006897/404160880500000 * rho^84 !=0,
mode 84=(0,6).                                         (27)
```

This is the first nonprimitive term.  It lies in the `Q(zeta_13)` main-effect
sector.  At order fourteen,

```text
[s^14]X_0
 =-18258883976232522940057
   /20343399512726430000000000000 * rho^78 !=0,
mode 78=(1,0),                                         (28)
```

which lies in the other `Q(zeta_7)` main-effect sector.

Consequently:

```text
no C91-equivariant conductor-preserving full-germ factorization
can pass through the doubly-centred interaction alone.              (29)
```

Indeed, such a source has only primitive mixed characters, while `(27)` and
`(28)` are nonzero in both missing margin summands.  This does not contradict
the marked mode-one isomorphism `(11)`: `(11)` deliberately projects away
all other Keller modes.

## 6. The missing object is now sharply typed

A full D5 construction would need a **gauge-linked Rees connection**: a
formal LRC lift `Y_j(s)` satisfying

```text
Y_j(s)=zeta_91^j Y_0(zeta_91^(-j)s),                  (30)
```

so that its order-`m` coefficient lies in mode `1-m`, together with

1. both ANOVA main effects and the primitive interaction;
2. synchronized inertia/root/active-action gauges;
3. an amplitude calibration compatible with `(10)`; and
4. a common ancestry/owner/clock interpretation.

THM-3437's corrected Tate module has a formal-power-series shape, but it
supplies no identification of its parameter with `t`, no `C91`
linearization, and no LRC map.  Its Prüfer direct limit is a different object
again.  Thus `(30)` names a missing sidecar; it is not a construction.

The cheapest hostile tests for any future proposal are now fixed:

```text
target-only inertia reversal       tests gauge independence;
the order-eight (0,6) coefficient  tests centred full-germ factorization;
the order-fourteen (1,0) term      tests restoration of both margins.      (31)
```

## 7. Connection and loss ledger

| field | exact content |
|---|---|
| source | THM-2512 rational doubly-centred interaction, projected to the marked primitive line |
| target | THM-3443 Keller mode-one module `E[[t]]e_1` |
| map | mixed Fourier extraction followed by multiplication by `u_K(t)` |
| preserved | linked deck character, `t`-adic valuation, zero/nonzero, replica/nonreplica, and rational interaction with Galois sidecar |
| destroyed | ANOVA main effects, positivity, owner/deep ancestry, cut labels, clocks, physical time, and all non-mode-one Keller components |
| required markings | inertia, branch zero, Puiseux root, factor roots, action convention, `rho`, and amplitude calibration |
| first fixed-line hostile | order two, mode `(6,12)` |
| first centred-interaction hostile | order eight, mode `(0,6)` |
| other-margin hostile | order fourteen, mode `(1,0)` |

The result is a marked carrier theorem plus an obstruction theorem.  It is
not a semantic D5 word-current/JC-flux map, not a bispectrum identity, not a
physical current, and not an LRC(14) row exclusion.  It has no `JC(2)`
consequence.

## 8. Exact companion

The deterministic companion verifies `(2)--(7)` by exact reduction modulo
`Phi_91` and a `72 x 72` rational rank computation; directly replays both cut
factorizations `(3),(5)` on a nonreplica control; solves `(24)` through order
fourteen with exact fractions; verifies `(26)--(28)` and the two margin
addresses; and freezes the torsion and gauge type boundaries.  Run

```text
python -B 04-computation/d5_marked_carrier_full_germ_obstruction_thm3450.py
python -B -O 04-computation/d5_marked_carrier_full_germ_obstruction_thm3450.py
```

Normal and optimized runs reproduce the pinned transcript byte-for-byte after
LF normalization.  The structural proof, exact recurrence, and type/gauge
boundaries have passed independent audit.  This completes the proof.  QED.
