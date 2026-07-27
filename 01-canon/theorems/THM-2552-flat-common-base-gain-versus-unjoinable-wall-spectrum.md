---
id: THM-2552
title: "Flat common-base deck gain, free-wall saturation, and the coherent-face obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The closest proved common-future groupoid carries only a deck coboundary,
  so every closed gain is zero.  The opposite relaxation, which composes the
  individually positive THM-2540 walls independently, has full F13 gain after
  two layers.  A sharp three-edge probability model has identical one- and
  two-edge marginals but prescribed nonzero face curvature, proving that the
  missing semantic 2-cell is irreducibly joint data.  No lawful live-row
  seven-layer automaton is constructed; no row is removed; LRC(14) is OPEN.
source: codex-2026-07-27-holotopy
depends_on:
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
  - THM-2540-weighted-live-event-kakeya-flux-and-transverse-gain-boundary-refinement
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
  - THM-2548-seven-step-c91-transfer-and-full-norm-separation
related:
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2466-delayed-word-simultaneous-drift-service-retention
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
  - THM-2549-future-pullback-target-neutrality-and-cemetery-hall-boundary
script: 04-computation/lrc14_common_base_gain_coherent_face_thm2552.py
output: 05-knowledge/results/lrc14_common_base_gain_coherent_face_thm2552.out
script_sha256: 1ebe5a8d30d57ef1b04a78f9f03998078e6f1b1ee261c8d7ce9d026215cb393b
output_sha256: 07c57743c408b789885d7fa73e1b2b2573236054d2bb671609added4bfc2605a
hash_basis: LF-normalized bytes
---

# THM-2552 -- the current carrier is flat or unjoinably free

THM-2548 isolates a precise invoice.  Seven semantic/chart edges with
corrections `c_0,...,c_6 in F_13` trivialize the root holonomy only when

```text
c_0+...+c_6=-7a,                 a!=0.                    (1)
```

The repo has many positive local walls and one genuine common-future
groupoid.  Neither supplies (1).  The reason is now exact: the common-future
gain is a coboundary, while the local walls lose the joint face on which
their corrections would have to compose.

## 1. The lawful common-future deck correction is flat

THM-2518's physical carrier is the equivalence-relation groupoid

```text
R_L={(x,x'):T^Lx=T^Lx'},          T(x)=13x mod 1.          (2)
```

On one oriented fibre write its inverse branches as `x_r`,
`r in Z/13^L Z`.  The **pure deck projection** to the first root residue is

```text
c(x_r -> x_s)=s-r mod 13.                                 (3)
```

For every composable closed walk

```text
x_(r_0)->x_(r_1)->...->x_(r_m)=x_(r_0),                   (4)
```

the gain telescopes:

```text
sum_(i=0)^(m-1)c(x_(r_i)->x_(r_(i+1)))=0 mod 13.           (5)
```

The same holds after every lawful vertex gauge

```text
c_phi(u->v)=c(u->v)+phi(v)-phi(u).                         (6)
```

Hence any neutral closed-section gain set built only from these arrows is

```text
Z^0 subseteq {0}.                                         (7)
```

Equality means that a neutral closed section actually exists.  Since
`-7a!=0`, the deck coboundary itself can never pay (1).

This is not yet a positive THM-2548 cut.  THM-2518 does not identify its
arrows with the seven chart layers, does not designate the one target-active
role, and does not give the later semantic root used by THM-2545.  The result
is a flatness theorem for the closest lawful carrier, not an arrival theorem.
Indeed, on seven artificially layered copies the affine correction

```text
c_a(x_r->x_s)=s-r-a                                       (7a)
```

has closed gain `-7a`.  The extra `-a` is exactly the missing chart/clock
one-cocycle: it is not supplied by THM-2518's pure deck arrow.  Thus THM-2518
alone neither constructs nor rules out a correctly clock-twisted
identification with THM-2548.

## 2. The positive wall atlas saturates when composition is forgotten

THM-2540 separately constructs a positive word/owner-marked spatial wall
for every

```text
tau in F_13^*.                                            (8)
```

These walls are not proved to share one temporal ancestry fibre or one
compatible gauge.  If that missing typing is discarded and one *relaxes*
each layer to choose an arbitrary wall correction, then

```text
D_k=F_13^*                                                (9)
```

and already

```text
D_0+D_1=F_13.                                             (10)
```

Indeed, for any `x in F_13`, choose `u` different from `0,x`; then
`x=u+(x-u)` is a sum of two nonzero residues.  Equivalently (10) is the
sharp Cauchy--Davenport bound `min(13,12+12-1)=13`.  Therefore the relaxed
seven-layer neutral spectrum is all of `F_13` and always contains `-7a`.

This gives the exact present dichotomy:

```text
lawful common fibre  -> deck correction is a coboundary -> closed gain 0;
all positive walls   -> no common fibre product         -> free gain F_13.
                                                                  (11)
```

The desired object lies strictly between them: a chart-dependent semantic
correction on one common ancestry fibre.

## 3. Pairwise data cannot reconstruct the missing face

The failure is not cured by recording every individual edge and every pair
of edges.  For `t in F_13`, let `mu_t` be the uniform probability measure on

```text
Omega_t={(c_01,c_12,c_20) in F_13^3:
         c_01+c_12+c_20=t}.                               (12)
```

Each set has `13^2=169` atoms.  For every one-coordinate or two-coordinate
projection, the pushforward of `mu_t` is uniform and independent of `t`.
Thus `mu_0` and `mu_1` have identical:

- edge supports and edge spectra;
- pairwise co-supports and every pairwise moment; and
- all statistics obtained from proper-face marginals.

But the closed-face curvature is identically `0` under `mu_0` and identically
`1` under `mu_1`.  In Fourier form, with `zeta` a thirteenth root,

```text
E_(mu_t) zeta^(alpha c_01+beta c_12+gamma c_20)
 =zeta^(gamma t) 1_(alpha=beta=gamma).                     (13)
```

All nonconstant proper-face characters vanish; only the genuinely
three-edge diagonal character sees `t`.  This is a finite, exact
Massey-type warning, not a claim about a topological Massey product.

Equation (12) is minimal in arity: with only two corrections, their joint
law is itself the full face.  Three composable edges are the first place
where all proper marginals can agree while curvature changes.

THM-2545's aligned/swap transportation hostile is the independent arrival
analogue once a later active root has already been supplied.

## 4. Exact map/loss audit

The current canon occupies the following sides of the missing fibre product.

| input | proved carrier | retained | first missing datum |
|---|---|---|---|
| THM-2461/2466 | source atom plus delayed word on a supplied common root base | owner, word, drift, service | no root transition or map from active role to correction |
| THM-2471/2478 | finite ancestry stalk plus future owner-word graft | old root/deep data and literal future event | graft is root-constant/target-neutral; rebasing loses ancestry residue |
| THM-2518 | actual common-future inverse-branch groupoid | deck arrows, old owner/deep arm, common future word | correction is a coboundary; role-wise carries remain |
| THM-2537/2540 | positive wall for each nonzero displacement | word, delayed owner, carrier re-entry | empty same-horizon head; walls are not temporally composable |
| THM-2545 | exact Hall criterion once `(sigma,h,b)` exists | word-stratified diagonal arrival | construction of the genuinely later active root `b` |
| THM-2549 | common-base selected head and strictly future owner-word algebra | positive chronology on all 165 typed rows | future algebra is root-constant/target-neutral; Hall table stops at the cemetery |

THM-2550 is a typed-row nondegeneracy candidate, not a scalar-cover carrier.
Neither it nor THM-2549 changes (11).

## 5. Cheapest decisive live test

A useful future artifact must output exact rational atoms with at least

```text
layer, word, owner, source_atom, target_atom,
source_root, target_root, owner_sheet, deep_sheet,
role_carry, correction, active_flag, mass.                 (14)
```

On those *joint* rows, not on separately integrated marginals:

1. run a seven-layer bitset DP on `(typed_vertex,gain in F_13)` using only
   neutral edges;
2. in a separate all-edge min-plus pass, test the diagonal state at gain
   `-7a` and its minimum active count;
3. on forced active survivors, form THM-2545's exact transportation table
   and run weighted Hall/max-flow.

The DP costs `O(7*13*|E|)` before Hall.  The blocking task is the lawful
generation of (14).  No theorem currently emits that chart x semantic x
ancestry fibre product.

No row is excluded.  LRC(14) remains open.

## 6. Exact companion

The dependency-free companion verifies the deck/gauge telescoping controls,
the full two-wall sumset, every proper marginal of `Omega_0,Omega_1`, and all
`2*13^3` Fourier characters in (13).  Run

```bash
python3 04-computation/lrc14_common_base_gain_coherent_face_thm2552.py
python3 -O 04-computation/lrc14_common_base_gain_coherent_face_thm2552.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_common_base_gain_coherent_face_thm2552.out
```

after LF normalization.  Every check raises explicitly under optimized
Python.

## 7. Independent hostile audit

An independent audit rederived the deck coboundary and gauge invariance,
the affine `-a` control, the two-wall sumset, all proper marginals of the
coherent-face bank, and its complete Fourier law.  It separately checked the
THM-2549 row and the separation between the neutral bitset and all-edge
min-plus passes.  The companion contains no Python `assert`; `35,170`
explicit checks execute under both ordinary and optimized Python, and both
transcripts byte-match the stored output after LF normalization.  The audit
reproduced the frontmatter hashes and found no remaining mathematical or
typing defect.
