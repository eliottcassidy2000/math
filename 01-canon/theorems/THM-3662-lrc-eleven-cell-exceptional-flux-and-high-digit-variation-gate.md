---
id: THM-3662
title: "LRC eleven-cell exceptional flux and high-digit variation gate"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  The 156-cell carry response has 17 zero rows, 128 nonzero rows on the
  generic reversal-even line, and exactly 11 rows with an exceptional
  reversal-odd component.  Every weighted transverse response is one
  explicit eleven-coefficient flux; after reversal symmetry it has only four
  independent interior coefficients.  Natural carry multiplicities annihilate
  the full response, and positive reversal-symmetric generic weights can also
  cancel.  Translated fluxes detect exactly high-digit variation, but no
  physical current-to-weight table is constructed.
source: kps-s191 / THM-3660 response-flux continuation, 2026-08-21
audit: >
  PASS -- agent Euler independently reconstructed the 17/128/11 response
  partition, the ten transverse projective classes, both rank-11 coordinate
  matrices, the rank-4 reversal-odd projection, every projector and flux
  identity, both cancellation hostiles, the noncocycle witness, and the two
  rank-156 translate spans with their row-constant annihilator.  Normal and
  optimized transcripts and all pinned hashes agreed.  It found no scope or
  positivity overclaim.
depends_on:
  - THM-3660-lrc-exceptional-leakage-functional-and-fourteen-edge-boundary
  - THM-3661-lrc-exceptional-detector-simple-spectrum-convolution-rigidity
related:
  - THM-3659-lrc-mod169-carry-response-cochain-obstruction
  - THM-2334-relation-residue-current-and-character-twist-pushforward
script: 04-computation/lrc_eleven_cell_exceptional_flux_thm3662.py
output: 05-knowledge/results/lrc_eleven_cell_exceptional_flux_thm3662.out
script_sha256: be144898d386b703e7eb9e754b0cd99e729b06aabb2f823c86b10e0c208f0f5e
output_sha256: 7e03fb4c31a9c7ced2d1504ececa27d393639a3dce1630becab4795663561804
semantic_sha256: 5606124428d10919903c0d7cff433e2550a5e14010c971e4cdcf261a6ff28342
hash_basis: raw LF bytes
---

# THM-3662 -- the transverse carry response is one eleven-cell flux

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  Taking the correct reversal eigenspace turns THM-3659's apparently
dense rank-two response into a sparse signed boundary.  This theorem records
the exact linear functional that any future physical weight table must pass.

## 1. Exact response partition

Let `e:F13^2->E` be the correction map, and on the carry-reachable chart set

```text
R(t0,t1)=e(t0,t1+1)-e(t0,t1),
0<=t0<=11, 0<=t1<=12.                                (1)
```

Among the 156 rows of (1), exact reconstruction gives

```text
17 zero rows,
128 nonzero rows on the generic line L_+=span(1,b),
11 rows transverse to L_+.                            (2)
```

The three parts take respectively `1,44,10` distinct vector values, for the
55 values of THM-3659.  The 11 transverse cells are

```text
(0,10),(0,11),(0,12),
(3,8),(3,9),
(6,4),(6,5),(6,6),(6,7),
(9,2),(9,3).                                         (3)
```

Their ten projective classes have multiplicities

```text
(2,1,1,1,1,1,1,1,1,1),                              (4)
```

and the repeated exact row occurs at `(0,10)` and `(6,7)`.  If each response
coordinate is arranged as a `12 by 13` matrix, both matrices have rank 11.

## 2. Reversal-odd projection is the signed exceptional boundary

THM-3657's reversal action has eigenlines

```text
L_+=span(1,b),     L_-=span(1,-b),                   (5)
```

where

```text
b=371578917865089240854253.
```

Let `Pi_-=(I-J)/2` be the odd projector.  With THM-3660's detector
`g=1_(X_+)-1_(X_-)` and

```text
c=A/(2b)=28939885940104996879767,                    (6)
```

one has on all 169 addresses and all 156 response cells

```text
Pi_- e(t)=-c g(t)(1,-b),
Pi_- R(t)=-c D_1g(t)(1,-b).                          (7)
```

The `12 by 13` scalar matrix `D_1g` has rank four, even though each
unprojected coordinate matrix has rank 11.  The useful low rank appears only
after retaining the reversal eigenspace.

## 3. The eleven-coefficient flux

For an arbitrary weight table `w=(w_(i,j))` define

```text
Phi(w)=sum_(i=0)^11 sum_(j=0)^12 w_(i,j)D_1g(i,j).   (8)
```

Then the entire transverse component of the weighted response is

```text
Pi_- sum_t w_t R(t)=-c Phi(w)(1,-b).                 (9)
```

Explicitly,

```text
Phi(w)=
 w_(0,10)-2w_(0,11)+w_(0,12)
-w_(3,8)+w_(3,9)
+w_(6,4)-w_(6,5)-w_(6,6)+w_(6,7)
+w_(9,2)-w_(9,3).                                   (10)
```

Thus the other 145 cells cannot help cancel or create a source correction
transverse to `L_+`.  A future current-to-digit construction needs only the
11 weights in (10) for this quotient test, while retaining its source,
owner, phase, chronology, and gauge sidecars elsewhere.

On the reversal-stable interior chart `1<=t0<=11`, impose

```text
w_t=w_(rho t),       rho(t0,t1)=(12-t0,11-t1).       (11)
```

Then (10) reduces to

```text
Phi_int(w)=
2[-w_(3,8)+w_(3,9)+w_(6,4)-w_(6,5)].                (12)
```

Reversal symmetry reduces the test to four independent weights but does not
force it nonzero.

## 4. Three exact hostile boundaries

### 4.1 Natural carry counting cancels everything

Every split sum in row `t0` has `13(12-t0)` carrying preimages.  These weights
are constant in `t1`, so vertical telescoping gives the vector identity

```text
sum_t 13(12-t0)R(t)=(0,0).                           (13)
```

This annihilates the full response, not only its odd projection.  Raw pair
counts therefore cannot prove exceptional exposure.

### 4.2 Positivity plus reversal symmetry is insufficient

The two generic cells

```text
t=(1,0),       rho t=(11,11)                         (14)
```

satisfy

```text
R(1,0)=
(289369708799585224299099,644718915560530562860929),

R(11,11)=-R(1,0).                                    (15)
```

The positive two-atom weight supported equally on (14) is reversal-symmetric,
avoids all 11 exceptional-flux cells, and has zero full response.

### 4.3 The projected carry response is not a 2-cocycle

Put

```text
gamma(x,y)=kappa(x,y)D_1g(x+_sp y).                  (16)
```

Although the carry bit `kappa` is the usual extension cocycle, (16) is not a
split-group 2-cocycle.  For

```text
x=(0,1), y=(1,0), z=(12,9),                          (17)
```

the standard trivial-action differential

```text
gamma(y,z)-gamma(x+y,z)+gamma(x,y+z)-gamma(x,y)
```

equals `-1`.  Thus no unmodified cohomological carry shortcut can replace
the address-dependent response.

## 5. Translated fluxes detect exactly high-digit variation

Extend a `12 by 13` carry weight table by zero on the missing column `t0=12`.
THM-3661 proves that translations of `D_1g` span exactly

```text
W={f:sum_(t1)f(t0,t1)=0 for every t0}.               (18)
```

The standard pairing is nondegenerate because `p` does not divide 169.
Consequently

```text
<w,tau_s D_1g>=0 for every address shift s

iff

w is constant in t1 on each fixed t0 row.             (19)
```

Thus a table has nonzero high-digit variation exactly when some translated
exceptional flux detects it.  The fixed flux (10) may still vanish, and a
detecting translation need not preserve the fixed exceptional set `X`.
These are now the two precise remaining geometric obligations.

## 6. Exact verification and scope

Reproduce with

```bash
python3 -B 04-computation/lrc_eleven_cell_exceptional_flux_thm3662.py
python3 -B -O 04-computation/lrc_eleven_cell_exceptional_flux_thm3662.py
```

The assertion-free companion source-pins THM-3660 and THM-3657, reconstructs
the full response, verifies (2)--(17), both matrix ranks, every projector
identity, the natural-weight cancellation, hostile pair, and noncocycle
witness.  Normal and optimized streams are byte-identical.  Its response
digest equals THM-3659's independent table digest

```text
0e48e8384bcb5d98f52eba69f3e3bf1a01d739558d5bf53d90c94767aa455939.
```

This is a static finite-field response and flux theorem.  It constructs no
THM-2334 current-to-digit weight table, proves no high-digit variation for a
covering row, supplies no physical address translation or positivity,
restores no translated detector to fixed `X`, and implies neither exceptional
entry nor LRC(14).  **QED.**
