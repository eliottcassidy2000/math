---
id: THM-3647
title: "LRC single point-reversal branch-orbit spectral projector"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE AUDIT.
  The point-order-reversal fixed subalgebra of THM-3625's thirteen branch
  addresses is exactly the two-channel algebra F I + F Pi_W.  Each of its
  seven branch-orbit packets is an invertible cyclic generator, so any one
  affinely recovers Pi_W; adjoining THM-3636's distinct native arc reversal
  then recovers Pi_K.  The fixed branch r=6 alone suffices.  This is static
  finite-field algebra, not chronology, current, characteristic-zero
  transport, row exclusion, or LRC(14).
source: kps-s189 / THM-3636 sparse branch-packet continuation, 2026-08-21
depends_on:
  - THM-3636-lrc-arc-reversal-and-middle-address-quotient-intertwiner
related:
  - THM-3625-lrc-point-by-branch-split-four-character-address-algebra
  - THM-3534-r5-middle-response-relative-cospan-and-twisted-h1-collapse
script: 04-computation/lrc_single_reversal_paired_branch_projector_thm3647.py
output: 05-knowledge/results/lrc_single_reversal_paired_branch_projector_thm3647.out
script_sha256: 2e42c26d43e5d94b5a5c33c6b425dc64f72b7cd38e0759d8dba61aab3fb4c11b
output_sha256: 300a43fcc6cbc679d0ec91f571d19115c292555d1ae1fa92db64e65678002950
semantic_sha256: 1f2ab658418c85b68b39f1f590c64ce690890df169c53e9292b91d5858b07e87
hash_basis: LF-normalized bytes for parents; raw LF bytes for this package
---

# THM-3647 -- one point-reversal branch orbit recovers the LRC projector

**PROVED + FINITE-EXACT + VERIFIED-EXACT; PENDING INDEPENDENT HOSTILE
AUDIT.**  This compresses the static projector mechanism of THM-3636 from
the full branch-address algebra to any one point-reversal orbit packet.  Two
different reversals are load-bearing and must not be conflated.

## 1. The two involutions

Work over the pinned field

```text
F=F_p,                   p=755373809845391722745761.    (1)
```

In the six-point order

```text
(0,0),(1,0),(1,6),(3,6),(3,12),(2,12),                (2)
```

let `A_r`, `r in F_13`, be THM-3625's branch-address operators.  Point-order
reversal is

```text
J=(0 5)(1 4)(2 3).                                     (3)
```

THM-3625 proves, and the companion reconstructs directly, the covariance

```text
J A_r J=A_(-r-1).                                      (4)
```

The native directed-arc reversal used by THM-3636 is instead

```text
S=(0 1)(2 3)(4 5).                                     (5)
```

The operators `J` and `S` are not interchangeable: `J` pairs branch digits,
whereas `S` splits the two doubled address-character spaces.

## 2. Seven orbit packets collapse to two channels

The involution in `(4)` has the seven orbits

```text
{0,12},{1,11},{2,10},{3,9},{4,8},{5,7},{6}.            (6)
```

For representatives `0<=r<=6`, put

```text
B_r=A_r+A_(-r-1)       for r<6,
B_6=A_6.                                                (7)
```

The singleton convention in `(7)` is important.  Since the orbits partition
all thirteen branches and `sum_r A_r=I`, one has

```text
sum_(r=0)^6 B_r=I.                                     (8)
```

Exact reconstruction gives

```text
rank B_r=6                         for every r,
dim span{B_0,...,B_6}=2.                              (9)
```

Let

```text
Pi_W=diag(1,1,0,0,1,1),             Pi_M=I-Pi_W.      (10)
```

For every orbit packet there are distinct nonzero scalars
`alpha_r,beta_r` such that

```text
B_r=alpha_r Pi_W+beta_r Pi_M.                          (11)
```

Their exact values are:

| `r` | mate | `alpha_r` on `W` | `beta_r` on the middle plane |
|---:|---:|---:|---:|
| 0 | 12 | 583299261868940099523963 | 460700142682906771021811 |
| 1 | 11 | 43106534571494424015806 | 107225529790717415656375 |
| 2 | 10 | 353459988660952602113073 | 107225529790717415656375 |
| 3 | 9 | 32271302543699805871617 | 189845688026952442643144 |
| 4 | 8 | 440916124876065854346257 | 745044336390295863842285 |
| 5 | 7 | 327919687265179724083408 | 122334475151651243469948 |
| 6 | 6 | 485148529749842658283160 | 533745727702934015947346 |

Consequently the point-reversal fixed algebra is exactly

```text
mathcal A^J=F I direct_sum F Pi_W.                     (12)
```

Indeed the Reynolds image of the spanning set `{A_r}` is spanned by the
packets `(7)`, while `(9)--(11)` identify that image with the right side of
`(12)`.

## 3. Every packet is a cyclic spectral detector

Equation `(11)` and `alpha_r!=beta_r` give the quadratic identity

```text
(B_r-alpha_r I)(B_r-beta_r I)=0.                       (13)
```

The two factors have ranks `2` and `4`, respectively.  Thus the two
eigenspaces are exactly the endpoint four-space `W` and the middle two-plane.
Lagrange interpolation gives the load-bearing affine formula

```text
Pi_W=(B_r-beta_r I)/(alpha_r-beta_r)                  (14)
```

for **each** `r=0,...,6`.  In particular, the fixed address `A_6=B_6` alone
recovers `Pi_W`; no branch summation is required in that orbit.

The exact affine coefficients and all seven identities are pinned in the
companion transcript.  The ordered packet digest is

```text
2740c4f8ef562c6393fd101528a5916b5a033e0240fcf564f15728bbf4d2b034. (15)
```

## 4. Adding native arc reversal recovers rigidity

THM-3636 proves

```text
Pi_K=Pi_W (I+S)/2,
im Pi_K=K2sharp.                                       (16)
```

Combining `(14)` and `(16)` gives seven sparse formulas

```text
Pi_K=
  (B_r-beta_r I)/(alpha_r-beta_r) * (I+S)/2.           (17)
```

Although the packets were formed using `J`, equation `(11)` implies that
each packet also commutes with `S`.  The two operations therefore separate
cleanly: point reversal forgets the distinction within a branch orbit and
detects endpoint versus middle support; arc reversal selects the even line
inside each endpoint character doublet.

This is a strict compression of static observables.  THM-3625 needed the
four-dimensional algebra to expose all four characters; THM-3647 shows that
the coarser endpoint/middle predicate needs only one two-spectrum packet.

## 5. Reproduction and strict boundary

Reproduce with

```bash
python3 04-computation/lrc_single_reversal_paired_branch_projector_thm3647.py
python3 -O 04-computation/lrc_single_reversal_paired_branch_projector_thm3647.py
```

The assertion-free companion source-pins THM-3615, THM-3625, and THM-3636;
rebuilds the branch-resolved tensor; checks all thirteen covariance equations,
all seven spectral decompositions, quadratic identities, affine projectors,
the orbit sum and two-dimensional span; and verifies the inherited `Pi_K`
image.  Normal and optimized streams must reproduce the stored transcript
after LF normalization.

The packet `B_r` is a sum of static branch-conditioned operators, not a
proved temporal transition or physical choice.  Even the singleton `A_6`
is an address operator in the pinned coefficient model, not an LRC child
map.  The theorem supplies no chronology, current, characteristic-zero
transport, row exclusion, or proof of LRC(14).  Its concrete next test is
whether a lawful digit/entry operation realizes one of the packet actions in
the physical carrier without introducing a noncanonical section.  **QED.**
