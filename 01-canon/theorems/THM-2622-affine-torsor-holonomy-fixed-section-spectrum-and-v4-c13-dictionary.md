---
id: THM-2622
title: "Affine-torsor holonomy fixed-section spectrum and V4/C13 dictionary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For a finite abelian deck torsor with cyclic affine holonomy x->Ax+c,
  compatible parallel sections are exactly the fixed points: none if
  c is not in im(I-A), otherwise one coset of ker(I-A).  Thus the C13
  spectrum is 0/13 for fixed-action transport and 1 for every nonidentity
  generator twist.  For V4=F2^2, AGL(2,2)=S4 gives the exact section-count
  spectrum 0^9,1^8,2^6,4^1.  For quartic monodromy G<=S4, its image in
  GL(2,2)=S3=S4/V4 is the cubic-resolvent monodromy (full S3 in the generic
  full-S4 case).  Double transpositions are invisible pure
  translations, while transpositions and four-cycles have the same
  order-two linear shadow but different affine cocycles.  The D4 subgroup
  has spectrum 0^5,2^2,4^1.  This unifies the LRC C13 ancestry obstruction
  with the quartic V4 origin loss, but constructs neither a physical LRC
  connector nor a Keller exclusion.
source: kind-pasteur-2026-07-28-affine-holotopy-spectrum
depends_on:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
related:
  - THM-2612-d4-deck-pole-tax-and-depressed-resolvent-gcd-gate
script: 04-computation/affine_torsor_holonomy_spectrum_thm2622.py
output: 05-knowledge/results/affine_torsor_holonomy_spectrum_thm2622.out
script_sha256: 6001c091451dba07ab92f940b3adbea867c6a3326d17294e81664527547857e0
output_sha256: e78374167ca727e5c9d1e33c09d2f51774e7271f9a9612401683b70e672a53a5
hash_basis: LF-normalized bytes
---

# THM-2622 -- affine holonomy counts the surviving sections

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The repeated holotopy failure in LRC and the quartic resolvent has one finite
normal form.  A quotient remembers the linear action on its deck group but
forgets the translational cocycle.  A section depends on both.  This theorem
makes that statement exact and shows that the two live examples are not an
analogy:

```text
C_13 fixed generator:       0 or 13 sections;

V_4 semidirect S_3=S_4:     0, 1, 2, or 4 sections.       (1)
```

The missing datum is the affine term `c`.  In the quartic case it is exactly
what distinguishes an invisible double transposition from the identity, and
a four-cycle from a transposition sharing the same resolvent shadow.

## 1. The affine fixed-section theorem

Let `V` be a finite abelian group and let

```text
P_0 --phi_0--> P_1 --phi_1--> ... --phi_(n-1)--> P_0      (2)
```

be a cyclic family of `V`-torsors.  Suppose each transition is affine: after
choosing origins it has the form

```text
phi_i(x)=A_i x+c_i,                 A_i in Aut(V).        (3)
```

Composition around the cycle gives one affine holonomy

```text
H(x)=Ax+c,                                                (4)

A=A_(n-1)...A_0,

c=c_(n-1)+A_(n-1)c_(n-2)+...+A_(n-1)...A_1 c_0.          (5)
```

A parallel section is a tuple `(x_0,...,x_(n-1))` satisfying
`x_(i+1)=phi_i(x_i)`, including the closing edge.  Every entry is determined
by `x_0`, and closure is exactly

```text
(I-A)x_0=c.                                              (6)
```

Therefore the section set is empty when `c` is not in `im(I-A)`.  If one
solution `x_*` exists, all solutions are

```text
x_*+ker(I-A).                                            (7)
```

In particular, when `V` is a vector space over `F_p`, the number of sections
is either zero or

```text
p^(dim ker(I-A)).                                        (8)
```

Changing the chosen origins changes (4) by affine conjugacy.  It bijects the
fixed points, so the existence and number in (7)--(8) are intrinsic.  This
is the finite holotopy invariant of the cyclic local system.

When every `A_i=I`, equation (4) is translation by `c`.  A free torsor has a
fixed point only for `c=0`; then all `|V|` points are fixed.  This recovers
THM-2611's fixed-action holonomy theorem, now as one row of the affine
spectrum.

## 2. The `C_13` spectrum and its lawful boundary

Take `V=F_13`.  Every affine holonomy is

```text
H(x)=lambda x+c,                lambda in F_13^*.         (9)
```

Equation (6) gives the complete classification:

```text
lambda=1,c=0:       13 sections;
lambda=1,c!=0:       0 sections;
lambda!=1:           1 section.                          (10)
```

Across all `13*12=156` affine maps, the histogram is

```text
0^12,               1^143,               13^1.          (11)
```

The `lambda=1` row is precisely the physical fixed-generator situation in
THM-2611: zero kernel holonomy gives thirteen references and nonzero
holonomy gives none.  The `143` unique-section cases do not repair LRC.
They change the generator of the `C_13` deck, hence reindex the root character
and target action whose physical identification is the point at issue.
Allowing semilinear rebasing would solve a different problem.

## 3. The quartic `V_4` torsor is the affine four-point action

Identify the reconstruction deck in THM-2598 with

```text
V_4=F_2^2.                                               (12)
```

Its three nonzero elements are the three double transpositions.  Conjugation
permutes them, so

```text
GL(2,2) isomorphic S_3.                                  (13)
```

The affine action on the four points of `V_4` is faithful and has order
`4*6=24`; hence

```text
AGL(2,2)=V_4 semidirect GL(2,2) isomorphic S_4.           (14)
```

Under (14), the linear part is exactly the action on the three quartic
pairings.  For any quartic monodromy subgroup `G<=S_4`, its image in this
`S_3` quotient is the cubic-resolvent monodromy; it is the full `S_3` only in
the generic/full-`S_4` case.  The affine translation is the lost root-origin
sidecar.  Applying (7) gives the complete table:

| quartic class | number | linear shadow | affine condition | root sections |
|:---|---:|:---|:---|---:|
| identity | 1 | `A=I` | `c=0` | 4 |
| transposition | 6 | `ord(A)=2` | `c in im(I-A)` | 2 |
| double transposition | 3 | `A=I` | `c!=0` | 0 |
| 3-cycle | 8 | `ord(A)=3` | always solvable | 1 |
| 4-cycle | 6 | `ord(A)=2` | `c notin im(I-A)` | 0 |

Thus the exact fixed-section histogram is

```text
0^9,                 1^8,                 2^6, 4^1.      (15)
```

Two losses in THM-2598 now have a single explanation.  A double
transposition is a nonzero pure translation, so the resolvent sees the
identity linear action while the quartic has no fixed root.  A transposition
and a four-cycle both have order-two linear shadow; the former has two fixed
roots and the latter none because their affine cocycles lie on opposite sides
of `im(I-A)`.  Linear resolvent monodromy alone cannot distinguish them.

## 4. The `D_4` boundary

A Sylow-two subgroup of (14) is obtained from all four translations and one
order-two linear involution.  Its eight affine maps have section spectrum

```text
0^5,                 2^2,                 4^1,            (16)
```

and cycle atlas

```text
identity^1, transposition^2, double-transposition^3,
four-cycle^2.                                               (17)
```

This is the holotopy boundary behind THM-2612.  The linear `C_2` quotient
cannot distinguish the two transpositions from the two four-cycles in its
nontrivial coset, and it erases all three nonidentity translations in the
kernel.  THM-2612's pole and gcd invoices are additional geometric data that
may constrain this affine cocycle; (16) does not by itself prove such a pole,
exclude `D_4`, or close a Keller branch.

## 5. Exact evidence and stopping boundary

Run

```text
python 04-computation/affine_torsor_holonomy_spectrum_thm2622.py
python -O 04-computation/affine_torsor_holonomy_spectrum_thm2622.py
```

The dependency-free companion exhausts all `156` affine maps of `F_13`, all
`24` maps of `AGL(2,2)`, the selected eight-element `D_4` subgroup, all `576`
affine conjugacy pairs, and all `576` two-edge cycle systems.  It checks the
fixed-point counts, linear orders, full permutation cycle atlas, gauge
invariance, and the equality between compatible sections and holonomy fixed
points with optimized-mode-safe guards.  Two independent audits rederived the
general coset theorem, every class in the `S_4` and `D_4` tables, the LRC
fixed-generator boundary, and the quartic-subgroup scope.  Normal and
optimized runs both byte-match the stored transcript.

The theorem classifies an already supplied affine deck local system.  It does
not construct THM-2611's missing physical ancestry carrier, identify the
THM-2613 local shift with the THM-2585 next-target state, prove that a quartic
resolvent comes from a Keller source, exclude any LRC row, prove JC(2), or
prove LRC(14).  QED.
