---
id: THM-3599
title: "Universal Danielewski affine-linear rational Darboux classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On every smooth squarefree Danielewski surface c^N e=Sigma(b) with N>=2
  and deg Sigma>=2, an affine-linear target delta+lambda b+mu c+nu e has a
  rational nonzero-constant-bracket mate if and only if nu=0 and the target
  is nonconstant.  The nu!=0 time form is an interior Newton differential,
  except for one logarithmic conic boundary with nonzero residues.  No
  polynomial Darboux pair or counterexample to JC(2) is constructed.
source: root / planar-Jacobian affine time-form session, 2026-08-21
audit: >
  An independent hostile audit rederived generic-fibre regularity, every
  Newton mask and margin, all toric faces, the resonant double-root
  resolution, affine-arm units, the unique logarithmic conic residues, both
  rational-mate formulas, and the exponent/degree boundaries.  Normal and
  optimized companions are byte-identical to the stored 2,611-gate output.
depends_on:
  - THM-3596-a13-invoice-paid-mixed-coordinate-toric-time-form-nonentry
related:
  - THM-3595-danielewski-separated-transverse-time-form-rational-nonentry
  - THM-3598-danielewski-rational-exact-polar-graph-family-and-classification
script: 04-computation/jc2_universal_danielewski_affine_linear_rational_classification_thm3599.py
output: 05-knowledge/results/jc2_universal_danielewski_affine_linear_rational_classification_thm3599.out
script_sha256: 621314d2b76809d3297fc63e433d800923818c7a90d562061aa89ab59c1da430
output_sha256: 37fef42a4577b01af4d1b8b636a7390de3bf94c8858f546aeae900d33aab9a2b
hash_basis: raw LF bytes
---

# THM-3599 -- Universal Danielewski affine-linear rational Darboux classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
gives a complete rational-mate classification in the
affine-linear target cell.  The hypotheses `N>=2`, `deg Sigma>=2`, and
squarefreeness are explicit and load-bearing.  It proves no case of `JC(2)`.

## 1. Statement and the easy half

Work over `C`.  Let `N,D>=2`, let `Sigma in C[b]` be squarefree of degree
`D`, and put

```text
Y: c^N e=Sigma(b),
{b,c}=c^N,       {c,e}=-Sigma'(b),
{b,e}=-N c^(N-1)e.                                      (1)
```

For

```text
Q=delta+lambda b+mu c+nu e,                             (2)
```

there are `P in C(Y)` and `kappa in C*` with `{P,Q}=kappa` if and only if
`nu=0` and `Q` is nonconstant.  Scaling `P` reduces `kappa` to one.  When
`nu=0`, exact mates are

```text
lambda!=0:          P=c^(1-N)/[lambda(N-1)],
lambda=0, mu!=0:   P=b/(mu c^N).                        (3)
```

Direct substitution in `(1)` gives `{P,Q}=1`.  The additive constant
`delta` is irrelevant.  It remains to exclude every `nu!=0` row.

## 2. Generic fibre and forced time form

Translate `b` so that `Sigma(0)!=0`; the resulting constant in `(2)` is
absorbed into `delta`.  On `c!=0`, eliminate `e` and suppress `delta`:

```text
F(b,c)=lambda b+mu c+nu Sigma(b)c^-N,
{R,S}=c^N(R_bS_c-R_cS_b).                               (4)
```

On the generic fibre `F=w`, the forced Hamiltonian time form is

```text
eta=db/(c^N F_c)=-dc/(c^N F_b),        eta(D_F)=1.      (5)
```

Equivalently, with

```text
G_w=c^N(w-lambda b-mu c)-nu Sigma(b),                  (6)
```

one has `eta=-db/(G_w)_c=dc/(G_w)_b`.  Thus `{P,Q}=kappa`
would imply `dP=kappa eta` on the generic fibre.

The generic fibre is geometrically integral.  For any root `beta` of
`Sigma`, the arm `L_beta={b=beta,c=0}` maps isomorphically to the `Q`-line:

```text
Q|_(L_beta)=delta+lambda beta+nu e.                     (7)
```

The arm valuation is trivial on `C(Q)` and has residue field exactly
`C(Q)`.  Its restriction to a finite algebraic intermediate field is still
trivial, while residue embeds that field into `C(Q)`.  Hence `C(Q)` is
algebraically closed in `C(Y)`; characteristic zero makes the extension
regular.

## 3. The Newton point

The Laurent residue attached to a lattice point `(i,j)` is

```text
b^(i-1)c^(j-1) db/F_c.
```

Therefore `(5)` corresponds to

```text
p=(1,1-N).                                               (8)
```

Put

```text
A=(0,-N), B=(D,-N), O=(0,0), L=(1,0), U=(0,1).
```

The four coefficient masks have polygons

| `(lambda,mu)` | Newton polygon |
|---|---|
| `(0,0)` | `conv(A,B,O)` |
| `(!=0,0)` | `conv(A,B,L,O)` |
| `(0,!=0)` | `conv(A,B,U)` |
| `(!=0,!=0)` | `conv(A,B,L,U)` |

For the pure-`e` triangle, the slanted margin at `p` is

```text
ND-D-N=(N-1)(D-1)-1.                                   (9)
```

It is positive for `N,D>=2` except at `(N,D)=(2,2)`, where it is zero.
The remaining strict margins in the `lambda`-only polygon include

```text
(D-1)(N-1)>0,             N-1>0,                       (10)
```

and the `mu`-only slanted margin is

```text
N(D-1)-1>0.                                             (11)
```

The both-nonzero polygon contains the `lambda`-only polygon.  Hence `p` is
strictly interior in every row except

```text
lambda=mu=0,                  (N,D)=(2,2).              (12)
```

## 4. Every boundary face, including the resonance

All toric faces are nondegenerate except for one explicitly resolvable
resonance.

1. The bottom face is `Sigma(b)`, nondegenerate because `Sigma` is
   squarefree.
2. The left face is, up to a monomial,

   ```text
   nu Sigma(0)+mu c^(N+1)-w c^N.
   ```

   A repeated torus root would force an algebraic equation making the
   transcendental `w` constant.
3. All remaining faces are binomials unless

   ```text
   lambda mu!=0,                  D=N+1.                (13)
   ```

In `(13)`, with `s=lc(Sigma)`, the exceptional face is

```text
h(t)=nu s t^(N+1)+lambda t+mu.                          (14)
```

It can have a repeated root for special coefficients, but this does not
create a pole.  Homogenize `(6)` to degree `N+1` and use
`t=B/C,z=Z/C`.  Locally,

```text
g(t,z)=wz-lambda t-mu
       -nu sum_(j=0)^(N+1) a_j t^j z^(N+1-j).          (15)
```

At a repeated root `t_0` of `(14)`,

```text
g_z(t_0,0)=w-nu a_N t_0^N !=0,
h''(t_0)=N(N+1)nu s t_0^(N-1) !=0.                     (16)
```

Thus the root is exactly double.  With `u=t-t_0`,

```text
z=u^2 unit,       g_t=u unit,
eta=-z^(N-2) dz/g_t=u^(2(N-2)) unit du.                (17)
```

So the resonant form is holomorphic, and is a unit when `N=2`.  The
subleading `wz` term is the missing sidecar that resolves the degenerate
face.

At each affine arm, `(6)` gives

```text
(G_w)_b(beta,0)=-nu Sigma'(beta)!=0,
eta=dc/(G_w)_b,                                         (18)
```

so the form is a unit there too.  Toric adjunction now makes `eta` a
nonzero holomorphic differential on the smooth projective normalization
whenever `p` is interior.

## 5. The unique logarithmic boundary

In `(12)`, the generic fibre is the conic

```text
w c^2=nu Sigma(b).                                      (19)
```

Let `s=lc(Sigma)`.  At its two infinity points put `t=b/c=+/-rho`, where
`rho^2=w/(nu s)`.  Since `eta=dc/[-nu Sigma'(b)]`,

```text
res_(+rho) eta= 1/(2nu s rho)=rho/(2w),
res_(-rho) eta=-1/(2nu s rho).                          (20)
```

Both residues are nonzero.  Thus `(19)` is logarithmic rather than
holomorphic, but still not rationally exact.

In the interior cases, a nonzero holomorphic differential on a complete
smooth curve cannot equal `dP`: a pole of `P` would give a pole of `dP`,
while a pole-free rational function is constant.  In `(19)`, exact rational
differentials have zero residues.  Both cases contradict `(5)`, proving the
`nu!=0` half.

## 6. Sharp boundaries and ledger

1. `N>=2` is essential.  At `N=1`, `Q=nu e` has mate
   `P=-b/(nu e)`.
2. `D>=2` is essential.  If `Sigma=s(b-beta)`, then `Q=nu e` has mate
   `P=-c/(nu s)`.
3. Squarefreeness is the smoothness and bottom-face boundary; no singular
   extension is claimed.
4. The exceptional conic is not a positive row: its nonzero residues are
   exactly the obstruction.

```text
source       affine-linear Q on the smooth Danielewski surface
target       complete normalization of the generic Q-fibre
map          Laurent time form eta and Newton lattice point p
preserved    constant bracket iff eta has a rational primitive
destroyed    affine arm and infinity charts before compactification
sidecar      coefficient mask, resonant face, and logarithmic residues
cheap test   interior margin; if zero, compute the two infinity residues
```

## 7. Exact verification contract

The companion checks:

- the Poisson contraction and both explicit `nu=0` mates;
- every mask for `2<=N,D<=20`, its hull, and the unique boundary row;
- the symbolic margins `(9)--(11)`;
- the resonant derivative and local-order identities `(16)--(17)`;
- the opposite nonzero residues `(20)`;
- the `N=1` and `D=1` hostiles.

The grid is an exact control of the universal inequalities, not an
extrapolation.  Reproduce with

```bash
python3 04-computation/jc2_universal_danielewski_affine_linear_rational_classification_thm3599.py
python3 -O 04-computation/jc2_universal_danielewski_affine_linear_rational_classification_thm3599.py
```
