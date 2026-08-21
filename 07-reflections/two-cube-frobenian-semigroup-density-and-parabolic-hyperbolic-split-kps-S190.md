# Two-cube support density and the parabolic/hyperbolic split

Status: **FINITE-EXACT PROBE + CONJECTURAL ASYMPTOTIC SYNTHESIS**
(2026-08-21).  This is not a theorem and creates no proof dependency.
Its exact input is THM-3645; its asymptotic and transfer statements are
explicitly conjectural.

## 1. The local screen is three correlated prime semigroups

For a primitive parity-correct slope put

```text
U=n-m,                 V=2m-n,                 T=2m+n. (1)
```

THM-3645 proves that the slope conic is soluble modulo every prime exactly
when the exceptional `p=3` condition holds and

```text
p|nV  implies (6/p)=1,
p|T, p!=3 implies (-2/p)=1.                            (2)
```

After `p=3`, one has `3|T`, so define the multiplicative semigroups

```text
S_6   ={a>=1: p|a implies (6/p)=1},
S_-2  ={a>=1: p|a implies (-2/p)=1}.                   (3)
```

Then the complete mod-prime gate is the correlated three-form condition

```text
n in S_6,                 V in S_6,
(2n+V)/3 in S_-2,                                     (4)
```

together with `0<V<n`, `gcd(n,V)=1`, `n+V=0 mod 4`, and
`n=V!=0 mod 3`.  Thus the old prime-by-prime scan is not intrinsically a
168-coordinate object.  It is a two-variable lattice problem for three
Frobenian multiplicative semigroups.

## 2. Exact census through denominator 49,999

Let `A(N)` count slopes satisfying the complete mod-prime gate with odd
denominator `n<=N`.  The vectorized exact companion gives:

| `N` | primitive parity candidates | `p=3` eligible | `A(N)` | `A(N)(log N)^(3/2)/N^2` |
|---:|---:|---:|---:|---:|
| 401 | 8,195 | 2,053 | 104 | 0.00949108 |
| 997 | 50,507 | 12,654 | 481 | 0.00877964 |
| 1,999 | 202,713 | 50,741 | 1,668 | 0.00874632 |
| 4,999 | 1,266,692 | 316,812 | 8,941 | 0.00889302 |
| 9,999 | 5,065,838 | 1,266,211 | 31,859 | 0.00890688 |
| 19,999 | 20,264,681 | 5,066,584 | 114,230 | 0.00890108 |
| 49,999 | 126,652,918 | 31,664,297 | 633,416 | 0.00901764 |

The first two rows independently recover the maintained THM-3640/3645
counts.  The last column is not used in any exact decision.  Its stability is
the new signal.

The exact script is

```text
04-computation/berggren_two_cube_support_density_kps_s190.py
SHA256 d929fe021f371526a74282be89a0ff2b7e8ad85766f0703503aaad43e92e9f3f

05-knowledge/results/berggren_two_cube_support_density_kps_s190.out
SHA256 bfe5927eba9ea6a62116a481265f39e31ed69eb32ee5e5befe218accd0e4bba3
```

and prints only the integer ledger and its semantic digest.  The displayed
logarithmic normalization is a diagnostic computed from that ledger.

## 3. Frobenian `3/2`-law conjecture

Each allowed prime set in `(3)` has Dirichlet density `1/2`.  The Dirichlet
series of integers supported on one such set therefore has a square-root
singularity at `s=1`; by the usual Landau--Selberg--Delange heuristic, one
form contributes a factor `(log N)^(-1/2)`.  There are three genuinely
different forms in `(4)`.  This predicts

```text
A(N) ~ kappa N^2/(log N)^(3/2),                        (5)
```

with a positive singular-series constant apparently near

```text
kappa approximately 0.009.                             (6)
```

Equations `(5)--(6)` are **CONJECTURAL**.  The linear relation among the three
forms, the coprimality condition, and the exceptional primes `2,3` all enter
the constant.  Multiplying three one-variable constants would ignore exactly
the correlations that matter.

A plausible rigorous first target is only the upper bound

```text
A(N) << N^2/(log N)^(3/2).                             (7)
```

This should be approachable as a half-dimensional Frobenian sieve on each of
the three forms.  A full asymptotic needs a two-variable Selberg--Delange or
comparable correlation theorem and a convergent local singular series.  A
lower bound may meet a genuine parity/correlation barrier; the finite table
does not decide that issue.

## 4. Parabolic screens versus hyperbolic existence

The user's Berggren-tree branch supplies a useful dynamical distinction.  On
the Pythagorean `U`-spine the generating matrix is unipotent; `N=U-I` has
`N^3=0`, so orbit coordinates grow quadratically.  The fixed-screen hostile
from THM-3645 is analogous.  With

```text
M=product_(p<=997)p,
U_s=507+2Ms,          V_s=5,
m_s=512+2Ms,          n_s=1019+4Ms,                   (8)
```

translation in `s` preserves every residue seen by the fixed screen.  This is
a parabolic/affine orbit: it preserves a finite congruence packet while new
prime support enters at larger scale.  Its first member already fails at the
new support prime `1019`.  Thus unipotent residue preservation is a powerful
hostile generator, not a global-solubility mechanism.

The successful positive two-cube rays of THM-3375/3640 use the opposite
dynamics.  Once a slope survives complete generalized-Pell class and compiler
address tests, multiplication by a norm-one unit is hyperbolic: coordinates
grow exponentially while a positive cone is invariant.  The proof architecture
is therefore naturally two-stage:

```text
parabolic/Frobenian stage: classify or count viable slope support,
hyperbolic/Pell stage:     certify an infinite positive ray.             (9)
```

Confusing the stages causes both common errors: a fixed finite screen is not
global, and a locally viable slope is not an admissible Pell orbit.

## 5. Crisp next targets

1. Derive the Euler product predicted for `kappa`, including the coupled
   local factors for `(n,V,2n+V)` and the primitive/parity conditions.
2. Prove `(7)` by a three-form Frobenian upper-bound sieve.
3. Extend the exact census by residue cone and fixed-`V` slice to identify
   secondary terms or exceptional concentrations.
4. Measure the conversion from local survivors to complete Pell/compiler
   survivors beyond denominator `401`; the existing `104 -> 42` conversion
   is one bounded datum, not a density.
5. Search for a Vieta or tree operation on `(U,V)` that preserves the full
   support semigroups and the compiler address.  Preserving only finitely many
   residues, as in `(8)`, is an automatic hostile control.

The main transferable point is that the sparse set of two-distinct-cube
slopes appears to be governed first by a *fractional sieve dimension* and only
then by Pell dynamics.  That is a sharper organizing principle than either a
raw prime screen or a raw tree analogy alone.
