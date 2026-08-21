# Two-cube support density and the parabolic/hyperbolic split

Status: **FINITE-EXACT PROBE + RIGOROUS SIEVE-DIMENSION DERIVATION +
CONDITIONAL ASYMPTOTIC SYNTHESIS** (2026-08-21).  This reflection is not a
theorem and creates no proof dependency.  Its exact input is THM-3645.  The
upper-bound mechanism in Section 3 is rigorous; the full asymptotic still
requires the two stated correlation-theorem interfaces.

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

## 3. The sieve dimension is rigorously `3/2`

For `p>=5`, put `chi_6=(6/p)=chi_24` and
`chi_-2=(-2/p)=chi_-8`.  The bad residue set in `F_p^2` is the union of

```text
n=0 and V=0  if chi_24(p)=-1,
2n+V=0       if chi_-8(p)=-1.                          (5)
```

These three lines are distinct.  If `rho(p)` is their union size, then

```text
(chi_24,chi_-8)   (+,+)  (+,-)  (-,+)  (-,-)
rho(p)                 0       p    2p-1    3p-2.      (6)
```

Writing `e_D(p)=1_(chi_D(p)=-1)`, one has

```text
rho(p)/p^2=(2e_24(p)+e_-8(p))/p+O(1/p^2),
2e_24+e_-8=3/2-chi_24-chi_-8/2.                        (7)
```

The prime number theorem in the fixed residue classes modulo `24` therefore
gives

```text
sum_(p<=z) rho(p) log(p)/p^2=(3/2)log z+O(1).          (8)
```

This is exactly the coupled sieve-dimension axiom with dimension `3/2`.
It must be applied once to the union `(5)`, not as three formally independent
half-dimensional sieves.

For squarefree `d`, CRT gives at most `3^omega(d)d` bad residue classes.
Elementary counting in `0<V<n<=N` has total remainder

```text
sum_(d<=D)|r_d| << ND(log D)^2+D^2(log D)^2.           (9)
```

Taking `D=N^(1/3)` and then a fixed-level fundamental-lemma sieve proves the
rigorous upper bound

```text
A(N) << N^2/(log N)^(3/2).                            (10)
```

Dropping primitivity and the conditions at `2,3` only enlarges the set, so
they create no gap in this upper bound.  A maintained theorem package and a
primary-source citation audit remain to be written before `(10)` is used as a
canon dependency.

## 4. The coupled constant and the honest asymptotic boundary

Selberg--Delange gives, for the one-form support indicator attached to a
quadratic character `chi_D`, the constant

```text
C_D=1/sqrt(pi) * [ L(1,chi_D)
      product_(chi_D(p)=-1)(1-p^-2)
      product_(p|D)(1-p^-1) ]^(1/2).                  (11)
```

Here

```text
C_24  approximately 0.30845624,
C_-8  approximately 0.40503888,
L(1,chi_24)=log(5+2sqrt(6))/sqrt(6),
L(1,chi_-8)=pi/(2sqrt(2)).                             (12)
```

After dividing the actual local density by the three marginal densities,
only the equal-sign Frobenius classes contribute.  The absolutely convergent
correction is

```text
mathfrak C=
 product_(p>=5,chi_24=chi_-8=+1)(1-p^-2)
 * product_(p>=5,chi_24=chi_-8=-1)
       (1-2/p)/(1-1/p)^2
 approximately 0.95871001.                            (13)
```

The triangular area supplies `1/2`; the coupled local condition at `p=3`
supplies another `1/2`; the apparent mod-`24` cone factor is already coupled
to the third form's mod-`8` support and must not be multiplied again.  The
predicted constant is therefore

```text
kappa=(1/4) C_24^2 C_-8 mathfrak C
      approximately 0.0092365781.                     (14)
```

At `N=49999`, the exact normalized count is `0.0090176366`, or
`0.9762963` times `(14)`.  The agreement is strong evidence but not a proof
of

```text
A(N) ~ kappa N^2/(log N)^(3/2).                       (15)
```

The elementary sieve proving `(10)` cannot supply the matching lower bound.
A full asymptotic appears to fit the geometry of linear-correlation theorems
for small-mean multiplicative functions: the two residue cones give the
pairwise independent triples

```text
(24a+5, 24b+23, 16a+8b+11),
(24a+23,24b+5,  16a+8b+17).                           (16)
```

Two interfaces remain explicit: verify the fixed-character support
indicators uniformly in the required multiplicative-function class,
including exceptional-character terms; and add primitivity by truncated
Moebius inversion, using `(10)` to bound the common-prime tail.  Until those
are written and the primary theorem is checked, `(15)` and `(14)` are
**CONDITIONAL / CONJECTURAL**, not proved.

## 5. Parabolic screens versus hyperbolic existence

The user's Berggren-tree branch supplies a useful dynamical distinction.  On
the Pythagorean `U`-spine the generating matrix is unipotent; `N=U-I` has
`N^3=0`, so orbit coordinates grow quadratically.  The fixed-screen hostile
from THM-3645 is analogous.  With

```text
M=product_(p<=997)p,
U_s=507+2Ms,          V_s=5,
m_s=512+2Ms,          n_s=1019+4Ms,                  (17)
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
hyperbolic/Pell stage:     certify an infinite positive ray.            (18)
```

Confusing the stages causes both common errors: a fixed finite screen is not
global, and a locally viable slope is not an admissible Pell orbit.

## 6. Crisp next targets

1. Package `(5)--(10)` as a maintained upper-bound theorem with a
   primary-source citation and an independently audited remainder estimate.
2. Verify the multiplicative-class and exceptional-character interfaces for
   the two support indicators; then derive `(15)` and `(14)` from the precise
   linear-correlation theorem.
3. Extend the exact census by residue cone and fixed-`V` slice to identify
   secondary terms or exceptional concentrations.
4. Measure the conversion from local survivors to complete Pell/compiler
   survivors beyond denominator `401`; the existing `104 -> 42` conversion
   is one bounded datum, not a density.
5. Search for a Vieta or tree operation on `(U,V)` that preserves the full
   support semigroups and the compiler address.  Preserving only finitely many
   residues, as in `(17)`, is an automatic hostile control.

The main transferable point is that the sparse set of two-distinct-cube
slopes appears to be governed first by a *fractional sieve dimension* and only
then by Pell dynamics.  That is a sharper organizing principle than either a
raw prime screen or a raw tree analogy alone.
