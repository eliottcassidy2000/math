---
id: THM-4396
title: "LRC14 finite dual and exact-pair hybrid certificate"
status: >
  PROVED ANALYTICALLY RELATIVE TO THM-4392 + RIGOROUS-INTERVAL VERIFIED +
  INDEPENDENTLY AUDITED. Two-coordinate Fejer smoothing plus the exact
  ordered-pair primal remainder certifies mu(F_(11,13,17))<6/77 from eleven
  resonance sites. The positive-part quotient is sharp after forgetting the
  third incidence and is strictly slack on (1,5,11) for every coordinate pair
  and every pair of finite Fejer degrees. Universal nonresonance, seam entry,
  synchronization, and LRC(14) remain OPEN.
source: root + lrc_defect3_extra_cleanroom + independent referee / LRC14 continuation session, 2026-09-03
depends_on:
  - THM-4392-lrc14-raw-carrier-box-spline-fourier-poisson-duality
related:
  - THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence
  - THM-4398-lrc14-one-zero-relation-residue-dichotomy-and-small-norm-atlas
primary_script: 04-computation/lrc14_finite_dual_exact_pair_hybrid_certificate_thm4396.py
primary_output: 05-knowledge/results/lrc14_finite_dual_exact_pair_hybrid_certificate_thm4396.out
primary_script_sha256: 1b35c8a4e9ff27b6b5b8a399b9853f706755f9ce08e65d020cb91e6f83fa25d4
primary_output_sha256: 16146a7469eff6abef16b29d129b9338cd798193639ca2087f2e8e87817fab6b
independent_audit_script: 04-computation/lrc14_finite_dual_exact_pair_hybrid_certificate_thm4396_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_finite_dual_exact_pair_hybrid_certificate_thm4396_independent_audit.out
independent_audit_script_sha256: 7862e0b5f17698c2759529b7de0021e39b4f4a53b39eb8ccdbcabedf21d13044
independent_audit_output_sha256: 7b5e01b5cc3b6c52a4855197b0839357f063caf95c5808263b6572d562901384
hash_basis: raw LF bytes
audit: >
  PASS. The standalone referee imports no producer or repository mathematics,
  reconstructs the finite identity, pair geometry, interval arithmetic, raw
  controls, and all-degree equality obstruction, then compares with the
  producer only after freezing its own result. Normal, optimized, and
  hash-seeded runs agree. No mathematical, numerical, or status discrepancy
  remains.
---

# THM-4396 -- finite dual and exact-pair hybrid LRC certificate

**Status: PROVED ANALYTICALLY RELATIVE TO THM-4392 + RIGOROUS-INTERVAL
VERIFIED + INDEPENDENTLY AUDITED.**  The hybrid inequality is correct, the
pointwise quotient is sharp after the stated loss of third-sheet information,
and the concrete `(11,13,17)` certificate lies strictly below `6/77` with a
fully rational outward-rounded enclosure.  The strict equality-comb
obstruction also holds for every coordinate pair and every pair of finite
Fejer degrees.

This is only a local three-speed scale-three statement.  It proves no
universal local bound, synchronization, seam entry, or all-tail theorem.
**LRC(14) remains OPEN.**

## 1. Definitions and the finite hybrid term

On `T=R/Z`, with normalized Haar measure, set

```text
lambda=1/14,
h(x)=1_{||x||<lambda},
hhat(0)=1/7,
hhat(n)=sin(pi n/7)/(pi n) for n!=0.
```

For a speed `a` prime to three and a sheet `s in F_3`, write

```text
f_(a,s)(x)=h(a(x+s/3)).                                  (1)
```

Let `K_H` be the normalized Fejer kernel,

```text
tau_H(n)=(1-|n|/(H+1))_+,
h_H=K_H*h,
g_(a,s;H)(x)=h_H(a(x+s/3)).                              (2)
```

Positivity and unit mass of `K_H`, together with `0<=h<=1`, give

```text
0<=g_(a,s;H)<=1.                                        (3)
```

For distinct coordinate indices `i,j,k`, smooth only the first two chosen
coordinates and leave the third indicator exact:

```text
M_(i,j)|k(H,K)
 =sum_(pi in S_3) integral_T
   g_(w_i,pi(i);H) g_(w_j,pi(j);K) f_(w_k,pi(k)).        (4)
```

Both `g` factors have finite Fourier support.  After expanding those two
finite polynomials, the integral against the remaining exact indicator is
one of its ordinary Fourier coefficients.  Thus no infinite interchange is
needed, and character orthogonality yields exactly

```text
M_(i,j)|k(H,K)
 =sum A_w(n) tau_H(n_i) tau_K(n_j)
       product_l hhat(n_l),                              (5)

w.n=0,       |n_i|<=H,       |n_j|<=K,
n_k=-(w_i n_i+w_j n_j)/w_k in Z.                        (6)
```

Here the sum of the six sheet phases is

```text
A_w(n)=6  if w_1n_1,w_2n_2,w_3n_3 are equal mod 3,
       -3 if they are all distinct mod 3.                (7)
```

The resonance equation makes these the only two residue cases.

## 2. Positive-part remainder and quotient sharpness

For one permutation abbreviate its factors by `f_i,f_j,f_k,g_i,g_j` and put

```text
X=f_i f_j-g_i g_j.                                      (8)
```

The exact two-coordinate replacement identity is

```text
f_i f_j f_k=g_i g_j f_k+f_k X.                         (8a)
```

Thus all approximation error is isolated in the one scalar product `f_k X`;
there is no hidden Fourier tail in the identity itself.  Since
`0<=f_k<=1`, pointwise

```text
f_k X<=X_+.                                             (9)
```

Moreover `f_i f_j` is zero or one while `0<=g_i g_j<=1`, so

```text
X_+=f_i f_j(1-g_i g_j).                                (10)
```

Summing (9) over the six permutations proves

```text
mu(F_w)<=M_(i,j)|k(H,K)+R_(i,j)(H,K),                   (11)

R_(i,j)(H,K)
 =sum_(s,t in F_3, s!=t)
   integral_(D_(w_i,s) intersect D_(w_j,t))
      (1-g_(w_i,s;H)g_(w_j,t;K)) dx.                   (12)
```

The word “ordered” in (12) matters: each `(s,t)` determines the unique
remaining sheet assigned to coordinate `k`.

This positive-part envelope is exactly sharp in the quotient that forgets
the third incidence.  Once `X` is fixed and the only retained information is
`0<=z<=1`,

```text
sup_(0<=z<=1) zX=X_+.                                  (13)
```

Therefore a strictly smaller pointwise majorant cannot be obtained from the
same pair data alone.  This does not claim optimality among certificates that
retain any third-sheet sidecar, couple different pair choices, or use a
different positive spectral taper.

## 3. Exact finite spectrum for `(11,13,17)`

Take `w=(11,13,17)`, smooth coordinates `(11,13)` with degrees `(5,9)`, and
leave speed 17 exact.  Direct enumeration of (6) gives exactly

```text
(-5,-1, 4)  (-4, 6,-2)  (-3,-4, 5)  (-2, 3,-1)
(-1,-7, 6)  ( 0, 0, 0)  ( 1, 7,-6)  ( 2,-3, 1)
( 3, 4,-5)  ( 4,-6, 2)  ( 5, 1,-4).                    (14)
```

There are three sites of character weight `6` and eight of weight `-3`.
Exactly `+/-(-1,-7,6)` vanish: their second coordinate is a **nonzero**
multiple of seven.  The zero frequency survives because `hhat(0)=1/7`.
Consequently the eleven-site spectrum has nine live sinc products.

The independently certified interval is

```text
M=[0.012578340838550376,
   0.012578340838550377].                               (15)
```

## 4. Exact ordered-pair geometry and remainder

For speeds `(11,13)`, literal rational construction of the shifted danger
arcs shows that each of the six ordered distinct-sheet pairs has four
intersection pieces and exact measure `20/1001`.  Therefore

```text
number of pieces=24,
total ordered-pair measure=120/1001.                    (16)
```

On one rational piece `[l,r]`, expand the two finite Fejer polynomials.  For
integers `n,m`, set

```text
q=11n+13m,
p=11ns+13mt.
```

The real exponential integral needed for each term is exactly

```text
integral_l^r cos(2 pi q x+2 pi p/3) dx
 =(r-l)cos(2 pi p/3),                                  q=0,
 =[sin(2 pi q r+2 pi p/3)-sin(2 pi q l+2 pi p/3)]
   /(2 pi q),                                           q!=0.  (17)
```

There is no quadrature and no spectral tail in this pair calculation.  The
rational interval evaluation gives

```text
integral_pair g_11 g_13
 =[0.054902598885952923,
   0.054902598885952924],

R=120/1001-integral_pair g_11 g_13
 =[0.064977520994166956,
   0.064977520994166957].                               (18)
```

Combining (15) and (18), with outward rounding after the exact interval
addition, proves

```text
M+R=[0.077555861832717332,
     0.077555861832717333],

6/77-(M+R)
   =[0.000366216089360589,
     0.000366216089360590] >0.                          (19)
```

Thus this finite dual/exact-pair primal certificate rigorously proves

```text
mu(F_(11,13,17))<6/77.                                  (20)
```

## 5. Independent transcendental enclosure

All geometry and all Fejer weights are exact `Fraction`s.  Transcendental
values are enclosed without binary floating point:

1. Machin's identity
   `pi=16 atan(1/5)-4 atan(1/239)` is bounded by consecutive alternating
   rational partial sums (through indices 50 and 12).
2. That fine enclosure is rounded outward to a rational interval with
   denominator `10^24`; it remains strictly inside the classical convergent
   bounds `103993/33102 < pi < 104348/33215`.
3. Every trigonometric argument is reduced exactly modulo `2*pi` and by
   symmetry to `pi*x` with `0<=x<=1/2`.
4. Sine and cosine are evaluated by rational interval Taylor polynomials
   through powers 29 and 28.  Since the reduced squared argument is below
   three, the omitted alternating tails are bounded respectively by the
   next `z^31/31!` and `z^30/30!` terms.
5. The displayed decimals are then obtained by exact floor/ceiling at 18
   decimal places.

The intervals in (15), (18), and (19) are therefore mathematical enclosures,
not high-precision guesses.

## 6. Sheet-blind fallback and why it loses

Let

```text
epsilon_H=||h-h_H||_1.                                  (21)
```

On the danger interval `h-h_H>=0`, while outside it `h-h_H<=0`.  Since `h`
and `h_H` have the same integral,

```text
integral (h-h_H)_+=epsilon_H/2.                         (22)
```

Telescoping the two replacements in (4), bounding the remaining factors by
one, and summing six sheets gives the valid but sheet-blind estimate

```text
mu(F_w)<=M_(i,j)|k(H,K)+3(epsilon_H+epsilon_K).          (23)
```

Fourier integration over the danger interval gives the exact finite formula

```text
epsilon_H=2[1/7-1/49
 -2 sum_(n=1)^H tau_H(n) sin(pi n/7)^2/(pi^2 n^2)].      (24)
```

The same independent interval engine yields

```text
epsilon_5=[0.108950255634203925,0.108950255634203926],
epsilon_9=[0.075930178359372161,0.075930178359372162],

M+3(epsilon_5+epsilon_9)
 =[0.567219642819278636,
   0.567219642819278637].                               (25)
```

This is far above `6/77`.  The successful information is the exact ordered
pair domain in (12), not merely a higher Fourier cutoff.

## 7. Exact physical controls

A separate nearest-integer sweep on `y=3x` constructs physical raw carriers
without using the hybrid sum.  It gives

```text
w=(11,13,17):
  +/-(-5,-1,4), each of length 10/1547,
  mu(F_w)=20/1547;

w=(1,5,11):
  +/-(-1,-2,1), each of length 3/77,
  mu(F_w)=6/77.                                         (26)
```

The first verifies that the certificate is a genuine upper bound rather than
an evaluation of the target.  The second is the equality hostile used next.

## 8. Strict finite-degree obstruction at the equality comb

For any finite degree `H`, `g_(a,s;H)(x)>0` at every `x`.  Indeed, the Fejer
kernel is nonnegative and is positive away from finitely many points; its
integral over the translated danger interval of positive length is strictly
positive.

Fix the sheet assignment `(0,1,2)` for speeds `(1,5,11)`.  For each possible
smoothed coordinate pair, the following nonempty open interval lies inside
the omitted exact sheet and outside one of the two corresponding unsmoothed
danger indicators:

| smoothed speeds | exact speed | open interval | excluded pair sheet |
|---|---:|---:|---|
| `(1,5)` | 11 | `(67/462,73/462)` | `D_(1,0)` |
| `(1,11)` | 5 | `(1/14,17/210)` | `D_(1,0)` |
| `(5,11)` | 1 | `(0,11/210)` | `D_(5,1)` |

These inclusions and exclusions are checked directly against the complete
rational sheet-interval lists.  On any one of these intervals,

```text
f_k=1,       f_i f_j=0,       X=-g_i g_j<0.             (27)
```

The pointwise slack in (9) is

```text
X_+-f_kX = (1-f_k)X,  X>=0,
          =-f_kX,     X<0.                              (28)
```

It is therefore strictly positive throughout an open interval, for every
pair of finite degrees.  Integration proves the scoped statement

```text
M_(i,j)|k(H,K)+R_(i,j)(H,K)
 >mu(F_(1,5,11))=6/77                                  (29)
```

for all three coordinate pairs and all finite nonnegative integers `H,K`.
This obstructs only the stated pair-Fejer/positive-part quotient from proving
the equality case.  It does not obstruct a certificate retaining third-sheet
incidence, another positive taper, a coupled-pair argument, or a raw primal
proof.

## 9. Reproducibility and firewall

The proof verifier uses only the Python standard library.  All theorem
comparisons use exact rational arithmetic or rigorously outward-rounded
rational intervals, and every check raises explicitly under `python -O`.

```powershell
python -B 04-computation/lrc14_finite_dual_exact_pair_hybrid_certificate_thm4396.py
python -B -O 04-computation/lrc14_finite_dual_exact_pair_hybrid_certificate_thm4396.py
python -B 04-computation/lrc14_finite_dual_exact_pair_hybrid_certificate_thm4396_independent_audit.py
python -B -O 04-computation/lrc14_finite_dual_exact_pair_hybrid_certificate_thm4396_independent_audit.py
$env:PYTHONHASHSEED='4392'
python -B -O 04-computation/lrc14_finite_dual_exact_pair_hybrid_certificate_thm4396_independent_audit.py
```

The normal, optimized, and hash-seeded streams are byte-identical to the
frozen canonical outputs.  The referee output records 7,706 live checks and
a semantic hash of the exact interval endpoints and finite objects.  The four
raw-LF hashes are frozen in the theorem metadata.

The certificate is finite and local.  No claim in this report closes the
universal three-speed problem or supplies any missing global LRC14 step.
**LRC(14) remains OPEN.**

## 10. Post-clean-room comparison with the primary verifier

Only after the derivation, verifier, frozen output, and three replay modes
above were complete did the referee inspect the primary verifier and output.
Its replay also passes byte-for-byte in normal, optimized, and fixed-hash-seed
modes.  Every primary numerical enclosure agrees endpoint-for-endpoint with
the independent computation:

```text
M                  [0.012578340838550376,0.012578340838550377]
pair product       [0.054902598885952923,0.054902598885952924]
R                  [0.064977520994166956,0.064977520994166957]
M+R                [0.077555861832717332,0.077555861832717333]
6/77-(M+R)         [0.000366216089360589,0.000366216089360590]
sheet-blind bound  [0.567219642819278636,0.567219642819278637]
```

The finite spectrum, `3/8` split of the two character weights, two nonzero
mod-seven sinc zeros, 24 ordered pair pieces, pair measure `120/1001`, raw
physical masses `20/1547` and `6/77`, and all three hostile open intervals
also agree exactly.

The implementations are genuinely different.  This verifier uses
`Fraction` interval endpoints, a longer Machin enclosure, direct full-mode
sums rather than conjugate pairing, and a literal nearest-integer physical
control in `y=3x`; the producer uses fixed-point interval arithmetic,
conjugate pairing, and a raw-carrier enumeration.  This verifier additionally
checks four pieces and measure `20/1001` for each ordered sheet pair and emits
the two individual `epsilon` intervals.  The differing check counts and
semantic hashes therefore are expected.

The referee's sole wording recommendation has been incorporated in the
canonical primary output: the scope is “every coordinate pair and every pair
of finite Fejer degrees,” not arbitrary spectral cutoffs.  No discrepancy
remains.  In particular, the negative result is confined to this
pair-Fejer/positive-part quotient and must not be generalized to all hybrid
or primal certificates.
