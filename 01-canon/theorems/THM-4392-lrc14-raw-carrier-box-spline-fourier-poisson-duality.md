# The raw LRC carrier and the old Fourier resonance lattice are Poisson dual

**Status.**  The box-spline/Poisson identity below is **PROVED ANALYTICALLY
RELATIVE TO THM-1092 + VERIFIED-EXACT + INDEPENDENTLY AUDITED**, using
the absolute-convergence theorem
[`THM-1092-fejer-regularized-kfold-resonance-identity.md`](../../01-canon/theorems/THM-1092-fejer-regularized-kfold-resonance-identity.md)
to justify the dual series.  The finite character table and five rational
physical controls are **VERIFIED-EXACT**.  The displayed Fejer truncations are
only **NUMERIC DIAGNOSTICS**.  The raw-carrier theorem
[`THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence.md`](../../01-canon/theorems/THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence.md)
is canonical; the proof below repeats the address calculation so every
normalization in the duality remains explicit.
**LRC(14) remains OPEN.**

## Outcome

For a primitive distinct positive odd three-unit triple
`w=(w_1,w_2,w_3)`, both of the following exact expressions compute the same
scale-three failure comb:

```text
PRIMAL / COMPONENT FORM

sum_{C in Lambda_w, every C_i nonzero mod 3} L_w(C),

Lambda_w={C in Z^3:C.w=0};                              (1)
```

and

```text
DUAL / FOURIER FORM

sum_{n in Lambda_w} A_w(n) product_i hhat(n_i),

h=1_{||t||<1/14},
hhat(0)=1/7,
hhat(n)=sin(pi n/7)/(pi n) for n!=0,                    (2)

A_w(n)=sum_{pi in S_3}
  exp((2 pi i/3) sum_i w_i n_i pi(i))
 = 6   if the three w_i n_i are equal mod 3,
 =-3   if they are all distinct mod 3.                  (3)
```

The equality is a genuine planar box-spline Poisson formula.  It is not a
termwise identification of the two coordinate deletions:

- the raw `mod 3` gate selects two primal cosets in
  `Lambda_w/3 Lambda_w`;
- its finite Fourier transform is the signed sheet factor `A_w=6 or -3`;
- the box-spline Fourier transform, sampled on the one-third dual lattice,
  is `27 product hhat(n_i)`;
- the nonzero `mod 7` zeros come from `sin(pi n_i/7)`, hence from the interval
  width after the scale-three change of variables.

Thus the two primes have different jobs.  The mod-three deletion becomes a
dual **character weight**; it does not become the mod-seven deletion.  The
mod-seven deletion is an analytic zero of the transformed spline.  Also,
`7|0` but `hhat(0)=1/7`, so only **nonzero** multiples of seven are killed.

The bridge also gives an exact obstruction to bare shortest-vector control.
The two triples

```text
w =(11,13,17),       w'=(101,103,107)                   (4)
```

have all of the following identical data:

```text
w mod 3 = w' mod 3 = (2,1,2),
lambda_1(Lambda)^2=14,
complete Euclidean minimizing orbit = +/-(-2,3,-1),
m_7(Lambda)=6,
complete m_7 minimizing orbit       = +/-(-2,3,-1).
```

Nevertheless

```text
mu(F_w)  =20/1547,
mu(F_w') =36488/7791847,
difference=14198572/1721998187,
ratio     =5565605/2015962 = 2.7607688041... .           (5)
```

So even the scalar Euclidean minimum, the multiplicative minimum, the full
minimizer orbit, and its mod-three sheet sign do not determine the comb sum.
This is a local exact counterpart to the nonuniformity diagnosed in
[`THM-1085-multmin-nonuniform-sumfree-refuted.md`](../../01-canon/theorems/THM-1085-multmin-nonuniform-sumfree-refuted.md).

## 1. Inheritance and loss ledger

The closest proved mechanism is
[`THM-1075-kfold-resonance-lattice.md`](../../01-canon/theorems/THM-1075-kfold-resonance-lattice.md),
with its analytic repair in THM-1092: integrating a product of danger
indicators selects the integer resonance lattice.  The corrected near miss is
[`THM-1080-shortvector-multiplicative-minimum.md`](../../01-canon/theorems/THM-1080-shortvector-multiplicative-minimum.md)
as explicitly scoped by THM-1085: `m_7` is a leading scale, not a uniform
sharp bound, because near-minimal populations and signs remain.

The canonical hostile is `MISTAKE-226/MISTAKE-235`: a shared integer kernel is
only syntax unless the map preserves the weights, constraints, and target
predicate.  Here those obligations are discharged explicitly:

| field | raw side | Fourier side | connecting map / loss |
|---|---|---|---|
| lattice role | component carrier `C` | frequency `n` | planar duality `n=w cross eta`; same embedded integer set, different semantics |
| local weight | nonnegative compact roof `L_w(C)` | signed `A_w(n) product hhat(n_i)` | Fourier transform of the box spline |
| arithmetic gate | `C_i != 0 mod 3` | `A_w=6,-3` plus sinc zeros | finite Fourier transform of two mod-three cosets; no termwise match |
| deletion modulus | three | seven for nonzero coordinates | three is sampling/sheet data; seven is interval-width data |
| retained sidecar | raw scale and all three coordinates | residue character, support, signs, all frequencies | shortest length alone loses all of these |
| preserved target | exact local Haar measure | exact local Haar measure | equality (1)=(2), not an LRC entry theorem |

The least-used sidecar is the finite character of the two retained raw cosets.
It is precisely what turns a misleading shared-kernel observation into a typed
identity.

## 2. The local object and its direct Fourier expansion

Put `lambda=1/14`.  For a speed `a` prime to three and a sheet
`j in Z/3Z`, define

```text
D_(a,j)={x in R/Z: ||a(x+j/3)||<lambda}.
```

The local scale-three failure comb is, up to endpoints,

```text
F_w = disjoint union over pi in S_3 of
      intersection_i D_(w_i,pi(i)).                     (6)
```

The six pieces are disjoint almost everywhere because the three sheet
intervals for one three-unit speed are separated by `1/3`, whereas their
total pairwise radius is `1/7`.

Expand

```text
h(w_i(x+pi(i)/3))
 =sum_{n_i in Z} hhat(n_i)
   exp(2 pi i n_i w_i x)
   exp(2 pi i n_i w_i pi(i)/3).
```

Character orthogonality in `x` imposes `w.n=0`.  Summing the six sheets gives
exactly (2)-(3).  This is the shifted three-factor version of THM-1075.
Absolute convergence is legitimate: THM-1092 supplies it on full support;
one-zero-coordinate strata are rank-one pair sums of order `1/k^2`; two
zero coordinates force the zero vector.  Since `|A_w|<=6`, the sheet factor
does not alter convergence.

For `n in Lambda_w`, put `r_i=w_i n_i mod 3`.  Their sum is zero.  Over
`F_3`, a zero-sum triple is either all equal or all distinct.  Directly
summing the six third roots of unity gives (3).  Equivalently, if

```text
Lambda_eq={n in Lambda_w:
           w_1 n_1=w_2 n_2=w_3 n_3 mod 3},
```

an index-three sublattice, then the exact signed splitting is

```text
mu(F_w)
 =9 sum_{n in Lambda_eq} product_i hhat(n_i)
  -3 sum_{n in Lambda_w} product_i hhat(n_i).            (7)
```

The second sum is the unshifted THM-1075 triple intersection.  The first is
the congruence sidecar forced by the sheet assignment.  Equation (7) makes
the missing signed information visible; applying absolute values before this
split discards the mechanism.

## 3. The raw carrier is a two-coset box-spline sample

Set `r=3/14=3 lambda` and pass to `y=3x mod 1`.  For the unique nearest
integers `N_i` write

```text
e_i=w_i y-N_i,       |e_i|<r,
C=w cross N=-w cross e.                                 (8)
```

For primitive `w`, the cross-product matrix has Smith form
`diag(1,1,0)`: its entries have gcd one and its two-by-two minors have gcd
one.  Hence

```text
Z^3/Zw -> Lambda_w,       [N] -> w cross N              (9)
```

is an integral bijection. This independently recovers THM-4386's raw address.

The three error intervals for a fixed carrier have common length

```text
L_w(C)=[min(
  2r/w_1, 2r/w_2, 2r/w_3,
  r/w_1+r/w_2-|C_3|/(w_1w_2),
  r/w_1+r/w_3-|C_2|/(w_1w_3),
  r/w_2+r/w_3-|C_1|/(w_2w_3))]_+.                       (10)
```

The owner differences are the three carrier coordinates times invertible
speed factors.  Thus distinct owners are equivalent to every `C_i` being
nonzero modulo three, proving (1).

Now let `H=w^perp`, with Euclidean area, and push forward the cube
`[-r,r]^3` by

```text
T(e)=-w cross e.
```

Define its planar box-spline density `B_w` by

```text
integral_H phi(C) B_w(C) dC
 =integral_{[-r,r]^3} phi(T(e)) de.                      (11)
```

The two nonzero singular values of `T` are both `||w||`.  A fibre segment
over `C` has arclength `||w|| L_w(C)`.  Coarea therefore gives the exact
normalization

```text
B_w(C)=L_w(C)/||w||.                                    (12)
```

This is the centered three-direction box spline with directions
`-w cross e_i`.  Its planar Fourier transform is

```text
Bhat_w(xi)
 =product_i H_r((w cross xi)_i),

H_r(s)=integral_{-r}^r exp(-2 pi i s u)du.              (13)
```

## 4. Poisson summation and the exact `27/9` factor

The lattice `Lambda_w` has covolume `||w||` in `H`; its dual is

```text
Lambda_w^*=projection_H(Z^3).
```

The quarter-turn-and-scale map

```text
J_w:Lambda_w^* -> Lambda_w,       eta -> w cross eta    (14)
```

is a bijection.  Indeed, if `eta=projection_H(z)`, then
`J_w eta=w cross z`; the Smith calculation in (9) says these fill
`Lambda_w`, and `J_w` is injective on `H`.

Modulo three, the raw gate in (1) has exactly two classes.  Since each
`w_i` is a unit and its own inverse in `F_3`, they are

```text
+/-u,       u=(w_1^(-1),w_2^(-1),w_3^(-1))=(w_1,w_2,w_3) mod 3.
```

Choose a lift `c_0 in Lambda_w` of `u`.  Such a lift exists because
`w.u` is divisible by three and `w` is primitive.  Then the raw sum is

```text
||w|| [sum_{C in c_0+3Lambda_w} B_w(C)
      +sum_{C in-c_0+3Lambda_w} B_w(C)].                (15)
```

The covolume of `3Lambda_w` is `9||w||`.  Poisson summation in `H`, followed
by `xi=eta/3` and (14), turns (15) into

```text
(1/9) sum_{eta in Lambda_w^*}
 [exp(2 pi i eta.c_0/3)+exp(-2 pi i eta.c_0/3)]
 Bhat_w(eta/3).                                         (16)
```

If `n=J_w eta`, then `w cross (eta/3)=n/3`.  Since `r=3/14`, including at
zero,

```text
H_r(n_i/3)=3 hhat(n_i).
```

Consequently

```text
Bhat_w(eta/3)=27 product_i hhat(n_i).                   (17)
```

The dimension-two sampling covolume contributes `1/9`, while the three
one-dimensional transforms contribute `27`; the remaining factor is exactly
three.

It remains to identify the finite character.  Write
`eta=projection_H(z)`, set `x_i=w_i z_i mod 3`, and put
`W=w_1w_2w_3`.  A direct cross-product calculation gives

```text
(w_1n_1,w_2n_2,w_3n_3)
 =W(x_3-x_2, x_1-x_3, x_2-x_1) mod 3,                  (18)

eta.c_0 = x_1+x_2+x_3 mod 3.                            (19)
```

The three entries in (18) are equal exactly when (19) is zero; otherwise
they are all distinct.  The transform of the two-point mask `{+u,-u}` is

```text
2                    when eta.c_0=0 mod 3,
-1                   otherwise.
```

Multiplication by the residual `27/9=3` gives precisely `A_w=6,-3`.
Substitution in (16) proves the equality of (1) and (2), including every
constant.

For analytic hygiene, one may apply Poisson first to product-Fejer
regularizations.  THM-1092's absolute convergence then permits the limit and
gives the same raw series.  Equivalently, the compact box spline is continuous
and the sampled dual coefficients are absolutely summable after splitting by
coordinate support.

## 5. Why the natural same-vector identification does not match deletions

Two exact carriers in the obstruction family make the failure visible.

First,

```text
n=(-2,3,-1)
```

lies in both resonance lattices in (4), is nonzero and seven-free in every
coordinate, and has sheet weight `A_w(n)=A_w'(n)=-3`.  Thus it contributes on
the Fourier side.  But its middle coordinate is divisible by three, so the
same integer vector is deleted from the raw carrier sum.

Conversely, for `w'=(101,103,107)`,

```text
C=(-35,-1,34),       L_w'(C)=10/11021>0
```

is a live raw carrier: all its coordinates are nonzero modulo three.  If the
same vector were read as a Fourier frequency, its first sinc factor would
vanish because `7|-35`.

So Poisson duality preserves the **total weighted sum**, while the natural
identity on the shared embedded lattice preserves neither deletion predicate
term by term. It also preserves neither support order nor positivity. This is
the exact loss ledger required by MISTAKE-226/MISTAKE-235; no stronger claim
about arbitrary set-theoretic bijections is needed.

## 6. Sharp shortest-vector obstruction

Both triples in (4) satisfy the same primitive relation

```text
-2w_1+3w_2-w_3=0.                                      (20)
```

The verifier exhausts every possible vector below the discovered Euclidean
and multiplicative minima.  For both lattices it proves

```text
shortest squared Euclidean norm =14,
all minimizers                 =+/-(-2,3,-1),
m_7                            =6,
all m_7 minimizers             =+/-(-2,3,-1).           (21)
```

Even the shortest Fourier orbit has the same coefficient product and the
same sheet factor because the speed residue profiles agree.  Yet the first
raw sum consists of only

```text
+/-(-5,-1,4),       length 10/1547 each,
```

while the second has ten nonzero raw carriers and the different total in
(5).  The populations of seven-free full-support frequencies with
`product |n_i|<=M` are:

| `M` | first lattice | second lattice |
|---:|---:|---:|
| 6 | 2 | 2 |
| 12 | 2 | 2 |
| 24 | 4 | 2 |
| 48 | 6 | 4 |
| 96 | 12 | 4 |

These counts explain the first place where the two dual sums can diverge,
but they are not asserted to be a sufficient statistic.  The proved
obstruction is stronger and simpler: the complete minimum orbit does not
determine the answer.

This sharpens the status split across the older theorems:

- **THM-1075 + THM-1092 (PROVED):** the resonance lattice and its absolutely
  convergent sinc sum are exact;
- **THM-1080 (SCOPED/VERIFIED):** `m_7` is a meaningful leading invariant and
  detects Schur triples at value one;
- **THM-1085 (PROVED convergence + measured variability):** a bound based on
  `m_7` alone loses near-minimal counts and cancellation;
- **THM-4392 (PROVED exact local hostile):** even identical full minimizing
  vector data and character sign can give different local comb measures.

## 7. What the bridge makes available—and what remains open

The exact local problem can now be attacked in either of two equivalent
languages:

```text
primal:  discrepancy of a compact anisotropic box spline on two cosets of
         3Lambda_w;

dual:    a signed sinc series on Lambda_w with the index-three character
         A_w and nonzero mod-seven coordinate zeros.
```

A useful future bound must retain more than a first minimum.  Natural
sidecars are the successive-minima/near-minimal population, coordinate-axis
incidence, the index-three character, zero-coordinate support strata, and
the signs of the sinc product.  A hybrid estimate that Poisson-splits low
dual frequencies and controls the remaining primal box-spline discrepancy is
a **CONJECTURAL PROGRAM**, not a result here.

Nothing in the identity forces an arbitrary body-safe component to meet the
complement of `F_w`, synchronizes multiple tail triples, or supplies the
universal bound `mu(F_w)<=6/77`.  It is a local exact re-expression, not a
seam-entry or all-tail theorem.  **LRC(14) remains OPEN.**

## 8. Exact and diagnostic verification

The primary verifier imports no mathematical implementation from the
repository and performs 1,119 explicit optimization-safe checks:

- all 216 combinations of unit speed residues and dual residue lifts verify
  the finite Fourier transform `3*(2 or -1)=6 or -3`;
- the finite character counts are `72` cases of weight `6` and `144` of
  weight `-3` (three representatives per dual character);
- five raw carrier sums agree exactly with a separate implementation of the
  original six shifted danger-comb intersections;
- the complete Euclidean and `m_7` minima in (21) are certified by exhaustive
  lower boxes;
- the two opposed termwise-deletion witnesses are checked exactly;
- product-Fejer sums at cutoffs `24,48,96,192` approach the exact rational
  totals for three controls.  These last floating-point rows are diagnostics,
  not proof of the identity.

An independent clean-room verifier rederives the finite transform,
normalizations, physical controls, shortest-orbit hostile, and the precise
natural-identity obstruction in 1,492 further checks.  It also caught and
removed a stronger unproved reading: arbitrary bijections or blockwise
correspondences are not excluded.

Reproduction:

```powershell
python -B 04-computation/lrc14_raw_carrier_boxspline_fourier_duality_thm4392.py
python -O -B 04-computation/lrc14_raw_carrier_boxspline_fourier_duality_thm4392.py
$env:PYTHONHASHSEED='4386'
python -B 04-computation/lrc14_raw_carrier_boxspline_fourier_duality_thm4392.py
python -B 04-computation/lrc14_raw_carrier_boxspline_fourier_duality_independent_referee_thm4392.py
python -O -B 04-computation/lrc14_raw_carrier_boxspline_fourier_duality_independent_referee_thm4392.py
```

Normal, optimized, and seeded outputs are byte-identical.

```text
primary script       eac06db776d31ec09986d99b17144c3b3d17d46ebf43bdf70ed8edbe450fd427
primary output       9fc7a4202b74479682321c1737c79a7deacfb85dea920744e6f61cfbd14ee799
independent script   0e473ce08e89d7a11f30de8abd431edffcb9e0d89b0174fdb476a9df8d447f36
independent output   6b15cc9916f080504096e07ffca67b97a7311b594f1f68d751a2273756c4752c
```

Normal, optimized, and fixed-hash-seed replays match the two frozen canonical
outputs byte for byte.
