# Exact closure of the width-three degree-six mixed Catalan cell

**Status: VERIFIED-EXACT + FINITE-EXACT research artifact, 2026-08-18.**
This note does not reserve a theorem identifier and does not edit THM-3557.
It closes the sole open `N=3,D=6` affine coefficient branch left by
[THM-3557](../01-canon/theorems/THM-3557-low-width-mixed-catalan-thickening-no-go.md).
Combined with that theorem, there is no characteristic-zero solution of

```text
P=v^2+a(v)w+c(v)w^2+e(v)w^3,
Q=v^3-v+b(v)w+d(v)w^2+f(v)w^3,
Jac(P,Q)=1
```

when every coefficient polynomial has `v`-degree at most six.  This remains
a sharply restricted ansatz.  It does not prove `JC(2)`, does not address
width at least four or coefficient degree at least seven, and makes no claim
about a projective coefficient-space closure.

Companion:
[`jc2_catalan_mixed_thickening_degree6_exact.py`](../04-computation/jc2_catalan_mixed_thickening_degree6_exact.py),
with stored output
[`jc2_catalan_mixed_thickening_degree6_exact.out`](../05-knowledge/results/jc2_catalan_mixed_thickening_degree6_exact.out).

```text
source_sha256 = a3f94ab9eb4157bec7effc15ac1f1a8ab0842814c3a4f4ae4f6705a97fe4346f
output_sha256 = 02d23d80e4c50edc1323105a9d6fcc2cf5c77d729731a4be354c415251c31616
hash_basis    = LF-normalized bytes
```

## 1. Inheritance pass and concept board

- **Closest proved mechanism:** THM-3557 derives the exact `w`-graded
  recurrence, proves width at most two impossible without a degree cap, and
  closes width three through coefficient degree five.
- **Canonical hostile:** the three-row Catalan truncation really satisfies
  `(E0,E1,E2)=(1,0,0)` and first leaks at `E3=-135/16`.  Any obstruction that
  rejects this prefix is mistyped.
- **Corrected near miss:** the cap-five ledger cannot simply be rerun at cap
  six.  The new common-power degrees `(4,6)` must be admitted explicitly.
- **Least-used sidecar:** every degree statement needs the affine leading-
  coefficient saturation.  Dropping it folds lower-degree boundary cells into
  the alleged degree-six component.

The working board had five live objects:

1. the six exact recurrence rows `E0,...,E5`;
2. the four `e/f` zero-pattern branches;
3. the UFD common-power law from the top weighted Wronskian;
4. the nonzero-leading-coefficient saturation;
5. the first lower-row leak after each common-power substitution.

The relevant method card was **Refine and saturate before transporting a
factor or shadow** from `META-PATTERNS.md`.  Here its concrete effect was to
retain `t*C*F*(C+t^2) != 0` until the leading E2/E3 equations reduced it
honestly to `t != 0`.  No new meta-pattern is proposed from one thread.

## 2. The complete degree-six ledger

Write

```text
A=deg a, B=deg b, Cc=deg c, Dd=deg d, E=deg e, Ff=deg f.
```

The bottom row still forces `B=A+1`, hence `0<=A<=5`.  The top row splits
into the four exhaustive zero patterns used by THM-3557.

### 2.1 Both `e` and `f` nonzero

Put `f=lambda e`, `t=b-lambda a`, and `s=d-lambda c`.  If `s` is nonzero,
the last transformed row gives

```text
2E=3 deg(s),
```

so cap six allows exactly

```text
(deg(s),E)=(0,0),(2,3),(4,6).
```

The constant pair has the old uniquely-highest-degree obstruction.  For
`(2,3)`, the necessary states after transformed E2 are

```text
(deg(t),deg(c)) =
(1,5),(3,3),(4,2),(5,-infinity),(5,0),(5,1),(6,1),
```

where `-infinity` denotes `c=0`; transformed E3 closes all of them.  For the
new pair `(4,6)`, transformed E2 leaves

```text
(4,5),(5,4),(6,-infinity),(6,0),(6,1),(6,2),(6,3),
```

and transformed E3 again has a unique exact top term in every state.

If `s=0`, then `e=mu*t^3`.  Cap six permits `deg(t)=1,2`; the first dies in
transformed E1.  The second forces `c` constant there, after which
`3eH` has degree eight in transformed E2 and its only competitor has degree
one.  Thus the whole `e,f!=0` branch is empty.

### 2.2 `e=0,f!=0`

If `c=0`, the exact one-sided relations give

```text
(A,Dd,Ff)=(1,1,0),(2,3,3),(3,5,6),
```

and the E3 leading coefficient is always `3`, so all three states close.
If `c,f` are nonzero, E4 gives `3Cc=2Ff`.  The constant pair and `(2,3)`
pair are the old exact branches; the replayed THM-3557 companion closes them
with the stored unit ideals and terminal `9/4` and `-9/4` contradictions.

The new pair `(Cc,Ff)=(4,6)` has the complete filtration

```text
after E1:
 (A,Dd)=(0,5),(1,5),(2,5),(3,-infinity),(3,0),...,(3,5),
after E2:
 (A,Dd)=(3,-infinity),(3,0),...,(3,5),
after E3:
 (A,Dd)=(3,5).
```

Therefore the sole degree-ledger survivor is exactly

```text
e=0,        deg(a,b,c,d,f)=(3,4,4,5,6).               (1)
```

### 2.3 `f=0,e!=0`

For `d=0`, the only states are

```text
(A,Cc,E)=(2,2,2),(3,4,5),
```

and E3 has leading coefficient `-7` in both.  For nonzero `d,e`, the
common-power pairs `(Dd,E)=(2,3),(4,6)` have no simultaneous E1/E2 degree
survivor.  The constant pair is the unchanged old `A=3` branch and closes by
the replayed degree-seven obstruction.

### 2.4 `e=f=0`

This is the cap-free width-two obstruction in THM-3557.  Hence `(1)` is not a
selected search target; it is the unique output of the exhaustive ledger.

## 3. Exact common-power and E0/E1 parametrization

In branch `(1)`, E4 is

```text
3c'f-2cf'=0.                                           (2)
```

Unique factorization gives

```text
c=C h^2,                  f=F h^3,                    (3)
```

with `C,F` nonzero constants and `h` quadratic.  Absorb the leading
coefficient of `h` into `C,F` and write

```text
h=v^2+u v+r.
```

Every degree-`(3,4)` solution of E0 has the exact form

```text
g=t v^2+A v+B,
a=1+2v g,
b=3v/2+(3v^2-1)g,                  t!=0.              (4)
```

Indeed E0 is then identically one; conversely E0 evaluated at `v=0` forces
`a(0)=1`, so `g=(a-1)/(2v)` and `(4)` follows.

Equation E1 solves `d`:

```text
d = [2c(3v^2-1)-(a'b-ab')]/(4v).                      (5)
```

Polynomiality of `(5)` is exactly the vanishing of its remainder

```text
2A-4B^2+4Cr^2-3=0,                                   (6)
```

and its degree-five leading coefficient is

```text
[v^5]d = (3/2)(C+t^2).                                (7)
```

Thus the honest prescribed-degree chart is saturated by

```text
Delta=t*C*F*(C+t^2) != 0.                             (8)
```

This records every nonzero leading coefficient: those of `a,b`, `c`, `f`,
and `d`, respectively.  The monic normalization has already retained the
quadratic leading coefficient of `h`.

## 4. The transparent E2/E3 contradiction

The highest coefficients of E2 and E3 are three times

```text
L2=-3Ct+2F+t^3,
L3= 2Ft-C^2-Ct^2.                                     (9)
```

They satisfy the exact eliminant

```text
t L2-L3=(C-t^2)^2.                                    (10)
```

Hence `C=t^2`, and then `L2=0` gives `F=t^3`.  Under `(9)--(10)`, the
saturation `(8)` becomes `2t^8!=0`, so no leading cell was lost.

Set `k=t h=t v^2+U v+R`.  Then `c=k^2,f=k^3`.  Equation `(6)` becomes

```text
A=2B^2+3/2-2R^2.                                      (11)
```

After substituting `(11)`, the next two coefficients are

```text
[v^5]E2 = 9t(A-U)^2,
[v^6]E3 = 3t^2(A-U)^2.                                (12)
```

Since `t!=0`, they force `U=A`.  The next common coefficient factors as

```text
(B-R)(3(B-R)+2t)=0.                                   (13)
```

There are two branches.

1. If `B=R`, then `(11)` gives `A=U=3/2`, while

   ```text
   [v^2]E2=-6t^2,
   ```

   contradicting `t!=0`.

2. If `R=B+2t/3`, put `X=48Bt+16t^2`.  Two consecutive equations are

   ```text
   [v^2]E2=(2t^2/9)(X-81),
   [v^3]E3=(4t^3/27)(2X-135).
   ```

   They demand `X=81` and `2X=135`, but `162!=135` in characteristic zero.

This is the terminal contradiction.  It explains why the degree-six cell
fails: the top equations align the two square/cube bases almost completely,
and the only two possible constant offsets produce incompatible lower leaks.

## 5. Gröbner/Nullstellensatz audit

The script independently forms every coefficient of E2 and E3 after the
exact substitutions `(10)--(12)`.  It adjoins `zeta*t-1` and computes the two
branch ideals obtained by adding the factors in `(13)`:

```text
I_A = (coefficients(E2,E3), B-R,             zeta*t-1),
I_B = (coefficients(E2,E3), 3(B-R)+2t,       zeta*t-1).
```

After the displayed substitutions these are ideals in
`QQ[zeta,B,t]`.  Both Gröbner bases are exactly `[1]`.  Therefore both
saturated affine branch varieties are empty over an algebraic closure of
`QQ`; the unit certificates persist over every characteristic-zero field.
The same calculations give `[1]` over `GF(5)` and `GF(7)` as controls, not as
a modular-to-characteristic-zero inference.

This branchwise ideal is the smallest honest algebraic certificate.  A raw
untriangularized coefficient Gröbner basis is not used and no runtime claim
about it is part of the result.

## 6. Controls, truth gates, and reproduction

The focused companion has `64` explicit runtime gates and contains no Python
`assert`, so optimization cannot erase a check.  Its controls are:

- the LF-normalized hashes of the THM-3557 script/output;
- ordinary and optimized replays of all inherited `86` gates, each matching
  the stored output byte-for-byte after line-ending normalization;
- the identity-map recurrence control `(E0,E1)=(1,0)`;
- the genuine Catalan prefix
  `(1,0,0,-135/16,-405/64,-729/128)`;
- a degree-`(3,4,4,5,6)` hostile tuple with `t=1`, `B=R=0`,
  `A=U=3/2`, which satisfies `E0,E1,E4,E5=(1,0,0,0)` but exposes exactly
  `[v^2]E2=-6` and `[v^3]E3=-4`;
- the `t=0` degree-drop hostile, confirming why saturation is mandatory.

Reproduce from the repository root:

```text
python 04-computation/jc2_catalan_mixed_thickening_degree6_exact.py
python -O 04-computation/jc2_catalan_mixed_thickening_degree6_exact.py
```

Both transcripts are LF-byte-identical and equal the stored output.

## 7. Connections and honest remaining frontier

The exact connection contract is:

```text
source: E4 weighted Wronskian on (c,f)
target: common quadratic base h and its saturated coefficient chart
map: (c,f) -> (C,F,h) with c=C h^2, f=F h^3
preserved: E4=0 and the exact degrees (4,6)
destroyed: every lower Jacobian row
sidecar: E0/E1 parametrization plus Delta!=0
decisive test: the first nonzero E2/E3 coefficient after each factor split.
```

Direct progress is the exact emptiness of `N=3,D=6`.  The reusable niche is
the branch compiler pattern `weighted Wronskian -> common powers -> leading
saturation -> first two leak rows -> branchwise unit ideals`.  The cubic-root
cover interpretation remains a wildcard lens: this calculation closes the
first even `r^6` cell but does not prove a general parity or ramification-line
obstruction.

The honest next cells are width three with coefficient degree at least seven
and width at least four.  No survivor from the degree-six ledger remains.
