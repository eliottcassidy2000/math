# The U_full H guard is an exact sheet–chamber skew product

**Status: PROVED ELEMENTARY FACTORIZATION + VERIFIED-EXACT SIDECAR; NOT
INDEPENDENTLY AUDITED.**  This constructs one actual, character-independent
address coordinate for the H-guard factor of the refined U_full bank.  It
does not construct the common-ancestry support relation, the full endpoint
coupling, a physical current, a row exclusion, or LRC(14).

## 1. Inheritance pass and the narrowed question

Promoted THM-3479 refines the U_full endpoint bank by characters

```text
ell=(tau,-alpha,-beta,0,0,0,0,alpha,beta) in F_13^9.   (1)
```

The frozen role classes

```text
q_H =(1,0,1),
q_q5=(1,0,0)                                           (2)
```

differ only in the Fourier coordinate dual to `tau`.  The minimum joint-
address gate proves that one character-independent mixed coordinate is
necessary and sufficient in the abstract `2x2` model, but it does not map
actual U_full cells to that coordinate.

The closest positive mechanism is THM-3505's owner/root shear: a translated
root coordinate becomes a projective slope pencil once the two labels are
retained before Fourier transformation.  The canonical hostile is the
checkerboard marginal kernel.  The corrected near miss is that equality of
circle points or components does not define ancestry.  The least-used actual
field is the H-guard translation already visible in `fast_build_set`.

For U_full the guard speed is one.  The `tau` dependence of the guard-safe
factor is exactly deletion of

```text
[-13-7tau,13-7tau) modulo 91.                           (3)
```

All other factors are independent of `tau` when `alpha,beta` are fixed.
Thus (3), unlike a synthetic parity label, is an actual pre-merge endpoint
coordinate.

## 2. Thirty-nine common guard atoms

Write a normalized guard-circle point uniquely as

```text
s=7a+r modulo 91,
a in F_13,        0<=r<7.                              (4)
```

The union of the boundaries of all thirteen translates in (3) splits each
sheet into the three half-open chambers

```text
C_L=[0,1),       C_M=[1,6),       C_R=[6,7).           (5)
```

The half-open convention is load-bearing for endpoint sums: `r=1` belongs to
the middle chamber and `r=6` to the right chamber.  Direct reduction of (3)
gives the forbidden sheet arcs

```text
D_L={12,0,1},
D_M={11,12,0,1},
D_R={11,12,0}.                                         (6)
```

If `g_C(u)=1_(u notin D_C)`, then the guard-safe indicator at the point
`(a,r)` is

```text
Guard_tau(a,r)=g_C(a+tau),       r in C.               (7)
```

Equations (4)--(7) give one common, `ell`-independent partition

```text
F_13 x {C_L,C_M,C_R},                                  (8)
```

with `39` atoms.  The tau twist translates only the sheet coordinate.  This
is the first requested address field obtained from the actual U_full factor
geometry rather than calibrated after seeing the bridge.

## 3. Every single-sheet Fourier mode survives

Let `xi` be a primitive thirteenth root and use transform convention

```text
G_C(k)=sum_(u in F_13) g_C(u) xi^(-ku).                (9)
```

Changing variables `u=a+tau` in (7) gives the exact covariance

```text
sum_tau Guard_tau(a,r)xi^(-k tau)
  =xi^(ka)G_C(k).                                      (10)
```

At `k=0`, the three values are the safe counts

```text
(10,9,10).                                             (11)
```

For `k!=0`, the full-circle sum vanishes and `G_C(k)` is the negative
Fourier sum of one consecutive cyclic arc of length three or four.  Hence it
is, up to a root-of-unity factor,

```text
-(1-xi^(h k))/(1-xi^k),        h in {3,4}.             (12)
```

Since thirteen is prime and `0<h<13`, neither numerator nor denominator in
(12) vanishes.  Thus all `39` chamber-frequency kernels are nonzero over
`Q(xi)`.  The exact split-field audit checks all `507=3*13*13` addressed
instances in `F_547`.

## 4. A left/right pair factors through relative drift

For a left point `(a,r)` and right point `(b,r')`, put

```text
d=b-a in F_13.                                         (13)
```

The common tau twist of the two guard factors obeys

```text
sum_tau g_C(a+tau)g_D(b+tau)xi^(-k tau)
 =xi^(ka) K_(C,D,d)(k),                                (14)

K_(C,D,d)(k)
 =sum_u g_C(u)g_D(u+d)xi^(-ku).                        (15)
```

So the pair needs only

```text
common sheet a
+ relative drift d
+ left/right residual chambers.                       (16)
```

The drift is genuinely joint: it cannot be computed from the two sheet
marginals.  Exact enumeration in `F_547` gives

```text
all pair kernels:       1521/1521 nonzero,
primitive k!=0 kernels: 1404/1404 nonzero,
k=1 kernels:             117/117 nonzero.              (17)
```

Thus no chamber or drift is spectrally silent at the frequency distinguishing
(2).  This is a capacity theorem for the guard factor, not a claim that every
bucket has nonzero physical endpoint weight.

## 5. The drift hostile and why marginal order matters

Take both residual chambers equal to `C_M` and couple uniform sheet marginals
by `b=a+d`.  Every value of `d` has exactly the same uniform left and right
marginal.  Nevertheless the number of jointly safe sheets is

```text
(9,8,7,6,5,5,5,5,5,5,6,7,8)                          (18)
```

as `d=0,...,12`.  Equation (18) is a thirteen-state version of the minimum
checkerboard hostile: equal marginals, different joint response.

There is a complementary cancellation.  For every `k!=0`, summing the common
sheet `a` before applying its character gives

```text
sum_a xi^(ka)K_(C,D,d)(k)=0.                           (19)
```

Thus the order of operations is forced:

```text
retain common sheet and drift
 -> twist
 -> only then sum.                                     (20)
```

This is exactly the order missing from a scalar-only endpoint API.

## 6. What has moved, and what remains absent

The previous minimum record requested

```text
(base,root,owner_sheet,word_sheet,source_sheet,
 left_horizon,right_horizon,address).                  (21)
```

The present calculation constructs one honest piece of `address`:

```text
guard_sheet=a,       guard_chamber=C,
pair_guard_drift=d.                                    (22)
```

It also explains why the boundary-provenance field in the actual pre-merge
extractor must be retained: it decides the half-open chamber in (5).

It does **not** identify (22) with THM-2471's root, owner, word, or source
sheets.  It does not define which left and right atoms share a physical
ancestry base.  The independent Cartesian product of endpoint atoms carries
(22), but that product has no LRC support relation and therefore cannot be
declared a current.

The next exact engine can now be much smaller than an arbitrary pair census.
Bucket each endpoint contribution by the `39` labels `(a,C)`, then combine
the two bucket tables through the `117` joint types `(C,D,d)`.  The decisive
test is whether the frozen `q_H-q_q5` bridge decomposes into a nonzero typed
drift bucket after factor lineage and boundary provenance are retained.  A
successful decomposition would still need a theorem selecting the lawful
ancestry support inside those buckets.

## 7. Connection contract

| field | exact content |
|---|---|
| source | actual H-guard factor `[-13-7tau,13-7tau)` in the U_full endpoint engine |
| target | a common `F_13` sheet, three residual chambers, and pair drift |
| map | `s=7a+r`, then `d=a_R-a_L` |
| preserved | all tau characters, half-open boundary convention, joint guard overlap |
| destroyed by sheet marginal | every primitive tau mode |
| needed sidecar | factor lineage, boundary side, and a lawful common-ancestry support relation |
| hostile | uniform diagonal/shift couplings have identical margins but profile (18) |
| cheapest next test | 39-bin endpoint aggregation and 117 drift-bucket bridge decomposition |

No D5 coefficient map, K4 factor, grouped relation current, all-unit
projector, physical chronology, row exclusion, or LRC(14) consequence is
claimed.

## 8. Reproduction

```text
python3 04-computation/lrc_ufull_guard_sheet_joint_drift_address_probe_20260816.py
python3 -O 04-computation/lrc_ufull_guard_sheet_joint_drift_address_probe_20260816.py
```

The companion pins the THM-3479 endpoint-engine hash and guard convention,
checks the cyclic danger sets and endpoint boundaries, verifies (10) and
(14) at every address, chamber, drift, and frequency, and freezes both
hostiles.  Its semantic ledger is

```text
f345d40c9b589910d83d2fd490ca9376b3cbb86d7aa6da4825677d06075bda7a.
```

The LF SHA-256 hashes of the frozen script and output are, respectively,

```text
79b78637b2cc0ff54051fde02a6651ef10c8694a8d7a865ae403696370125179
39ced0ea361bb0268e607d42d030a61a22e0e2a571e66e585c196c352ac88f3c
```
