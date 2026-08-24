# THM-4009 Fourier/LP bottleneck audit — 2026-08-24

**Status: FINITE-EXACT bounded laboratory + PROVED information-loss
diagnosis.  LRC(14) remains OPEN.**  This note reserves no theorem ID and
changes no navigation surface.

Companion computation:

```text
04-computation/lrc14_graver_bouquet_fejer_owner_audit_20260824.py
05-knowledge/results/lrc14_graver_bouquet_fejer_owner_audit_20260824.out
```

Normal and optimized runs agree byte-for-byte.  All signs called rigorous
below are certified at 256-bit precision by Arb interval arithmetic.  The
polynomial/minimal-root sign for the legacy order-119 control also has a
separate exact algebraic isolation and the rational bound `pi<355/113`.

## 1. Inheritance and the live concept board

The closest proved mechanism is
[THM-4009, Euclidean covering-transference short-relation compression](../../01-canon/theorems/THM-4009-euclidean-covering-transference-short-relation-compression.md):
a hypothetical counterexample has a Graver row `a` with

```text
sum a_i^2 <= 195,       ||a||_1 <= 50,       |a_i| <= 13.
```

The canonical equality hostiles are the arithmetic progression
`(1,...,13)` and the Goddyn--Wong row `(1,...,11,13,24)`.  Both are weakly
lonely at `t=1/14`, with opposite-slope endpoint owners `(1,+)` and `(13,-)`,
but their open danger combs cover almost everywhere.

The corrected near misses are the max/mean and strict-boundary failures in
`MISTAKE-127`, `MISTAKE-129`, and `MISTAKE-130`, together with the typed-kernel
failure in `MISTAKE-235`.  Their common warning is load bearing here:

- a continuum average is not a maximum or an admissible grid witness;
- a uniform grid cannot see rational strict-boundary pinches;
- two sums indexed by integer kernels are not equivalent unless the map also
  preserves their weights, phases, and target predicate.

The least-used relevant sidecar is the **prefix ordering of a relation**, not
just its coefficient multiset or autocorrelation.  The live board was:

```text
short Graver row | prefix bouquet | Toeplitz inertia | endpoint owner
BONF5 support gate | Vaaler relation-free gate | AP/GW equality hostile
```

## 2. Exact short-bouquet experiment

Use the bounded strict-control family

```text
V_w={1,...,11,13,w},       12<=w<=200,       w!=13.
```

For `V_38`, distinctness rules out relation norm squares one and two, while
`1+2-3=0` gives norm square three.  Hence every Euclidean-shortest relation
has three coefficients of absolute value one and is exactly an additive
triple

```text
x+y=z.
```

Exact enumeration gives thirty such triples.  For a cap-respecting scalar
bank, `k(x+y-z)` has `l1=3k<=50` for `1<=k<=16`.

For one prefix ordering, use

```text
B_x=union_(1<=k<=16){0,kx,kz};
```

the other non-equivalent ordering uses `B_y`.  The full ordering-retaining
bouquet is

```text
B_xy=B_x union B_y
    =union_(1<=k<=16){0,kx,ky,kz}.
```

Let `N_V(t)` count open radius-`1/14` danger combs and let
`H_V=N_V-1`.  On a finite frequency set `B`, form the Hermitian moment matrix

```text
M_B(lambda,mu)=hat(H_V)(lambda-mu).
```

For every coefficient vector `c`,

```text
c* M_B c = integral H_V(t) |sum_(lambda in B)c_lambda e(lambda t)|^2 dt.
```

Thus a negative value proves a positive-measure strict lonely set.

### Exact outcome

**FINITE-EXACT.**  All sixty one-order matrices (two orderings for each of
thirty shortest relations) are positive definite.  An interval `LDL^T`
factorization certifies every pivot strictly above `2/5`; the smallest
certified pivot is approximately `0.4043275367`.

The full two-order matrices split sharply:

```text
28 positive definite (every LDL pivot >1/3),
 2 with an exact negative integer direction.
```

For `2+6=8`, the coefficient word frozen in the script has square norm `394`
and gives

```text
integral H_(V_38)|P|^2
 = -5.9113725705681529025... < 0,
normalized value = -0.0150034836816450581... .
```

The second negative triple is `2+4=6`.  Because both one-order principal
spaces are positive definite, the cross-order terms are logically
load-bearing.  This is a small but genuine positive result: a bounded
trigonometric certificate can exploit the sharpened relation atlas, but it
must retain ordering/phase information.

A centered Fejer kernel does not find this direction: every order from one
through seventeen at the exact center `1/12` is positive, with order seventeen
equal to `0.4352363002...`.  The legacy scalar bouquet under the old `356`
budget first turns the same centered Fejer test negative at order `119`.
Relation shortness therefore does not automatically lower the localization
degree.

The cap `50` is not a logical ceiling on harmonic multiples: once a relation
exists, all its scalar multiples exist.  The `k<=16` bank is instead the
finite, cap-respecting architecture audited here.  No conclusion relies on
silently treating `50` as an absolute Fourier cutoff.

## 3. What THM-4009 genuinely improves

**PROVED / inherited.**  In support two, THM-4009 leaves exactly forty-seven
coprime ratios `p<q` with `p^2+q^2<=195`, giving 3,666 unoriented
ratio/support packets and 7,332 oriented labelled coefficient assignments
modulo global sign; all have `p+q<=19`.  Across all supports it leaves 55,459
absolute coefficient histograms.  This is a dramatic finite-orbit reduction
and is suitable input for an exact SDP/LP traversal.

**PROVED support gate.**  The improvement does not directly close
[THM-2051, Fejer--BV BONF5](../../01-canon/theorems/THM-2051-fejer-bv-small-relation-alternative-for-lrc14.md):

- if the short row has support three through five, its height at most thirteen
  substantially sharpens that theorem's residual height `2^20`;
- if it has support at least six, it cannot occur in any Fourier term of a
  product involving at most five labelled coordinates, so degree-five BONF5
  does not see it directly;
- support two is already paid exactly in THM-2051's stronger triangle-packed
  form.

Nor does a relation automatically improve the relation-free hypothesis of
[THM-2085, the Vaaler/Selberg rank-seven gate](../../01-canon/theorems/THM-2085-explicit-height-57-rank-seven-selberg-gate.md).
Vaaler's tensor calculation is exact when low-height resonances are absent;
a newly forced resonance creates a finite correction term whose sign and
coordinate ownership still have to be controlled.

This agrees with the primary literature pins in
[CORE-PAPERS](../reference/CORE-PAPERS.md): Vaaler supplies interval
sandwiches, Tao supplies a general finite checking reduction, and Bedert's
Riesz products improve asymptotic maximum loneliness.  None supplies a
positive-measure lower bound or converts a bulk spectral statement into a
tight endpoint witness.

## 4. Exact surviving obstructions

### Equality remains invisible

**PROVED + FINITE-EXACT.**  For AP and Goddyn--Wong, every one of the thirty
full two-order matrices above is positive definite; every audited LDL pivot
is greater than `2/5`.  More generally, `H_V>=0` almost everywhere on either
row, so every absolutely continuous Toeplitz quadratic is nonnegative at
every degree.  Their valid lonely time is measure zero and requires the
signed endpoint-owner sidecar.

This is the same boundary isolated by
[THM-2121, effective Fejer localization](../../01-canon/theorems/THM-2121-effective-fejer-localization-for-binary-shifted-combs.md)
and
[THM-2143, Gibbs zero-count soundness](../../01-canon/theorems/THM-2143-artanh-riesz-gibbs-certificate-soundness-boundary.md):
finite Toeplitz/Gibbs data are complete for positive-measure strict loneliness,
not for an isolated tight point.

### Relation fibre and autocorrelation remain phase blind

**PROVED / inherited hostile.**  Section 5 of
[THM-4002, signed-endpoint cross-phase](../../01-canon/theorems/THM-4002-lrc14-signed-endpoint-cross-phase-and-fixed-scale-two-family.md)
exhibits

```text
(1,...,11,33,154),       (1,...,11,77,110).
```

They have the identical complete minimum-`l1` relation fibre: norm three,
thirty relations, all inside the AP11 body, but their relevant covariances
have opposite signs.  Since norm three is already far below the new cap, the
THM-4009 sharpening does not touch this hostile.  Full Fourier magnitude is
also phase blind under translation.  A successful two-component certificate
must retain a cross-correlation with the signed body endpoint word.

## 5. Connection contract and stopping boundary

```text
source:      THM-4009 short Graver relation
target:      finite Hermitian danger matrix
map:         ordered closed-walk prefixes -> frequency bouquet
preserved:   sign partition, harmonic scale, prefix order, cross-order phase
destroyed:   endpoint owner, physical component incidence, arrival, parity
sidecar:     signed endpoint word / owner-address cross-correlation
hostiles:    AP, Goddyn--Wong, THM-4002 identical-fibre opposite covariance
```

The cheapest decisive next test is a finite exact SDP atlas, not another
scalar Fejer scan:

1. intersect the forty-seven support-two THM-4009 ratios with the seventeen
   live THM-3910 two-component types;
2. for each compatible cell, take the union of **all** prefix orderings and
   retain the component/owner word;
3. use numerical SDP only for discovery, then freeze an integer coefficient
   vector and certify its sign with Arb, as above;
4. run AP, Goddyn--Wong, and both THM-4002 phase twins as mandatory hostiles.

A single pair of phase twins with identical proposed LP input but opposite
required conclusions decisively refutes that input type.  If the owner-word
augmented bank separates all hostiles, it yields a concrete finite certificate
compiler rather than another measure-only heuristic.

## 6. Reproduction and exact scope

```bash
python3 -B 04-computation/lrc14_graver_bouquet_fejer_owner_audit_20260824.py
python3 -B -O 04-computation/lrc14_graver_bouquet_fejer_owner_audit_20260824.py
```

Both streams byte-match the frozen output.  This audit claims **no global
Fourier closure**, no weak-boundary certificate, and no new LRC(14) scope.
Its positive result is the exact bounded `V_38` cross-order certificate; its
negative result is that shortness, one-order bouquets, scalar Fejer weights,
and phase-free relation fibres still lose indispensable information.
