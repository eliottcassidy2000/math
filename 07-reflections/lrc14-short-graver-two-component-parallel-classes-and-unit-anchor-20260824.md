# LRC(14): short Graver rows, two parallel classes, and the unit-anchor crossing lemma

**Session status (2026-08-24):** LRC(14) is **OPEN**.  This note proves a
narrow `t<U` crossing lemma, proves that graphic cut/cycle and ordinary
matroid augmentation cannot by themselves force a crossing row in the exact
rank-eleven two-component quotient, and audits the effect of the incoming
Euclidean short-relation theorem on THM-3910's 17 conditional certificate
types.  It closes none of the 17 types.  Statements using THM-4009 were
**CONDITIONAL on its concurrent promotion** when this note was written; the
unit-anchor and parallel-class statements do not depend on it.

## Inheritance pass and concept board

- **Closest proved mechanism:**
  [THM-3818](../01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md)
  identifies the bounded support-two row matroid as graphic.  In its exact
  `dim W=11` equality branch, the decoder graph has two components and every
  bounded crossing row lies outside `W=V_dec`.
- **Closest positive LRC sidecar:**
  [THM-3878](../01-canon/theorems/THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse.md)
  uses component gluing and harmonic absorption; THM-3910 leaves 16
  scale-one pair types plus scale-two `(1,9)` when `t>=U`.
- **Canonical hostile:** THM-3818's safe `2+11` row has decoder rank eleven,
  `W=V_dec`, no height-`Q` crossing row, and the internal pair relation
  `(3,-1)` of Euclidean square norm `10`.  It is the clean witness that a very
  short relation need not cross the component cut.
- **Corrected near miss:** MISTAKE-460 separates the arbitrary-common-scale
  finite decoder from the stronger table-free inert-scale decoder.  This
  session found and repaired the same stale filter in the two-component
  companion; see MISTAKE-486.
- **Least-used relevant sidecar:** THM-3995's owner-oriented wall numerators.
  Its scale-two parity escape is the endpoint analogue of the divisibility
  gate in the unit-anchor lemma below.

The live board was:

| object | faithful representation | operation | lost coordinate / test |
|---|---|---|---|
| rank-eleven code `W` | two projective parallel classes | matroid quotient | relative scale and crossing-edge existence |
| decoder graph | direct sum of two graphic matroids | cut/cycle augmentation | a bridge belongs to no cycle |
| short Graver row | Euclidean relation ball | intersect with `W` | internal versus crossing component charge |
| `t<U` scale quotient | three speeds `sd,sU,tp` | Dirichlet approximation | divisibility anchor `sd|tp` |
| owner/arrival | oriented endpoint word | compare wall images | algebraic relation has no semantic time |

## 1. Why every triple can saturate without a crossing row

Let the two decoder components be `I,J`, with primitive shapes `u,v`.  In the
exact equality branch,

```text
W=V_dec,
W^perp={alpha u direct-sum beta v:alpha,beta in Q}.
```

Choose the displayed basis of `W^perp`.  Its coordinate columns are

```text
i in I: (u_i,0),               j in J: (0,v_j).       (1)
```

Projectively these are exactly two parallel classes.  Every three columns
contain two from one class, so their required THM-2052 dependence can always
be the internal pair dependence on those two columns.  In the `11+2` case
all `C(13,3)=286` triples saturate in this way.  The pigeonhole theorem's
every-triple dependence therefore has an exact equality model with no
crossing row.

This also settles the ordinary matroid-augmentation hope.  A basis of the
rank-two column matroid is one label from each class; there are `11*2=22`
bases, and basis exchange is already complete inside the supplied ground
set.  Augmentation never manufactures a new admissible element.

On the row side, the decoder cycle matroid is

```text
M(G_I disjoint_union G_J)=M(G_I) direct_sum M(G_J).    (2)
```

If a crossing edge is independently supplied, it joins the components as a
bridge and is a coloop: rank rises from eleven to twelve, but the new edge is
in no cycle.  Thus cycle equations cannot determine its coefficient or force
its existence.  The component cut is invisible on the old edge ground set;
it becomes nonzero only after the bridge has already been supplied.

**PROVED stopping conclusion.** Pure cut space, cycle space, basis exchange,
and every-three-columns dependence do not force an admissible bounded
crossing row.  A successful augmentation must add a metric or semantic
predicate—coefficient height, a quotient shortest vector, an endpoint owner,
or an actual arrival event—not another unweighted graphic identity.

## 2. Exact effect of the Euclidean short-relation compression

The concurrently audited THM-4009 strengthens THM-3743, for a hypothetical
counterexample, to a Graver relation `a` satisfying

```text
||a||_2^2<=195,       ||a||_1<=50,       |a_i|<=13.   (3)
```

Let

```text
L_int=(ker_Z u direct-sum 0) direct_sum
      (0 direct-sum ker_Z v),
lambda_int^2=min{||z||_2^2:0!=z in L_int}.             (4)
```

Inside `W=V_dec`, (4) is the exact internal absorption threshold.

- If `lambda_int^2>195`, every row supplied by (3) crosses the two components,
  lies outside `W`, and raises the relation rank to twelve.  This enters the
  inherited finite terminal; it does not itself prove loneliness.
- If `lambda_int^2<=195`, the norm cap alone cannot force crossing.  The
  global shortest row may be internal or crossing, and a component/owner
  sidecar is still necessary.

For a two-vertex primitive component `(p,q)`, its internal lattice is
generated by `(q,-p)` and has minimum square norm `p^2+q^2`.  The exact
support-two Euclidean atlas contains 47 coprime `p<q`.  On THM-3910's 16
scale-one residual pairs, 14 have `p^2+q^2<=195`; the only tail outliers are

```text
(8,21): 505,                  (9,11): 202.             (5)
```

The scale-two `(1,9)` tail has norm square `82`.  Consequently the new norm
bound cannot force crossing for 14 scale-one types or the scale-two type.
Even (5) do not close: the eleven-body component may contain a shorter
internal row.  Exact safe controls with body `(1,3,...,3^10)` have 19 short
support-two rows, all body-internal, with minimum square norm `10`, for each
of the two tails in (5).

The half-lattice parity sidecar in THM-4009 suggests a narrower future probe,
but not yet a theorem.  The tail circuit has coefficient sum `q-p mod 2`.
It already absorbs an odd character for 12 of the 16 scale-one types.  Its
even-character residuals are

```text
(1,3), (1,9), (3,7), (9,11),                           (6)
```

and scale-two `(1,9)` is also even.  THM-4009 proves existence of some
odd-character relation but does not bound it by (3).  A character-sensitive
transference bound could therefore be useful only after also excluding an
odd internal body row; attaching parity silently to the shortest vector would
be invalid.

## 3. A proved `t<U` divisor-anchor crossing lemma

The remaining `t<U` branch admits one clean exact reduction that does not use
THM-4009.

### Lemma

Let `Q>=1`.  Suppose a speed row contains two same-component speeds

```text
A=s d,                   M=s U,             0<d<U,
```

and a speed `N=t p` in the other component.  Assume

```text
A divides M,             A divides N,
0<N<M,                   M/A<=Q^2.                    (7)
```

Then there is a genuine crossing relation supported on `{A,M,N}` with
coefficient height at most `Q`.

### Proof

Normalize

```text
m=M/A,                  n=N/A,
```

so `0<n<m<=Q^2`.  Among the `Q+1` fractional parts
`0,n/m,...,Qn/m`, two lie in the same one of `Q` equal subintervals.
Therefore there are integers

```text
1<=c<=Q,          0<=b<=c<=Q,          |cn-bm|<=m/Q<=Q. (8)
```

Put `a=cn-bm`.  Multiplication by `A` gives

```text
a A+b M-c N=0,                |a|,|b|,|c|<=Q.         (9)
```

The coefficient of `N` is nonzero and the body partial sum equals `cN>0`,
so (9) genuinely crosses the component cut.  This proves the lemma.

The divisibility `A|M` is essential: without it `m=M/A` need not be an
integer, so the coefficient `a=cn-bm` need not be integral.  In the original
coordinates, the two divisibilities in `(7)` are exactly `d|U` and `sd|tp`.

In the canonical specialization `Q=91^6`, the finite counterexample box is
`B=Q^2`, so the last inequality in (7) is automatic whenever the normalized
speeds occur in that box.  But THM-3818's exact rank-eleven equality branch
has no crossing support-at-most-three height-`Q` row.  Hence (7) is forbidden
there.

The most transparent consequence is:

```text
scale one + body coordinates 1,U + pt<U
    ==> a height-Q crossing row
    ==> not the W=V_dec equality branch.               (10)
```

Thus a scale-one equality row with `t<U` and `p=1` cannot have coordinate
`1` in its primitive eleven-body shape.  More generally, a divisor anchor
`d|U` is forbidden whenever `d|tp` and `tp<U`.

This is a real narrowing, not a closure.  Primitive component gcd one does
not imply that a coordinate equals one or divides the maximum, and for
`p>1` the interval `U/p<=t<U` remains outside (10).  In the scale-two
exception, `s=2`, `p=1`, and coprimality makes `t` odd, so the natural anchor
`A=2` does not divide `N=t`.  That exact divisibility failure is the algebraic
shadow of THM-3995's parity escape; it does not reproduce THM-3995's oriented
owner holes or prove arrival.

## 4. Hostile controls and the repaired all-scale replay

The new companion

[`lrc14_two_component_parallel_class_unit_anchor_audit_20260824.py`](../04-computation/lrc14_two_component_parallel_class_unit_anchor_audit_20260824.py)
checks:

- the 47 Euclidean support-two ratios and the exact split (5);
- all 286 triples and all 22 bases of the two-parallel-class quotient;
- the bridge/coloop rank jump `11 -> 12`;
- the unit-anchor lemma for every `2<=Q<=25`, every `2<=U<=Q^2`, and every
  `1<=N<U`: `1,074,060` exact cases, including `5,500` cases at `U=Q^2`;
- THM-3818's finite-box rank-eleven hostile, where the short norm-10 row is
  internal and the dominance inequality excludes every height-`Q` crossing
  support-at-most-three row; and
- the two safe residual-tail controls from (5).

Normal and optimized streams byte-match the frozen
[`output`](../05-knowledge/results/lrc14_two_component_parallel_class_unit_anchor_audit_20260824.out):

```text
script sha256 = 60a20b244e945fe8548631ce34398d12428ae8e585959d677c282b3ea0b5b703
output sha256 = 5f0bc5d628d0d0c84f9d110d4bc03e7e3c10946e13bc2c5d6c680ba40d239397
requirements  = 5,371,764
RESULT        = PASS; LRC14=OPEN
```

During this replay, the inherited two-component companion was found still to
use the table-free common-scale filter corrected elsewhere by MISTAKE-460.
On its canonical `2+11` hostile, the all-scale decoder restores the internal
pairs `(3,9)` and `(9,24)`, with reduced ratios `(1,3)` and `(3,8)`, changing
the edge count `25 -> 27`.  MISTAKE-486 records the repair.  The refrozen
primary companion has

```text
script sha256 = aae0fa24c0f1eb2eee04445a16a8a899aee52b663a8d7757d44c5284e3697697
output sha256 = e360e2e1a4805850132e25b559b14091e7b164349c74ee30c78fa405c275d742
semantic      = 0606c4dd0c129a3feef83ccb31d09428e5ed6e2a49741d15bc5db9b335b16e71
```

The two restored edges stay internal.  Component sizes, rank, dominance,
scale recovery, singleton fibre, and every THM-3818 implication are unchanged.

## 5. Connection contract and next decisive tests

The proved new connection is

```text
source:      a divisor-normalized t<U three-speed packet
target:      THM-2052's bounded relation code
map:         Dirichlet approximation n/m -> aA+bM-cN=0
preserved:   coefficient height and the actual component cut
destroyed:   endpoint owner, sign side, phase, first arrival, loneliness
sidecar:     a divisor anchor sd|tp and an owner-labelled endpoint word
test:        classify divisor anchors in the 16 scale-one component shapes.
```

The most useful next tasks are therefore:

1. Compute or bound the internal Euclidean minimum of the eleven-body
   component.  THM-4009 forces a crossing rank increment exactly after both
   internal components are shown to have minimum square norm at least `196`.
2. In `t<U`, enumerate divisor anchors `d|U` carried by the connected decoder
   shape and intersect them with `sd|tp`.  This is a projective-shape
   computation, not a scale-blind graph count.
3. Seek a character-sensitive transference bound.  Its first honest target is
   the even-character list (6), after excluding odd body-internal rows.
4. Apply owner incidence before scalarization: compare the actual signed body
   endpoint word with the pair-obstruction endpoint word.  A crossing relation
   without this cospan raises rank but gives no safe time.

The conclusion is sharper than “augmentation failed.”  The exact failure
mechanism is a two-parallel-class saturation whose missing edge would be a
coloop.  The first successful escape is not another cycle identity; it is a
metric crossing vector or a divisor/owner sidecar that lives transverse to
the direct sum.
