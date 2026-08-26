---
id: THM-4149
title: "Minimum-boundary completion of the first-window universal odd-tail LRC(14) transfer"
status: >
  PROVED ELEMENTARY BOUNDARY COMPLETION + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED; LRC(14) OPEN. THM-4148's first-window width transfer holds for
  every nonempty finite positive body, with no minimum-label hypothesis.
  The minimum-one and minimum-two boundary strata use six explicit labelled
  physical clocks. The complete eleven-body width census is therefore
  60,301,653,510. Arbitrary bodies outside the width gate, entry, and LRC(14)
  remain open.
source: codex-lrc14-planar-jc-breakthrough-20260825
depends_on:
  - THM-4136-fixed-body-universal-odd-tail-lrc14-completion
  - THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
related:
  - THM-689-dead-zone-lemma-k7-rigidity-and-consecutive-maximality
  - THM-4142-common-safe-arc-clock-pool-universal-odd-tail-lrc14-completion
script: 04-computation/lrc14_minimum_boundary_universal_odd_tail_transfer_thm4149.py
output: 05-knowledge/results/lrc14_minimum_boundary_universal_odd_tail_transfer_thm4149.out
independent_audit_script: 04-computation/lrc14_minimum_boundary_universal_odd_tail_transfer_thm4149_independent_audit.cpp
independent_audit_output: 05-knowledge/results/lrc14_minimum_boundary_universal_odd_tail_transfer_thm4149_independent_audit.out
script_sha256: cf146ae0b9a5fbbbdb6bf8597550c7923848cea1ae13ccba453b11c1b855e1c7
output_sha256: 84bca834c0ccbbb656c7d06ebc5f4fa30d4a9fb682d46464518dff400a0589d8
semantic_sha256: 90de52f9075568cead74de79315daf0494ad7f1f49e5ae7303eff0eae7881f6d
semantic_fnv64: 6784b7e0b01a759d
independent_audit_script_sha256: 283477210905e2ae4729ad978a9e67c7496581a413ac8a8fb6dc5b616a2df35b
independent_audit_output_sha256: b8457993e1ae0802aa6b9d65d36cdb6b0a501dfb6f9e9f8ad084ab3a1ea5f036
independent_semantic_fnv64: 6784b7e0b01a759d
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact Fraction arithmetic checks both small-minimum caps, all 68
  primitive residual ratios, the labelled 56+11+1 and 58+8+2 clock
  partitions, direct full-superbody clearances, the (1,13) first-window
  hostile component, the isolated 3/14 body clock, endpoint hostiles and
  rescues, and the complete enlarged family count.
independent_audit: >
  ACCEPT. A separate warning-clean standard-library C++ implementation uses
  normalized integer rationals, reconstructs the hostile cross-comb cell,
  independently exhausts both clock partitions and the census, and matches
  the semantic FNV. Optimized, unoptimized, and sanitizer runs agree.
---

# THM-4149 -- minimum-boundary first-window odd-tail transfer

**PROVED ELEMENTARY BOUNDARY COMPLETION + VERIFIED-EXACT + INDEPENDENTLY
AUDITED; LRC(14) REMAINS OPEN.**

THM-4148 proved its transfer under `min(H)>=3`. The width inequality itself
compresses the two omitted minima into tiny universes. The minimum-one
universe contains a hostile showing that the old endpoint argument genuinely
stops, but it also exposes an isolated safe clock that completes the theorem.

## 1. The strengthened theorem

Let `H` be any nonempty finite set of positive integers and put

```text
m=min H,                 M=max H,
W=13/(14M)-1/(14m).                                      (1)
```

Assume only

```text
W>=2/189,                                                   (2)
```

or equivalently

```text
27(13m-M)>=4mM.                                            (3)
```

> **Theorem.** For every pair of distinct positive odd integers `a,b`,
> there is an `x in R/Z` such that
>
> ```text
> min_(v in 2H union {a,b}) ||vx|| >= 1/14.                (4)
> ```

Thus THM-4148's conclusion holds with no minimum-label hypothesis and for
bodies of every cardinality.

## 2. The inherited nonresidual scales

Interchange the tails if necessary and write

```text
a=pt,                    b=qt,
t=gcd(a,b),              0<p<q,                           (5)
```

where `p,q,t` are odd and `gcd(p,q)=1`. The first safe window

```text
I=[1/(14m),13/(14M)]                                     (6)
```

is a closed body-safe interval of length `W`. Its two physical half-lifts

```text
x_s=(y+s)/2,              s in {0,1},                    (7)
```

keep `2H` safe. THM-4136's open cross-comb `C_(p,q)` consists of the quotient
phases `ty` for which both half-lifts are tail-bad, and every one of its
components has length

```text
beta(p,q)<=2/63,          beta(p,q)<=2/(7q) for q>=9.    (8)
```

The compact-inside-open argument in THM-4148 therefore applies without any
condition on `m`: if `t>=3`, then `tW>=2/63`; if `t=1` and `q>=27`, then
`W>=2/(7q)`. In either case the closed image of `I` cannot lie in one open
cross-comb component. Hence only

```text
t=1,                    q<=25                            (9)
```

remains. There are exactly 68 coprime odd pairs `0<p<q<=25`.

## 3. Minimum two: a common interval and two endpoints

Suppose `m=2`. Equation `(3)` gives

```text
702>=35M,                M<=20.                          (10)
```

Every such body is contained in `{2,...,20}`, and the common interval

```text
J_2=[1/28,13/280]                                         (11)
```

is body-safe: for `2<=h<=20` and `y in J_2`,

```text
1/14<=hy<=13/14.                                         (12)
```

The residual pairs in `(9)` split into three exact clock classes.

1. If `q<=23`, take `y=1/28` and the upper sheet `x=29/56`. For every odd
   `r<=23`,

   ```text
   ||rx||=(28-r)/56>=5/56>1/14.                          (13)
   ```

2. If `q=25` and `p>=5`, take the lower sheet `x=1/56`. Both tail gaps are
   at least `5/56`.

3. The remaining pairs are `(1,25)` and `(3,25)`. Take `y=13/280` and its
   upper sheet `x=293/560`. The gaps of tails `1,3,25` are respectively

   ```text
   267/560,             241/560,             45/560=9/112>1/14.  (14)
   ```

These classes contain `58`, `8`, and `2` primitive ratios. This proves `(4)`
when `m=2`.

## 4. Minimum one: an isolated clock repairs a real trap

Suppose `m=1`. Equation `(3)` gives

```text
351>=31M,                M<=11.                          (15)
```

It is therefore enough to keep the full body

```text
H_1={1,2,...,11}                                          (16)
```

safe. The 68 ratios again have three exact clock classes.

1. If `p>=3`, take `y=1/14` and the lower sheet `x=1/28`. Every odd
   `3<=r<=25` has

   ```text
   ||r/28||>=3/28>1/14.                                  (17)
   ```

2. If `p=1` and `q!=13`, take

   ```text
   y=6/77,              x=(y+1)/2=83/154.                (18)
   ```

   The body clearance at `y` is `6/77>1/14`, while the tail `1` has gap
   `71/154`. For odd `3<=q<=11`, the other gap is

   ```text
   (77-6q)/154>=11/154=1/14;                             (19)
   ```

   for odd `15<=q<=25`, it is `(6q-77)/154>=13/154`.

3. For the single pair `(p,q)=(1,13)`, take the isolated body-safe phase

   ```text
   y=3/14,              x=y/2=3/28.                      (20)
   ```

   Since `gcd(3,14)=1` and `1<=h<=11<14`, the residue `3h mod 14` is
   nonzero for every body label, so every member of `H_1` has gap at least
   `1/14`.
   The two tail gaps at `x` are `3/28` and `11/28`.

The three classes contain `56`, `11`, and `1` primitive ratios. This proves
`(4)` for `m=1`; THM-4148 already proves it for `m>=3`, completing the
theorem. **QED.**

The last clock is not cosmetic. For `(p,q)=(1,13)`, the first window

```text
[1/14,13/154]                                             (21)
```

lies strictly inside the open cross-comb component

```text
(6/91,8/91).                                              (22)
```

Thus both physical sheets are tail-bad everywhere on that window. At
`y=3/14`, body speeds `5` and `9` attain `1/14` with opposite local slopes;
every sufficiently close phase on either side is body-bad. The successful
clock is an isolated boundary point, exactly the kind of object lost by an
open-component-only search.

## 5. Complete eleven-body width census

For an eleven-element body with minimum `m`, equation `(3)` is equivalent to

```text
M<=U(m)=floor(351m/(4m+27)).                              (23)
```

After fixing the unique minimum, the other ten labels may be chosen freely
from `{m+1,...,U(m)}`. THM-4148 counted the minima `3<=m<=70`. The two new
rows are

```text
m=1:  U(1)=11,       binom(10,10)=1,
m=2:  U(2)=20,       binom(18,10)=43,758.                (24)
```

Therefore the complete count becomes

```text
60,301,609,751 + 1 + 43,758 = 60,301,653,510.            (25)
```

The minimum range is now exactly `1,...,70`, and the maximum label remains
`80`. Each of these bodies accepts every distinct positive odd tail pair
after doubling. This is a large universal family theorem, not LRC(14): it
does not show that an arbitrary primitive thirteen-speed row enters this
width stratum.

## 6. Typed connection and audit scope

```text
source:       a width-qualified positive body H and an arbitrary odd tail pair
target:       a labelled physical half-lift x_s over a body phase y
map:          y |-> (y+s)/2 after primitive gcd reduction of the tails
preserved:    the 1/14 clearance of every doubled body speed
destroyed:    the physical sheet label if the quotient phase is kept alone
sidecar:      the six small-minimum labelled clocks and the open cross-comb
decisive test: 68 primitive ratios, including the trapped (1,13) cell.       (26)
```

The primary exact audit uses `Fraction` arithmetic. It checks all 68 residual
ratios against their full superbodies, reconstructs the hostile component,
tests both neighbours of the isolated phase, and recomputes the census. The
independent C++ audit uses normalized integer rationals, reconstructs the
wall cell independently, and matches the semantic FNV under optimized,
unoptimized, and sanitizer builds. Neither computation is a substitute for
the infinite cross-comb proof in Section 2.
