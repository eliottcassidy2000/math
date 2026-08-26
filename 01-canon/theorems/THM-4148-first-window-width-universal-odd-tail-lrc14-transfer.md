---
id: THM-4148
title: "First-window width transfer for universal odd-tail LRC(14) completion"
status: >
  PROVED ELEMENTARY FIRST-WINDOW/CROSS-COMB TRANSFER + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED; LRC(14) OPEN. Every finite positive body H with
  min(H)>=3 and first-window width at least 2/189 has universal 1/14-safe
  completion after doubling by every pair of distinct positive odd tails.
  The 46-label block {15,...,60} therefore supplies 13,340,783,196
  eleven-body families, far beyond THM-4142's fixed pool. Forty-six is sharp
  among consecutive blocks for this first-window mechanism. Arbitrary bodies,
  entry, and LRC(14) remain open.
source: codex-frontier-synthesis-creative-20260825ah
depends_on:
  - THM-689-dead-zone-lemma-k7-rigidity-and-consecutive-maximality
  - THM-4136-fixed-body-universal-odd-tail-lrc14-completion
related:
  - THM-526-lrc-arc-width-lemma-large-stranger-discharge
  - THM-4129-universal-two-speed-completion-of-the-eleven-speed-lrc14-body
  - THM-4142-common-safe-arc-clock-pool-universal-odd-tail-lrc14-completion
script: 04-computation/lrc14_first_window_width_universal_odd_tail_transfer_thm4148.py
output: 05-knowledge/results/lrc14_first_window_width_universal_odd_tail_transfer_thm4148.out
independent_audit_script: 04-computation/lrc14_first_window_width_universal_odd_tail_transfer_thm4148_independent_audit.cpp
independent_audit_output: 05-knowledge/results/lrc14_first_window_width_universal_odd_tail_transfer_thm4148_independent_audit.out
script_sha256: ffee89d8ec8833a5533ade63edf1ee3d34ca89c0dfe77b7a2651575f09d2bdf4
output_sha256: ed11811dbb0dc7fef83e00431d30276318360938b58b459bf7637479079e9677
semantic_sha256: 6f4c7b6fd62bca62f0a8baccb63822e6f19df831d2ece22001c2412cfce44e2e
semantic_fnv64: 1f67e51174709013
independent_audit_script_sha256: e617521359ad7ea8ee02c646e78fec69cad8c07a430fe4e47ce0e8352e6880cc
independent_audit_output_sha256: 222c4cc324b87c6f4a2464db22d2fa369b39320dff53d0ff14cb4289dc53650d
independent_semantic_fnv64: 1f67e51174709013
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact Fraction arithmetic checks both analytic scale gates, the sharp
  residual clock, all nine maximal 46-label blocks, impossibility of 47
  labels, the complete {15,...,60} endpoint clearances, the moving-base
  hostile and rescue, and the eleven-subset count. Normal, optimized, and
  hash-seeded replays byte-match. The cross-comb topology is inherited from
  THM-4136's independently audited exact proof.
independent_audit: >
  ACCEPT. A separate warning-clean standard-library C++ implementation uses
  normalized integer rationals and direct modular gaps. It independently
  reproduces every gate, discriminant, endpoint, hostile, rescue, family
  count, and the shared semantic FNV. Optimized, unoptimized, and sanitizer
  controls agree.
---

# THM-4148 -- first-window width odd-tail transfer

**PROVED ELEMENTARY FIRST-WINDOW/CROSS-COMB TRANSFER + VERIFIED-EXACT +
INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

THM-4142 enlarged one fixed body by intersecting its arc and clock
certificates before selecting eleven speeds. The pool boundary is not a
boundary of the cross-comb mechanism. The underused coordinate is the
elementary first safe window, which depends only on the smallest and largest
body speeds and can be moved after seeing the tail ratio.

## 1. Theorem

Let `H` be a nonempty finite set of positive integers, put

```text
m=min H,                     M=max H,                    (1)
W=13/(14M)-1/(14m),                                      (2)
```

and suppose

```text
m>=3,                        W>=2/189.                   (3)
```

> **Theorem.** For every two distinct positive odd integers `a,b`, there is
> a phase `x in R/Z` such that
>
> ```text
> min_(v in 2H union {a,b}) ||vx|| >= 1/14.             (4)
> ```

Equivalently, the width condition is the integral inequality

```text
27(13m-M)>=4mM.                                         (5)
```

The theorem applies to bodies of any cardinality. The LRC(14) corollary uses
eleven-element subbodies.

## 2. The first window and its two physical sheets

Put

```text
I=[1/(14m),13/(14M)].                                   (6)
```

For `h in H` and `y in I`, no wrap occurs and

```text
1/14 <= hy <= 13/14.                                    (7)
```

Thus the closed interval `I` is `1/14`-safe for all of `H`. This is the
THM-689 first-window mechanism, included here directly so that both closed
endpoints are explicit. The two half-lifts

```text
x_0=y/2,                       x_1=(y+1)/2              (8)
```

satisfy `2h x_i=hy mod 1`, so both keep the entire doubled body `2H` safe.

Interchange the tails if necessary and write

```text
a=pt,                         b=qt,
t=gcd(a,b),                   0<p<q,                    (9)
```

where `p,q,t` are odd and `gcd(p,q)=1`. THM-4136 defines the open quotient
danger set `C_(p,q)` of phases `w=ty` for which both physical lifts in `(8)`
are tail-bad. Its cross-comb theorem gives

```text
0 notin C_(p,q),
every component has length beta(p,q)<=2/63,
beta(p,q)<=2/(7q).                                      (10)
```

The uniform and shrinking bounds in `(10)` will handle every nonresidual
scale without another body-specific clock census.

## 3. Compact-to-open width closure

Assume first, contrary to `(4)`, that both lifts are tail-bad for every
`y in I`. Then the circular image `tI` is a compact connected subset of the
open set `C_(p,q)`. If `tW>=1`, this image is the whole circle, contradicting
`0 notin C_(p,q)`. If `tW<1` but the image wraps, it contains zero and gives
the same contradiction. Otherwise it is a compact interval in one open
component. The endpoint types then force the strict inequality

```text
tW<beta(p,q):                                             (11)
```

a closed interval of equal length cannot lie inside an open component.

If `t>=3`, hypothesis `(3)` and the first bound in `(10)` give

```text
tW>=3(2/189)=2/63>=beta(p,q),                            (12)
```

contradicting `(11)`. It remains to take `t=1`. If `q>=27`, the second bound
in `(10)` gives

```text
W>=2/189=2/(7*27)>=2/(7q)>=beta(p,q),                    (13)
```

again contradicting `(11)`.

## 4. One endpoint clock closes every residual ratio

The only remaining case has `t=1` and odd `p<q<=25`. Take the lower closed
endpoint of `(6)` and its upper physical lift:

```text
y=1/(14m),                    x=(y+1)/2=(14m+1)/(28m).  (14)
```

For every odd `r<=25`, hypothesis `m>=3` gives

```text
0<r/(28m)<=25/84<1/2,
||rx||=1/2-r/(28m)>=17/84>1/14.                         (15)
```

Apply `(15)` to `r=p,q`. The body is safe because `2x=y mod 1`. This closes
the residual ratios and proves `(4)`. **QED.**

## 5. A 46-label superbody and 13,340,783,196 families

Take

```text
H={15,16,...,60}.                                       (16)
```

Its first window is

```text
I=[1/210,13/840],       W=3/280,
W-2/189=1/7560>0.                                       (17)
```

Hence for every distinct positive odd pair `a,b`, the full 48-speed set

```text
2{15,...,60} union {a,b}                                (18)
```

is `1/14`-safe. Safety is hereditary under deleting speeds. Every
eleven-subset `B` of `{15,...,60}` therefore gives a thirteen-speed LRC(14)
row `2B union {a,b}`, and there are

```text
binom(46,11)=13,340,783,196                              (19)
```

such bodies. These bodies are not restricted to THM-4142's 26-speed pool.

## 6. Consecutive-block sharpness for this mechanism

For a consecutive block `H={A,...,B}`, condition `(5)` is exact. A block of
at least 47 labels has `B>=A+46`. Increasing `B` only shrinks the first
window, so a necessary condition is

```text
4A^2-140A+1242<=0.                                      (20)
```

Its discriminant is `-272`, so `(20)` is impossible. For exactly 46 labels,
`B=A+45`, and `(5)` becomes

```text
4A^2-144A+1215<=0.                                      (21)
```

The roots are `13.5` and `22.5`. Thus the maximal consecutive blocks for
this first-window transfer are exactly

```text
{A,...,A+45},                    14<=A<=22.              (22)
```

This is sharpness of the stated min/max mechanism, not a claim that another
body certificate cannot handle 47 labels.

## 7. Moving-base hostile and information ledger

The residual endpoint clock is deliberately used only for `q<=25`. For the
body `(16)`, the lower base `y=1/210` has lifts

```text
x_0=1/420,                      x_1=211/420.             (23)
```

The tails `(1,211)` kill opposite lifts:

```text
||x_0||=1/420,        ||211x_1||=1/420,
211^2=1 mod 420.                                         (24)
```

So no fixed choice of the endpoint sheet proves the large-ratio lane. Moving
inside the same body interval to

```text
y=1/105,                       x=53/105                  (25)
```

rescues the complete 48-speed set, with tail gaps `52/105` and full-row
clearance exactly `1/7`. The movable compact interval in Section 3 is
load-bearing.

```text
source:       first body window I determined by (m,M)
target:       every doubled-body row 2H union {a,b}, a,b odd
map:          y -> {y/2,(y+1)/2}, then w=gcd(a,b)y
preserved:    all body gaps, closed endpoints, both physical sheets
destroyed:    internal body structure beyond min/max
sidecar:      open cross-comb C_(p,q), not one fixed lift
hostile:      (23)--(24) defeats a universal endpoint clock
decisive test: W versus 2/63 and 2/(7q), then clock (14). (26)
```

## 8. Audit and scope

The primary and independent artifacts exactly reproduce `(12)--(25)`, the
nine blocks in `(22)`, and count `(19)`. Replay with

```text
python3 -B 04-computation/lrc14_first_window_width_universal_odd_tail_transfer_thm4148.py
python3 -B -O 04-computation/lrc14_first_window_width_universal_odd_tail_transfer_thm4148.py
PYTHONHASHSEED=271828 python3 -B 04-computation/lrc14_first_window_width_universal_odd_tail_transfer_thm4148.py
clang++ -std=c++17 -O3 -Wall -Wextra -Werror \
  04-computation/lrc14_first_window_width_universal_odd_tail_transfer_thm4148_independent_audit.cpp \
  -o /tmp/thm4148-independent-audit
/tmp/thm4148-independent-audit
```

The theorem does not cover bodies failing `(3)`, mixed/even tails, physical
entry into this branch, or arbitrary LRC(14) instances. LRC(14) remains
open. **QED.**
