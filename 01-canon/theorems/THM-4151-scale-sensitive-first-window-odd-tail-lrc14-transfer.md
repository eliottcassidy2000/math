---
id: THM-4151
title: "Scale-sensitive anchored first-window transfer for universal odd-tail LRC(14) completion"
status: >
  PROVED ELEMENTARY ANCHORED-CARRIER TRANSFER + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED; LRC(14) OPEN. Every finite positive body H with
  m=min(H), M=max(H), and 16M<=156m+13 has universal 1/14-safe completion
  after doubling by every pair of distinct positive odd tails. The affine
  cap is sharp for this arbitrary-cardinality theorem: the first excluded
  row 2{2,...,21} union {1,25} is unsafe. Every finite additive pattern is
  certified after an explicit translation. Among eleven-subsets of [1,N],
  the anchored family has asymptotic density (35/39)^10; combining THM-4148
  adds the lone minimum-one body {1,...,11}. Entry, arbitrary thirteen-speed
  rows, and LRC(14) remain open.
source: codex-lrc14-planar-jc-creative-20260825
depends_on:
  - THM-689-dead-zone-lemma-k7-rigidity-and-consecutive-maximality
  - THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
related:
  - THM-4129-universal-two-speed-completion-of-the-eleven-speed-lrc14-body
  - THM-4136-fixed-body-universal-odd-tail-lrc14-completion
script: 04-computation/lrc14_scale_sensitive_anchored_first_window_transfer_thm4151.py
output: 05-knowledge/results/lrc14_scale_sensitive_anchored_first_window_transfer_thm4151.out
independent_audit_script: 04-computation/lrc14_scale_sensitive_anchored_first_window_transfer_thm4151_independent_audit.cpp
independent_audit_output: 05-knowledge/results/lrc14_scale_sensitive_anchored_first_window_transfer_thm4151_independent_audit.out
script_sha256: df324c947f753edf1d0c175e198c3e9c4da25957e60abbc495366f6d6f5a7ec5
output_sha256: bf4c413749623b28ea6884c66bf1acbe002e628b3912fa83f7824c6c653948d1
semantic_sha256: 84adc34f0dc8ffb415455827a227f86cab1bea3169d73d3912721bd815e1b514
semantic_fnv64: ff11c60c57deec04
independent_audit_script_sha256: 6bfc0d438e7c04d4a5b27ed9bde53fcd3bebd79c958cbba3b69080f4fca93de4
independent_audit_output_sha256: a3a2e244ca95886fa1e08ad3e4876a378027e9f127856fbde7cd14343ac1609d
independent_semantic_fnv64: ff11c60c57deec04
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact Fraction arithmetic checks the affine cap, the endpoint and
  odd-wall slack identities, 1,355 primitive strict-wall rows, 25 exact
  carrier failures, the literal unsafe cap hostile, translated-pattern
  minima, finite family counts, and the asymptotic density constant. Normal,
  optimized, and hash-seeded replays agree.
independent_audit: >
  ACCEPT. A separate warning-clean standard-library C++ implementation uses
  normalized integer rationals, explicit lower/upper lift walls, and direct
  physical safe-component reconstruction. It independently reproduces every
  gate, hostile, translation, family count, density numerator and denominator,
  and the shared semantic FNV. Optimized and unoptimized controls agree.
---

# THM-4151 -- scale-sensitive anchored odd-tail transfer

**PROVED ELEMENTARY ANCHORED-CARRIER TRANSFER + VERIFIED-EXACT +
INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

THM-4148 used only the length of the first body-safe interval. Its compact
interval could therefore be compared only with an arbitrarily placed open
cross-comb component. The missing coordinate is that this interval is not
arbitrarily placed: its lower endpoint is exactly `1/(14m)`. Odd tail walls
cannot approach that rational endpoint continuously. Their parity forces one
lattice unit of slack, and that one unit changes the nonlinear width condition
into an affine maximum cap.

## 1. The anchored theorem

Let `H` be a nonempty finite set of positive integers. Put

```text
m=min H,                         M=max H,
W=13/(14M)-1/(14m),
G_m=(4m-1)/(14m(12m+1)).                              (1)
```

> **Theorem.** If `W>=G_m`, then for every two distinct positive odd
> integers `a,b`, there is an `x in R/Z` such that
>
> ```text
> min_(v in 2H union {a,b}) ||vx|| >= 1/14.            (2)
> ```

The width hypothesis is equivalent to the affine integral gate

```text
(12m+1)(13m-M)>=M(4m-1)
                    iff 16M<=156m+13.                  (3)
```

Thus the exact maximum admitted by this theorem is

```text
U(m)=floor((156m+13)/16).                              (4)
```

The statement applies to bodies of every cardinality. Taking `|H|=11`
produces thirteen-speed LRC(14) rows in the indicated parity branch.

## 2. First window and Boolean sheet capacity

Put

```text
I=[y_0,y_1]=[1/(14m),13/(14M)].                        (5)
```

For every `h in H` and `y in I`, no wrap occurs and

```text
1/14<=hy<=13/14.                                      (6)
```

Hence `I` is safe for `H`, and both physical half-lifts

```text
x_0=y/2,                         x_1=(y+1)/2            (7)
```

are safe for `2H`.

For an odd tail `r`, let `K_r(y)` be the set of sheets in `(7)` on which
`r` is strictly bad. The two tail phases differ by `1/2`, so

```text
|K_r(y)|<=1.                                           (8)
```

If both physical sheets fail for the two tails, then

```text
K_a(y) union K_b(y)={0,1}.                             (9)
```

Thus equality holds in the Boolean join-rank bound
`|K_a union K_b|<=|K_a|+|K_b|`: each tail kills exactly one sheet and the
two killers are complementary. This recovers the two cross-comb orientations
without quotienting the physical tail labels.

Interchange the tails so that `a<b`. Condition `K_b(y)!=empty` is simply

```text
||by||<1/7.                                            (10)
```

Indeed, an even nearest integer in `(10)` means that `b` kills `x_0`, and
an odd nearest integer means that it kills `x_1`. Consequently the complete
`b`-carrier is the disjoint union of open teeth

```text
E_b=union_(k in Z) ((k-1/7)/b,(k+1/7)/b),              (11)
```

each of width `2/(7b)`.

## 3. The odd-wall slack

Suppose, contrary to `(2)`, that both lifts fail for every `y in I`. By
`(8)--(10)`, the connected interval `I` lies in one open tooth of `E_b`.
Every tooth wall has the form

```text
lambda=(7k+-1)/(7b).                                   (12)
```

At the anchored endpoint `y_0=1/(14m)`,

```text
|y_0-lambda|=|b-2m(7k+-1)|/(14mb)>=1/(14mb).           (13)
```

The numerator is an odd integer, hence is nonzero and has absolute value at
least one. In particular, the distance from `y_0` to the left wall of its
tooth is at least `1/(14mb)`. The room remaining to the right wall is at most

```text
2/(7b)-1/(14mb)=(4m-1)/(14mb).                         (14)
```

There are now two exhaustive cases.

If `b<=12m`, then oddness gives `b<=12m-1`. At the lower endpoint choose the
upper lift

```text
x=(1+y_0)/2=(14m+1)/(28m).                             (15)
```

For either odd tail `r<=b`,

```text
||rx||=1/2-r/(28m)>1/14.                               (16)
```

This contradicts total failure.

If `b>12m`, oddness gives `b>=12m+1`, so `(14)` is at most

```text
(4m-1)/(14m(12m+1))=G_m.                              (17)
```

But a closed interval of width `W>=G_m` cannot start inside and remain inside
an open tooth whose right room is at most `G_m`. Equality reaches a safe wall.
This again contradicts total failure and proves `(2)`. **QED.**

## 4. Exact sharpness and the minimum-one rescue

The constants in the carrier proof are exact. Take

```text
b=12m+1.                                               (18)
```

The tooth

```text
(6/(7b),8/(7b))                                        (19)
```

contains `y_0`; its left slack is exactly `1/(14mb)` and its right room is
exactly `G_m`. If `M=U(m)+1`, then

```text
13/(14M)<8/(7b),                                       (20)
```

so the whole first window lies strictly inside `(19)`. Tails `(1,b)` kill
opposite lifts throughout that window. Thus the affine cap cannot be enlarged
by even one using the anchored first-window carrier.

At `m=1`, the first excluded row is not globally unsafe. For

```text
H={1,...,11},                    (a,b)=(1,13),          (21)
```

both first-window lifts fail, but the physical phase `x=5/24` gives the full
row clearance `1/12`. This is the exceptional clock already available through
THM-4148.

At `m=2`, the cap is globally sharp for the theorem's arbitrary-cardinality
scope. Exact wall decomposition for

```text
H={2,...,21}                                            (22)
```

gives its complete quotient safe set

```text
[1/28,13/294] union [281/294,27/28].                   (23)
```

The first component lies strictly inside the tail-`25` carrier
`(6/175,8/175)` and the reflected component lies inside
`(167/175,169/175)`. Tail `1` kills the complementary physical sheet on each
component. Therefore

```text
2{2,...,21} union {1,25}                               (24)
```

has no `1/14`-safe phase. This is a 22-speed hostile, not a counterexample to
LRC(14).

## 5. Every finite additive pattern eventually enters

Let `S` be any finite subset of the nonnegative integers with `min S=0` and
span `D=max S`. Translate it to

```text
H=m+S.                                                 (25)
```

Then `min H=m`, `max H=m+D`, and `(3)` reduces to

```text
140m>=16D-13.                                          (26)
```

Hence every such pattern is universally completable by every odd tail pair
after any translation

```text
m>=max(1,ceil((16D-13)/140)).                           (27)
```

This sharply changes the consecutive-block picture from THM-4148's constant
width gate. A 47-label block has `D=46`, so its first admitted translate is

```text
{6,...,52}.                                            (28)
```

In fact `U(6)=59`, and the maximal consecutive block beginning at six has 54
labels, `{6,...,59}`. More generally the maximal length at minimum `m` is

```text
U(m)-m+1=floor((140m+13)/16)+1.                        (29)
```

It grows linearly rather than stopping at a universal block length.

## 6. A positive-density infinite eleven-body family

Let `A(N)` count eleven-element subsets of `{1,...,N}` satisfying the anchored
cap `(4)`. Grouping by the unique minimum gives the exact formula

```text
A(N)=sum_(m=1)^N binom(min(N,U(m))-m,10),               (30)
```

where the binomial coefficient is zero if its upper argument is below ten.
The first nonzero minimum is `m=2`. Exact controls include

```text
A(20)=75,582,                     A(80)=3,548,681,310,136,
A(200)=131,378,242,150,108,190.                         (31)
```

Put `c=39/4`. Since

```text
U(m)=floor(cm+13/16)=cm+O(1),                           (32)
```

and `binom(k,10)=k^10/10!+O(k^9)` uniformly for `0<=k<=N`, division by
`binom(N,11)` turns `(30)` into the Riemann sum

```text
lim_(N->infinity) A(N)/binom(N,11)
 =11 integral_0^1 (min(cx,1)-x)^10 dx
 =(1-1/c)^10
 =(35/39)^10
 =2,758,547,353,515,625 / 8,140,406,085,191,601
 ~=0.33887097580227127.                                (33)
```

Combining with THM-4148 adds exactly one eleven-body for every `N>=11`,
namely `{1,...,11}` at minimum one. Thus the combined family has count
`A(N)+1` and the same density `(33)`.

## 7. Information ledger and cross-frontier boundary

```text
source:       anchored first body window [1/(14m),13/(14M)]
target:       every row 2H union {a,b} with a,b distinct positive odd
map:          y -> {y/2,(y+1)/2}; tail r -> sheet-kill set K_r(y)
preserved:    closed body safety, physical sheet labels, strict wall parity
destroyed:    internal body structure beyond min/max
sidecar:      the larger actual tail b and its open 1/7 carrier E_b
hostile:      2{2,...,21} union {1,25} is literally unsafe
decisive test: one odd wall unit converts width into 16M<=156m+13. (34)
```

The Boolean join equality in `(9)` has a genuine analogue in planar-Jacobian
monodromy: cycle partitions replace sheet-kill sets, and transitivity means
their join is the top partition. What does **not** transfer is interval length.
The LRC carrier components are ordered open intervals with trivial two-sheet
cover, whereas a punctured complex Jacobian carrier can support nontrivial
sheet monodromy around a compact loop. Degree two alone is not a common proof
mechanism; the carrier topology is load-bearing.

## 8. Audit and scope

The primary and independent artifacts reconstruct `(3)--(33)`. They also
enumerate every strict lower/upper-lift wall for `1,355` primitive controls
with `m<=5,q<=51`, confirm the first carrier failure for `m<=25`, and rebuild
all physical safe components of the literal hostile `(24)`.

Replay with

```text
python3 -B 04-computation/lrc14_scale_sensitive_anchored_first_window_transfer_thm4151.py
python3 -B -O 04-computation/lrc14_scale_sensitive_anchored_first_window_transfer_thm4151.py
PYTHONHASHSEED=271828 python3 -B 04-computation/lrc14_scale_sensitive_anchored_first_window_transfer_thm4151.py
g++ -std=c++17 -O2 -Wall -Wextra -Werror \
  04-computation/lrc14_scale_sensitive_anchored_first_window_transfer_thm4151_independent_audit.cpp \
  -o /tmp/thm4151-independent-audit
/tmp/thm4151-independent-audit
```

This theorem does not handle mixed/even tails, prove physical entry into the
two-odd-tail branch, or treat arbitrary thirteen-speed rows. The unsafe row
`(24)` has 22 speeds and says nothing against LRC(14). The Lonely Runner
Conjecture for 14 runners remains open. **QED.**
