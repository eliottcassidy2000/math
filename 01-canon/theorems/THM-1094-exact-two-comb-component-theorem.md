---
id: THM-1094
title: Exact two-comb component theorem for the three-killer clustered stratum
status: PROVED — elementary all-scale mass/component tail plus a finite-exact rational endpoint certificate over all 66 ten-speed cores and 9,246,070 guarded pairs; zero failures; independently replayed core minimum, bank cardinality, and hardest finite-bank row
source: codex-2026-07-18-S73 frontier audit
depends_on: [THM-1051]
related: [THM-1042, THM-1061, THM-1071, THM-1081, MISTAKE-163]
script:
  - 04-computation/lrc14_r3_two_comb_component_exact_codex_S73.cpp
  - 04-computation/lrc14_r3_two_comb_extremal_replay_codex_S73.py
result:
  - 05-knowledge/results/lrc14_r3_two_comb_component_exact_codex_S73.out
  - 05-knowledge/results/lrc14_r3_two_comb_extremal_replay_codex_S73.out
---

# THM-1094 — the exact two-comb component theorem

For a positive integer `k`, write

```text
D_k = {t in [0,1] : ||kt|| < 1/14}
```

for its danger comb.  For a ten-element core
`P subset {1,...,12}`, put

```text
S(P) = [0,1] minus union_{p in P} D_p.
```

Endpoint conventions do not affect any interval length below.

## Theorem

For every ten-element `P subset {1,...,12}` and every pair of integers

```text
13 max(P) < k1 < k2,                                      (1)
```

the exact remainder

```text
S(P) minus (D_k1 union D_k2)                              (2)
```

has a component of length

```text
L(P;k1,k2) > 1/(3k2).                                    (3)
```

Consequently, if `k3>k2`, the thirteen-speed family
`P union {k1,k2,k3}` has a time at which every runner is at distance at
least `1/14` from the origin.  Thus the three-killer clustered stratum is
uniformly closed; the bounded horns in THM-1061 are independent checks, not
load-bearing parts of the proof.

## 1. The 66-core exact atlas

Every endpoint of `S(P)` belongs to the finite rational arrangement

```text
(14j-1)/(14p), (14j+1)/(14p),    p in P.                 (4)
```

Classifying the cells between consecutive endpoints gives `S(P)` exactly.
Let `ell(P)` be the length of its largest component.  Exact enumeration of
all `C(12,10)=66` cores gives

```text
ell(P) >= 1/112.                                         (5)
```

Equality occurs only for

```text
P0={1,2,3,4,7,8,9,10,11,12},                            (6)
```

whose two longest components are

```text
[43/112,11/28] and [17/28,69/112].                      (7)
```

The C++ certificate obtains the same atlas by exact tooth subtraction.  The
independent Python replay instead constructs the complete breakpoint set (4),
tests one rational midpoint in every cell, and reproduces (5)--(7).

## 2. A uniform analytic tail

Fix a largest core-safe interval `J` of length `ell=ell(P)`.  A tooth of
`D_k` meeting `J` has centre `j/k` in the interval obtained by enlarging `J`
by `1/(14k)` on both sides.  After multiplying by `k`, that centre interval
has length `k ell+1/7`.  An interval of real length `d` contains at most
`d+1` integers, so the number `N_k(J)` of teeth meeting `J` satisfies

```text
N_k(J) <= k ell + 8/7.                                   (8)
```

There are two useful occupied-length bounds.  Summing the teeth counted in
(8) gives the elementary estimate

```text
|J intersect D_k| <= ell/7 + 8/(49k).                   (9)
```

The exact periodic discrepancy is sharper:

```text
|J intersect D_k| <= ell/7 + 6/(49k).                  (10)
```

To prove (10), scale `J` by `k` and write its length as `m+s`, with `m` an
integer and `0<=s<1`.  The `m` full periods contribute `m/7`.  The residual
interval meets the circle danger arc, of length `1/7`, in at most
`min(s,1/7)`.  Its excess over the mean `s/7` is at most

```text
max_{0<=s<=1}(min(s,1/7)-s/7)=6/49.
```

Dividing by `k` proves (10).

After removing `D_k1` and `D_k2`, the surviving mass in `J` is therefore at
least

```text
A = 5ell/7 - 6/(49k1) - 6/(49k2),                      (11)
```

and the number of surviving components is at most

```text
C = N_k1(J)+N_k2(J)+1
  <= ell(k1+k2)+23/7.                                  (12)
```

Therefore one component has length at least `A/C`.  Direct multiplication
shows that `A/C>1/(3k2)` whenever

```text
56ell k2 - 49ell k1 - 18(k2/k1) - 179 > 0.            (13)
```

Put `r=k2/k1>1`.  Condition (13) is

```text
ell k1(56r-49) > 18r+179.                              (14)
```

If `ell k1 >= 197/7`, then the left side minus the right side is at least

```text
(197/7)(56r-49)-(18r+179) = 1558(r-1) > 0.             (15)
```

Thus every pair with `ell(P)k1>=197/7` is proved analytically.  In
particular, (5) makes every `k1>=3152` automatic.

For a smaller fixed `k1`, the coefficient
`56ell-18/k1` is positive: (1) gives `k1>=131`, while (5) gives
`56ell>=1/2`.  Solving (13) for `k2` shows that every

```text
k2 > k1(49ell k1+179)/(56ell k1-18)                    (16)
```

is also automatic.  Equations (15)--(16) leave a finite, explicit endpoint
bank.

## 3. The exact finite bank

For audit simplicity, the certificate deliberately scans a larger bank than
the sharp tail in Section 2 requires.  It uses the weaker tooth-sum bound (9),
whose sufficient condition is

```text
56ell k2-49ell k1-24(k2/k1)-185>0.                     (17)
```

For each core it scans the guarded superset

```text
13 max(P) < k1 <= floor(209/(7ell(P)))+1,              (18)
k1 < k2 <= floor(k1(49ell(P)k1+185)
                 /(56ell(P)k1-24))+1.                 (19)
```

The added `+1` layers are deliberately redundant.  Any `k1` omitted from
(18) satisfies the weaker tail because `ell k1>209/7`, and for a scanned
`k1` any `k2` omitted from (19) satisfies its exact solved form.  Hence
(18)--(19), together with the analytic tail, cover
every pair in (1) with no extrapolation or sampled-residue assumption.

For every row in the bank the program:

1. constructs `S(P)` by subtracting rational teeth from `[0,1]`;
2. subtracts the `k1` and `k2` teeth in exact endpoint order;
3. computes the longest component as a rational endpoint difference; and
4. rejects the row if the exact integer comparison `3k2 L<=1` holds.

All endpoint comparisons, tail cutoffs, and acceptance decisions use signed
`__int128` cross-products.  There is no floating-point arithmetic.  The
complete result is

```text
cores                 66
guarded finite pairs  9,246,070
failures              0
```

The hardest row *inside this finite bank* is

```text
P  = {1,2,3,5,6,7,8,9,10,11},
(k1,k2) = (153,159),
L  = 158/56763,
3k2 L = 158/119,
R := 1/(3k2 L) = 119/158.                               (20)
```

The two longest components in (20) are the reflected pair

```text
[505/2142,177/742], [565/742,1637/2142].                (21)
```

The first interval is literally the gap between the right endpoint of a
`k1=153` tooth and the left endpoint of a `k2=159` tooth.  This is the exact
endpoint-interleaving phenomenon that the old generic `7/18` scaling picture
missed.  We do not call (20) the extremum of the infinite tail: the analytic
argument proves the required strict inequality there but does not optimize
its margin.

Combining Sections 2 and 3 proves (3).

## 4. A computation-light proof of the closure consequence

The strong component bound (3) uses the 9,246,070-row bank.  Uniform `r=3`
closure itself needs only the weaker

```text
L(P;k1,k2) > 1/(7k2),                                  (22)
```

and (22) has a much shorter proof using (10) and the 66-core atlas only.

With `a=k1`, `b=k2`, and `r=b/a`, the mass lower bound (11) and component
cap (12) give (22) whenever

```text
ell a(28r-7) > 6r+29.                                  (23)
```

If `ell a>=5/3`, the difference in (23) is at least

```text
(5/3)(28r-7)-(6r+29)=(122/3)(r-1)>0.                  (24)
```

The exact core atlas gives the stronger fact

```text
ell(P)(13 max(P)+1) >= 1727/1008 > 5/3                (25)
```

for every core except the unique minimum core `P0` in (6).  Thus (24)
handles all those 65 cores immediately.

For `P0`, `ell=1/112`, `max(P0)=12`, and hence `a>=157`.  If `a>=187`,
then `ell a>5/3` and (24) again applies.  If `a<187`, expression (23) is
already positive for `b>=187` (its minimum on this boundary is the exact
value `23/248` at `(a,b)=(186,187)`).  Indeed, for `ell=1/112` it increases
with `b` because `1/4-6/a>0`; at `b=187` it decreases with `a` on
`157<=a<=186` because `-1/16+1122/a^2<0`.  It remains only

```text
157 <= a < b <= 187.                                   (26)
```

For every speed in (26), (8) sharpens to `N_k(J)<=2`, so after two removals
there are at most five components.  The sharp mass (11) is greater than
`5/(7b)` exactly when

```text
5ab-96b-656a>0.                                        (27)
```

The left side is increasing in both variables throughout (26), since its
unit increments are `5a-96>0` and `5b-656>0`; at the smallest pair
`(157,158)` it equals `5870>0`.  Therefore `A>5/(7b)`, one of
at most five components has length greater than `1/(7b)`, and (22) follows
without a pair scan.

Now let `I` be the component from (22), of length `L`, and let `k3>k2`.
The sharp periodic bound (10) gives

```text
|I intersect D_k3| <= L/7+6/(49k3) < L,               (28)
```

because `L>1/(7k2)>1/(7k3)`.  The third comb cannot cover `I`, so a point of
`I` is safe for all thirteen speeds.  Thus the closure consequence of
THM-1094 depends only on the 66-core atlas and elementary inequalities;
the 9,246,070-row certificate proves the stronger structural statement (3)
and independently stress-tests the endpoint carrier.

For comparison, the older coarse estimate

```text
L/7 + 2/(7k3).
```

would require (3).  Bound (28) exposes the natural final-comb target
`1/(7k2)`.  THM-1097 subsequently uses exactly this sharper target to close
the four-killer (`r=4`) clustered stratum uniformly.  The first still-open
iteration is the four-removal bridge for `r=5`.
No covering assumption is used.

## 5. Toothpick self-similarity: what is true

The full two-comb endpoint word has no naive homogeneous self-similarity, but
two exact local sectors do recur.

- If `k2 ell>=13/7`, then `J` contains a complete `k2`-safe gap of length
  `6/(7k2)`.  At most one `k1` tooth can meet that gap: two such teeth would
  force `k2/k1<=1`.  If `k2/k1<4/3`, deleting the one tooth leaves a side
  component longer than `1/(3k2)`.
- If `k1 ell>=13/7` and `k2/k1>=13/6`, a complete `k1`-safe gap contains a
  complete `k2`-safe gap, again longer than `1/(3k2)`.

By (5), `k_i>=208` supplies the relevant `13/7` hypothesis.  These near-ratio
and far-ratio sectors are not needed by the certificate, but they identify
the genuine toothpick recursion: complete safe gaps nest, whereas the middle
ratio endpoint word is the finite irreducible piece.  The finite bank is not
an arbitrary height cutoff; it is a proof-obligation bank left after the
mass/component tail, with these local recurrences explaining much of its
geometry.

## 6. Carrier and Tournament Analysis audit

The predicate-preserving carrier is not a tournament on runners.  It is the
ordered list of core-safe components together with rational boundary
coordinates and endpoint owners.  The two killer combs refine that list by
wall-crossing events.  This carrier preserves exactly whether a positive
length safe component remains and its metric length.

Several candidate vertex sets were challenged: runners, individual teeth,
residue classes, core components, exposed boundary events, and proof
obligations.  On exposed endpoints the natural pair observable is left/right
order, with the endpoint coordinate as gauge and exact equality as the tie.
Its tournament fingerprint is always transitive: score multiset
`{0,1,...,m-1}`, no directed cycles, singleton SCCs, and one Hamiltonian path
(the sorted endpoint word).  An edge flips only at an exact endpoint
coincidence.  Thus the tournament contains no information beyond the order;
without the coordinate and owner sidecars it destroys the component lengths
needed in (3).  A tournament on runners is still coarser and also loses which
teeth meet which core component.

The useful Kakeya-needle picture is therefore one-dimensional and metric:
the periodic teeth are parallel needles piercing a fixed core-safe interval;
(8) controls their incidence count and (9) their total occupied length.
Pairwise orientation alone cannot certify the uncovered interval.  This is
why the exact referee uses the interval/endpoint carrier and records that a
plain Tournament Analysis quotient is not faithful here.

## 7. Reproducibility and independent replay

The C++ output is deterministic across an unoptimized build, an optimized
build, and an AddressSanitizer+UndefinedBehaviorSanitizer build.  The three
runs produce the same frozen output hash:

```text
O0 output SHA-256       10153a80b79da94d179a9f1b0f83402fa41023fdb20922a9de192d3b7d13c550
O2 output SHA-256       10153a80b79da94d179a9f1b0f83402fa41023fdb20922a9de192d3b7d13c550
ASan+UBSan output SHA-256 10153a80b79da94d179a9f1b0f83402fa41023fdb20922a9de192d3b7d13c550
```

The independent `fractions.Fraction` replay verifies the core minimum (5),
the scaled second-core minimum (25), the guarded-bank cardinality, and every
rational in (20)--(21) by a different
core-atlas construction.

```text
C++ source SHA-256      eb59b13c29170f62dab28e783ffeca5bf497fc438e13295d6ed282dfd37951e8
C++ frozen output       10153a80b79da94d179a9f1b0f83402fa41023fdb20922a9de192d3b7d13c550
Python source SHA-256   70b56d8845a50141da182fe84c009a26fbbfbcc35b5d1f2bd71acc0defb6be6d
Python frozen output    43066508d9bf39c148d703912eb1681a6cdef5642d3f400bae07f274d8e2e435
```

Audit commands:

```bash
clang++ -std=c++20 -O0 -pthread \
  04-computation/lrc14_r3_two_comb_component_exact_codex_S73.cpp -o /tmp/thm1094-o0
clang++ -std=c++20 -O2 -pthread \
  04-computation/lrc14_r3_two_comb_component_exact_codex_S73.cpp -o /tmp/thm1094-o2
clang++ -std=c++20 -O1 -g -pthread -fsanitize=address,undefined \
  -fno-omit-frame-pointer \
  04-computation/lrc14_r3_two_comb_component_exact_codex_S73.cpp -o /tmp/thm1094-san
python3 04-computation/lrc14_r3_two_comb_extremal_replay_codex_S73.py
```

## Scope

THM-1094 repairs exactly the missing all-scale implication in THM-1061 and
leaves MISTAKE-163 intact as the audit of the invalid sampled proof.  It does
not close the analogous three-removal bridge for the four-killer stratum and
does not by itself prove global LRC(14).
