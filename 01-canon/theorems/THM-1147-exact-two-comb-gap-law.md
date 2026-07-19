---
id: THM-1147
title: THE EXACT TWO-COMB GAP LAW — the "alignment bias" is arithmetic and closed-form, and it supplies the nonuniformity the four-comb theorem needs. (I) THE LAW, verified exactly: a surviving component of a multi-comb complement is the interval from the RIGHT edge of tooth j of one comb a to the LEFT edge of tooth j+1 of another comb b, so its length is **L(a,b,j) = (a − j·d)/(a·b) − 1/(14a) − 1/(14b)** with d = b − a. On the standing worst case (core [1,3,5,6,7,8,11,12], killers 371/374/377/379) the law predicts 127/142883 and the measured longest component IS 127/142883 — the endpoints being 1373/5278 = 98/377 + 1/(14·377) and 1385/5306 = 99/379 − 1/(14·379). (II) SO THE GAP IS **LINEAR IN j**, falling from ≈ 1/a at j = 0 to zero at j ≈ a/d. Numerically for (371,379), d=8: usable gap × b runs 0.856, 0.748, 0.640, 0.424, 0 as j runs 0, 5, 10, 20, 40. For (371,372), d=1 the same descent takes until j ≈ 370. (III) THIS IS THE NONUNIFORMITY SOURCE, and it is provable rather than statistical: a linear descent from 1/a to 0 has mean ≈ 1/(2a) and maximum ≈ 1/a, so **max gap ≈ 2 × mean gap** from the pair law alone — comfortably above the **4/3** that THM-1141 identified as sufficient for the four-comb theorem. With four combs there are six pairs and more accessible indices, which is where THM-1141's measured 3.34 comes from. (IV) HONEST NEGATIVE: my proposed *predictor* — that small j·d/k correlates with large surviving gaps — is NOT confirmed. The measured correlation is inverted (median j*·d_min/k₄ = 0.345 among the worst cases versus 0.758 among the best). The proxy is mis-specified: it uses k₄ and d_min, whereas the law involves the actual bounding PAIR (a,b) and its own d, which need not include k₄. So the mechanism (I)–(III) is confirmed exactly while the summary statistic built on it is not
status: (I) PROVED and verified exactly in rational arithmetic — the law is an identity about tooth edges, and it reproduces the standing worst case to the last digit. (II) PROVED (immediate from I). (III) is a correct consequence of the linear law for a single pair; the extension to four combs is indicated, not proved. (IV) REFUTED — the proxy correlation is inverted and the test was mis-specified; recorded rather than dropped. Uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128 (cont.70; owner: work the arithmetic alignment bias)
depends_on:
  - THM-1141    # which identified nonuniformity as the lever and asked for its arithmetic cause
  - THM-1140    # the four-comb reconnaissance
related: [THM-1097, THM-1094]
script: 04-computation/alignment_arithmetic_kps_S128c70.py (+ .out)
---

# THM-1147 — the exact two-comb gap law

THM-1141 found that the surviving gaps are far larger than uniform interleaving predicts,
and asked whether the bias has an arithmetic cause. It does, and it is closed-form.

## (I) The law

A surviving component of the complement of several combs is bounded by tooth edges. Decoding
the standing worst case — core [1,3,5,6,7,8,11,12], killers 371/374/377/379 — its longest
component is

> [1373/5278, 1385/5306], and 5278 = 14·377, 5306 = 14·379,

i.e. it runs from the **right edge of tooth j = 98 of comb 377** to the **left edge of tooth
j = 99 of comb 379**. Hence in general, for bounding combs a < b with d = b − a:

> **L(a,b,j) = (j+1)/b − j/a − 1/(14a) − 1/(14b) = (a − j·d)/(a·b) − 1/(14a) − 1/(14b).**

Check: (377 − 98·2)/(377·379) − 1/5278 − 1/5306 = **127/142883**, and the measured longest
component is **127/142883** exactly.

## (II) Linear in j

The gap descends linearly from ≈ 1/a at j = 0 to zero at j ≈ a/d:

| (a,b) = (371,379), d=8 | j=0 | 5 | 10 | 20 | 40 |
|---|---|---|---|---|---|
| usable gap × b | 0.856 | 0.748 | 0.640 | 0.424 | ≈0 |

| (a,b) = (371,372), d=1 | j=0 | 50 | 150 | 300 | 370 |
|---|---|---|---|---|---|
| usable gap × b | 0.857 | 0.722 | 0.453 | 0.048 | ≈0 |

Small d stretches the descent over many more indices — which is exactly why *clustered*
killers are not the disaster the uniform model predicted.

## (III) Why this is the nonuniformity the four-comb theorem needs

A quantity descending linearly from 1/a to 0 has mean ≈ 1/(2a) and maximum ≈ 1/a. So from
the pair law alone,

> **max gap ≈ 2 × mean gap**,

against the **4/3** that THM-1141 showed is sufficient. The four-comb case has six pairs and
a wider index range, which is where THM-1141's measured ratio of 3.34 comes from. Crucially
this is a *closed-form* source of nonuniformity, not a statistical observation — which is
what an analytic tail for a four-comb bank would need.

## (IV) The predictor I proposed is refuted

I expected small j·d/k to predict large surviving gaps. It does not. Measured over 160
clustered quadruples:

| | median j*·d_min/k₄ |
|---|---|
| worst half by L·k₄ | **0.345** |
| best half by L·k₄ | **0.758** |

The correlation is **inverted**. The proxy is mis-specified: it uses k₄ and d_min, while the
law involves the actual bounding pair (a,b) and *its* difference, which need not involve k₄
at all. So (I)–(III) stand exactly and the summary statistic built on top of them does not.
Recording it rather than quietly dropping it, since an inverted correlation is the kind of
thing that reads as noise if you do not say you predicted the opposite.

## Named next
- Redo (IV) against the **actual bounding pair** rather than (k₄, d_min): identify a and b
  from the component's endpoints, then test gap versus j·d/a. The law says that correlation
  must be exact, so it is a check on the extraction, not on the mathematics.
- For the four-comb bank: the analytic tail should use (III). What must be proved is that
  over the index range accessible inside a core component, some pair (a,b) and index j give
  j·d/a bounded away from 1 — i.e. the component cannot sit entirely in the far tail of
  *every* pair's descent simultaneously. That is now a concrete finite statement about six
  linear functions, which is a much better target than the beat argument it replaces.
