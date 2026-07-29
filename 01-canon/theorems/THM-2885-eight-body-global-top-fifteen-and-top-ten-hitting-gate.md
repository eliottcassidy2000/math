---
id: THM-2885
title: Eight-body global top-fifteen seal and top-ten hitting gate
status: >
  PROVED + FINITE-EXACT + VERIFIED.
  For every one of the 3003 eight-speed bodies in {1,...,14}, the exact
  global top fifteen external tooth coverages are tail-sealed and their
  ranks 11 through 15 have sum strictly below the body's good-set measure.
  Hence every five-speed external set whose scalar coverages can exhaust
  the body meets its selected global top ten.  This is a hitting reduction,
  not a ranked-apex closure or a proof of the j=5 rung.
source: root/lrc-rank-impossible-overlap-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
related:
  - THM-2883-ranked-apex-hitting-closure-all-thm741-roots
  - THM-741-near-AP-four-slot-closure-all-2002-bodies-in-1-14
verification:
  - 04-computation/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.py
  - 05-knowledge/results/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.out
---

# THM-2885 — eight-body global top-fifteen seal and top-ten hitting gate

## 1. Statement

Let

```text
E in C({1,...,14},8),
G_E = {t in R/Z : every e in E is at distance at least 1/14 from 0},
m_E = |G_E|,
r_E = number of interval components of G_E,
D_w = {t in R/Z : wt is at distance less than 1/14 from 0},
c_E(w) = |G_E intersect D_w|                         (1)
```

for integer `w>=15`.  Order all pairs `(c_E(w),w)` by decreasing coverage,
breaking a tie by increasing speed, and write their coverage coordinates as

```text
q_1(E)>=q_2(E)>=... .                                (2)
```

For every one of the `C(14,8)=3003` bodies, the global first fifteen terms
exist as an exactly certified finite list and satisfy

```text
q_11(E)+q_12(E)+q_13(E)+q_14(E)+q_15(E) < m_E.       (3)
```

Let `H_E` be the ten speeds at ranks `1,...,10`.  If `Q` is any set of five
distinct integers at least `15` and

```text
sum_{w in Q} c_E(w) >= m_E,                           (4)
```

then

```text
Q intersect H_E is nonempty.                          (5)
```

Equivalently, if `Q` avoids `H_E`, the union bound and `(3)` give

```text
|G_E intersect union_{w in Q} D_w|
 <= sum_{w in Q} c_E(w)
 <= q_11(E)+...+q_15(E)
 < m_E,                                               (6)
```

so `E union Q` has positive lonely-time measure.  Thus only extensions
meeting one of ten body-dependent apex speeds can remain after this gate.

## 2. Why the top fifteen are global

The verifier reconstructs `G_E` exactly and scans every integer speed from
`15` through `600`.  Let `q_15^(600)` be rank fifteen in that finite scan.
Every body satisfies

```text
q_15^(600)>m_E/7.                                     (7)
```

The strict discrepancy consequence of THM-735 (via its THM-731/732
dependency) is

```text
c_E(w)<m_E/7+(99/70)r_E/(7w).                         (8)
```

Set

```text
T_E=(99/70)r_E/[7(q_15^(600)-m_E/7)].                 (9)
```

The verifier also scans `601<=w<ceil(T_E)`.  For every unscanned
`w>=ceil(T_E)`, `(8)` and `(9)` imply

```text
c_E(w)<m_E/7+(99/70)r_E/(7w)<=q_15^(600).             (10)
```

The strict first inequality is essential when `T_E` is an integer: equality
in the rational majorant still cannot admit a new rank-fifteen value.
Consequently sorting the finite scanned ledger produces the true global top
fifteen in `(2)`.

The largest threshold is

```text
T_E=1002456/953,  ceil(T_E)=1052
```

at

```text
E=(2,8,9,10,11,12,13,14),
```

and the largest speed appearing in any retained global top fifteen is only
`243`.

## 3. Exact census and hostile boundary

For comparison, the raw scalar five-cover gate

```text
m_E-(q_1+...+q_5)>0                                  (11)
```

holds on only `16` of the `3003` bodies and is nonpositive on the other
`2987`.  The top-ten hitting reformulation therefore retains information
that the raw five largest coverages discard.

In contrast, `(3)` holds on all `3003` bodies.  Its smallest margin is

```text
m_E-(q_11+...+q_15)=61/12936                          (12)
```

at

```text
E=(1,3,5,7,8,9,11,13),
```

where the rank-`11,...,15` speeds, in decreasing-coverage order, are

```text
36, 42, 84, 25, 72.                                  (13)
```

The positivity in `(12)` is strict, so `(6)` also handles every equality
case in the necessary scalar condition `(4)`.

## 4. Verification audit

The exact verifier is

```text
04-computation/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.py
SHA-256 dff97f67b1104c25589802a6a2f216b6e7bfedd58eebfa1bcce615d59c1e872f.
```

Its stored output is

```text
05-knowledge/results/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.out
SHA-256 21a89f15fb144c406936ff62eaf039c0643e36d82ed99a9d28495181fa13e402.
```

The verifier hash-pins the inherited exact interval kernel, checks the full
body universe and ordering, and uses rational arithmetic throughout.  Its
batched integer primitive has explicit signed-`int64` overflow guards.
For every body, both speed `15` and the last scanned boundary speed are
independently recomputed with a scalar rational primitive: `6006`
vector/scalar controls in total.

The canonical serialization of all `3003` body rows, including `m_E`, `r_E`,
`T_E`, the tail boundary, all fifteen speed/coverage pairs, and the margin
in `(3)`, has digest

```text
f1c478fc386171e4bfbdab52100808d7a37853ae755a49c9e45c6412708bb20a. (14)
```

Ordinary `python3` and optimized `python3 -O` replays are byte-for-byte
identical and both match the stored output.

## 5. Exact scope

This theorem proves a global finite hitting gate.  It does **not** prove that
the ten apex branches close, does not certify a proposed recursive carrier
scheme, and does not prove the eight-body/five-slot rung or LRC(14).  Those
are separate obligations; treating `(5)` as an equivalence would lose the
overlap geometry among the five dangerous combs.
