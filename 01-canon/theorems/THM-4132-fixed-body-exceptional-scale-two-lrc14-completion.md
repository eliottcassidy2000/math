---
id: THM-4132
title: "Fixed-body completion of the exceptional scale-two LRC(14) cell"
status: >
  PROVED ELEMENTARY COMPACT-TO-OPEN COMPLETION + VERIFIED-EXACT. Let
  U=(1,4,6,8,10,12,14,15,16,18,22). For every positive odd integer t, the
  thirteen distinct speeds 2U union {t,9t} have a common phase of clearance
  at least 1/14. This closes the U-body slice of THM-3878's exceptional
  scale-two (2,1,9) cell at every relative scale. It proves no arbitrary-body
  scale-two closure, physical entry, or LRC(14).
source: codex-frontier-synthesis-creative-20260825p
depends_on:
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
  - THM-4129-universal-two-speed-completion-of-the-eleven-speed-lrc14-body
related:
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
  - THM-4002-lrc14-signed-endpoint-cross-phase-and-fixed-scale-two-family
script: 04-computation/lrc14_fixed_body_scale_two_completion_thm4132.py
output: 05-knowledge/results/lrc14_fixed_body_scale_two_completion_thm4132.out
script_sha256: 2a928cc1d622fa3b5a5d351edf3594058b2631cf6c8224b0a203d1f7258ac70f
output_sha256: a4084955d9adbde07e2459757dd3c305268fe84a99c436475f0d9379faa428cd
semantic_sha256: aa6314ca8c25b5ad512bef0874dfb7262bd65cab8b200abba32d1cfe6b7710a6
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone Fraction implementation reconstructs the closed body
  interval and its endpoint owners, derives the two-sheet quotient from its
  four strict danger inequalities, proves the sharp first-scale surplus,
  checks the t=1 clock and a fixed-sheet hostile, and constructs literal
  verified lifts for all 1,000 odd scales from 3 through 2,001. Normal,
  optimized, and hash-seeded replays agree with the frozen output.
---

# THM-4132 -- fixed-body exceptional scale-two completion

**PROVED ELEMENTARY COMPACT-TO-OPEN COMPLETION + VERIFIED-EXACT.**

THM-4129 closes every two-speed completion of one eleven-speed body but
correctly stops before scaling the body separately from two odd tails.
THM-3910 supplies exactly the missing coordinate: the two physical lifts of
one body-safe quotient phase. Joining those results closes the exceptional
fixed-body scale-two cell at every scale.

## 1. Statement

Put

```text
U=(1,4,6,8,10,12,14,15,16,18,22),       delta=1/14.       (1)
```

> **Theorem.** For every positive odd integer `t`, the thirteen distinct
> positive speeds
>
> ```text
> 2U union {t,9t}                                           (2)
> ```
>
> have a phase `x in R/Z` such that
>
> ```text
> min_(v in 2U union {t,9t}) ||vx|| >= 1/14.               (3)
> ```

The row is primitive because `t` is odd, and parity makes the two tails
disjoint from the eleven even body speeds.

## 2. Closed body arc and two physical lifts

THM-4129 proves that

```text
J=[33/70,27/56] subset G_U,             |J|=3/280,          (4)
```

where `G_U` is the closed `delta`-safe set of `U`. Runner `15` uniquely owns
the left equality and runner `4` uniquely owns the right equality.

For `y in J`, both

```text
x_0=y/2,                   x_1=(y+1)/2                     (5)
```

keep every runner of `2U` safe, because `2u x_i` is congruent to `uy`
modulo one. Put `w=ty` and `z=w/2`. Since `t` is odd, the two tail-phase
packets at `(5)` are `(z,9z)` and `(z+1/2,9(z+1/2))`.

Let `D_a={z:||az||<delta}`. Both lifts are tail-bad exactly when

```text
w in C := 2((D_1 union D_9) intersect
             ((D_1 union D_9)-1/2)).                       (6)
```

THM-3910 computes, and the primary audit reconstructs directly from the
strict inequalities, that

```text
C=(2/21,8/63) union (55/63,19/21),       beta(C)=2/63.      (7)
```

The components in `(7)` are open. This strictness is load-bearing.

## 3. Compact-to-open closure

Suppose first that `t>=3`. If `(2)` were a counterexample, every `y in J`
would make both lifts in `(5)` tail-bad. Hence the image `tJ` would be a
compact connected subset of the open set `C`. If it wraps around the circle,
this is impossible immediately. Otherwise it must fit strictly inside one
component, forcing

```text
t|J| < beta(C).                                             (8)
```

Already at the first admissible scale,

```text
3|J|-beta(C)=9/280-2/63=1/2520>0.                          (9)
```

Thus `(8)` fails for every odd `t>=3`. Equality would also close: a compact
arc cannot lie inside an open arc of equal length.

The only smaller positive odd scale is `t=1`. Here

```text
x=1/13,             min_(v in 2U union {1,9})||vx||=1/13, (10)
```

which is strictly stronger than `(3)`. This proves the theorem. **QED.**

## 4. Hostile and information ledger

No fixed physical sheet proves the result. At

```text
w=1/9,                   z=w/2=1/18,                      (11)
```

the first lift is killed by tail `1`, with gap `1/18`, while the second lift
is killed by tail `9`, with gap zero. Thus the interval and the choice between
the two lifts are genuine proof coordinates.

The connection contract is

```text
source:       THM-4129's closed U-safe interval J
target:       the separately scaled row 2U union {t,9t}
map:          y -> (y/2,(y+1)/2), followed by w=ty
preserved:    all eleven body clearances and one common physical time
destroyed:    which physical sheet makes the two tails safe
sidecar:      THM-3910's open two-sheet quotient obstruction C
decisive test: whether the compact arc tJ fits in one component of C.       (12)
```

THM-4129 could not be applied by common dilation: halving the body in `(2)`
makes the odd tails nonintegral. The new result closes precisely that missing
fixed-body cell rather than repairing arbitrary bodies.

## 5. Scope and replay

This theorem closes the `U`-body slice of the exceptional
`(s,p,q)=(2,1,9)` type, including the formerly conditional `t>=22` region and
all lower relative scales. It does not close an arbitrary eleven-speed body,
the general scale-two certificate type, physical entry into the `11+2`
branch, even-`t` bookkeeping outside the primitive scale-two normalization,
or LRC(14).

The exact audit additionally constructs and directly checks a literal safe
lift for every odd `3<=t<=2001`; this is a control, not the source of the
all-`t` proof. Run

```text
python3 -B 04-computation/lrc14_fixed_body_scale_two_completion_thm4132.py
python3 -B -O 04-computation/lrc14_fixed_body_scale_two_completion_thm4132.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc14_fixed_body_scale_two_completion_thm4132.py
```

All three streams agree with the frozen output. **QED.**
