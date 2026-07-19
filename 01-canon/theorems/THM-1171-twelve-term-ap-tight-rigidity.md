---
id: THM-1171
title: Twelve-term arithmetic-progressions are tight only on the homogeneous ray
status: PROVED — for C={a,a+d,...,a+11d}, M(C)=1/13 iff a=d.  A nontrivial primitive step has an explicit common-phase witness of clearance at least 1/3; primitive step one has an explicit symmetric-endpoint witness, strict unless the block is {1,...,12}.  This classifies only APs: it neither extracts an AP from an arbitrary tight twelve-set nor collapses the Cover14 crown
source: codex-2026-07-18 incoming-frontier synthesis; independent short proof concurrent with boxeph-S118
depends_on:
  - dilation invariance of M
  - the elementary circle-pigeonhole equality M({1,...,12})=1/13
related: [THM-1149, THM-1150, HYP-4382, HYP-7708, HYP-7710]
script: 04-computation/lrc13_twelve_term_ap_tight_rigidity_codex_20260718.py
output: 05-knowledge/results/lrc13_twelve_term_ap_tight_rigidity_codex_20260718.out
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCAPCentering.lean
  - 04-computation/lean/TournamentH7/TournamentH7/LRCMod13Blocking.lean
---

# THM-1171 — exact tight rigidity inside the twelve-term AP locus

For a finite positive speed set `C`, write

```text
M(C)=max_(t in R/Z) min_(c in C) ||ct||.
```

Let

```text
C(a,d)={a+kd: 0<=k<=11},                 a,d positive integers.       (1)
```

The earlier mod-13 argument proved the necessary congruence `a=d (mod 13)`
and the consecutive calculation closed the primitive-step-one face.  Neither
condition is needed for the remaining spread face: after gcd normalization,
all twelve AP terms have one common phase at a suitable rational time.

> **Theorem (twelve-term AP tight rigidity).**  For positive integers `a,d`,
>
> ```text
> M(C(a,d))=1/13    iff    a=d.                         (2)
> ```
>
> More quantitatively, put
>
> ```text
> g=gcd(a,d),       A=a/g,       D=d/g.                 (3)
> ```
>
> If `D>=2`, then
>
> ```text
> M(C(a,d)) >= floor(D/2)/D >= 1/3.                    (4)
> ```
>
> If `D=1`, then
>
> ```text
> M(C(a,d)) >= A/(2A+11),                              (5)
> ```
>
> and the right side is strictly greater than `1/13` exactly when `A>1`.

## 1. The primitive-step collapse

Multiplying every speed by `g` does not change `M`: the circle map `t -> gt`
is surjective.  It therefore suffices to study

```text
C(A,D)={A+kD:0<=k<=11},        gcd(A,D)=1.              (6)
```

Assume `D>=2` and put `r=floor(D/2)`.  Since `A` is a unit modulo `D`, choose
an integer `j` satisfying

```text
Aj=r mod D.                                             (7)
```

At the rational time `t=j/D`, every AP term has exactly the same phase:

```text
(A+kD)t = Aj/D+kj = r/D mod 1.                          (8)
```

Because `0<r<=D/2`, its circle distance is `r/D`.  For `D=2` this is `1/2`;
for odd `D>=3` it is `(D-1)/(2D)>=1/3`; and for even `D` it is again `1/2`.
This proves (4).  In the original scale the same clock is `t=j/d`.

This is the whole spread argument.  In particular the former `D>=17`
residual is not merely non-tight: every such primitive-step AP is at least
`1/3`-lonely, independently of its offset and of the number of displayed AP
terms.

## 2. Primitive step one

Now let `D=1`.  At

```text
t=1/(2A+11),                                            (9)
```

the twelve phases are

```text
A/(2A+11), (A+1)/(2A+11), ..., (A+11)/(2A+11).         (10)
```

The first and last add to one, so the whole interval (10) is symmetric about
`1/2`.  Every phase therefore has circle distance at least

```text
A/(2A+11).                                              (11)
```

Moreover

```text
A/(2A+11)>1/13  iff  13A>2A+11  iff  A>1.              (12)
```

Thus a tight AP must have `D=1` and `A=1`, which by (3) is exactly `a=d=g`.

Conversely, `a=d=g` gives `C(a,d)=g*{1,...,12}`.  Dilation invariance reduces
to `{1,...,12}`.  The time `1/13` gives clearance `1/13`.  For the reverse
inequality, place the thirteen circle points

```text
0,t,2t,...,12t.                                        (13)
```

in cyclic order.  One of their thirteen cyclic gaps has length at most
`1/13`; its two endpoint indices differ in absolute value by some
`q in {1,...,12}`, and hence `||qt||<=1/13`.  Therefore
`M({1,...,12})<=1/13`, proving equality and (2).

Only the lower bound (5), not an exact formula for shifted consecutive
blocks, is used here.  Thus no unproved upper-bound formula is hidden in this
classification.

## 3. Consequence for the LRC14 crown interface

THM-1149's tight-deletion branch and HYP-4382 had used the open statement

```text
M(W)=1/13  =>  W=c*{1,...,12}.                          (14)
```

Within that route, (14) now factors honestly as

```text
tight twelve-set  =>  twelve-term AP        [OPEN]
tight twelve-term AP  =>  c*{1,...,12}      [THM-1171]. (15)
```

The second arrow, including every shifted and spread AP, is closed.  The first
arrow is still the genuine AP-extraction/equality-classification problem.
Separately, CrownCollapse14 must first produce a tight deletion from an
all-loose Cover14 crown.  THM-1171 does neither of those jobs and therefore
does not by itself prove CrownCollapse14, PrimitiveCarrierINV, or LRC(14).

Boxeph-S118 independently proves the same classification with the stronger
centering witness `q=2a+11d` (and gives nearly one-half clearance away from
the homogeneous ray).  Its kernel-pure `LRCAPCentering` module formalizes the
reusable centered-band arithmetic.  The common-phase proof above is shorter
on the primitive-step `D>=2` branch and serves as an independent proof and
exact replay; neither argument supplies the still-open first arrow in (15).

## 4. Tournament and carrier audit

On the twelve runner vertices, orient `i -> j` when `i<j`, using `d>0` as the
speed-order gauge.  This is the transitive tournament with score histogram

```text
0,1,...,11; directed cycles=0; SCC sizes=1^12;
Hamiltonian paths=1, with no ties.                       (16)
```

It cannot distinguish `D=1` from `D=41` and carries none of the witness in
(8).  Candidate vertices challenged here are runners, adjacent AP gaps, AP
indices, residues modulo `13`, primitive-step residues modulo `D`, and proof
obligations.

The proof-preserving quotient on the spread branch is instead the single
primitive residue class

```text
(A mod D, r=floor(D/2), j=A^(-1)r mod D),               (17)
```

with the explicit clock `j/D`.  It discards the AP index `k` precisely because
`kD*(j/D)=kj` is integral.  Thus it preserves the LRC predicate actually used:
all twelve runners share a phase at distance at least `1/3`.  On the `D=1`
branch, the faithful sidecar is the endpoint-owner pair `k=0,11` at the
symmetric clock (9).  The natural carrier is therefore a common-phase
certificate, not a runner tournament.

## 5. Exact replay

The dependency-free Python replay checks every pair `1<=a,d<=200`, constructs
the modular common-phase clock on every `D>=2` row, checks all twelve phases
with `Fraction`, verifies the consecutive endpoint witness on every `D=1`
row, and confirms that only the 200 homogeneous rows evade a strict
`>1/13` certificate.  It also computes the tournament fingerprint and the
representative spread rows in exact arithmetic.  Ordinary and
`PYTHONOPTIMIZE=1` runs are required to be byte-identical.

The cited Lean modules kernel-check the centered-block witness and the
shifted-consecutive endpoint branch.  The particularly short common-phase
normalization in section 1 remains a paper/Python producer until its gcd
split is exposed as a reusable Lean theorem.

```text
04-computation/lrc13_twelve_term_ap_tight_rigidity_codex_20260718.py
SHA-256 9f97c7910ba56ee4c2c21ea8c0fe37d0f967753525058eb1ccf4253db9735501

05-knowledge/results/lrc13_twelve_term_ap_tight_rigidity_codex_20260718.out
SHA-256 e025fba764e43289bd8288713f1ffb7e3b6d2cee5deb61cd309640308fec3fe9
```
