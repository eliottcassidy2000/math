---
id: THM-2707
title: "Full physical lift fibre, common simplex, and fixed-packet SCC"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  Above the THM-2694 terminal private vertex, all 12*13^5 nonzero
  THM-2657 physical lifts have been scanned.  Exactly 3,042 retain the full
  fixed six-label packet skeleton on the entire inherited open cylinder.
  The complete orbit has 3,346 such packet addresses in eleven residue
  parts; their restricted supports have a common Delta^3345 nerve and their
  lift graph is the complete directed eleven-partite graph, one SCC with
  10,177,910 edges.  Explicit three- and eleven-state physical support cycles
  exist.  This is fixed-skeleton support composition, not a target action,
  frozen following atom, semantic endpoint current, row exclusion, or
  LRC(14) conclusion.
source: codex/thm2694-full-lift-fibre-scout-2026-07-28
depends_on:
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2694-mixed-dilation-slope-seven-present-unit-long-word-and-first-gap
related:
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2687-slope-seven-global-configuration-switching-positive-thirteenfold-no-go
script: 04-computation/lrc14_full_physical_lift_fibre_thm2707.py
output: 05-knowledge/results/lrc14_full_physical_lift_fibre_thm2707.out
script_sha256: f05a07b2fb22cb5b39ed7d14e66d26154ecc50fc214861dc6576c3bcfaed2412
output_sha256: c5cf7eaef9393c2c551bb2bc8d1c01ff40f0bea763f473c6da11a9fb03d05173
hash_basis: LF-normalized bytes
---

# THM-2707 -- full physical lift fibre, common simplex, and fixed-packet SCC

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2694 proved that the twelve least positive representatives of the next
slope step all leave its present factor.  It explicitly left open the other
`13^5-1` representatives in each THM-2657 lift fibre.  They behave in the
opposite way: thousands return to the same positive packet cylinder, and
some compose into literal closed support cycles.

## 1. Inherited cylinder and packet skeleton

Retain the exact THM-2694 data

```text
p=13,                       R=p^6=4826809,
x=649039434905733/1304692766858936,
z={13x}=46873542509301/100360982066072,

I=(960117507257/1930018885886,
   324519717452867/652346383429468),
|I|=1/652346383429468.                                  (1)
```

For `n in Z/RZ`, write

```text
q_n(x')={13x'}+7n/R.                                    (2)
```

The six labels

```text
(rail,sector,edge,kappa,h,shallow)=(0,0,0,1,6,1)        (3)
```

are fixed.  The predecessor carry and private root are **not** frozen; they
covary as

```text
carry(n)=2+7n,                 root(n)=6+n              (mod 13). (4)
```

The retained typed packet is exactly the one in THM-2694: literal rail,
present factor, delayed Boolean prefix, predecessor carry cell, future digit
and half-digit, private deep half-tooth, and primitive coefficient-unit flag.
It is important that `(3)` is a fixed skeleton rather than an unqualified
fixed eight-label configuration.

## 2. Complete packet census and one common cylinder

Scan all `n in Z/RZ`.  Exactly `3346` addresses retain the packet at the
midpoint.  Their residue parts have sizes

```text
r:  0   1   2   3   4   5   6   9  10  11  12
N: 304 305 304 305 304 305 304 301 304 305 305.         (5)
```

The active residue set is therefore

```text
A={0,1,2,3,4,5,6,9,10,11,12}.                           (6)
```

This is stronger than a midpoint census.  For every one of the `3346`
addresses, every factor listed after `(4)` contains the **whole** open
cylinder `I`.  The exact residual wall slacks after subtracting the radius
of `I`, in the natural rail/present grid coordinate, have positive minima

```text
rail left      82675509840/2197,
rail right     12211300179540/371293,
present left   10051258140/2197,
present right  9961158781035/371293.                    (7)
```

The delayed coordinate is independent of `n`, the carry fractional part
does not cross an integer wall on `I`, and the private tooth translates with
`root(n)`.  The companion separately checks the future digit `h=6` and
half-digit `kappa=1` at both endpoints.  Thus, after restriction to the
inherited D component, all `3346` packet supports equal `I`.  Their common-
base support nerve is the full simplex

```text
Delta^3345.                                               (8)
```

This is the common higher simplex left missing by THM-2657, but only for the
fixed packet skeleton `(3)`.

## 3. Every physical lift from the terminal vertex

The THM-2694 terminal private address is

```text
n_0=110,                    n_0 mod13=6.                  (9)
```

A THM-2657 nonzero physical lift has numerator `k in Z/RZ` with `13` not
dividing `k`.  Since

```text
7^(-1)=4137265 mod R,                                      (10)
```

the lift sends `(9)` to

```text
n'=n_0+7^(-1)k mod R.                                    (11)
```

Its root increment is

```text
delta=n'-n_0=2k mod13,
k=7delta mod13.                                          (12)
```

Consequently `(11)` bijects all

```text
12*13^5=4455516                                           (13)
```

nonzero lifts with the twelve target residue cosets different from `6`.
Exactly `3042` target addresses survive the whole typed packet on `I`.
Their counts by nonzero root increment are

```text
delta:  3   4   5   6   7   8   9  10  11  12
count: 301 304 305 305 304 305 304 305 304 305.          (14)
```

There are no survivors for `delta=1,2`, because those target residues are
the two private root/unit holes `7,8`.  The shortest positive numerator is

```text
k=2,   delta=4,   n'=3447831.                            (15)
```

The signed numerator `k=-1` also survives.  The stored `3062`-line transcript
records, for every one of the `3042` survivors, its `k`, `delta`, target
address, residue, carry, root, common interval, and exact number of surviving
next lifts.

The old boundary remains sharp: the twelve least positive section values
`k=7delta`, `1<=delta<=12`, land at `n=111,...,122`, and every one is still
present-free on all of `I`.  They were a section of the lift torsor, not an
exhaustion of it.

The source packet together with the `3042` nonzero-lift target packets spans
a common-intersection cone `Delta^3042`.  This is not called the simplicial
star of `(9)`, because the full support nerve `(8)` also contains same-residue
vertices which are not joined to `(9)` by nonzero quotient lifts.

## 4. Exact lift graph and closed support words

Install a directed edge from packet address `a` to packet address `b` when
the unique translation numerator

```text
k=7(b-a) mod R                                             (16)
```

is a nonzero THM-2657 lift.  Equation `(16)` is admissible exactly when
`a mod13 != b mod13`.  Hence the graph on the `3346` packet addresses is the
complete directed eleven-partite graph with part sizes `(5)`.  It has

```text
3346^2-sum_r N_r^2=10177910                               (17)
```

directed edges, one strongly connected component, no terminal vertex, and
directed diameter at most two.  From a vertex in residue part `r`, the exact
number of next full-packet lifts is `3346-N_r`, hence one of

```text
3041, 3042, 3045.                                        (18)
```

In particular, every one of the `3042` terminal survivors continues.

The smallest explicit nonbacktracking closed word uses signed lift
numerators

```text
(-1,-1,2)                                                 (19)
```

and visits addresses/residues

```text
110 -> 689654 -> 1379198 -> 110,
  6 ->      4 ->       2 ->   6.                         (20)
```

There is also a closed word visiting all eleven active residue parts.  Its
cumulative lift numerators and successive steps are

```text
K=(0,2,4,5,8,10,11,12,9,6,3,0),
dK=(2,2,1,3,2,1,1,-3,-3,-3,-3).                         (21)
```

Every vertex of `(20)` and `(21)` carries the same open cylinder `I`, and
the cumulative numerator returns exactly to zero.  Thus either word can be
repeated indefinitely as a fixed-skeleton **support** word.

## 5. Proof mechanism

Multiplying `(2)` by `R` kills the `7n/R` displacement modulo one, so the
delayed Boolean coordinate and its future digit/half-digit are independent
of `n`.  The integer part changes by `7n`, giving the carry in `(4)`.  The
high-speed probe changes by `2n/13`, giving the private root in `(4)`.  This
proves the covariant part of the packet algebraically.

For the two low-speed factors, scale the exact rational rail and present
endpoints by the common denominator of `Tz`.  Advancing `n` then adds one
integer step on a finite cyclic grid of length `R`.  A complete scan gives
`(5)`.  For each retained address, pull the containing rail, present, delayed,
carry, and private-tooth inequalities back through `x' -> q_n(x')`; their
intersection contains `(1)`, with the explicit slacks `(7)`.  This proves
the common-cylinder claim without sampling.

Finally, multiplication by seven is invertible modulo `R`, proving the
affine lift bijection `(11)`.  Reducing `(16)` modulo thirteen proves the
complete multipartite graph description.  Equations `(19)`--`(21)` are then
direct exact substitutions, not quotient-only cycles.

## 6. Sharp scope and surviving obstruction

This theorem repairs one assumption from THM-2694: checking the twelve least
positive lifts cannot establish maximality in the full physical lift torsor.
It does **not** contradict that theorem's stated fixed-section maximality.

Nor does a common support SCC split the cyclic extension of THM-2657.  The
edges in `(19)` and `(21)` are tailored representatives; they do not define
a homomorphic `C_13` section.  More importantly, arbitrary `k/R` shifts move
the low-speed factors by fine phases rather than by a lawful global target
action.  Positive overlap therefore supplies composition in the support
nerve, but not the signed/complex endpoint-current covariance needed by the
LRC proof.

The frozen weighted `following` THM-2680 atom is deliberately excluded.  An
exact hostile check finds it absent already at terminal `n=110` and at ten
of the eleven distinct vertices in `(21)`; only the cumulative `K=10` vertex
hits it at the midpoint.  Thus “full typed packet” here means precisely the
seven retained factors listed in Section 1, not the old whole following atom.

No configuration-switch theorem, canonical target action, semantic source
or endpoint current, all-row statement, row exclusion, or proof of LRC(14)
follows.  The next sharp problem is to decide whether one closed support word
admits a common endpoint/current phase sidecar, or whether the nonsplit
odometer class forces cancellation around every such word.

## 7. Exact reproduction

Run

```bash
python3 04-computation/lrc14_full_physical_lift_fibre_thm2707.py
python3 -O 04-computation/lrc14_full_physical_lift_fibre_thm2707.py
```

Both executions byte-match
`05-knowledge/results/lrc14_full_physical_lift_fibre_thm2707.out`.  The
companion contains no optimized-away assertions.  It scans all `13^6`
addresses, thereby enumerating all `12*13^5` nonzero quotient lifts; verifies
whole-open-cylinder containment at all `3346` packet addresses; checks the
complete multipartite edge invoice, the canonical-section hostile, and both
closed cycles; and prints all `3042` survivor rows.
