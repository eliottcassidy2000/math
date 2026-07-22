---
id: THM-2072
title: "Fixed owner-clock banks are not uniform, and an antipodal safe pair closes a dyadic seam"
status: >
  PROVED. No finite THM-2066 owner-word clock bank fixed independently of
  the quotient core can close every primitive divisor-complete eleven-core:
  a single core speed divisible by every bank modulus makes all sampled safe
  packets empty and every owner constraint vacuous. Independently, if the
  closed weak-safe core contains a phase together with its half-translate,
  then no odd tail is dangerous everywhere, so the strict dyadic two-tail
  seam is impossible. This is a certificate-strategy theorem, not LRC(14).
source: codex-2026-07-21-LRC-fixed-bank-audit
depends_on:
  - THM-2061
  - THM-2066
related:
  - THM-2060
  - THM-2068
  - THM-645
---

# THM-2072 -- fixed-bank no-go and antipodal safe-pair certificate

Put `delta=1/14`. For a finite set `C` of positive integers, write

```text
G_C={theta in R/Z: ||c theta||>=delta for every c in C}.  (1)
```

For a positive clock `N`, retain THM-2066's weak-safe packet and residue-pair
relation

```text
A_N(C)={r mod N:14|cr|_N>=N for every c in C},
R_N(C)={(u,v) in E_N(C)^2:
          omega_v(r)=1-omega_u(r) for every r in A_N(C)}. (2)
```

Here `E_N(C)` is the set of eligible odd classes modulo `2N`, and `omega` is
the parity word on the labelled packet `A_N(C)`.

## 1. A fixed finite owner-clock bank cannot be uniform

Let `F` be any nonempty finite set of positive clocks. There is a primitive,
divisor-complete eleven-core `C_F` for which

```text
A_N(C_F)=empty,       R_N(C_F)=O_N^2       for every N in F, (3)
```

where `O_N` is the set of odd residue classes modulo `2N`. Consequently the
generalized-CRT bank relation `R_F(C_F)` of THM-2066 is nonempty. In
particular, no finite bank chosen independently of `C` can prove
`R_F(C)=empty` for every primitive eleven-core satisfying the divisor pins
through `14`.

### Proof

Set

```text
B=lcm(F union {2,3,...,14}),       C_F={1,2,...,10,B}.     (4)
```

Since `B` is a multiple of `lcm(2,...,14)=360360`, the entries in `C_F` are
eleven distinct positive integers. The core is primitive because it contains
`1`, and it is divisor-complete through `14` because every integer from `2`
through `14` divides `B`.

Fix `N in F`. For every residue `r mod N`, the core entry `B` gives

```text
|Br|_N=0,
```

so the inequality `14|Br|_N>=N` in (2) fails. Thus `A_N(C_F)` is empty.
Eligibility is then a universal condition over the empty packet, so every odd
class is eligible. The owner word has empty domain, and the complementary-bit
condition in (2) is also vacuous. Hence `R_N(C_F)=O_N^2`.

Let `L_F=lcm_(N in F)(2N)`. The two actual odd integers `1` and `3` give a
compatible residue pair modulo `L_F`; their reduction lies in `O_N^2` for
every `N in F`. Therefore `(1,3) mod L_F` belongs to `R_F(C_F)`, proving that
the bank relation is nonempty. This residue witness is not asserted to be an
LRC counterexample: it proves that this fixed certificate bank is blind. QED.

The same argument gives the more general diagnostic:

```text
if C contains a common multiple of all N in F,
then every A_N(C) is empty and R_F(C) contains all globally odd pairs. (5)
```

This does not conflict with THM-2066 or the finite compression reserved in
THM-2068, whose cores satisfy `max(C)<=24`. Formula (4) has
`max(C)>=360360`. It rules out only an extension by one bank fixed once and
for all; clocks chosen adaptively from the actual core remain available.

## 2. An antipodal safe pair closes the seam

Suppose

```text
there exists theta with theta in G_C and theta+1/2 in G_C. (6)
```

Then for no distinct positive odd `x,y` can

```text
S=2C union {x,y}
```

be a strict LRC(14) counterexample.

### Proof

For every odd integer `z` and every real `theta`,

```text
||z(theta+1/2)||=||z theta+1/2||=1/2-||z theta||.         (7)
```

Thus the two distances in (7) sum to `1/2`. They cannot both be strictly
less than `1/7`; in fact they cannot both be at most `1/7` either.

By THM-2061, a strict counterexample on the dyadic seam would require each
odd tail `z in {x,y}` to satisfy

```text
||z phi||<1/7       for every phi in G_C.                 (8)
```

Applying (8) at the two phases in (6) contradicts (7), already for one tail.
Hence the seam is impossible. QED.

The endpoint convention is load-bearing in its benign direction: `G_C` is
the **closed** weak-safe core from THM-2061, while tail danger in a strict
counterexample is **open**. The argument is stronger than needed because it
also excludes simultaneous weak tail inequalities at the two phases.

Condition (6) has the exact one-phase band form

```text
||c theta||>=1/14                  for every even c in C,
1/14<=||c theta||<=3/7             for every odd  c in C. (9)
```

Indeed a half-shift fixes every even phase, whereas for odd `c` it replaces
the distance `d` by `1/2-d`. Thus (9) is an adaptive continuous certificate
whose feasibility can be attacked without choosing tails or owner words.

## 3. What survives and what must change

The source object in THM-2066 is the continuous safe set `G_C` with its two
dyadic lifts. The fixed-bank map samples it on the rational packets
`A_N(C)` and retains exact tail ownership there. It preserves the lift killed
at every sampled phase, but it destroys the entire carrier when a core speed
is zero modulo the sampling clock. Empty owner words then turn both
eligibility and complementarity into vacuous truths. Any adaptive-clock
strategy therefore needs, at minimum, the side condition
`A_N(C) != empty`; merely adding finitely many more predetermined clocks
cannot repair the loss. The cheapest hostile test for a proposed bank is the
common-multiple core (4).

Two precise replacement targets remain open:

1. **Adaptive owner-clock lemma.** For every primitive divisor-complete
   eleven-core `C`, find a clock `N=N(C)` with `A_N(C)` nonempty and
   `R_N(C)=empty` (or an adaptive finite bank with empty CRT join).
2. **Antipodal safe-pair lemma.** Determine which primitive
   divisor-complete eleven-cores satisfy (9), and close those that do by
   Section 2 before any tail search.

These targets are complementary: the first retains labelled rational
ownership, while the second bypasses ownership through a two-point geometric
obstruction.

## 4. Assumption challenge and tournament relevance

Runners are not faithful vertices here. Owner words are faithful only over a
nonempty packet; on the cores (4), every residue class carries the same empty
word, so a tournament on tails or residue classes would orient pure vacuity.
For Section 2 the useful vertices are the antipodal phase pairs
`{theta,theta+1/2}`. Their relation to an odd tail is a symmetric exclusion:
the tail cannot be dangerous at both endpoints. There is no intrinsic
orientation and hence no tournament theorem to exploit. The faithful objects
are the labelled complement graph plus the nonemptiness sidecar in the
rational lane, and the antipodal-pair incidence system in the continuous
lane.
