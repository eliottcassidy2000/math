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
  seam is impossible. An exact quarter-anchor recursion lifts a safe phase
  of the rescaled multiple-of-four layer, and proves the certificate whenever
  that layer forms a controlled fan. This is a certificate-strategy theorem,
  not LRC(14).
source: codex-2026-07-21-LRC-fixed-bank-audit
depends_on:
  - THM-2061
  - THM-2066
related:
  - THM-2060
  - THM-2068
  - THM-2073
  - THM-2075
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

Moreover, the blind core (4) is not geometrically hard: it already has an
explicit antipodal safe pair. Namely, because `56|B`, put

```text
theta_B=15/56+1/(14B).                                  (4a)
```

At `15/56`, the distances for speeds `1,...,10`, in units of `1/56`, are

```text
15,26,11,4,19,22,7,8,23,18.                            (4b)
```

The perturbation in (4a) changes each of these distances by at most
`10/(14B)<1/56`. The sole lower-bound equality in (4b), at speed `4`, moves
in the safe direction. Every other distance stays above `4/56=1/14`, and
the odd-speed distances stay below `24/56=3/7` (the closest upper case is
speed `9`, at `23/56`). Finally

```text
||B theta_B||=1/14.
```

Thus the band criterion (9) below holds for `C_F`, so Section 2 closes this
core continuously. The example therefore separates **sensor blindness**
from mathematical hardness: the fixed rational packets vanish even though
a two-point certificate is explicit.

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

### Exact quarter-anchor toothpick recursion

Assume `C` has at least one member divisible by `4`, and form its rescaled
four-layer

```text
D={c/4:c in C and 4|c}.
```

Let `t in [0,1]` be a chosen real representative satisfying

```text
t in G_D,
ct<=5/7   for every odd c in C,
ct<=12/7  for every c in C with c=2 mod 4.               (10)
```

Then

```text
theta=1/4+t/4                                           (11)
```

satisfies (9), so `{theta,theta+1/2}` is an antipodal safe pair and the
strict dyadic seam over `C` is impossible.

Indeed, if `c=4d`, then

```text
||c theta||=||d+dt||=||dt||>=1/14
```

by `t in G_D`. If `c` is odd, its phase at `1/4` is a quarter or
three-quarters. Its nonnegative displacement `ct/4` is at most `5/28`, so
the resulting distance lies in the closed band `[1/14,3/7]`. If
`c=2 mod 4`, its phase starts at one-half and its displacement is at most
`3/7`; hence its distance is `1/2-ct/4>=1/14`. These residue cases exhaust
`C`, proving (9). All equalities are allowed, so the closed-core endpoint
convention is preserved. QED.

This is a literal recursive toothpick: the `4`-divisible layer is replaced
by the smaller speed set `D`, while the other two residue layers become
linear drift budgets on the selected phase `t`. The recursion retains the
actual representative and the products `ct`; reducing `t` modulo a coarse
clock would lose the bounds in (10).

### Depth-two H-drift corollary

Suppose

```text
C=4Q union {h_0,2h_1},       h_0,h_1 odd,               (12)
```

and `G_Q` is nonempty. By symmetry it meets `[0,1/2]`; define its first safe
time

```text
tau(Q)=min(G_Q intersect [0,1/2])>0.                    (13)
```

Then the seam over `C` is impossible whenever

```text
tau(Q)<=min(5/(7h_0),6/(7h_1)).                         (14)
```

Indeed the four-layer in the recursion is exactly `D=Q`. Taking
`t=tau(Q)`, the two conditions in (14) become

```text
h_0t<=5/7,          (2h_1)t<=12/7,
```

which are precisely (10). Equivalently, define the dimensionless drift

```text
H_Q(h_0,h_1)=tau(Q) max(7h_0/5,7h_1/6).                 (15)
```

The certificate is `H_Q<=1`, including equality. Therefore any hypothetical
strict seam with the depth-two THM-2073 normal form must satisfy the sharp
necessary condition

```text
H_(Q_2)(h_0,h_1)>1.                                    (16)
```

This is the functional form of the first nontrivial guard drift. It retains
the depth-two carrier through `tau(Q_2)` and the first two ordered guards; a
guard multiset, component count, or unlabelled wall word cannot recover it.

### Quarter-anchor bounded-fan corollary

Let `c_0` be the smallest member of `C` divisible by `4`, and assume

```text
c_0<=c<=7c_0 for every c in C divisible by 4,
c<=5c_0/2 for every odd c in C,
c<=6c_0   for every c in C with c=2 mod 4.               (17)
```

Then the recursion applies with

```text
t=2/(7c_0),       theta=1/4+1/(14c_0).                  (18)
```

For `c=4d`, condition (17) gives

```text
1/14<=dt=c/(14c_0)<=1/2,
```

so `t in G_D`. The last two inequalities in (17) are exactly the two drift
budgets in (10). This proves the corollary. In particular, (17) holds when
`c_0` is the unique multiple of `4` and `c_0=max(C)`. The bounded fan spreads
the `4`-divisible teeth across `[1/14,1/2]`, while the quarter anchor leaves
wide safe bands for both remaining residue layers. QED.

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
3. **Tower drift lemma.** On every depth-at-least-two THM-2073 residual,
   contradict (16) from divisor completeness, the safe-child addresses, or
   the metric guard bounds; equivalently, force one sufficiently early safe
   phase of `Q_2`.

These targets are complementary: the first retains labelled rational
ownership, the second bypasses ownership through a two-point geometric
obstruction, and the third transports that obstruction down the dyadic
safe-child tower.

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
lane. On a depth-two tower a sufficient drift sidecar is the ordered
tuple `(Q_2,tau(Q_2),h_0,h_1)` from (15), not a tournament on its speeds or
guards.
