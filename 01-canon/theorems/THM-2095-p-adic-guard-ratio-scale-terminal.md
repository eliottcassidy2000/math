---
id: THM-2095
title: "The bounded guard-ratio branch has a finite commensurability scale"
status: >
  PROVED from settled LRC through six speeds and an exact finite valuation-
  profile wheel. If THM-2087 Branch I has h=s d and q=r d, then every
  p-power dividing d must satisfy a six-comb gcd-deck capacity inequality at
  every intermediate p-adic layer. Hereditary terminal primitivity forces at
  least two unit valuations. The resulting exact table gives
  d|2^8*3^4*5^2*7^2*11*17*23*29 in the general odd-guard containment lemma.
  In the live terminal branch d is odd and 7 does not divide h, so actually
  d|3^4*5^2*11*17*23*29=252576225 and h,q<=14396844825.
  The ratio q=h is separately excluded by the mixed-overlap scalar wall.
  This bounds the common scale and the marked commensurate pair, not the
  other six terminal speeds, so it is not LRC(14).
source: codex-2026-07-22-LRC14-p-adic-guard-ratio
depends_on:
  - THM-2073
  - THM-2080
  - THM-2081
  - THM-2086
  - THM-2087
related:
  - THM-2069
  - THM-2082
  - THM-2094
  - THM-2097
script: 04-computation/lrc14_p_adic_guard_ratio_terminal_referee_codex_20260722.py
output: 05-knowledge/results/lrc14_p_adic_guard_ratio_terminal_referee_codex_20260722.out
---

# THM-2095 -- p-adic terminal for the guard-ratio scale

Retain a rank-seven terminal containment

```text
G_Q subset E_h,
Q={q_*,q_1,...,q_6} hereditarily primitive,             (1)
```

in the notation of THM-2081. Suppose THM-2087 Branch I supplies coprime
positive integers `r,s<=57` and `d>=1` with

```text
h=s d,                    q_*=r d.                      (2)
```

The last guard `h` is odd, so `s` and `d` are odd. THM-2086 also proves on
every live rank-seven terminal that

```text
7 does not divide h.                                    (3)
```

We first prove a slightly more general prime-power lemma before imposing
these two eliminations.

## 1. A safe invariant subdeck at every p-adic layer

Let `p^e|d`, and fix `1<=f<=e`. Put `m=p^f` and let

```text
R_f={w in Q:m|w}.                                       (4)
```

The distinguished speed `q_*` lies in `R_f`. Hereditary primitivity says at
least two terminal speeds are not divisible by `p`: if at most one were a
`p`-unit, deleting that one would leave a row of gcd divisible by `p`.
Consequently

```text
|R_f|<=5.                                               (5)
```

Apply settled LRC to the at most six positive integers

```text
{h/m} union {w/m:w in R_f}.                             (6)
```

Duplicates, if any, only lower the number of distinct speeds. There is a
phase `u` at which every integer in (6) has circle distance at least `1/7`.
Consider the complete `m`-deck

```text
t_j=(u+j)/m,                 0<=j<m.                    (7)
```

The guard and every speed in `R_f` have the same circle value at all these
lifts, so the guard is weakly safe and those terminal speeds are even safer
than the required `1/14`. Under containment (1), the speeds outside `R_f`
must therefore cover all `m` deck points by their strict danger combs.

If `w` has `v=nu_p(w)<f`, its orbit on the deck has length `p^(f-v)` and
each orbit point occurs `p^v` times. An open circle interval of length `1/7`
contains at most

```text
floor(p^(f-v)/7)+1                                     (8)
```

points of that uniform orbit. Hence the exact necessary capacity inequality
is

```text
p^f <= sum_(i:nu_p(q_i)<f)
  p^nu_p(q_i) * (floor(p^(f-nu_p(q_i))/7)+1).           (9)
```

The six indices in (9) exclude `q_*`, which is invariant. The inequality is
valid with room to spare at endpoints: (8) is an upper bound for an open
tooth, so endpoint coincidences can only lower actual capacity.

## 2. The finite valuation-profile wheel

Cap the six valuations outside `q_*` at `e` and sort them:

```text
0=v_1=v_2<=v_3<=...<=v_6<=e.                           (10)
```

The two leading zeros are forced by hereditary primitivity. In this notation
the complete necessary system is

```text
p^f <= C_f(p;v)
     =sum_(i:v_i<f) p^v_i(floor(p^(f-v_i)/7)+1),
1<=f<=e.                                                (11)
```

This is a finite wheel: only the four entries `v_3,...,v_6` vary in
`{0,...,e}`, with repetitions allowed.

There is already a uniform prime cutoff. If `k` valuations are below `f`,
then `k<=6`, and

```text
C_f<=k p^f/7+sum_(i:v_i<f)p^v_i
   <=6p^f/7+6p^(f-1).                                  (12)
```

Thus (11) forces `p<=42`. At `f=1`, all contributing valuations are zero,
so the sharper exact test is

```text
p<=k(floor(p/7)+1),             2<=k<=6.               (13)
```

Testing the primes through `41` leaves only

```text
p in {2,3,5,7,11,17,23,29}.                            (14)
```

In particular the superficially small primes `13` and `19` are excluded:
six danger teeth cover at most `12<13` and `18<19` points on their prime
decks.

Exhausting the four capped entries in (10) gives the exact terminal table

```text
p       maximum e      one last-level profile
2           8          (0,0,1,2,2,5)
3           4          (0,0,0,1,2,2)
5           2          (0,0,0,0,0,1)
7           2          (0,0,0,0,0,1)
11          1          (0,0,0,0,0,0)
17          1          (0,0,0,0,0,0)
23          1          (0,0,0,0,0,0)
29          1          (0,0,0,0,0,0).                 (15)
```

The companion checks (11) with integer arithmetic. The first infeasible
levels are respectively `9,5,3,3,2,2,2,2`. This proves all higher levels
infeasible as well: any profile at a higher exponent, capped at a lower
level, would satisfy the already impossible lower system.

It follows that the general containment lemma forces

```text
d divides 2^8*3^4*5^2*7^2*11*17*23*29.                (16)
```

This is a genuine prime-power result, not merely a radical bound. The
valuation wheel is the p-adic lift of THM-2069's low-weight evaluation code,
with deck capacity retained.

## 3. Specialization to the live odd guard

Equation (2) and oddness make `d` odd. Equation (3) gives `7 not|d`.
It also gives `7 not|s`.
Therefore (16) sharpens to

```text
d divides D_*=3^4*5^2*11*17*23*29
             =252576225.                               (17)
```

Since `r,s<=57`, the marked commensurate pair is absolutely bounded:

```text
h<=57D_*=14396844825,
q_*<=57D_*=14396844825.                                (18)
```

The surviving marked-pair ledger is already small enough to name exactly.
The number `D_*` has

```text
tau(D_*)=(4+1)(2+1)2^4=240
```

positive divisors. There are exactly `1165` coprime pairs

```text
1<=r,s<=57,       s odd,       7 not|s,
(r,s)!=(1,1).                                           (18a)
```

Thus only

```text
240*1165=279600                                         (18b)
```

triples `(r,s,d)` survive before the other terminal labels are attached.
This count keeps whether `7|r`, hence whether the marked speed `q_*` belongs
to THM-2094's remaining one-to-three seven-carrier packet.

THM-2094 independently reduces the live terminal residue profile to one,
two, or three speeds divisible by seven. The p-adic proof above needs only
the earlier `7 not|h` conclusion; the stronger fourfold moment obstruction
remains an additional filter on the six unbounded terminal coordinates.

## 4. The diagonal ratio is impossible

The ratio `(r,s)=(1,1)` would give `q_*=h`. Its mixed overlap is

```text
measure(E_h intersect D_(q_*))=1/7.                    (19)
```

Among the other six distinct speeds, THM-2086's odd-guard spectrum permits
at most one exception of each type `q=6h` and `q=h/11`; every other overlap
is at least `1/35`. Hence

```text
sum_(q in Q) measure(E_h intersect D_q)
 >=1/7+4/35+1/42+2/77
 =709/2310
 =2/7+7/330.                                           (20)
```

But THM-2081 containment forces this sum to be at most `2/7`. Thus

```text
q_*!=h.                                                 (21)
```

This exclusion is analytic and separate from the p-adic deck: the diagonal
ratio has `d=h` but fails before any valuation classification is needed.

## 5. Frontier effect and scope

The formerly unbounded part of THM-2087 Branch I is no longer its common
scale. There are only finitely many triples

```text
(r,s,d),      r,s<=57, gcd(r,s)=1, s,d odd,
7 not|s,      d|D_*,        (r,s)!=(1,1).              (22)
```

The remaining six terminal speeds can still be unbounded. A pair relation
can be reused as a zero-private-coefficient relation in every incident
triple, so THM-2087's complete cut does not automatically parameterize those
six speeds. The next lawful split is on the induced short-relation graph after
deleting all guard-commensurate vertices, with THM-2091/2094 energy and
residue sidecars retained. Claiming that the old cut propagates through the
contaminated vertex would be circular.

THM-2097 subsequently closes that infinitary residual by a different route:
its mixed-threshold two-torus cell makes every fixed rank-two terminal
coefficient template finite. The present theorem remains a uniform arithmetic
sharpening which bounds `d,h,q_*` independently of the coefficient template;
it should feed the finite enumeration rather than duplicate THM-2097's
geodesic terminal.

## 6. Assumption challenge and Tournament Analysis

The challenged assumption is that only the prime support of `d` matters.
Every intermediate power `p^f` supplies a different invariant deck, and the
six noninvariant valuations must cover it. Forgetting those layers loses the
finite exponent bounds in (15). Conversely, orienting speeds by valuation
keeps only their order and loses the capacity summands in (11).

The useful vertices are the p-adic layers `f=1,...,e`, or equivalently the
six valuation columns against those layers, not the runners themselves.
Orienting layers by depth gives a transitive tournament with score histogram
`(0,1,...,e-1)`, no directed cycles, singleton SCCs, and one Hamiltonian path.
That fingerprint contains no capacity information. The faithful carrier is
the valuation Ferrers diagram decorated by the exact deck loads `C_f`.
QED.

## 7. Exact referee

The companion exhausts every capped profile through the first infeasible
level for each prime allowed by (13), verifies the monotone cutoff and the
constants (16)--(18), and checks the exact mixed-overlap margin `7/330` in
the diagonal case. Runtime checks remain active under optimization; normal
and `python -O` outputs match the stored transcript and end in `PASS`.
