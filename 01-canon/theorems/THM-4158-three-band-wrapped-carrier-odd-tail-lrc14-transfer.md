---
id: THM-4158
title: "Three-band wrapped-carrier odd-tail LRC(14) transfer"
status: >
  PROVED ELEMENTARY COMPLETE WRAPPED-CARRIER ALPHABET + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED; LRC(14) OPEN. For every m>=2, every finite body
  drawn from the exact three-band alphabet P_m accepts every pair of distinct
  positive odd tails after doubling; at m=1 the exact alphabet has one extra
  band and the same transfer holds. Among eleven-subsets of [1,N], a
  canonical minimum-indexed certified subfamily has asymptotic density
  0.3421696706978653..., strictly above THM-4151's first-band density. The
  m=7 minimum-indexed seam has 81,572,506,886,508 bodies; its
  38,620,298,376-body divisor-complete subfamily lies outside the stated
  THM-4148/4151 gates and is disjoint from THM-4156's explicit family. Entry,
  mixed/even tails, arbitrary thirteen-speed rows, and LRC(14) remain open.
source: codex-frontier-synthesis-creative-20260825as
supplementary_source: codex-open-frontiers-overnight-20260825
depends_on:
  - THM-4151-scale-sensitive-first-window-odd-tail-lrc14-transfer
related:
  - THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4154-mod-six-fixed-clock-and-haar-pool-inheritance-correction
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4160-anchored-haar-deletion-cover-and-content-tower
script: 04-computation/lrc14_three_band_wrapped_carrier_odd_tail_transfer_thm4158.py
output: 05-knowledge/results/lrc14_three_band_wrapped_carrier_odd_tail_transfer_thm4158.out
independent_audit_script: 04-computation/lrc14_three_band_wrapped_carrier_odd_tail_transfer_thm4158_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_three_band_wrapped_carrier_odd_tail_transfer_thm4158_independent_audit.out
cross_language_audit_script: 04-computation/lrc14_three_band_wrapped_carrier_odd_tail_transfer_thm4158_independent_cpp_audit.cpp
cross_language_audit_output: 05-knowledge/results/lrc14_three_band_wrapped_carrier_odd_tail_transfer_thm4158_independent_cpp_audit.out
script_sha256: 8595784844856f61326728d80f53980d239e14e1400574bea374b382794e37ab
output_sha256: 3eb3ae21b6b644b1d5211a2b890f6043d6bf462e3707b406d721535bd8aed753
semantic_sha256: 7b71523951db221ce08ff16c2c87629e1ff896b052f7a1324032bd4bca7314eb
semantic_fnv64: fbb8b1eeaef5dd11
independent_audit_script_sha256: 4c0da48086a3bab5bcf947b5b113ba71c18047060cb4db90094a4df894680b5d
independent_audit_output_sha256: 21da776146c19d39ba5a15f477b58d417b4f57a956d60aed3c8657fde7aedffc
cross_language_audit_script_sha256: fa657f2aef9a093ca0b3fffd385a6f3683df49a3899df78cee66dfe0e5a0c9d3
cross_language_audit_output_sha256: 7577edb6de80f4b6f7e4b32c1c70cd539252ecdc981aaf4f839cda65dd59286f
cross_language_semantic_fnv64: 435c75c9215201cb
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact Fraction arithmetic reconstructs the complete pool for
  m=1,...,1000, checks the odd-wall and sharp endpoint identities through
  m=300, directly resolves 7,030 strict-wall tail rows, verifies both m=1
  censuses, the m=7 families and hostiles, and reproduces the finite counts.
  Every mandatory check survives optimized execution; normal, optimized, and
  hash-seeded replays byte-match the frozen output.
independent_audit: >
  ACCEPT. A no-import implementation reconstructs the complete pool for
  m=1,...,250 from endpoint inequalities, independently computes the full
  m=7 safe set, all four safe walls, absence of isolated safe atoms, and Haar
  measure, then reproduces the m>=2 finite census and rational density.
  Normal, optimized, and hash-seeded outputs byte-match.
  A separate analytic audit verified the infinite odd-tail proof and its
  compact/open endpoint case.
cross_language_audit: >
  ACCEPT. Warning-clean standard-library C++ independently checks 15,456,403
  speed gates and 127,500 odd-tail pairs, reconstructs both finite censuses
  and the m=7 seams, and recomputes the density with a separate BigUInt
  implementation. O0 and O3 outputs agree byte-for-byte.
---

# THM-4158 -- three-band wrapped-carrier odd-tail transfer

**PROVED ELEMENTARY COMPLETE WRAPPED-CARRIER ALPHABET + VERIFIED-EXACT +
INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

THM-4151 retained only the first unwrapped safe band above the anchored phase
`1/(14m)`. The missing coordinate is the winding index of the same carrier.
Keeping the next two windings enlarges the alphabet while leaving the
odd-wall proof unchanged: the transfer sees the carrier, not which safe band
a body label occupies.

## 1. The exact wrapped alphabet

Fix an integer `m>=1` and put

```text
J_m=[1/(14m),8/(7(12m+1))],
G_m=|J_m|=(4m-1)/(14m(12m+1)).                         (1)
```

For `m>=2`, define the three closed integer bands

```text
P_m=
 [m,floor(13(12m+1)/16)]
 union [15m,floor(27(12m+1)/16)]
 union [29m,floor(41(12m+1)/16)].                     (2)
```

Intervals in `(2)` mean all integers between their endpoints.

At the boundary `m=1`, use instead the exact four-band alphabet

```text
P_1=[1,10] union [15,21] union [29,33] union [43,44]. (2a)
```

> **Theorem.** For every finite `H subseteq P_m` and every two distinct
> positive odd integers `a,b`, there exists `x in R/Z` such that
>
> ```text
> min_(v in 2H union {a,b}) ||vx|| >= 1/14.            (3)
> ```

The alphabet is complete for this carrier:

```text
P_m={h in Z_(>0): ||hy||>=1/14 for every y in J_m}.    (4)
```

Thus `(2)` is not a sampled pool or a sufficient truncation. It is the exact
set of all positive integer labels safe throughout `J_m`.

## 2. Completeness of the three bands

A positive integer `h` is safe throughout `J_m` if and only if the connected
interval `hJ_m` lies in one closed component

```text
[k+1/14,k+13/14]
```

of the safe set on the real line. Since `h>0` and `J_m` is positive, `k>=0`.
The two endpoint inequalities are exactly

```text
m(14k+1) <= h
             <= floor((12m+1)(14k+13)/16).             (5)
```

The rows `k=0,1,2` are precisely the three bands in `(2)`. Before taking the
floor, upper endpoint minus lower endpoint has numerator

```text
m(140-56k)+14k+13.                                     (6)
```

For `m>=2`, the numerator at `k=3` is `55-28m<0`; every later row is still
smaller. At `m=1`, row `k=3` is `[43,44]`, while row `k=4` has numerator
`-15` and later rows decrease by `42`. Thus `(2)` and `(2a)` are the exact
fixed-carrier alphabets, including both closed endpoints.

For every `h in H` and `y in J_m`, equation `(4)` says that `hy` is safe.
Consequently both physical half-lifts

```text
x_0=y/2,                         x_1=(1+y)/2            (7)
```

keep the doubled body safe, because `2h x_i=hy mod 1`.

## 3. One odd tail kills at most one sheet

For an odd positive integer `r`, the two values `r x_0,r x_1` differ by
`1/2` modulo one. Hence `r` is strictly bad, with gap `<1/14`, on at most
one sheet. Moreover, it is bad on one of the sheets exactly when

```text
||ry||<1/7.                                            (8)
```

Interchange the actual tails so that `a<b`. Suppose for contradiction that
both lifts fail for every `y in J_m`. Each tail can kill at most one lift, so
each tail must kill one and the killers must be complementary. In particular,
`b` is bad on one sheet at every `y`, and `(8)` places the connected interval
`J_m` inside one open `b`-tooth

```text
((7k-1)/(7b),(7k+1)/(7b))                              (9)
```

for some integer `k`.

## 4. Endpoint clock and odd-wall closure

Write `y_0=1/(14m)` for the left endpoint. If `b<=12m`, oddness gives
`b<=12m-1`. Choose the upper lift at `y_0`:

```text
x=(1+y_0)/2.
```

For either tail `r in {a,b}`,

```text
||rx||=1/2-r/(28m)
       >=1/14+1/(28m)>1/14.                            (10)
```

The body is safe by `(7)`, contradicting total failure.

It remains that `b>=12m+1`. Let `lambda=(7k-1)/(7b)` be the left wall of
the tooth in `(9)`. Containment of the closed interval in the open tooth
forces `lambda<y_0`, while

```text
y_0-lambda=[b-2m(7k-1)]/(14mb)>=1/(14mb).              (11)
```

The numerator is a positive odd integer. Since a tooth has width `2/(7b)`,
the room from `y_0` to its right wall is at most

```text
2/(7b)-1/(14mb)
 =(4m-1)/(14mb)
 <=(4m-1)/(14m(12m+1))
 =|J_m|.                                                (12)
```

A closed interval of length `|J_m|` cannot lie inside that open tooth; at
equality its right endpoint reaches the excluded wall. This contradicts
`(9)` and proves `(3)`. No gcd, primitive-ratio, or residual finite-tail split
is needed. **QED.**

## 5. The carrier endpoint is sharp for this mechanism

Take the tails

```text
(a,b)=(1,12m+1).                                       (13)
```

Throughout `J_m` with its right endpoint removed, tail `1` kills the lower
lift and tail `12m+1` kills the upper lift. At the full right endpoint,

```text
y=8/(7(12m+1)),
```

the latter gap is exactly `1/14`, so the closed endpoint rescues the pair.
Thus any strict shortening of this carrier at the right loses the universal
two-sheet certificate for `(13)`. This is sharpness of the anchored carrier
mechanism, not a claim that the resulting speed row is globally unsafe.

## 6. Exact finite census and asymptotic density

Let

```text
p_m(N)=|P_m intersect [1,N]|.
```

Define the canonical minimum-indexed subfamily by requiring an eleven-element
body `H subset [1,N]`, with `m=min(H)>=2`, to satisfy `H subset P_m`. For
fixed `m`, its other ten labels are chosen from
`(P_m intersect [1,N])\{m}`. Hence this subfamily has exact count

```text
F(N)=sum_(m=2)^N binom(p_m(N)-1,10),                    (14)
```

with the binomial interpreted as zero when fewer than ten choices exist.
This is a collision-free certified subfamily, not the total union over all
possible carrier parameters. Exact arithmetic gives

| `N` | `F(N)` |
|---:|---:|
| 20 | 75,582 |
| 40 | 812,850,987 |
| 80 | 3,595,550,244,611 |
| 120 | 397,529,462,747,261 |
| 160 | 10,616,582,432,233,990 |
| 200 | 132,777,517,674,540,845 |

The boundary carrier certifies all `binom(24,11)=2,496,144` eleven-subsets
of `P_1`. To extend the collision-free minimum-indexed census, require the
body to contain `1` and put

```text
E_1(N)=binom(|P_1 intersect [1,N]|-1,10),
A(N)=F(N)+E_1(N).                                      (14a)
```

Thus `E_1(N)` stabilizes at `binom(23,10)=1,144,066` once `N>=44`:

| `N` | `E_1(N)` | `A(N)` |
|---:|---:|---:|
| 20 | 3,003 | 78,585 |
| 40 | 352,716 | 813,203,703 |
| 80 | 1,144,066 | 3,595,551,388,677 |
| 120 | 1,144,066 | 397,529,463,891,327 |
| 160 | 1,144,066 | 10,616,582,433,378,056 |
| 200 | 1,144,066 | 132,777,517,675,684,911 |

This finite boundary slice does not change the limiting density.

For `m/N -> u`, the normalized pool size tends to

```text
f(u)=
 (63/4)u,       0<=u<=4/123,
 1-15u,         4/123<=u<=1/29,
 14u,           1/29<=u<=4/81,
 1-(25/4)u,     4/81<=u<=1/15,
 (35/4)u,       1/15<=u<=4/39,
 1-u,           4/39<=u<=1.                            (15)
```

The Riemann-sum limit of `(14)`, normalized by `binom(N,11)`, is therefore

```text
11 integral_0^1 f(u)^10 du
=848953086913769850118498851618778832628468542103282298683365532079
 /2481088067163593416217816176836483026480276818419826456353950662656
=0.3421696706978653....                                 (16)
```

The first band alone is exactly THM-4151 and has density
`(35/39)^10=0.338870975802271291...`; the two wrapped bands add
`0.003298694895594002...` to the certified density.

## 7. A divisor-complete m=7 family

At `m=7`, the complete alphabet is

```text
P_7=[7,69] union [105,143] union [203,217],
|P_7|=117.                                              (17)
```

Every eleven-subset of `P_7` containing `7` is a distinct minimum-indexed
certificate, giving the larger unstructured seam

```text
binom(116,10)=81,572,506,886,508.                       (17a)
```

No primitivity, divisor coverage, or first-window failure is asserted for
all bodies in `(17a)`. For that stronger structured conclusion, require the
four anchors

```text
A={7,120,126,143}.                                      (18)
```

Their gcd is one. The anchor `120` owns moduli `2,3,4,5,6,8,10,12`, `126`
owns `7,9,14`, and `143` owns `11,13`. Thus every eleven-body

```text
H=A union K,                K subset P_7\A, |K|=7       (19)
```

is primitive and contains a multiple of every modulus `2,...,14`. There are

```text
binom(113,7)=38,620,298,376                              (20)
```

such bodies. Because every one has minimum `7` and contains `143`, all fail
the stated THM-4148 and THM-4151 gates. Every one contains `7`, which is
absent from THM-4156's explicit pool, so `(19)` is disjoint from that named
family. The old fixed clock is genuinely unavailable because

```text
||2*120*(1/12)||=0.                                     (21)
```

This comparison is deliberately mechanism-specific. The complete safe set
of the full alphabet `(17)` is exactly

```text
G_(P_7)=[1/98,13/966] union [953/966,97/98],
mu(G_(P_7))=22/3381<4/63.                               (22)
```

Thus THM-4150's Haar threshold does not subsume the full `P_7` carrier.
Deleting labels enlarges safe sets, however, so `(22)` does not prove that
each of the `38,620,298,376` subbodies lies outside the abstract THM-4150
gate. Their individual Haar classification remains open.

## 8. Reproduction and boundary control

The primary exact script checks the algebraic identities, direct pool
reconstruction, parity-wall bounds, hostile endpoint, finite rows, census,
and density. Its mandatory checks use explicit exceptions and therefore
remain live under `python -O`. The independent no-import script rebuilds the
pool from endpoint inequalities and the `m=7` safe set from its full rational
wall arrangement.

The `m=1` computation is not merely a sampled control: Section 2 proves the
fourth band is exact, and Sections 3--5 use only `m>=1`. Hence

```text
P_1=[1,10] union [15,21] union [29,33] union [43,44].   (23)
```

is the proved boundary corollary used in `(14a)`. The Python interval audit
also tests every safe wall of the full `P_7` arrangement and rules out
isolated safe atoms. A warning-clean C++ referee imports neither Python
implementation and independently rebuilds the carrier, direct tail pairs,
finite seams, and density.

```text
python3 -B 04-computation/lrc14_three_band_wrapped_carrier_odd_tail_transfer_thm4158.py
python3 -B -O 04-computation/lrc14_three_band_wrapped_carrier_odd_tail_transfer_thm4158.py
python3 -B 04-computation/lrc14_three_band_wrapped_carrier_odd_tail_transfer_thm4158_independent_audit.py

g++ -std=c++17 -O3 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_three_band_wrapped_carrier_odd_tail_transfer_thm4158_independent_cpp_audit.cpp \
  -o thm4158-independent-cpp-audit
./thm4158-independent-cpp-audit
```

The theorem handles the isolated dyadic two-odd-tail form. It does not prove
entry into that form, handle mixed/even tails, close arbitrary thirteen-speed
rows, or settle LRC(14). **QED.**
