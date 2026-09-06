# The scale-three consumer needs both parity and body addresses

**Status: PROVED for the adjacent-tail family; FINITE-EXACT for the stated
heads; REFUTED for the proposed parity-free ceiling and universal ten-body
floor. Independent audit: [analytic and literal referee](overnight3_20260906_lrc_consumers_audit.md).**
These are obstructions to sufficient certificates, not to LRC(14).
The body extremizer is recovered prior canon, explicitly credited below.
No external novelty or general mixed-parity ceiling is claimed.

## 1. Inheritance and the actual consumer

The incoming [universal local certificate](lrc14_global_slope_empty_core_certificate_sep06.md)
and its [Haar consumer](lrc14_universal_haar_consumer_empty_core_certificate_sep06.md)
prove, in their declared odd-tail domain,

```text
mu(G_C)>=6/77  ==> G_(3C union T) is nonempty,
G_C={y in R/Z: ||cy||>=1/14 for every c in C}.
```

Here `T` has three distinct positive odd ternary units. Their proof keeps
the common branch label of all three lifts. The closest hostile is
[THM-4004 / three-detuned divisor combs, Section 3](../../01-canon/theorems/THM-4004-lrc14-three-detuned-divisor-comb-profile.md):
`C={1,...,10}`, `T=(1,10,11)` can spoil all lifts of a valid body phase,
although another phase wins. The corrected near miss is treating body
nonemptiness as a quantitative measure floor. The recovered sidecar is
[THM-530 / G_P-intersection global witness floor, Section A](../../01-canon/theorems/THM-530-lrc-gp-intersection-global-witness-floor.md)'s
finite small-part extremizer, including its actual component addresses.
Only that finite exact section is inherited here; no later covering-route
claim in THM-530 is needed.

The live concepts are owner-labeled sheet danger, relation parity,
arithmetic sampling of a tent, complete safe components, and weak endpoints.
Both tempting new inputs were given literal hostile probes before use.

## 2. Parity is load-bearing

For three distinct positive integers prime to three define the physical set

```text
F_T={x: for every j in {0,1,2}, some w in T satisfies
             ||w(x+j/3)||<1/14}.
```

Each tail can spoil at most one sheet. Consequently all three sheets fail
exactly when the tails occupy a permutation of the three sheet labels.
This statement uses the ternary-unit condition, but uses no parity.

The smallest primitive positive mixed-parity counterexample by maximum
speed is

```text
T=(2,5,7),              mu(F_T)=22/245 > 6/77.
```

This minimality is a **FINITE-EXACT** statement: the complete primitive
ternary-unit universe with maximum below seven has no violation. The
larger head through sixty contains 8,664 primitive triples and 136 violations.
Its maximum is `6/55`, uniquely at `(1,10,11)`. That maximum is only a
finite statement for arbitrary triples.

The first failed implication is importing an odd-tail bound after dropping
oddness. The relation `(1,1,-1)` can vanish on additive tails `a+b=c` and
has norm three; it cannot vanish on three odd speeds. This supplies an
infinite mechanism, proved sharply on one family in Section 3.

Deleting only additive triples does not repair the claim:

```text
T=(2,11,20),            mu(F_T)=11/140 > 6/77.
```

This arithmetic progression has relation `(1,-2,1)` of norm four. It is
the sole nonadditive violation in the height-sixty head, which does not
classify the unbounded nonadditive domain. A parity-free analogue needs
its own short-relation analysis, including norm four.

## 3. A sharp all-height additive family

**PROVED.** Let `b>=4` and `b=1 mod3`. Put `M(b)=mu(F_(1,b,b+1))`.
Then

```text
M(4)=1/28;
6/77 < M(b) <= 6/55       for every b>=7, b=1 mod3;
M(b)=6/55                iff b=10;
lim M(b)=9/98            as b -> infinity in this progression.
```

The family is primitive, distinct, and consists of ternary units. One of
the two larger speeds is even, so every member lies outside the odd theorem.

Here is a direct physical proof, independent of the odd carrier formula.
The set `F_T` is invariant under translation by `1/3` and reflection.
On `0<x<1/14`, speed one spoils sheet zero. The other tails must occupy
the other sheets. Since `b=1 mod3`, their bad phases have common centers
`t=k/3`, `k>=1`, `3` not dividing `k`. The two centers use the same integer
translate: their difference is `x`, whose magnitude is too small to allow
a different translate. The corresponding interval for `x` is

```text
((t-1/14)/b, min(1/14,(t+1/14)/(b+1))).
```

It is nonempty precisely when `t<(b+1)/14`. Its length is

```text
((2b+1)/14-t)/(b(b+1))       if t<=b/14,
((b+1)/14-t)/b              if b/14<t<(b+1)/14.
```

The center intervals are disjoint. Multiplying their sum by six accounts
for reflection and the three translates. Strict endpoints have zero
measure and cannot change this formula.

Write `b=42h+r`, `r in {1,4,...,40}`. Summing integers except multiples
of three gives

```text
M(b)=(1134h^2+B_r h+C_r)/(7b(b+1)).
```

| r | B_r | C_r |
|---:|---:|---:|
| 1 | 54 | 0 |
| 4 | 204 | 5 |
| 7 | 396 | 31 |
| 10 | 630 | 84 |
| 13 | 738 | 120 |
| 16 | 846 | 156 |
| 19 | 1080 | 253 |
| 22 | 1188 | 307 |
| 25 | 1422 | 444 |
| 28 | 1530 | 516 |
| 31 | 1638 | 588 |
| 34 | 1872 | 769 |
| 37 | 2064 | 935 |
| 40 | 2214 | 1080 |

For completeness, put `Q=1134h^2+B_r h+C_r`.
The lower-bound difference is `77Q-42b(b+1)`. Every coefficient is
strictly positive, using `h=u+1` for residues one and four (the condition
`b>=7`), and using `h>=0` for the other residues. The upper difference is
`42b(b+1)-55Q`; all its coefficients are nonnegative. Its only admissible
zero is `r=10,h=0`, where it equals `11718h^2+2394h`.
The displayed table proves these finite polynomial assertions by direct
expansion; the executable certificate also checks them symbolically.
The ratio of leading coefficients is `1134/(7*42^2)=9/98`.

The [independent general additive analysis](overnight3_20260906_lrc_additive_audit.md)
gives the broader mechanism: for primitive ternary-unit additive triples,
the physical mass differs from `9/98` by at most `6/(7c)`. Thus every such
triple with `c>=62` exceeds `6/77`. The exact adjacent formula above is
stronger on its stated one-parameter subclass. Neither statement bounds
every mixed-parity triple from above by `6/55`.

## 4. The forgotten ten-body extremizer already defeats the new floor

Let

```text
C_*=(1,2,3,5,7,8,9,11,12,13).
```

THM-530 Section A already recorded the exact identity

```text
mu(G_(C_*))=14249/252252 < 6/77.
```

This session independently reconstructs it and connects it to the new
consumer. It refutes the universal ten-speed floor proposed in that
consumer's next-probe section. It does **not** establish a universal
lower bound of `14249/252252` for arbitrary ten-speed bodies; THM-530's
small-part universe is `C subset {1,...,13}`.

Among ten-subsets of `{1,...,12}` the exact minimum is `16277/194040`,
above `6/77`. At maximum thirteen there are twelve sets below `6/77`.
Thus thirteen is the first possible maximum speed for a positive integral
ten-body counterexample. This follows from exhaustive finite enumeration,
with no primitivity filter needed.

The positive-length components of `G_(C_*)` in the lower half-circle are

```text
[15/154,13/126], [29/182,9/56], [29/168,27/154],
[43/182,27/112], [29/112,41/154], [29/98,55/182].
```

Include their reflected components and the six isolated safe points
`r/14`, `r in {1,3,5,9,11,13}`. These zero-length components matter for
weak safety even though they do not affect Haar mass.

The sufficient gate fails on an actually safe odd-tail row:

```text
mu(G_(3C_* union {1,5,11}))=25331/756756 > 0,
x=1/14 is weakly safe.
```

Its tail failure-set mass is `6/77`, above the body's mass. This is not
a new LRC family closure: the fixed denominator fourteen already certifies
the displayed row. It proves the intended logical boundary—failure of the
scalar Haar comparison does not imply failure of the literal completion.

## 5. Restored coordinates and decisive next questions

The source is `G_C` together with the full three-sheet tail failure set.
The map to the completion keeps `y=3x` and the common lift label. Scalar
measure preserves a sufficient noncontainment test but loses component
addresses and isolated endpoints. Retaining those addresses yields the
exact test `G_C \ m_3(F_T) != empty`. The scalar ten-body floor has now
stopped for a recorded reason; another arbitrary-body floor census would
not repair the lost information.

For parity, the source is the literal mixed-parity comb. Passing to the
odd-only relation catalogue loses norm-three relations and changes the
norm-four boundary. The strongest surviving statement here is the proved
odd-triple theorem plus the exact additive family and the two small
mixed-parity hostiles. The open next question is the correct parity-free
short-relation bound, or a direct body-to-comb noncontainment argument.
Neither is supplied by this report.

## 6. Reproduction and finite scope

Run from the repository root:

```text
python3 -B 04-computation/overnight3_20260906_lrc_consumers.py
```

The [program](../../04-computation/overnight3_20260906_lrc_consumers.py),
[frozen output](overnight3_20260906_lrc_consumers.out), and
[exact certificate](overnight3_20260906_lrc_consumers_certificates.json)
declare all universes. The physical engine intersects the six literal
sheet assignments on an exact common denominator; the independent referee
instead partitions at rational breakpoints and tests the actual predicate.
The adjacent-family proof is unbounded because of the exact residue
polynomials, not because of the height-sixty head. Body intersections
preserve weak endpoints. Normal and optimized output must agree before
checkpointing. No Lean proof or arbitrary ten-body floor is claimed.
