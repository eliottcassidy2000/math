---
id: THM-4129
title: "Universal two-speed completion of the eleven-speed LRC(14) body"
status: >
  PROVED ELEMENTARY UNIVERSAL COMPLETION + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. Let U=(1,4,6,8,10,12,14,15,16,18,22). For every two distinct
  positive integers a,b outside U, the thirteen-speed set U union {a,b} has
  a common phase of clearance at least 1/14. The proof combines a closed
  body-safe interval of length 3/280, a sharp periodic danger-comb discrepancy
  bound, nineteen exact high-tail measures, eleven low-tail body intervals,
  and twenty-six residual pairs covered by three named clocks. It closes
  every genuine scale-one U-body slice of all 5,855 THM-3878 ratio types, but
  no arbitrary-body type or scale-two slice. LRC(14) remains open.
source: codex-frontier-synthesis-creative-20260825j
depends_on: []
related:
  - THM-764-covering-small-period-signed-pair-deck-and-q25-refutation
  - THM-3878-lrc14-eleven-plus-two-harmonic-absorption-seam-collapse
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
  - THM-4117-physical-eleven-plus-two-primitive-support-obstruction
  - THM-4119-infinite-supplier-free-eleven-plus-two-residue-family
  - THM-4121-sharp-projective-tail-multiplier-residue-compiler
  - THM-4125-arbitrary-multitail-projective-packing-and-sharp-density
script: 04-computation/lrc14_universal_two_speed_completion_thm4129.py
output: 05-knowledge/results/lrc14_universal_two_speed_completion_thm4129.out
independent_audit_script: 04-computation/lrc14_universal_two_speed_completion_thm4129_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_universal_two_speed_completion_thm4129_independent_audit.out
script_sha256: f027c93324096de71a63f6cefcf7d9f75a6ef36b0b878c8b4248d6871d783cfa
output_sha256: 9a31e2ab7d319e63df067cfbfd05059a52b8e032b2fa1ae5bf7de93342d087b8
semantic_sha256: a717248ff237c9cbadf1a15dbaeb2902c0b7b7da10ebca0e29d797ea09bbdbe5
independent_audit_script_sha256: f977c65658b2dad6374196866c6f6211a1b833d3d24c9b852cf6fa959ea2bbe7
independent_audit_output_sha256: a7e4df345ae9c50ce15ccef4a4e9dcfe31902df499e2733d27991d0d7ddfa456
independent_semantic_sha256: a717248ff237c9cbadf1a15dbaeb2902c0b7b7da10ebca0e29d797ea09bbdbe5
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone Fraction implementation proves the body and low-tail
  closed intervals by endpoint concavity, constructs strict danger-comb
  components directly, checks the exact nineteen-entry high-tail table and
  twenty-six low residuals, prints a literal witness for all 119,316 pairs
  of missing speeds through 500, rebuilds all 5,855 scale-one ratio types by
  trial factorization, and checks six dilation controls. Normal, optimized,
  hash-seeded, and frozen outputs byte-match.
independent_audit: >
  ACCEPT. A clean-room implementation imports no primary code. It computes
  danger measures from a periodic cumulative antiderivative, finds the 171
  small high-tail witnesses on an exact common wall grid, verifies the low
  compact-to-open and three-clock routes, and independently rebuilds the
  5,855-type atlas with a smallest-prime-factor sieve. Normal, optimized,
  hash-seeded, and frozen outputs byte-match; the semantic ledgers agree.
---

# THM-4129 -- universal two-speed completion of the eleven-speed LRC(14) body

**PROVED ELEMENTARY UNIVERSAL COMPLETION + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.**

THM-4125's seven-tail transition makes every nonzero parameter class reach
the natural value `1/19` at the fixed phase `9/19`. The load-bearing feature
is not the residue count itself: the eleven-speed body has a whole closed
`1/14`-safe interval around that phase. Moving from the single phase to that
interval replaces projective residue blockers by two short periodic danger
combs. Their union never covers the interval. A finite low-speed repair then
removes the requirement that the two adjoined speeds be tails at all.

## 1. Universal completion theorem

Put

```text
U=(1,4,6,8,10,12,14,15,16,18,22),       delta=1/14.    (1)
```

> **Theorem.** Let `a,b` be distinct positive integers with
> `a,b notin U`. Then the thirteen distinct positive speeds
>
> ```text
> S_(a,b)=U union {a,b}                                  (2)
> ```
>
> have a phase `x in R/Z` such that
>
> ```text
> min_(v in S_(a,b)) ||v x|| >= 1/14.                   (3)
> ```

The order of `a,b` is immaterial; below take `a<b`.

## 2. The closed body window

Define

```text
J=[33/70,27/56],                    |J|=3/280.           (4)
```

For each `u in U`, the interval `uJ` lies in one interval between consecutive
integers. On such an interval distance to the nearest integer is concave, so
its minimum occurs at an endpoint. The endpoint table is

| `u` | `||33u/70||` | `||27u/56||` |
|---:|---:|---:|
| 1 | `33/70` | `27/56` |
| 4 | `4/35` | `1/14` |
| 6 | `6/35` | `3/28` |
| 8 | `8/35` | `1/7` |
| 10 | `2/7` | `5/28` |
| 12 | `12/35` | `3/14` |
| 14 | `2/5` | `1/4` |
| 15 | `1/14` | `13/56` |
| 16 | `16/35` | `2/7` |
| 18 | `17/35` | `9/28` |
| 22 | `13/35` | `11/28` |

Hence every member of `U` is weakly `delta`-safe throughout `J`. The two
endpoints are sharp: runner `15` owns the left equality and runner `4` owns
the right equality. The centre inherited from THM-4121/4125 satisfies

```text
9/19 in int(J),             min_(u in U)||9u/19||=2/19. (5)
```

## 3. A periodic danger-comb lemma

For an integer `v>=1` and a real interval `I`, put

```text
D_v(I)={x in I:||v x||<1/14},            ell=|I|.       (6)
```

The inequality in `(6)` is strict because equality is a valid LRC witness.

> **Lemma (one-comb discrepancy).** For every `v` and `I`,
>
> ```text
> measure(D_v(I)) <= ell/7+6/(49v).                     (7)
> ```

### Proof

After the change of variable `y=vx`, the danger set is the periodic arc
`||y||<1/14`, of length `d=1/7` in every unit period. Write
`v ell=N+r`, where `N` is an integer and `0<=r<1`. The first `N` unit
intervals contribute exactly `Nd`; the remaining interval meets a circular
arc of length `d` in at most `min(r,d)`. Therefore its excess above density
`d` is at most

```text
min(r,d)-rd <= d(1-d)=6/49.
```

Divide by `v` to obtain `(7)`. **QED.**

There is also a closed exact formula. For `y>=0`, set

```text
A(y)=floor(y)/7+psi({y}),
psi(r)=r                 for 0<=r<=1/14,
       =1/14             for 1/14<=r<=13/14,
       =r-6/7            for 13/14<=r<1.               (8)
```

Then for `I=[L,R]`,

```text
measure(D_v(I))=(A(vR)-A(vL))/v.                        (9)
```

Formula `(9)` supplies an independent exact route to every finite value used
below.

## 4. Two arbitrary speeds above the body

Let `m(v)=measure(D_v(J))`. Exact substitution in `(9)` gives

| `v` | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `m(v)` | `1/161` | 0 | `1/200` | 0 | `5/1512` | 0 | `3/1624` | 0 | `1/1736` | 0 |

| `v` | 33 | 34 | 35 | 36 | 37 | 38 | 39 | 40 | 41 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `m(v)` | 0 | `3/2380` | 0 | `1/360` | 0 | `1/266` | 0 | `1/280` | 0 |

The first four entries, together with `(7)`, prove the uniform ceiling

```text
m(v)<=1/161                       for every v>=23,       (10)
```

because for `v>=27`,

```text
|J|/7+6/(49v) <= 107/17640 < 1/161.                    (11)
```

Now suppose `23<=a<b`. If `b>=42`, equations `(7)` and `(10)` give

```text
m(a)+m(b)
 <=1/161+|J|/7+6/(49*42)
 =3363/315560
 =|J|-9/157780
 <|J|.                                                  (12)
```

Thus `D_a(J) union D_b(J)` cannot cover `J`; a point outside it proves `(3)`.

It remains to check `23<=a<b<=41`. Among the `171` pairs in the displayed
table, the sum `m(a)+m(b)` is below `|J|` except at `(a,b)=(23,25)`.
Every nonexceptional pair is again closed by the union bound. For the unique
exception, the strict interior phase

```text
x=381/805 in int(J),
min_(v in U union {23,25})||vx||=16/161>1/14            (13)
```

is a direct witness. This proves the theorem whenever both adjoined speeds
exceed `22`.

## 5. Low-speed completion

Suppose now that `a<=22`. Because `a` is missing from `U`,

```text
a in M={2,3,5,7,9,11,13,17,19,20,21}.                  (14)
```

For each `a in M`, the twelve-speed body `U union {a}` has the following
closed `delta`-safe interval. The endpoint-owner column lists the unique
left and right equality owners.

| `a` | safe interval `I_a` | `|I_a|` | endpoint owners | first integer `B_a` with `B_a|I_a|>=1/7` |
|---:|---:|---:|---:|---:|
| 2 | `[29/196,13/84]` | `1/147` | `(14,6)` | 21 |
| 3,5,7,9,11,13 | `J` | `3/280` | `(15,4)` | 14 |
| 17 | `[113/238,27/56]` | `1/136` | `(17,4)` | 20 |
| 19 | `[29/84,69/196]` | `1/147` | `(6,14)` | 21 |
| 20 | `J` | `3/280` | `(15,4)` | 14 |
| 21 | `[29/196,13/84]` | `1/147` | `(14,6)` | 21 |

As in Section 2, every speed times its displayed interval stays between the
same consecutive integers at both endpoints. Endpoint distance is at least
`1/14`, so concavity proves uniform closed safety.

If `b|I_a|>=1/7` and the completed row had no witness in `I_a`, the compact
connected interval `I_a` would lie inside one component of the **open** set

```text
D_b={x:||bx||<1/14}.                                    (15)
```

Every component of `(15)` has length exactly `1/(7b)<=|I_a|`. A closed
interval cannot fit inside an open interval of no greater length. This is a
contradiction. Notice that equality `b|I_a|=1/7` is already enough; this is
where the closed target and strict danger convention is load-bearing.

Only twenty-six pairs fail that component-length gate. They split among
three explicit clocks. Put

```text
P_13=binom({2,3,5,7,9,11},2)
     union {(2,17),(2,19),(2,20),(17,19),(19,20)},
P_20={(2,13)},
P_19={(3,13),(5,13),(7,13),(9,13),(11,13)}.             (16)
```

The sets have sizes `20,1,5`. Exact residue arithmetic gives

```text
(a,b) in P_13:  x=1/13,    clearance=1/13,
(a,b) in P_20:  x=7/20,    clearance=1/10,
(a,b) in P_19:  x=9/19,    clearance=2/19.              (17)
```

All three values exceed `1/14`. Equations `(14)--(17)` close every remaining
low pair and complete the proof. **QED.**

## 6. Inheritance and exact frontier effect

The closest proved mechanism is THM-4125, with THM-4121 its two-tail
specialization. The canonical hostile is a parameter divisible by `19`: it
kills the displayed phase `9/19` outright, so a residue-only continuation
cannot prove the universal result. THM-2050's local-germ blindness is the
broader correction behind the same warning. The least-used relevant sidecar
is THM-3910's compact-to-open body-component response; Section 5 specializes
and sharpens it using literal intervals of this body. THM-764's signed deck
is the discrete analogue of the strict danger-comb calculation in Section 3.

The connection is

```text
source:       the THM-4121/4125 phase 9/19 and fixed body U
target:       the closed body-safe interval J and two strict danger combs
map:          enlarge the phase to its connected 1/14-safe body component
preserved:    every body clearance and the common-time requirement
destroyed:    the projective residue label of each adjoined speed
sidecar:      absolute tail sizes, then the eleven low-speed intervals
decisive test: whether D_a(J) union D_b(J) covers J.      (18)
```

Consequences are exact.

1. **THM-4121 is completed globally.** Every `U union {t,ct}` with
   `c>1,t>22` is lonely, including all residue classes that fail its fixed
   phase.
2. **THM-4119 is completed globally.** Every `U union {t,3t}`, `t>22`, is
   lonely, not merely its density-`14/19` phase subfamily. THM-4117's physical
   `t=2^45` row is included, although that row was already known safe.
3. **Every two-tail projection of THM-4125 is closed at the LRC(14)
   threshold.** For any two distinct multipliers in its bank, delete the
   other tails and apply the theorem to their two actual speeds. This says
   nothing about the larger bank at the stronger fixed threshold `1/14`.
4. **All `5,855` scale-one ratio types have their `U`-body slice closed.** A
   genuine thirteen-speed THM-3878 row `U union t{p,q}` is an instance of
   `(2)`. In particular, this closes the `U`-body cells of all sixteen live
   scale-one types listed by THM-3910. It does **not** remove an arbitrary-body
   type from that ledger. The scale-two `(s,p,q)=(2,1,9)` slice has body `2U`
   with odd tails and is not a common dilation of `(2)`; it is untouched.

Because `1 in S_(a,b)`, the row is primitive. For every common dilation
`k>=1`, if `x` proves `(3)` then

```text
||(kv)(x/k)||=||vx||,             nu(kS_(a,b))=S_(a,b). (19)
```

Thus the conclusion is invariant under common dilation and primitive
normalization. It does not cover translation, separate component dilation,
negative speeds, repeated-speed bookkeeping, arbitrary eleven-speed bodies,
physical entry into THM-3818, or LRC(14) in general.

## 7. Exact replay

The primary audit constructs every strict danger component by clipped
rational intervals. Besides all proof tables, it prints a literal witness for
all `119,316` pairs of missing positive speeds at most `500`, and rebuilds the
full `5,855`-type THM-3878 scale-one atlas by trial factorization. The
independent audit instead uses the cumulative formula `(8)--(9)`, scans the
exact common wall grid for all `171` small high-tail pairs, and rebuilds the
ratio atlas by a smallest-prime-factor sieve. Both include the named clocks,
endpoint owners, strict exceptional witness, and six dilation controls. Run

```text
python3 -B 04-computation/lrc14_universal_two_speed_completion_thm4129.py
python3 -B -O 04-computation/lrc14_universal_two_speed_completion_thm4129.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc14_universal_two_speed_completion_thm4129.py
python3 -B 04-computation/lrc14_universal_two_speed_completion_thm4129_independent_audit.py
python3 -B -O 04-computation/lrc14_universal_two_speed_completion_thm4129_independent_audit.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc14_universal_two_speed_completion_thm4129_independent_audit.py
```

All six streams byte-match their frozen outputs. The independent semantic
ledgers agree. **QED.**
