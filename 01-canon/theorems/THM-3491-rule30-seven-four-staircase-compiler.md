---
id: THM-3491
title: "Rule 30 seven-four staircase entropy and thirteen-four compiler"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The physical zero-tail
  Rule 30 staircase language has an exact height-21 overlap-graph certificate
  proving eta<=log_2(7/4)<13/16.  In the charged word-RAM model of THM-3480,
  this improves the nonrectangular macroblock query tariff from 7/2 to 13/4.
  No Rule 30 prize or literature novelty is claimed.
source: root-rule30-next-targets-20260816
audit: >
  An independent hostile audit rebuilt the physical factor language with a
  separate tuple transducer, rederived the positive-vector entropy bound and
  all charged word-budget inequalities, and checked the finite-only lower
  signal.  Ordinary and optimized replays match the stored transcript
  byte-for-byte: ACCEPT.
depends_on:
  - THM-3480-rule30-staircase-transducer-entropy-and-nonrectangular-macroblock-compiler
related:
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
script: 04-computation/rule30_staircase_entropy_seven_four_thm3491.py
output: 05-knowledge/results/rule30_staircase_entropy_seven_four_thm3491.out
script_sha256: d107f8603d379b9cf9f196b55369612f8db8763993c038d0c7b674f68636c58f
output_sha256: 66e7d0c0c2abbf108c79aa8ccb7a728a3401434c7aa24460039a10b0de90d570
hash_basis: raw bytes
---

# THM-3491 -- Rule 30 seven-four staircase entropy and thirteen-four compiler

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3480 replaced a rectangular Rule 30 cone boundary by its physically
reachable ternary staircase and proved

```text
eta<=log_2(117/64)<7/8.                               (1)
```

The present theorem pushes the exact factor language eight levels farther.
A positive integer Collatz certificate on the height-21 overlap graph proves

```text
boxed: eta<=log_2(7/4)<13/16.                         (2)
```

The new exponent permits a taller one-word macrostep.  In the same fully
charged uniform word-RAM model, the main lookup tariff becomes

```text
boxed: Q(n)=(13/4+o(1)) n^2/log_2(n)^2.              (3)
```

This is a better upper-bound compiler, not a computational lower bound.  It
does not prove Rule 30 Prize 3, either of the other Rule 30 prizes, or
literature novelty for the staircase language.

## 1. Inheritance and typed scope

Use THM-3480's exact three-state right-to-left Mealy quotient

```text
A=00,       B=01,       C=1*,                        (4)
```

and its height-`h` serial cascade.  The physical zero right tail initializes
the cascade at `A^h`; write

```text
S_h=Orb_(Phi_0,Phi_1)(A^h),       N_h=|S_h|.         (5)
```

Height concatenation embeds `S_(h+k)` in `S_h x S_k`, so
`log_2 N_h` is subadditive and

```text
eta=lim_(h->infinity) log_2(N_h)/h                   (6)
```

exists.  Every contiguous length-`r` factor of a state in `S_h` lies in
`S_r`.  This factorial-language implication is the only passage from a
physical staircase to the finite overlap graph below.

The inheritance controls are deliberately retained:

- the closest proved mechanism is THM-3480's `R_12/R_13` Collatz certificate;
- the hostile is that every state through height 13 is Myhill--Nerode
  distinguishable and every fixed-state input/output chunk map is a
  permutation;
- the corrected near miss is that a path in a factor SFT need not be a
  globally reachable staircase; and
- the required sidecars are the zero right tail, scan orientation, ordered
  height coordinate, and physical state cocycle.

## 2. The exact height-21 overlap graph

Put `R_h=S_h` as a factorial language of ternary words.  Exact breadth-first
orbit construction extends THM-3480's count table to

```text
h       14     15     16      17       18       19       20       21
N_h  18619  32885  57741  100901   175680   304714   526563   906525. (7)
```

Build `G_21` with

```text
vertices: R_20,
edge u->v: a word r in R_21 has prefix u and suffix v. (8)
```

Thus `G_21` has exactly

```text
526563 vertices and 906525 edges.                    (9)
```

The adjacency matrix `A_21` is zero-one.  Indeed a length-21 word is uniquely
determined by its length-20 prefix and final ternary symbol, so two different
words cannot produce the same ordered prefix/suffix pair.

For a reproducible ordering, encode a state `q=(q_0,...,q_(h-1))` by two
binary masks

```text
E_i=[q_i!=A],       C_i=[q_i=C],
code_h(q)=sum_i E_i 2^i+2^h sum_i C_i 2^i,           (10)
```

and sort by increasing `code_20`.

## 3. The exact `4 A v <= 7 v` certificate

Starting with the all-one vector on the ordered vertices, synchronously
iterate the whole vector by

```text
v_i^(r+1)=max(v_i^r, ceil(4 sum_(i->j)v_j^r/7)).     (11)
```

The first fixed-point check succeeds on synchronous iteration `5789`
(`5788` strict vector updates).
The resulting vector satisfies

```text
min(v)=11,
max(v)=2303,
sum(v)=253229011,
-6<=4(A_21 v)_i-7v_i<=0.                            (12)
```

In particular,

```text
boxed: 4 A_21 v <= 7v.                               (13)
```

The canonical certificate payload consists of the ascending pairs

```text
(code_20(q) as unsigned 64-bit little endian,
 v_q       as unsigned 64-bit little endian).        (14)
```

Its SHA-256 digest is

```text
cf7adb12a00c294013135ebd19e606c0ff553eb2f27d65f300fbfbae45dc284a. (15)
```

All certificate weights and row sums are far below signed 64-bit overflow.
The companion reconstructs the physical orbit, graph, vector, and digest; the
digest is a commitment, not a substitute for checking (13).

### 3.1 Universal entropy consequence

For `h>=20`, a physical word in `S_h` injects into its path of consecutive
length-20 factors in `G_21`.  With `n=h-20`, positivity and (13) give

```text
N_h
 <=1^T A_21^n 1
 <=(sum(v)/min(v))(7/4)^n
 <=318(7/4)^h.                                       (16)
```

The last prefactor is an exact cleared-integer comparison:

```text
253229011*4^20
 =278428242084716609536
 <279113347509046779498
 =318*11*7^20.                                       (17)
```

Taking logarithms, dividing by `h`, and using (6) proves the first inequality
in (2).  The clean dyadic coarsening is

```text
7^16=33232930569601
    <35184372088832=2^45,                            (18)
```

so `(7/4)^16<2^13` and `log_2(7/4)<13/16`.

Equation (16) remains an upper envelope.  No path-surjectivity, entropy
equality, or lower estimate is inferred from `G_21`.

## 4. The charged `13/4` compiler

Retain exactly THM-3480's model: a uniform random-access word-RAM with

```text
w=ceil(log_2(n+2)),                                  (19)
```

unit-cost Boolean word operations, shifts, addition, and addressed word
loads/stores; a constant program; no advice or oracle; and all preprocessing
charged.

For `h=16m>=32`, equations (16) and (18), together with `318<2^9`, give

```text
N_h<2^(13m+9).                                       (20)
```

Put

```text
d=ceil(log_2 w),       L=w-8d-8,
m=floor((L-9)/26),     h=16m,
b=floor((L-9)/2).                                    (21)
```

For all sufficiently large inputs, `m>=2`.  The reachable-state ID and fresh
input chunk fit in one word because

```text
13m+9+b<=L.                                          (22)
```

The deterministic direct ternary marker array remains charged and lower
order.  The exact block inequality

```text
3^16=43046721<67108864=2^26                         (23)
```

gives

```text
3^h<2^(26m)<=2^(L-9).                               (24)
```

The diagonal table has fewer than

```text
N_h 2^b<2^L                                         (25)
```

entries.  Consequently the direct marker, dense state construction, and
complete chunk table have the same lower-order charged bounds as in THM-3480;
in particular they use `o(n/w)` words and operations before the main scan.
The active packed row still costs `O(n/w)` words.

At each macro time, scan the explicitly zero-padded active row from right to
left, starting at `A^h`.  The exact staircase table advances it by `h` rows
while emitting `b` new bottom-row bits per call.  Since

```text
h=(8/13+o(1))w,
b=(1/2+o(1))w,
hb=(4/13+o(1))w^2,                                  (26)
```

the sum of active row widths over `n/h` macrosteps gives

```text
Q(n)
 =n^2/(hb)+O(n/b+n/h)
 =(13/4+o(1))n^2/w^2.                               (27)
```

Table construction and the fewer-than-`h` final ordinary packed updates do
not change the bounds

```text
T(n)=O(n^2/log^2 n),       S(n)=O(n/log n).          (28)
```

Thus (27) improves THM-3480's exact lookup constant from `7/2` to `13/4`, a
factor `14/13`.  It does not change the asymptotic exponent.

## 5. A finite-only lower-bound signal

There is one exact positive signal on the opposite side of the entropy
question.  In top-to-bottom height order, put

```text
U=CABC,       V=BCBC.                                (29)
```

For every `1<=m<=5`, exact orbit membership shows

```text
{U,V}^m subset S_(4m).                               (30)
```

The companion checks all

```text
2+4+8+16+32=62                                      (31)
```

words.  If (30) held for every `m`, it would prove `eta>=1/4`.  No induction,
synchronizing encoder, or all-height reachability proof is known, however.
Accordingly (30) is **FINITE-EXACT ONLY** and supplies no lower entropy bound.
Local compatibility in a factor graph is not enough to promote it.

## 6. Preservation, loss, and prize boundary

| source -> target | preserves | destroys / admits | required sidecar |
|---|---|---|---|
| physical orbit -> `R_21` | every zero-tail height-21 state | its driving word | zero-tail calibration |
| staircase -> `G_21` path | every consecutive length-21 factor | global reachability; admits spurious paths | exact state cocycle |
| `A_21` -> Collatz weight | a certified exponential upper rate | exact physical count and lower entropy | positive vector and ordering |
| state/chunk -> table entry | exact bottom chunk and successor state | intermediate spacetime rows | scan direction, `h,b`, padding |
| finite block test -> `{U,V}^m` | all 62 witnesses through `m=5` | every larger height | an all-height encoder or induction |

The full-permutation hostile remains decisive.  For each fixed staircase
state, the `b`-bit input-to-output chunk map is a permutation, so this compiler
compresses only the reachable diagonal boundary, never the `b` fresh data
bits.  The improved entropy upper bound therefore supplies neither a
single-seed lower bound nor a shortcut below (28).

Prize 3 would require a lower bound in a precisely fixed uniform computation
model, or a faster single-seed algorithm.  An improved table-lookup upper
constant proves neither.  Center nonperiodicity and limiting balance are also
untouched.

## 7. Exact verification

Run

```bash
python3 04-computation/rule30_staircase_entropy_seven_four_thm3491.py
python3 -O 04-computation/rule30_staircase_entropy_seven_four_thm3491.py
```

The companion has no assertion-dependent gates.  It checks:

1. the two-mask update against an independent serial ternary Mealy update;
2. all reachable counts in (7), the `G_21` prefix/suffix construction, and
   zero-one adjacency;
3. the deterministic recurrence (11), every row of (13), the statistics
   (12), and certificate digest (15);
4. every cleared integer inequality in (17)--(25), including floor-sensitive
   word-budget samples; and
5. all 62 finite block-code memberships in (30).

The finite computations certify the displayed graph and vector.  The
factor-language injection, positive-vector estimate, and charged compiler
analysis prove the all-height and all-`n` consequences.
