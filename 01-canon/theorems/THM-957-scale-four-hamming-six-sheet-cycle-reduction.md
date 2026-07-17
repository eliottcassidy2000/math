---
id: THM-957
title: Scale-four Hamming-six sheet-cycle and two-triangle reduction
status: PROVED STRUCTURAL + FINITE-EXACT — the primitive proper AP-centred c=4 Hamming-six face is empty: the 64 all-order-four presentations give 256 unit contexts, and their complete no-height-cutoff metric recursion visits 166,976,181 candidate edges and finds zero covers
source: codex-2026-07-17-S60 scale-four sheet-cycle audit
depends_on: [THM-765, THM-810, THM-815, THM-823, THM-857, THM-859, THM-860, THM-861, THM-862]
related: [THM-856, HYP-6820]
verification:
  - 04-computation/lrc13_scale_four_hamming_six_sheet_cycle_codex_S60.py
  - 05-knowledge/results/lrc13_scale_four_hamming_six_sheet_cycle_codex_S60.out
  - 04-computation/lrc13_scale_four_hamming_six_terminal_scout_codex_S60.cpp
  - 04-computation/lrc13_scale_four_hamming_six_terminal_combine_codex_S60.py
  - 05-knowledge/results/lrc13_scale_four_hamming_six_terminal_census_codex_S60.out
---

# THM-957 — the scale-four sheet bank is a cycle plus two triangles

Let `R subset F_13^*` have six elements, put `P=[12] minus R`, and consider
the primitive proper AP-centred packet

```text
A=4P union {w_r:r in R},
w_r=4r (mod 13),                    w_r>0,       w_r!=4r. (1)
```

Assume `M(A)<=1/13`.  For each replacement put

```text
D_r=4/gcd(4,w_r) in {1,2,4},
e_r=(D_r w_r/4) mod D_r in (Z/D_r Z)^*.               (2)
```

THM-860 gives both the hereditary leave-one-out law

```text
lcm(D_s:s!=r)=4                                      (3)
```

and common-sheet coverage at all six replacement owners.  The purpose of
this theorem is to classify that finite sheet condition exactly before any
metric heights are chosen.

## Theorem

### A. Every surviving colour has effective order four

The complete scale-four common-sheet bank is

```text
effective orders              presentations          unit contexts
4^6                                  64                    256
all mixed order words                 0                      0. (4)
```

Each presentation has exactly four unit words.  In particular, a
hypothetical packet in (1) satisfying `M(A)<=1/13` has

```text
D_r=4 and w_r odd                         for every r in R. (5)
```

Its six metric rays are therefore

```text
w_r=u_r+52k_r,       k_r>=0,                         (6)
```

where `u_r in [1,51]` is the odd CRT representative

```text
u_r=4r (mod 13),                 u_r=e_r in {+1,-1} (mod 4).
                                                               (7)
```

No zero-height omission is needed in (6), since an odd `u_r` cannot equal
the even unchanged speed `4r`.

### B. The 64 supports are exactly the scale-two signed-cycle supports

For a six-set `R`, define a directed provider relation

```text
p -> o       iff       o/p in {2,-2}={2,11} in F_13. (8)
```

The 64 label sets surviving (4) are exactly the 64 supports from THM-860(D)
on which (8) is one directed Hamiltonian six-cycle.

There is more structure than that sparse cycle.  Colour every remaining
unordered pair by its ratio.  Put

```text
V={3,4,-4,-3}={3,4,9,10},
{p,o} a variable edge       iff p/o in V.              (9)
```

On every survivor:

```text
fixed-provider relation (8)       = directed C6,
variable graph (9)                = K3 disjoint union K3,
remaining unordered pairs         = 3K2.               (10)
```

Thus the edge-coloured `K6` splits into a directed six-cycle, two triangles,
and a perfect matching.  The graph of all pairs not in the directed cycle is
a triangular prism.  The coverage-bearing union of the cycle pairs and the
two triangles is `K6` minus a perfect matching, namely the octahedral graph.

### C. The unit fibre is an affine rank-four two-triangle code

Write

```text
x_r=+1 if e_r=1 (mod 4),            x_r=-1 if e_r=-1 (mod 4),
epsilon(a)=+1 for a in {3,4},       epsilon(a)=-1 for a in {9,10}.
                                                               (11)
```

At an owner `o`, let `p,q` be its two neighbours in the variable graph.
The exact four-sheet coverage condition at `o` is

```text
x_p x_q = -epsilon(p/o)epsilon(q/o).                  (12)
```

The owner is the third vertex of its variable triangle, so its equation is
carried by the opposite edge `{p,q}`.  The six owners index the six triangle
edges exactly once.  On each triangle the product of the three right-hand
sides in (12) is `+1`.  Hence each triangle system has rank two and one free
sign.  The disjoint union has rank four, two free signs, and exactly

```text
2^2=4                                                   (13)
```

unit words.  Equations (11)--(13) are an exact description of the literal
four-sheet fibre, not merely a necessary parity test.

### D. Sheet actions, tournament audit, and the information they destroy

Multiplying all labels by an element of `F_13^*`, while transporting their
attached units, preserves the sheet bank.  Global unit reflection

```text
J:x_r -> -x_r                                           (14)
```

also preserves it.  The 256 contexts have orbit distributions

```text
under F_13^*:                 20 orbits of size 12, 4 of size 4;
under F_13^* x <J>:           10 orbits of size 24, 2 of size 8. (15)
```

These are sheet actions, not metric quotients.  On `R={1,2,3,4,5,6}`, the
least packet for the valid word `(1,1,1,-1,1,1)` is

```text
Q=(3,17,21,25,28,32,33,36,37,40,44,48),
M(Q)=7/31,                   witnesses 8/31,23/31.      (16)
```

Its global unit reflection is

```text
JQ=(7,11,28,29,32,36,40,43,44,47,48,51),
M(JQ)=14/79,                witnesses 39/79,40/79.     (17)
```

Thus even the global flip in (14) changes literal metric geometry.

For Tournament Analysis, take (8) as the pairwise observable.  Missing pairs
are ties.  Their graph is the triangular prism from (10); choose its
lexicographically first Hamiltonian path and orient every tie by that total
path order.  Across all 64 supports, the completed tournaments have

```text
score multiset             (1,2,2,3,3,4),
directed triangles         6,
SCC sizes                  (6),
sparse-cycle edge flips    {2:8,3:52,4:4},
Hamiltonian-path counts    {29:32,31:20,37:12},
joint fingerprints         5.                           (18)
```

The completion is a useful audit but not the proof object.  It forgets which
ties are variable-triangle edges and which are zero matching edges, and it
does not see the four affine unit words.  The exact pair carrier is the
**edge-coloured** `K6`; the exact coverage carrier is better described by its
24 owner-sheet obligations and six provider incidence masks.

### E. Exact first metric layer

For a surviving support let `L_0(P)` be the longest strict-safe component of
its unscaled six-speed retained core.  Since

```text
L_0(4P)=L_0(P)/4,                                      (19)
```

the THM-815 root real cap with six replacements remaining is

```text
B_root=132/[13 L_0(4P)]=528/[13 L_0(P)].               (20)
```

Over the 64 roots it ranges exactly over

```text
528 <= B_root <= 1440.                                 (21)
```

Intersecting every ray (6) with its exact cap gives

```text
256 logical root contexts,
25,132 cap-admissible first edges,
60 <= first edges per context <=168.                   (22)
```

The first literal child geometry depends only on `(R,x_1)`.  There are

```text
12,566 geometric keys, each with exactly two logical lanes. (23)
```

The factor two has a structural explanation: fixing `x_1` fixes the first
label's unit, while the two-triangle code leaves exactly two compatible unit
words.  Those two words have different future ray systems, so (23) is an
exact geometry cache and not a quotient of terminal proof states.

### F. Complete exact terminal closure

Run the exact strict-safe-component recursion on all 256 contexts, preserving
the labelled step-52 rays and numerical insertion order.  Its node vector is

```text
depth                 0          1           2             3        4    5  6
nodes               256     25,132   2,577,024   163,876,444  496,938  643  0.
                                                               (23a)
```

The complete proof accounting is

```text
candidate edges       166,976,181
cap-dead nodes         (0,0,1,163,695,372,496,386,643,0)
full-safe-tooth        (0,0,0,156,889,649,0,0,0)
streaming-cap          (0,0,1,6,805,723,496,386,643,0)
covering terminals     0
loose terminals        0
maximum cap            7,665
maximum candidate      7,659.                           (23b)
```

Every cap-admissible branch is certified dead by depth five, so no branch
reaches a six-replacement terminal.  Since THM-815's cap contains every
hypothetical continuation with maximum at most `1/13`, (23a)--(23b) prove

```text
M(A)>1/13                                                 (23c)
```

for every primitive proper packet in (1).  Thus the complete `c=4` face is
empty.

## Proof

### 1. Exhaustion of effective-order states

At scale four the possible effective orders and units are exactly

```text
(D,e)=(1,0),(2,1),(4,1),(4,3).                         (24)
```

Equation (3) is equivalent to requiring at least two order-four colours.
Consequently each label set has

```text
sum_(j=2)^6 binom(6,j) 2^j 2^(6-j)=3,648              (25)
```

order/unit state words.  There are `binom(12,6)=924` label sets, hence

```text
924*3,648=3,370,752                                    (26)
```

candidate contexts.

For a state `(D,e)`, provider `r`, and owner `o`, let `u` be the CRT
representative

```text
u=Dr (mod 13),                         u=e (mod D),
```

and define its literal four-sheet mask by

```text
E_(r,D,e)(o)={ell in Z/4Z:
 -D < [u(o^(-1)+13ell)]_(13D) <= D}.                  (27)
```

The verifier packs the four bits for each of the six owners into one 24-bit
integer and tests whether the six provider vectors have bitwise union
`2^24-1`.  Equations (24)--(27) describe every possible common-sheet context,
so this is an exhaustive finite proof.  The result is exactly (4), proving
Part A.

### 2. The local mask grammar and the edge-coloured carrier

In an owner-normalized sheet gauge, direct CRT reduction of (27) gives the
following support table:

```text
D=1:   owner ratio 1 supplies all four sheets; every other ratio supplies 0;
D=2:   ratios 1,6,7 supply two sheets; every other ratio supplies 0;
D=4:   ratio 1 supplies the self sheet;
       ratios 6,7 supply one fixed sheet;
       ratios 3,4 and 9,10 supply the two unit-sensitive sheets;
       ratios 2,5,8,11,12 supply 0.                    (28)
```

The all-order-four output of Part A therefore has, at each owner, one self
provider, one incoming provider from (8), and two providers from (9).  The
remaining two directed incidences supply nothing.  The exact support scan
identifies the 64 label sets with THM-860's 64 signed-cycle supports and
checks (10) on each.  This proves Part B.

For the two variable providers in (28), the sheet chosen is encoded by
`epsilon(p/o)x_p`.  They fill the two remaining sheets exactly when their
encoded signs differ, which is (12).  The two-triangle structure then gives
the rank and consistency calculation in (13).  The verifier independently
checks the direct-mask/parity equivalence on all

```text
64*2^6=4,096                                             (29)
```

support/unit pairs.  This proves Part C.

### 3. Actions, metric guardrail, and first-layer replay

Literal substitution in (27) proves the sheet actions; exact orbit partition
gives (15).  The piecewise-linear breakpoint set

```text
{a/(2v), a/(u+v), a/|u-v|}
```

for speeds `u,v` computes (16)--(17) exactly and proves that (14) is not a
metric action.  Completing the sparse relation (8) by the declared tie path
and exhaustive six-vertex enumeration gives (18).

Finally, exact closed-danger reconstruction on the 64 retained cores gives
(21).  Enumerating every arithmetic progression (6) up to the integer part
of (20) gives (22).  Grouping first children by `(R,x_1)` gives (23).

### 4. Terminal recursion and shortcut audit

At a prefix with longest component `I`, `s` later combs, and next candidate
`x`, the engine uses the exact THM-815 cap

```text
x <= floor(22s/[13(13-2s)|I|]).                        (29a)
```

It materializes literal intersections of rational open intervals.  Beginning
at metric depth two, a retained component containing one complete safe tooth
of the candidate certifies a surviving interval and makes the later cap
strictly smaller than the current speed.  The streaming version stops as soon
as one newly formed component makes the cap smaller than the least member of
every remaining labelled ray.  These are noncoverage certificates, not state
quotients.  Their counts add exactly to the dead vector in (23b).

Four contiguous 64-context shards cover indices `[0,256)` exactly.  Their row
summaries reproduce (23a)--(23b), and

```text
candidate edges = sum_(d=1)^6 nodes[d].                (29b)
```

For referee hardening, context 255 was rerun through depth six with both
shortcuts disabled.  Its gated and ungated runs have identical node vector,
dead vector, candidate count, maxima, and zero-cover verdict; the ungated
shortcut counters are zero.  This verifies one full tree independently, not
the entire bank without shortcuts.  AddressSanitizer plus UndefinedBehaviorSanitizer
runs on contexts `0`, `134` (the maximum-cap row), and `255` are clean.  The
optimized source also compiles without warnings under
`-Wall -Wextra -pedantic`.

This proves Part F and completes the theorem. ∎

## Faithful carrier and scope guardrail

This computation challenges three tempting vertex choices:

1. **Runners alone** retain the label support but lose the unit code and all
   owner-sheet obligations.
2. **A completed tournament on runners** retains one arbitrary orientation
   per pair but erases the essential distinction between variable triangle
   ties and zero matching ties.
3. **Literal owner-sheet obligations as 24 vertices** retain common-sheet
   coverage exactly; provider states are incidence hyperedges.  After Part A,
   this hypergraph contracts faithfully to the directed cycle plus signed
   two-triangle code.

For metric continuation the vertices must change again: the faithful state is

```text
literal strict-safe components
 + six labelled step-52 progressions
 + the last inserted speed.                            (30)
```

The sheet carrier proves that only 256 scale-four contexts can matter; the
component/progression carrier then closes every one.  This proves only the
primitive proper AP-centred common-scale-four Hamming-six face.  It does not
close any `c>=5` face, the ramified Hamming-five bank,
non-AP-centred/deep-sheet branches, or global `n=12` sporadic emptiness.

THM-958, THM-960, THM-962, THM-963, THM-969, THM-970, and THM-974
subsequently answer `c=5,6,7,8,9,10,11`.  The next scale frontier is therefore
`c=12`.  The successful later proofs
confirm that the right first question is whether the owner-sheet obligation
hypergraph contracts to a small signed or algebraic carrier before any metric
recursion is launched.

## Verification

```bash
python3 04-computation/lrc13_scale_four_hamming_six_sheet_cycle_codex_S60.py |
  cmp - 05-knowledge/results/lrc13_scale_four_hamming_six_sheet_cycle_codex_S60.out
python3 -O 04-computation/lrc13_scale_four_hamming_six_sheet_cycle_codex_S60.py |
  cmp - 05-knowledge/results/lrc13_scale_four_hamming_six_sheet_cycle_codex_S60.out
clang++ -std=c++20 -O3 -Wall -Wextra -pedantic \
  04-computation/lrc13_scale_four_hamming_six_terminal_scout_codex_S60.cpp \
  -o /tmp/lrc13-c4-S60
/tmp/lrc13-c4-S60 --context-start 0 --context-limit 64 --depth 6 > /tmp/c4-0-64.out
/tmp/lrc13-c4-S60 --context-start 64 --context-limit 64 --depth 6 > /tmp/c4-64-128.out
/tmp/lrc13-c4-S60 --context-start 128 --context-limit 64 --depth 6 > /tmp/c4-128-192.out
/tmp/lrc13-c4-S60 --context-start 192 --context-limit 64 --depth 6 > /tmp/c4-192-256.out
/tmp/lrc13-c4-S60 --context-start 255 --context-limit 1 --depth 6 \
  --no-early-cap-gate > /tmp/c4-ungated-255.out
python3 04-computation/lrc13_scale_four_hamming_six_terminal_combine_codex_S60.py \
  /tmp/c4-0-64.out /tmp/c4-64-128.out \
  /tmp/c4-128-192.out /tmp/c4-192-256.out \
  --ungated /tmp/c4-ungated-255.out |
  cmp - 05-knowledge/results/lrc13_scale_four_hamming_six_terminal_census_codex_S60.out
clang++ -std=c++20 -O1 -g -Wall -Wextra -pedantic \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  04-computation/lrc13_scale_four_hamming_six_terminal_scout_codex_S60.cpp \
  -o /tmp/lrc13-c4-S60-sanitize
/tmp/lrc13-c4-S60-sanitize --context-start 0 --context-limit 1 --depth 6
/tmp/lrc13-c4-S60-sanitize --context-start 134 --context-limit 1 --depth 6
/tmp/lrc13-c4-S60-sanitize --context-start 255 --context-limit 1 --depth 6
```

Frozen integrity data:

```text
sheet source SHA-256
               6f84e5c1f640c82f2caf4249fef5a106388b907f70b81dfd914db7c898a91df1
sheet output SHA-256
               6510faa7aa48c3c00ef2895a92b9e417616ee62b50f6194acfbd37799b7cdf2f
terminal source SHA-256
               f67334d7db973d7e68240dbe29c4be149191cb4407cf1f58e70c2ff5e3ea199d
combiner SHA-256
               2cfa26540768b1b6abfb87c08cfede9a9d505079eac7a4649917d1615d2de01d
terminal output SHA-256
               c4551b10bd27e23ea7ec2a3977b770ff10d36dc6555b9b63461dff0f6949f410
context payload SHA-256
               7a9f8134aa9edad191b74a824ac24aa96b564a8c11a167fde5af4e26882f3346
parity payload FNV-64
               a039f50b7391414d
terminal row payload SHA-256
               09a46adda0ffce3647075c563ab83e7616a9e892a22d6f172ace591e5ed5b11e
ungated context-255 SHA-256
               ee3a4241f81a1b97cd450bb35c1dc333c6707fb8e17cb25a11703f5e9b6048e7. (31)
```
