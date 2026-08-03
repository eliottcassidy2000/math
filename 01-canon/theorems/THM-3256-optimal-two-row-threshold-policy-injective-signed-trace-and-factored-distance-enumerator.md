---
id: THM-3256
title: "Optimal two-row threshold policy, injective signed trace, and factored distance enumerator"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED in the fixed
  4,319-state THM-3238 physical bank.
  THM-3244's lawful rows 2 and 10 admit a memoryless axis-threshold selector
  of minimum worst-case depth 5 and minimum leaf count 15, with one tree
  attaining both optima.  In sharp contrast, every Q-monotone compiler that
  retains signed coordinate edits has 4,319 distinct continuation words:
  the signed Parikh vector is Q-n and reconstructs the starting state.
  Forgetting signs still leaves at least 1,536 Parikh classes, although the
  commutative distance enumerator has one short factored closed form.  This
  compiles one face's adaptive selector and separates enumeration complexity
  from continuation-state complexity; it proves no other face or FC(3).
source: root/policy-compiler-beach/2026-08-03
audit: >
  The exact companion pins and replays THM-3244; reconstructs the complete
  physical bank, two lawful direction masks and 407 exclusive chart states;
  exhausts all 16 primitive multiplicity thresholds by bit-set dynamic
  programming; proves the depth and leaf minima; and verifies the recovered
  tree on all 4,318 nonreset states.  A deterministic least-coordinate
  control supplies exact continuation-class and switch histograms.  The
  policy-independent signed and unsigned Parikh statements are checked
  against the complete bank, and the coordinatewise distance factors are
  multiplied independently.  Normal and optimized runs byte-match the
  frozen transcript; the source has no assertion node or floating literal.
  A separate tuple/frozenset classifier DP reproduced both policy optima,
  and an independent literal controller reproduced every continuation,
  Parikh, collision, switch and distance count.
depends_on:
  - THM-3244-unique-reset-exposure-deletion-graph-nonmorse-boundary
related:
  - THM-3238-complete-physical-product-gamma-bank-unique-reset-stitch
  - THM-3248-q4-paired-owner-stirling-compiler
script: 04-computation/fc3_two_row_decision_policy_compiler_beach_20260803.py
output: 05-knowledge/results/fc3_two_row_decision_policy_compiler_beach_20260803.out
script_sha256: 8910a0bc1a6bc10ac27742781142c28d89209b57bc0f6f5e461da314d0e0e197
output_sha256: b4171403580c05a34e33f5cc3b0d4b67402dea81986c0459a7f3b29afc19529a
hash_basis: LF-normalized bytes
---

# THM-3256 -- the chart is small, the signed operation trace is not

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED in the fixed
THM-3238 bank.**

THM-3244 proves that lawful rows 2 and 10 cover every nonreset state by a
strict Q-directed one-pole ascent, while no one row and no fixed positive
blend of those two rows does so.  This theorem asks two different
complexity questions:

1. how small can the adaptive row selector be; and
2. how much state survives when the selected physical edits are retained?

The answers lie at opposite extremes.

## 1. Universe and policy class

Represent a physical state by its pole-multiplicity vector

~~~text
n=(n1,...,n8),       0<=n<=(4,3,2,2,2,1,1,1).          (1)
~~~

Delete the empty physical multiset, leaving the 4,319 states of THM-3238.
Let Q be its unique reset.  THM-3244's direction masks define the sets C2
and C10 on which rows 2 and 10 respectively admit a strict Q-directed
one-pole ascent.  Their complete census is

~~~text
|C2\C10|=304,       |C10\C2|=103,
|C2 intersect C10|=3911.                                (2)
~~~

Only the 407 exclusive states constrain a deterministic chart classifier.
The declared policy class consists of binary trees whose internal tests are

~~~text
n_j <= h,             0<=h<capacity_j.                  (3)
~~~

There are exactly 16 nontrivial primitive tests in (3).  Leaves are labelled
2 or 10 and must classify the exclusive states correctly.  On overlap states
either label is lawful.

## 2. Simultaneously optimal threshold tree

There exists a tree with

~~~text
worst-case depth 5,             number of leaves 15.     (4)
~~~

No tree in the declared class has depth at most 4, and no such tree at any
depth has at most 14 leaves.  The same tree therefore attains both separate
optima.

One optimal policy is:

~~~text
if n8=0:
    if n7=0:
        if n3=0: choose row 10
        elif n6=0: choose row 10 if n4<=1 else row 2
        else: choose row 10 if n5=0 else row 2
    else: choose row 2
else:
    if n5<=1:
        if n1<=1 and n2=0: choose row 2 if n4<=1 else row 10
        else: choose row 10
    else:
        if n2<=1: choose row 10 if n4<=1 else row 2
        elif n2<=2: choose row 2 if n3<=1 else row 10
        else: choose row 2.                              (5)
~~~

On the complete nonreset bank, (5) assigns

~~~text
row 2:  1955 states,
row 10: 2363 states,
unavailable choices: 0.                                  (6)
~~~

Reapplying (5) after every edit therefore gives a fixed, memoryless
implementation of THM-3244's adaptive cover.  Each selected mask contains a
Q-directed coordinate, so choosing any such coordinate reaches Q in exact
edit distance.

### Exact minimality recurrence

For a classified subset S, let L_d(S) be the least number of leaves in a
threshold tree of depth at most d.  A pure subset costs one; a mixed subset
at depth zero is impossible; otherwise

~~~text
L_d(S)=min_(j,h) [
  L_(d-1)(S intersect {n_j<=h})
 +L_(d-1)(S intersect {n_j>h})].                         (7)
~~~

The exact bit-set recurrence over all 16 tests returns

~~~text
L_d(all)=infinity for d=0,1,2,3,4,
L_5(all)=15.                                             (8)
~~~

The well-founded unbounded recurrence also returns 15.  Equations (7)--(8)
prove (4) inside the declared policy class; they do not lower-bound arbitrary
arithmetic circuits or lookup encodings.

## 3. One reproducible operation compiler

For a deterministic control, after (5) chooses a row, choose the least-index
Q-directed coordinate in that row's exact mask.  Every one of the 4,318
nonreset routes reaches Q in exact edit distance.  Its chart-switch census is

| switches | states |
|---:|---:|
| 0 | 553 |
| 1 | 3158 |
| 2 | 387 |
| 3 | 220 |

The three-switch routes do not contradict THM-3244's sharp global two-switch
existence theorem: the present memoryless controller minimizes neither future
switches nor complete routes.

Minimize this no-input transducer by continuation equivalence, including Q's
empty word.  The exact response-state ladder is

| output retained | continuation classes |
|---|---:|
| chart label only | 296 |
| unsigned coordinate only | 4002 |
| chart plus unsigned coordinate | 4318 |
| signed coordinate edit | 4319 |

The chart-plus-coordinate words have exactly one collision, between the
physical states

~~~text
(2,3,3) and (1,1,2,3,3).                                (9)
~~~

Thus a 15-leaf row selector does not imply a small physical operation-word
automaton.

## 4. Policy-independent trace obstruction

The noncompression result does not depend on policy (5) or its tie-break.
Along any Q-monotone one-pole route from n to Q, coordinate j is edited
exactly |n_j-q_j| times.  If insertion and deletion signs are retained, its
signed action count is q_j-n_j.  Therefore

~~~text
Parikh(unsigned coordinate word)=|n-Q|,
Parikh(signed edit word)=Q-n.                            (10)
~~~

The second vector reconstructs n.  Hence every deterministic compiler
outputting a complete signed Q-monotone edit word has exactly

~~~text
4319 distinct continuation words, including Q's empty word.               (11)
~~~

No two physical states can be merged, regardless of how rows and lawful
directions are selected.

After forgetting signs, the absolute-deviation vector still has exactly

~~~text
4*4*3*2^5=1536                                           (12)
~~~

values on this bank.  Thus every unsigned-coordinate continuation quotient
has at least 1,536 classes.

## 5. Closed-form enumeration survives abelianization

Although the signed operation carrier is injective, the commutative distance
enumerator factors coordinatewise:

~~~text
(1+2z+z^2+z^3)(1+z+z^2+z^3)(1+z+z^2)
  *(1+2z)^2*(1+z)^3-z^8                                 (13)

=1+11z+55z^2+169z^3+365z^4+598z^5+775z^6+810z^7
 +685z^8+467z^9+250z^10+101z^11+28z^12+4z^13.
~~~

Each factor is one coordinate's absolute-deviation enumerator; z^8 removes
the forbidden empty physical state.  The coefficients sum to 4,319 and give
the exact Q-distance histogram.

Equation (13) separates two notions often conflated under closed form:

~~~text
commutative counting:       short factored polynomial;
ordered signed continuation:one state per physical start.                (14)
~~~

## 6. Typed Q4 connection and scope

The connection to THM-3248 is an operation-response contract:

| field | FC two-row compiler | Q4 compiler |
|---|---|---|
| source | multiplicity state plus two direction masks | walk resolvent W=N/D |
| compressed owner | row label | paired variables Y,Z,C |
| preserved | existence of a lawful strict ascent | exact diagonal after retaining N |
| destroyed | margin, direction and signed displacement | labelled walk placement if N is dropped |
| required sidecar | direction mask or signed edit | full walk numerator |
| hostile | Parikh injectivity | denominator-only value 1 versus true value 5 |

The common principle is not a shared tournament.  It is that a small
owner-selector can compile the control decision while the operation trace
remains large.  Exact continuation needs the response sidecar that the owner
quotient forgets.

This theorem concerns one support-(1,3), bank-I2 face.  It supplies no
face-independent symbolic derivation of (5), no 22-row optimum, no global
two-switch memory controller, no transport to the other 229 maintained
faces, and no proof of FC(3) or SFC(3).

## 7. Exact reproduction

Run

~~~text
python3 04-computation/fc3_two_row_decision_policy_compiler_beach_20260803.py
python3 -O 04-computation/fc3_two_row_decision_policy_compiler_beach_20260803.py
~~~

and compare LF-normalized bytes with the declared transcript.  The companion
pins and replays THM-3244, performs every finite optimization and route check
exactly, and contains no floating point, randomness, or
optimization-sensitive assertion.

QED.
