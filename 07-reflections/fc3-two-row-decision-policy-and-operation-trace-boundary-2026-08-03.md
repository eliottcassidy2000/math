# FC3 two-row decision policy and the operation-trace boundary

**Date:** 2026-08-03

**Status:** SYNTHESIS SUPPORT for hostile-audited THM-3256; the theorem file is the truth source

**Exact companion:** `04-computation/fc3_two_row_decision_policy_compiler_beach_20260803.py`

**Frozen transcript:** `05-knowledge/results/fc3_two_row_decision_policy_compiler_beach_20260803.out`

Canonical statement:
[THM-3256](../01-canon/theorems/THM-3256-optimal-two-row-threshold-policy-injective-signed-trace-and-factored-distance-enumerator.md).

## Inheritance pass

- **Closest proved mechanism:** THM-3244, `01-canon/theorems/THM-3244-unique-reset-exposure-deletion-graph-nonmorse-boundary.md`, gives the exact row-2/row-10 one-pole cover and its min-plus switch dynamic program on the full 4,319-state THM-3238 bank.
- **Canonical hostile example:** THM-3244's two complementary trap intervals show that no fixed positive blend of rows 2 and 10 replaces the selector.
- **Corrected near miss:** the globally exposing THM-3238 functional `H` has 32 nonreset one-pole local maxima.  Global exposure is not a local deletion policy.
- **Least-used relevant sidecar:** the two exact direction masks in the THM-3244 certificate, not merely the two chart-cover bits.
- **Concurrent results not used:** THM-3249 and THM-3254 are now proved and hostile-audited.  Neither is a dependency of this artifact.

The live concept board was:

1. the two-row cover;
2. an axis-threshold decision policy on pole multiplicities;
3. the min-plus switch recurrence;
4. an operation-response sequence compiler; and
5. the signed edit Parikh vector as a loss sidecar.

The relevant method cards were **Correct object before technique**, **Operation-response**, and **Controlled forgetting owes a sidecar**.

## 1. The exact chart policy is small

Write a physical state as its multiplicity vector

```text
n=(n1,...,n8),  0<=n<=(4,3,2,2,2,1,1,1).
```

Only 407 states constrain a deterministic row choice: the 304 states covered only by row 2 and the 103 states covered only by row 10.  At each of the other 3,911 nonreset states either row is lawful.  There are exactly 16 nontrivial primitive axis tests `nj<=t` on this integer box.

Exhaustive exact decision-tree dynamic programming gives both sharp statements:

```text
minimum worst-case depth = 5,
minimum number of leaves = 15.
```

The same tree attains both minima.  In particular no depth-at-most-four axis-threshold policy exists, and no axis-threshold policy of any depth has at most fourteen leaves.  One optimal policy is:

```text
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
        else: choose row 2
```

On all 4,318 nonreset states this assigns row 2 to 1,955 states and row 10 to 2,363 states, with zero unavailable choices.  Thus the state-dependent cover of THM-3244 can be compiled into a fixed, memoryless, five-comparison chart rule.  Reapplying the rule after every edit still selects a lawful strict ascent, so any allowed direction reaches the reset in exact edit distance.

This is a real compression of the **chart selector**.  It is not a scalar height, a fixed positive blend, or a closed form for the actual edit word.

## 2. Exact minimality mechanism

For a classified subset `S`, let `L_d(S)` be the fewest leaves in an axis-threshold tree of depth at most `d`.  A pure subset has cost one.  A mixed subset at depth zero is impossible.  Otherwise

```text
L_d(S)=min_(j,t) [L_(d-1)(S intersect {nj<=t})
                  +L_(d-1)(S intersect {nj>t})].
```

The exact bit-set recurrence exhausts all 16 primitive thresholds.  It returns impossible for depths zero through four and 15 at depth five.  The unbounded version of the same well-founded recurrence returns 15.  This proves the two optimality assertions inside the declared axis-threshold policy class; it is not a lower bound for arbitrary arithmetic circuits or arbitrary lookup encodings.

## 3. A deterministic operation compiler and its response-state growth

For a reproducible control, after the tree chooses a row, choose the least-index certified `Q`-directed coordinate.  The resulting 4,318 routes all reach `Q` in exact edit distance.  Their chart-switch histogram is

```text
switches  states
0         553
1         3158
2         387
3         220.
```

Thus a compact memoryless chart rule need not realize THM-3244's globally optimized two-switch bound.  The distinction is useful: THM-3244 minimizes separately over complete future routes, whereas the tree makes one local chart decision from the current multiplicities.

Minimizing a finite no-input sequence transducer by continuation equivalence means merging two states exactly when their complete remaining output words agree.  For this deterministic control, including the reset's empty word, the exact counts are

```text
output retained                    continuation classes
chart label only                   296
unsigned coordinate only          4002
chart plus unsigned coordinate     4318
signed edit                        4319
```

The chart-plus-coordinate words have only one collision:

```text
(2,3,3)  and  (1,1,2,3,3).
```

So the 15-leaf chart policy does not imply a small operation-word automaton.  Almost all compression disappears when the actual direction is restored.

## 4. Policy-independent Parikh obstruction

The strongest obstruction is elementary and does not depend on the least-coordinate tie-break.

Let `Q=(q1,...,q8)` be the reset multiplicity vector.  Along any `Q`-directed one-pole route from `n` to `Q`, coordinate `j` is edited exactly `|nj-qj|` times.  If insertion and deletion signs are retained, its signed action count is exactly `qj-nj`.  Therefore:

```text
Parikh(unsigned coordinate word) = (|n1-q1|,...,|n8-q8|),
Parikh(signed edit word)          = Q-n.
```

The signed vector reconstructs `n`.  Hence **every deterministic compiler that outputs the full signed monotone edit word has 4,319 distinct continuation words on this bank**, regardless of how it chooses rows or directions.  No continuation quotient can merge two physical states.

If signs are forgotten, the exact lower bound is still

```text
4*4*3*2^5 = 1536
```

distinct absolute-deviation vectors.  This is the clean information boundary: the selector bit is compressible, while the physical operation word carries an injective displacement sidecar.

There is nevertheless a compact **commutative** closed form.  The exact
distance enumerator of the nonempty physical bank, including the reset at
distance zero, is

```text
(1+2z+z^2+z^3)(1+z+z^2+z^3)(1+z+z^2)(1+2z)^2(1+z)^3-z^8

= 1+11z+55z^2+169z^3+365z^4+598z^5+775z^6+810z^7
  +685z^8+467z^9+250z^10+101z^11+28z^12+4z^13.
```

Each factor is one coordinate's absolute-deviation enumerator; `z^8` removes
the forbidden empty physical state.  Thus efficient factored enumeration and
small continuation state are different notions.  Abelianizing the operation
word gives a short product formula, while retaining signed order distinguishes
all 4,319 states.

## 5. Typed connection to the Q4 sequence compiler

The comparison with THM-3248, `01-canon/theorems/THM-3248-q4-paired-owner-stirling-compiler.md`, is typed rather than syntactic.

| Field | FC two-row policy | Q4 compiler |
|---|---|---|
| source | physical multiplicity state plus two response direction masks | quotient walk resolvent `W=N/D` |
| target | chart label, then a monotone edit word | finite-difference/Stirling contraction |
| map | optimal threshold tree, followed by a direction rule | paired coordinates `Y,Z,C`, then diagonal extraction |
| preserved predicate | selected row admits a strict lawful `Q`-directed ascent | exact path-cover coefficient |
| destroyed information | margins, alternative directions, signed displacement | labelled walk numerator if only `D` is retained |
| required sidecar | direction mask or signed edit word | full numerator `N` |
| hostile test | depth four / fourteen leaves; Parikh injectivity | determinant-only value `1` versus true value `5` at `r=1` |

The shared mechanism is not “both use two owners.”  It is sharper:

> a small owner/selector quotient can compile the control decision while the operation trace remains large; exact computation requires the response sidecar that the quotient forgets.

For Q4 the load-bearing sidecar is the walk numerator.  For FC it is the direction/sign data.  This is evidence for the operation-response method card and a counterindication to treating a compact selector as a compact sequence generator.

## 6. What this changes, and what remains open

### Direct progress

- The THM-3244 adaptive row selector has a fixed, exact, memoryless implementation with optimal axis-threshold complexity: depth five and fifteen leaves.
- The physical signed edit trace has a policy-independent injectivity obstruction, so a full continuation automaton cannot compress this 4,319-state bank at all.

### Indirect / niche progress

- The response-state ladder `296 -> 4002 -> 4318 -> 4319` quantifies exactly how each restored operation coordinate destroys chart-word compression.
- The tree provides a cheap implementation layer for future experiments; evaluating the 22 enormous response coordinates is unnecessary merely to choose between rows 2 and 10.

### Claims narrowed

- “Finite-state selector” is true but weak: every finite bank is finite-state, and the substantive invariant is the optimal policy representation.
- “Compact chart policy gives a compact closed-form route sequence” is false on the signed operation carrier.
- The three-switch canonical control does not refute THM-3244's sharp two-switch route optimum; the two statements optimize different objects.

### Honest frontier

1. Derive the 15-leaf tree symbolically from response-sign inequalities rather than from the finite census.
2. Determine the smallest decision policy on all 22 lawful rows, not just the distinguished pair.
3. Add one bit of previous-chart memory and ask for the smallest threshold controller that realizes the global two-switch optimum.
4. Test whether the other 229 maintained FC faces admit uniformly bounded chart trees, or whether policy depth grows with the face.
5. Aggregate operation words by commutative statistics when a sequence question does not need the full signed trace; the Parikh obstruction says exactly which information that aggregation retains.

No statement here proves the three-variable Factorial Conjecture or transfers to another face without a new exact carrier map.
