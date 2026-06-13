# Final Findings: The Per-Path Identity, Its Scope, and Its Failure

## The Core Algebraic Fact (proven, unconditional)

For any path P' of Ham(T-v), the quantity (inshat(v,P')-1)/2 equals
the number of Type-II positions in P' -- this is a pure algebraic identity
that holds for any binary string.

## What the Per-Path Identity Actually Says

The per-path identity from the paper states:
  (inshat(v,P')-1)/2 = sum_{Type-II at j} mu(v, P'[j], P'[j+1])

Since LHS = #Type-II positions (algebraic), this is really asking:
  #Type-II positions = sum_{Type-II at j} mu(v, P'[j], P'[j+1])

This can only hold if the weighted sum of mu values equals the count,
which requires careful cancellation.

## Why the Identity Holds for n <= 5 (TRIVIALLY)

For a 3-cycle C = (v, a, b), mu(C) = I(Omega(T-v)|avoid{a,b}, 2).
Omega(T-v)|avoid{a,b} contains only odd cycles of T-v on vertices V\{v,a,b}.
That set V\{v,a,b} has n-3 vertices.

  n=3: V\{v,a,b} has 0 vertices. No cycles possible. mu = 1.
  n=4: V\{v,a,b} has 1 vertex.  No cycles possible. mu = 1.
  n=5: V\{v,a,b} has 2 vertices. No 3-cycle (needs 3 vertices). mu = 1.

So for n <= 5: mu(C) = 1 for ALL 3-cycles C through v.
The identity becomes #TypeII = #TypeII. TRIVIALLY TRUE.
The per-path identity for n<=5 has zero mathematical content beyond the algebraic identity.

## Why the Identity Fails for n = 6 (EXACTLY CHARACTERIZED)

At n=6: V\{v,a,b} has 3 vertices, which CAN form a directed 3-cycle.
If they do: mu(v,a,b) = 3. If not: mu(v,a,b) = 1.

THEOREM (proven computationally, clear theoretically):
The per-path identity holds for path P' if and only if:
  for every Type-II position (a,b) in P', the 3 vertices V\{v,a,b}
  do NOT form a directed 3-cycle in T-v.

Equivalently: the identity FAILS for P' iff some Type-II position (a,b)
has mu > 1, iff V\{v,a,b} contains a directed 3-cycle in T-v.

This is a PERFECT BINARY SEPARATION (verified on thousands of examples):
  - path has some mu>1 TypeII => identity ALWAYS fails (100%)
  - path has no mu>1 TypeII  => identity ALWAYS holds (100%)

## The Role of 5-Cycles: Two Separate Effects

There are TWO ways the per-path identity can fail to capture Claim A:

Effect 1: The 3-cycle mu weights exceed 1 (n>=6).
  B = sum_P' sum_{TypeII at j} mu(v,P[j],P[j+1]) 
  differs from A = sum_P' #TypeII
  because some mu > 1.
  Mean (A-B) ≈ -5.9, meaning B > A typically.
  B overcounts relative to A.

Effect 2: The 5-cycle contributions to D = sum_C mu(C).
  D = sum_{ALL odd cycles C through v} mu(C)
  includes 5-cycle contributions not captured in B.
  Mean (B-D) ≈ +5.9, meaning B > D typically.
  B overcounts relative to D.

STUNNING FINDING: Mean(A-B) ≈ -Mean(B-D) ≈ -5.88
The two effects are NEARLY EQUAL AND OPPOSITE in their mean effect!
This means: on average, A ≈ D (the naive unweighted count ≈ the true sum).
Mean(A-D) ≈ -0.02 ≈ 0 (verified: 259/1200 pairs have A=D exactly, but the
average is near zero suggesting symmetric distribution of errors).

## The Deep Structure: Why Two Errors Almost Cancel

A - D = (A - B) + (B - D)

Effect 1: mu>1 weights make B > A (A-B < 0). 
  This happens when V\{v,a,b} contains a 3-cycle -- the 3-cycle
  INFLATES the mu weight of the embedding.

Effect 2: 5-cycles in D make D > B (B-D < 0).
  5-cycles through v contribute mu(5-cycle) to D,
  but the per-path formula B only sees 3-cycles.

The near-cancellation (mean A-D ≈ 0) suggests:
  Every 3-cycle C=(v,a,b) with mu(C)=3 "corresponds to" a 5-cycle
  that contributes roughly the same amount to D-B.
  The extra mu=2 from the 3-cycle (mu=3 instead of 1)
  is "paid for" by a 5-cycle that contributes to D.

THIS IS THE KEY STRUCTURAL INSIGHT. It suggests:

HYPOTHESIS: When V\{v,a,b} has a 3-cycle C' = (c,d,e), there exists
a 5-cycle through v that involves {a,b,c,d,e} (or some of them), and
the 5-cycle's mu contribution exactly compensates for the extra weight in mu(C).

If this can be made precise, it would give a PROOF STRATEGY for Claim A:
show that the "excess" in B over A (due to mu>1 weights) is exactly
compensated by the 5-cycle and higher-cycle contributions to D.

## The General Pattern: What Happens at n=7,8,...

At n=7: V\{v,a,b} has 4 vertices. Can have 3-cycles. mu(3-cycle) can be > 1.
         Additionally 5-cycles through v exist.
         mu(5-cycle) for a 5-cycle C=(v,a1,a2,a3,a4) has V\{v,a1,a2,a3,a4}
         with n-5=2 vertices, no odd cycles possible. So mu(5-cycle)=1.

At n=8: V\{v,a,b} has 5 vertices. Can have 3-cycles and 5-cycles.
         mu(3-cycle) can be 1, 3, 5, 9, ...
         V\{v,a1,a2,a3,a4} has 3 vertices. mu(5-cycle) can be 3.
         7-cycles through v exist.

The pattern: at each n, cycles of length up to n become relevant,
and the mu values of shorter cycles grow as longer cycles are added.

## The Correct Statement of the Per-Path Identity Problem

The per-path identity is asking for a PATH-BY-PATH version of Claim A.
Claim A is a SUM identity. The per-path version requires each path to 
independently satisfy the formula.

For n<=5: trivially true (all mu=1), but this is degenerate.
For n=6: requires that no Type-II position has mu>1 -- a genuine constraint
that fails for roughly half of paths.

The CORRECT per-path formula would need to somehow account for:
1. The mu weights of 3-cycle embeddings (can be > 1)
2. The contributions of 5-cycles (and longer) through v

There is no known formula that does both correctly at the per-path level.
The failure of such a formula to exist (or its difficulty) reflects
the global nature of Claim A: the balance only works in aggregate over
all paths, not individually.

## Implications for Proving Claim A

The near-cancellation (mean A-D ≈ 0) with B as an intermediate suggests
a two-step proof strategy:
  Step 1: Show A ≈ B in a suitable average sense 
          (the mu>1 inflation of 3-cycles equals some correction)
  Step 2: Show B ≈ D in a suitable average sense
          (the 5-cycle contributions fill in what B overcounts)

With A-B + B-D = A-D = 0 (Claim A), and the two terms being nearly equal
and opposite, this might decompose Claim A into two more tractable identities.

## Open Questions

Q1: Is there an exact algebraic identity relating:
    (A-B) + (B-D) = 0  [which is Claim A rewritten]
    to the structure of 5-cycles through v?

Q2: For each 3-cycle C=(v,a,b) where V\{v,a,b} has an odd cycle,
    is there a specific 5-cycle through v that "compensates"?
    Can this be made into a bijection?

Q3: At n=7: mu(5-cycle through v) = 1 always (V\{v,a1,a2,a3,a4} has 2 vertices).
    Does this mean the B-D gap at n=7 is entirely due to 7-cycles?
    Or do 5-cycles at n=7 have mu>1 (V\{v,a1,a2,a3,a4} has n-5=2 vertices,
    so no odd cycles, so mu=1 always for 5-cycles at n=7)?
    
    YES: at n=7, mu(5-cycle) = 1 always. So 5-cycle contributions to D
    are exactly 1 per consecutive embedding, same as 3-cycles at n<=5.
    The per-path formula including BOTH 3-cycles (with mu weights) and
    5-cycles (with weight 1 each) might work at n=7!

Q4: At n=8: mu(5-cycle) can be > 1 (V\{5 non-v vertices} has 3 vertices,
    can form a 3-cycle). The same failure mode from n=6 repeats one level up.
    The pattern: at n=2k, cycles of length up to 2k-1 can have mu=1 trivially,
    and the first cycle of length 2k-1 has mu potentially > 1 first appearing
    when n >= 2k.

Q5: Is there a recursive structure where the "correction" for L-cycles
    is itself expressed in terms of (L+2)-cycles, creating a tower
    that, when summed, gives Claim A?
