# A rank-staircase obstruction to preserving already-correct targets

## Exact five-state example

Let the quotient graph be a star with center H and leaves a,b,c,e. Distinct colors on the same spoke are distinct parallel edges. Let

- A1 = {Ha,Hb,Hc,He};
- A2 = {Ha,Hb,He};
- A3 = {Hb,Hc};
- A4 = {Hc,Ha}.

The individual graphic ranks are 4,3,2,2, meeting the staircase 4,3,2,1. A rainbow basis must occupy every spoke exactly once. Write a state as its leaf choices in color order A1,A2,A3,A4. The complete state space is

    (e,a,b,c)
    (a,e,b,c)
    (c,e,b,a)
    (b,e,c,a)
    (e,b,c,a).

In the displayed order, consecutive states differ in exactly two coordinates; every nonconsecutive pair differs in at least three coordinates. Thus the width-two reconfiguration graph is exactly a path on five vertices. The source (e,a,b,c) and target (e,b,c,a) have distance exactly four.

Their A1 element He is already correct. Keeping it fixed leaves precisely the source and target, which differ in three coordinates. Every allowed path must therefore temporarily remove an already-correct target element. This is stronger than a poor choice of installation order: no target-preserving reconfiguration exists for this pair.

Equivalently, once A1=He is fixed, the remaining three colors form the cyclic choices {a,b}, {b,c}, {c,a}; only the two cyclic assignments are bases. The spare He option of A2 and the additional A1 choices provide the four-move repair.

## What this establishes

The rank-staircase theorem cannot generally be strengthened to a strategy that permanently preserves every initially correct target coordinate. This explains why its recursive buffer construction may disturb a correct earlier color. The example does not give an exponential diameter lower bound: its diameter is only four. It also does not contradict special connected-contact theorems with stronger hypotheses.

## Integer-difference realization: bounded controls, not a conclusion

We tested the precise analogous obstruction in actual graceful interval trees on labels 0,...,m, for m<=8:

1. Freeze one arbitrary color edge and check the induced width-two graph. Every such nonempty restricted fiber was connected.
2. Freeze all but any three colors and check the induced width-two graph. Every such nonempty restricted fiber was connected. The m=8 run checked 32,192 frozen groups having at least two completions, among 4,020 total graceful trees.

These tests cover arbitrary frozen sets of the stated sizes, and hence cover upper-tail fibers with those further already-correct coordinates frozen. They do not prove that an arithmetic analogue is impossible. In particular, intermediate numbers of free/frozen colors outside these cases were not exhaustively tested, nor were larger label intervals.

A structural necessary condition narrows a future three-free-color example. After contracting the fixed forest, there are four components. If even one free color has graphic rank at least three, the earlier three-bad-color rank theorem gives width-two connectivity (with no connected colors needed). Thus a disconnected three-free-color example must have rank at most two in every free color. This rules out merely matching the generic cyclic pattern while leaving an additional rank-three color contact available.

The actual integer-difference target-preservation obstruction remains unresolved here. No claim about global one-color connectivity or fixed-tree-shape graceful existence follows.

## Exact controls

`check_target_monotonicity.py` lists the complete five-state abstract graph and searches one-fixed-edge interval fibers through m=8. `check_three_free.py` partitions the same interval tree sets by frozen coordinates and checks connectivity with three free colors. Both use the co-located standard-library `exact_utils.py`; receipts are `controls.jsonl` and `three_free_controls.jsonl`.
