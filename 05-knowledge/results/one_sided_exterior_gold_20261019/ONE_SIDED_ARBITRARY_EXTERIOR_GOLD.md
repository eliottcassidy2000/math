# Strict Gold for all one-sided exteriors of the ten-vertex core

19 October 2026. Status: PROVED finite reduction with FINITE-EXACT complete profile certificates. Root reviewed and approved. This is a theorem for one specified core and arbitrary finite exterior sizes, not a general seven-thin Gold theorem or a morphology exclusion.

## 1. Precisely quantified result

Let P be any finite poset of incomparability degree at most 7. Suppose its full critical-root core at the oriented incomparable roots (0,3) is the labeled induced poset

    Q=(0,0,0,6,0,19,2,18,86,151).

Here the full core consists of the roots together with every vertex incomparable with at least one root; the tuple lists strict predecessor bitmasks on labels 0 through 9. Critical orientation means D(0) is contained in D(3) and U(3) is contained in U(0). The exterior splits as A below both roots and B above both roots, with A<B.

THEOREM.

1. If every vertex of B is comparable with every vertex of Q, then ask 8<9 first and 3<8 second, on both first branches. This is strict Gold for P. Every normalized Gold margin is at least

       5031/38734 > 1/8.

   The lower exterior A is arbitrary, including its size, internal order and allowed incidences with Q.

2. If every vertex of A is comparable with every vertex of Q, then ask 2<4 first and 1<4 second, on both first branches. This is strict Gold for P. Every normalized Gold margin is at least

       95/796 > 1/9.

   The upper exterior B is arbitrary under the same ambient hypotheses.

The constants are sharp for the respective SPECIFIED comparison trees, within their stated one-sided exterior classes. No optimality among other trees is claimed.

The passive opposite side is allowed to be large and non-chain. In the first case comparability implies Q<B, so that whole side is an ordinal terminal block. In the second case A<Q is an ordinal initial block. Its extension factor cancels from every core comparison probability. Merely having no active vertex on one side is exactly this passive condition; merely keeping that side small would not suffice.

## 2. Reduction to 46 and 48 actual models

Use the exact [boundary chain compression theorem](../boundary_chain_compression_20261019/BOUNDARY_CHAIN_COMPRESSION.md). The core incomparability degrees and exterior slacks are

    d_Q=(7,3,6,7,5,6,7,6,5,4),
    delta=(0,4,1,0,2,1,0,1,2,3).

A lower profile alpha is nonnegative, antitone on Q, zero at the roots and at every core vertex above a root, and bounded by delta. Its positive support has size at most 7. In this core these constraints are exactly

    alpha=(0,a,b,0,c,0,0,d,0,0),
    0<=a<=4, 0<=b<=1, 0<=c<=2,
    0<=d<=min(1,a,c),

with all variables integral. There are 46 such profiles: a=0 contributes 6, and each a=1,2,3,4 contributes 10. The zero value at label 8 also follows from 6<8 and delta(6)=0. The support condition is automatic because only four coordinates can be positive. Each model adds a lower chain of length at most 4, so has at most 14 vertices.

The upper profiles are exactly

    beta=(0,0,0,0,c,e,0,f,g,h),
    c in {0,1},
    c<=e<=1, c<=f<=1, c<=g<=2, f<=h<=3,

again integral. There are 48 profiles: c=0 contributes 42, and c=1 contributes 6. The constraint c<=1 follows from beta(4)<=beta(5) and delta(5)=1; the other inequalities are the relevant Q comparabilities. Positive support is at most five. These models add an upper chain of length at most 3, so have at most 13 vertices.

These are the COMPLETE profile universes, not sampled extensions or cutoffs on the original exterior size. Every listed profile gives an actual seven-thin chain model with the same full core. Conversely every exterior extension compresses to one listed profile. The verifier independently generates profiles from the general monotonicity/slack rules and checks exact equality with the two explicit parameter descriptions.

## 3. Correct mixture weights and transport of Gold margins

Consider the arbitrary-lower case. For alpha let m_alpha be the number of linear extensions of A inducing that profile, and let E_alpha be the total extension count of the corresponding core-plus-chain model. Then the ambient core law is the convex mixture of those model laws with weights

    w_alpha=m_alpha E_alpha / sum_gamma(m_gamma E_gamma).

The passive B contributes the same factor e(B) to numerator and denominator and cancels. In particular the weights are neither uniform over profiles nor m_alpha/e(A). The upper-side formula is identical with upper profiles. This exact compression preserves ALL joint core comparison counts simultaneously, so it applies to the four cells of a fixed two-question tree, not just its marginal probabilities.

If a model cell has count c and its first branch has count b, its normalized Gold margin is (E-b-c)/E. For a fixed tree each margin is affine in the normalized core law. A common lower bound over all 46 or 48 models therefore survives every actual exterior mixture. Thus the finite inequalities below establish the arbitrary-size theorem.

## 4. Exact finite certificates and sharp controls

For the lower tree (8<9, then 3<8), every listed profile satisfies all four inequalities

    38734 (E-b-c) >= 5031 E.

Equality is attained at

    alpha=(0,4,1,0,0,0,0,0,0,0).

Its total is 38734, with joint cells in order TT,TF,FT,FF equal to

    (11584,10535,13342,3273).

The four margins E-b-c are

    (5031,6080,8777,18846).

For the upper tree (2<4, then 1<4), every listed profile satisfies

    796 (E-b-c) >= 95 E.

Equality is attained at

    beta=(0,0,0,0,1,1,0,1,1,3).

Its total is 57312, its cells are

    (22730,5012,13738,15832),

and its margins are

    (6840,24558,14004,11910).

Both equality profiles are actual allowed ambient posets, establishing the claimed sharpness for these fixed trees.

The complete 94-row text table `one_sided_gold_profiles.tsv` records the profile, ambient order, total, four cell counts and minimum margin for every model. The standard-library verifier computes these counts directly by integer ideal dynamic programming with added comparisons, without using the boundary interleaving formula. In recurrence form, if s is the already-placed ideal, the count is the sum of the counts after adding each available unplaced vertex; the full ideal has count 1. Added comparisons may create a cycle, in which case no complete placement is available and the count is 0.

Run from this directory:

    python verify_one_sided_gold.py
    python -O verify_one_sided_gold.py
    python verify_one_sided_gold.py --table

The normal and optimized receipts agree on 665 active guards. The table SHA-256 is

    5f02d7ba29a5d6db397bc873bfa6862b1d55da6f72ec804d71aae46ed2eb0917.

The packet manifest fixes the script and receipt hashes too. The written compression proof justifies the unbounded exterior quantifiers; the complete exact finite certificate justifies these two specific fixed trees and margins.

## 5. Frontier

The arbitrary one-sided problem for this core is settled by the two fixed trees above. When BOTH exterior sides have active vertices, no one first core pair works universally: the separately certified 12- and 13-vertex attachments in [the adaptive selection note](../adaptive_gold_core_selection_20261019/ADAPTIVE_GOLD_CORE_SELECTION.md) already have disjoint balanced-core-pair lists. Their two-law segment has an adaptive strict-Gold repair, but no complete adaptive cover of all two-sided profiles has been proved here.

A global extension should retain the exact profile mixture normalization, the lower/upper multiplicity factorization, and the fixed-interface correlation constraints developed in the preceding packets. Passing from one-sided validity regions to all two-sided realizable parameters remains open.
