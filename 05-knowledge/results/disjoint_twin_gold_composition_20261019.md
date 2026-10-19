# Conditional Gold certificates and composition of disjoint twin gadgets

Date: 2026-10-19. Status: proved for arbitrary finite ambient posets. No independence of extension orders or comparison outcomes is assumed. This is a conditional and multistage strengthening of [the one-relation twin Gold theorem](one_relation_twin_gold_theorem_20261019.md), not a stronger universal two-comparison margin or a proof of the general Gold conjecture.

## 1. Gadget and invariant skeleton

A one-relation twin gadget is an ordered triple T=(u,v,t) in a finite poset satisfying

    u∥v, D(u)=D(v), U(u)=U(v)∪{t}, t∉U(v).

Thus u<t and v∥t. Its three marked extension orders have counts

    A: u<v<t, B: u<t<v, C: v<u<t.

The proof in the preceding theorem gives A=C and B≤A. Crucially, both its label swap and its three-cycle injection change only the labels u,v,t at their three already occupied positions. Every other label and position stays fixed.

Define the T-skeleton of a full extension to be the full word with the three labels in T replaced by identical blank marks. It retains every outside label and its position, as well as the three occupied marked positions. Both maps preserve this skeleton exactly.

**Invariant-weight lemma.** Give each extension any nonnegative weight depending only on its T-skeleton, and suppose the total weight E is positive. Then the weighted marked-order totals still satisfy

    A=C, B≤A, E=2A+B.

Indeed the bijection pairs equal-weight words, and the injection preserves the weight of every input word. Equivalently, the two relations hold separately on every individual skeleton, and nonnegative mixtures preserve them.

This allows any conditioning event or nonnegative reweighting that is invariant under the skeleton maps. In particular, a history of comparisons whose endpoints all lie outside T has this property. It does not permit arbitrary reweighting that distinguishes the three internal marked orders of T.

## 2. Conditional two-comparison certificate

Under every measure allowed by the lemma, compare u with v, then v with t. Its weighted joint table is

    [[A,B],
     [A,0]], with B≤A.

Every realized path obeys the same weak Gold inequalities:

    E−t1−t2 ∈ {0,A−B,B} ≥ 0.

Each final cell is at most A≤E/2. Each first-answer branch is at most2E/3. The determined second answer in the v<u branch is valid under [the repository's Gold convention](gold_terminal_context_20261019.md): it leaves t2=t1=A. If a branch is a chain, the counts are1 as required.

The order-dual gadget has the same skeleton invariance and the dual comparison rule from the preceding theorem. Upper- and lower-defect gadgets can therefore be mixed in the results below.

## 3. Disjoint gadgets compose despite correlations

Suppose P has r pairwise vertex-disjoint gadget triples T1,...,Tr, each satisfying its full ambient neighborhood condition. Start with uniform linear extensions, or more generally any nonnegative measure invariant under the skeleton maps for every triple.

Process unused gadgets in any order, possibly chosen adaptively from the observed history. For each selected gadget perform its two fixed comparisons. At the start of a gadget's block, every earlier designated comparison has endpoints outside that triple. The resulting conditional measure therefore still satisfies its invariant-weight lemma. The choice of which unused triple to process next adds no extra restriction beyond the history on which that choice was based.

Let E_j be the total surviving weight after j completed two-comparison blocks, with E_0 the initial total weight. For uniform linear extensions these are survivor counts and E_0=e(P). All subsequent counts are interpreted as weights when a weighted initial measure is used. For every realized path and every j,

    E_j ≤ E_(j−1)/2.

Within each block separately, if M_j is the count after its first comparison, then

    M_j+E_j ≤ E_(j−1),
    M_j ≤ 2E_(j−1)/3.

Induction yields the correlation-robust bound

    E_r ≤ E_0/2^r

after the2r designated comparisons. It also gives the total intermediate-count bound

    sum_(j=1..r) (M_j+E_j) ≤ 2E_0(1−2^(−r)).

These conclusions do not assert that the gadget outcome probabilities multiply. They follow by applying the same conditional certificate successively.

## 4. Permitted adaptive outside interleaving

Before a triple's block, arbitrary additional comparisons may be interleaved provided their outcomes are invariant under that triple's skeleton maps; a sufficient condition is that both endpoints lie outside it. To retain the theorem for every still-unused gadget, these extra comparisons must preserve the skeleton invariance of every such gadget. A simple safe policy uses only vertices outside all still-unused triples, including vertices in already completed gadgets.

Such observations only restrict the surviving extension set, and preserve the required conditional identities for future blocks. Extra outside comparisons between the two designated comparisons of the current block can only reduce its eventual second count, so its Gold inequality and half-survivor bound remain valid relative to its block-entry count.

Thus the same2^(−r) bound holds after completing the2r designated comparisons plus any such extra observations. Without extras, this is literally a2r-comparison bound. No assertion is made that arbitrary early comparisons inside an unused gadget preserve its later certificate.

## 5. Why two first comparisons are not a Gold composition

Take two independent three-point posets, each consisting of u_i<t_i and an isolated v_i, and take their disjoint union. There are180 extensions. Comparing only u1 with v1, then u2 with v2, has a realized heavy/heavy path with

    t0=180, t1=120, t2=80.

Thus t1+t2=200>180, or in normalized form2/3+4/9=10/9. Even genuine independence does not make these two first comparisons into a Gold strategy. The valid composition uses both comparisons of each gadget and guarantees a Gold inequality within each two-comparison block.

## 6. A genuinely correlated control

Use vertices(u1,v1,t1,u2,v2,t2)=(0,1,2,3,4,5) with relations

    u1<t1,t2; v1<t2;
    u2<t1,t2; v2<t1.

Both triples are valid one-relation twin gadgets. In its60 extensions, each event B_i=(u_i<t_i<v_i) has positive probability, but the two events cannot coexist: B1 would force u2,v2<t1<v1<t2, contradicting B2. The gadget outcomes are therefore dependent. The sequential conditional proof still applies.

## 7. Exact controls and scope

[The self-contained checker](../experiments/verify_disjoint_twin_gold_composition_20261019.py) generates every naturally labeled poset of sizes3 through6. It tests every upper- or lower-defect gadget on every individual outside skeleton, and every ordered pair of vertex-disjoint gadgets on all realized four-answer paths. It also checks the independent naive-composition failure and the correlated example above. The [normal](disjoint_twin_gold_composition_controls_20261019.json) and [optimized](disjoint_twin_gold_composition_controls_optimized_20261019.json) receipts agree.

Arbitrary finite size and arbitrary r follow from the invariant injections and induction, not from finite enumeration. Reproduce from the repository root with `python3 05-knowledge/experiments/verify_disjoint_twin_gold_composition_20261019.py`; the checker prints a fresh receipt without modifying files. No priority claim is made.
