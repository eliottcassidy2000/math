# A physical grammar for transporting exact-target repair directions

**PROVED universal response law and composition rule; FINITE-EXACT context
head and actual lifts. This is not an all-target source theorem.**

## 1. The context is a typed rule, not just a count

Use a context of order at least2 and the one-sided rooted amalgamation and quantities in
[CONTEXT](CONTEXT.md). For two actual rooted graphs B,B' sharing the current
count and the first root-child target, write their profile difference as

    Delta z_B=(0,v,w,v).

Put P=P_(A-r)(2)>0, δ=G_A-h_A-q_A+e_A,
α=G_A+e_A-2q_A, β=F_A-P-s_A. Subtracting the exact gluing formula gives

    Delta H=αv-δw,
    Delta s=s_Aw,       Delta q=q_Aw,
    Delta F=Pv+(F_A-P)w.                            (1)

Consequently both observations are neutral precisely when

    [ α  -δ ] [v] = [0],
    [ P   β ] [w]   [0].                           (2)

Because P>0, a nonzero neutral direction has w≠0. Such a direction exists
iff αβ+δP=0, in which case its unique slope is

    ρ=v/w=(P+s_A-F_A)/P.                            (3)

The output difference is then

    (0,s_Aw,q_Aw,s_Aw),

so the outgoing slope is **σ=s_A/q_A**. The context is a rule ρ→σ with
an actual graph witness and magnitude multiplier q_A. This is an exact
classification by response data for every finite one-sided context, not a
structural enumeration of all possible responses.

Physical q differences are even. Since q_A is odd, this constraint survives.
For a supplied actual pair, the gluing creates another actual pair and
keeps the named root. It does not automatically preserve strongness or
unit readiness. Those are separate predicates to verify when consumed.

## 2. Composable words and concrete lifts

If a context accepts the outgoing slope of the previous context, apply it
to both actual output graphs. Equation(1) proves by induction that the pair
shares both observations at every stage. The magnitudes multiply by the
successive q_A values. Retain the graph word, vertex relabellings and root;
a walk of rational slopes alone is not a source witness.

The exhaustive context scan through order6 yields eight slopes. Every one
has a matching actual repair pair in the inherited strong rooted bank
through7. The table uses one explicit context at each slope; all rows,
input graphs, output graphs and intermediate profiles are in the receipts.

| Input slope ρ | Output slope σ | Context profile(h,s,q,F) | Common output H | Common next H |
|---|---|---|---:|---:|
|0|0|(1,0,1,2)|13|49|
|1|1|(1,1,1,1)|49|97|
|1/2|0|(1,0,1,2)|9|27|
|2/3|1/3|(5,1,3,9)|605|1861|
|5/6|2/3|(5,2,3,6)|3203|7671|
|5/8|0|(9,0,3,18)|845|2807|
|7/10|6/11|(23,6,11,36)|5031|13481|
|5/7|4/9|(17,4,9,28)|3917|11417|

Equal context profiles in rows1 and3 do not identify contexts: P differs
(2 versus4). This is a minimal witness that the four-coordinate profile
alone is not the full response of a context on the other side of gluing.

The word5/6→2/3→1/3 is realized literally. At its second output the two
profiles are

    (49821,129,387,107361),
    (49821,111,333,107343),

with common next target157053. The initial pair had H121; neither the
scalar count nor a formal kernel vector predicted this continuation.

The scan uses all labelled rooted contexts of orders2..5 and all labelled
order6 contexts rooted at5; the latter includes every rooted isomorphism
class. There are38172 sampled representations,4238 singular response
matrices and exactly eight slopes. These are representation counts, not
A000568 counts. The outgoing slopes1/3,4/9,6/11 have no accepting context
of orders2..6. Their availability at larger sizes is **OPEN** here. The
singleton identity context accepts every slope but supplies no growth or
change; it is outside formula(1) and this nontrivial-context classification.

## 3. An all-size conversion of the H9 repair

The inherited one-arc H9 pair has profiles(9,3,3,9) and(9,4,5,10), slope1/2.
Take the transitive three-context a→r→b, a→b, with cross edges a,b→B-r.
It has P4 and profile(1,0,1,2), hence accepts1/2 and emits0. The two new
profiles are

    (15,0,3,30),       (15,0,5,30).

This operation is equivalently an ordinary R birth followed by an ordinal
source singleton. It transports the old off-root flip to a pure deletion
repair without changing either observation of the resulting pair.
Additional ordinal singletons double F on both sides, giving

    (15,0,3,30·2^a),   (15,0,5,30·2^a),       a≥0.

Both have the exact next target15+30·2^a. This is an infinite physical
conversion rule. It does not yield primitive chain entry: after the next R,
the chain contents are respectively3 and5. This hostile control prevents
confusing neutral transport with removal of every entry obstruction.

## 4. What this changes in the research plan

The finite slope graph organizes repairs by the contexts that can consume
them. The all-order matrix criterion permits new contexts to be added
without rescanning every historical repair. A source-free slope or a
source-rich slope with no outgoing accepting context is an explicit
different failure mode. An exact-surplus selector still has to furnish
the initial actual pair at the prescribed target, and a chosen terminal
rule must independently prove the required endpoint predicate.

Reproduce `direction_probe.py` then `transport_words.py` in this directory.
