# Exact low spectrum on the equal-response plane

Session 195, 2026-11-05. **PROVED elementary all-order classification; independently reviewed.** This supplements, rather than changes, `PROOF.md`. It concerns only its moving-root one-seam operation. The full `{7,21}` spectrum remains **OPEN**.

## The sharp statement

Use the rooted state `(h,e,R)` and operation `K(A,a;G,b)`, marked at the new b, from Theorem 1 of `PROOF.md`. The gadget G may have **arbitrary finite size**, and its choice can depend on the target and every preceding state. Include the empty word of operations.

**Theorem.** Suppose the physical initial rooted state satisfies

    R=h,    h>=5,    Delta=h-e>=2.

The set of Hamilton counts at most `3h+4` of all moving-root descendants is exactly

    {h, 2h-e, 3h}.                                     (1)

In particular neither `3h+2` nor `3h+4` is a descendant count.

A nontrivial strong tournament has positive endpoint counts at all vertices, so its state satisfies `h-e>=2`. The proof itself needs only the displayed numerical hypotheses, not strongness. The lower bound h>=5 is part of the theorem; no smaller-h classification is inferred by extrapolation.

This is stronger than the gcd obstruction. It permits primitive state content one. It does **not** apply to stationary-root moves, internal arc changes, replacing the graph by another witness, or the combined four-coordinate grammar.

## Proof

Write the fixed initial count as h throughout this proof, and put `C=3h+4`. For every moving-root step the inherited exact laws are

    H' = g H-s e+L R,
    e' = E H-delta e,
    R' = Q H+B R,

with `g=H(G)`, `s=s_b(G)`, `L=L_b(G)`, `Q=R_b(G)`, `B=P_(G-b)(2)`, and the other constants as in `PROOF.md`. Every step has `H'>=R` and `R'>=H+R>R`.

### 1. No nontransitive gadget can occur below the ceiling

A nontransitive gadget contains a directed triangle, so `g>=3`: the triangle's three Hamiltonian words inject into words of G by successively inserting missing vertices at deterministic legal positions. Also `L>=g` and `s<=g`.

If such a gadget is the **first** step, its count is

    H'=(g+L-s)h+s Delta.

The integer `g+L-s` is at least three. If it is at least four, then `H'>=4h>C` because h>=5. If it equals three, the inequalities force `g=L=s=3`, and then

    H'=3h+3 Delta>=3h+6>C.

After any first step, R is at least `2h`. It only increases afterward. Any later nontransitive gadget gives

    H' >= L R >= 3R >= 6h > C.

Hamilton counts never decrease, so a word whose terminal count is at most C cannot contain a nontransitive gadget anywhere.

### 2. Exact formulas for every transitive gadget

Let `G=TT_m` with source-first order `0,...,m-1`, and mark its vertex j. Its constants are

    g=1,  s=[j=0],  L=2^j,
    E=[j=m-1],  Q=2^(m-1-j),  B=2^(m-1).

The cover formulas follow by assigning the vertices preceding or following the marked vertex to the two ordered path slots; each allowed assignment has a unique increasing path in each slot.

If j=0, the Hamilton count changes by

    H' = H-e+R,    R'=2^(m-1)(H+R).

For m=1 the new endpoint count is `e'=H-e`; for m>=2 it is zero.

If j>=1, then

    H'=H+2^j R,
    R'=2^(m-1-j)H',
    e'=H if j=m-1, and 0 otherwise.

### 3. Exhaust every transitive word below C

A first **nonsource-marked** gadget has count `(1+2^j)h`. To stay below C it must have j=1, giving exactly `3h`.

- If m=2 its state is `(3h,h,3h)`. Every next source-marked gadget has count `5h>C`, and every next nonsource-marked gadget is larger still.
- If m>=3, its R is `3*2^(m-2)h>=6h`, so every further step has H greater than C.

Thus this branch contributes only `3h`.

A first **source-marked** gadget has count `h+Delta=2h-e`, which is below `3h`.

- If m>=2 its R is at least `4h>C`, so it has no further descendant below C.
- If m=1 its state is `(h+Delta,Delta,2h)`.

In the second case a next nonsource-marked gadget has count at least `5h+Delta>C`. A next source-marked gadget has count exactly `3h`.

- If this second gadget has m>=2, its R is at least `2(3h+Delta)>C`, so no continuation stays below C.
- If it is another singleton, the new state is `(3h,h,3h+Delta)`. Every next source-marked gadget has count `5h+Delta>C`; every nonsource-marked gadget is still larger.

The empty word contributes h. A singleton source-marked step realizes `2h-e`, and the transitive-edge sink-marked step realizes `3h`. Every remaining word has been excluded. This proves the equality in (1). QED.

## Primitive physical hostiles

A smallest-h example allowed by the theorem is the strong four-vertex graph

    rows=[12,9,2,4], root=2, (H,e,R)=(5,1,5).

Its exact descendant counts through19 are `{5,9,15}`; in particular17 and19
are excluded despite primitive state content. The following larger example
shows that the obstruction is not specific to a unit endpoint mark.


The five-vertex rooted graph

    rows=[26,16,3,6,12], root=2,
    (H,e,R)=(13,2,18)

is constructed and checked in `PROOF.md`. Apply the transitive-edge sink-marked operation once. The resulting strong graph is

    rows=[122,112,67,102,108,4,32], root=5,
    (H,e,R)=(49,13,49),    gcd(H,e,R)=1.

The theorem gives the exact descendant counts through151:

    {49,85,147}.

In particular 149 and151 are not moving-root descendants of this marked graph. They may of course be realized by independent tournament constructions. This is a primitive architecture obstruction, not a forbidden Hamiltonian-path value.

## Verification and limits

The [independent audit](../rooted-interfaces-audit/AUDIT.md) reviews the exact
low-spectrum proof separately from its finite controls.

`python3 verify_extensions.py` and its `python3 -O` counterpart check the primitive graph directly and run exact bounded moving-root reachability from every rooted labelled strong tournament of order at most five satisfying the theorem's hypotheses. There are304 such rooted graphs, reducing to four distinct states; the primitive49 state is the fifth. The70-gadget finite control bank is the five-gadget extension described in `FIVE_GADGET_EXTENSION.md`.

The controls verify both the ternary exclusions and the exact three-value sets. They are not the reason the theorem covers arbitrary gadget orders: that follows from the exhaustive structural argument above. Independent review separately checked the nontransitive and all transitive cases, including the m>=3, j=1 first-step case.
