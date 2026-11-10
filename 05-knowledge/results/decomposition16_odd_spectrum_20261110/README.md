# Five-part decompositions and the remaining odd-spectrum induction

2026-11-10. **PROVED scoped reductions and generators; FINITE-EXACT
controls and inherited audited censuses; the full spectrum remains OPEN.**

The target is still

    {H(T): T a finite tournament} = positive odd integers minus {7,21}.

This session makes the needed induction more precise, proves several
permanent obstructions to overly simple steps, and recovers relevant old
decompositions without identifying unlike objects by their cardinalities.

## A sharp candidate induction threshold

The main next statement is now:

> **OPEN:** Every realizable odd h>=95 has one-vertex extension pairs
> (h,3h+2) and (h,3h+4), allowing different source tournaments for the
> two targets.

Together with the established allowed counts through287 and multiplication
by3 via an ordered join with C3, this would prove the full target by strong
induction. Every odd N>=289 has a unique expression 3h+delta with h odd
and delta in {0,2,4}; then 95<=h<N. The implication is proved; the universal
branch statement is not. The existing audited census verifies both
branches for all465 odd h from95 through1023.

The threshold cannot be lowered for this particular rule. Exactly these
branch failures occur among realizable h<1024, and each fails at **every
order**:

| Branch | Permanently absent parent→child pairs |
|---|---|
|3h+2|3→11, 15→47, 35→107, 93→281|
|3h+4|1→7, 17→55, 35→109, 41→127, 87→265|

The extra edge13→27 is also permanently absent. These are missing
**extension pairs**, not new missing path counts. For example11,27,47,
55,107,109,127,265,281 themselves all occur.

The mechanism is short. If C is the SCC containing the added vertex v,
then for some integer r,

    (H(T-v),H(T)) = r*(H(C-v),H(C)),   r | gcd(H(T-v),H(T)).

All the stated pairs are coprime. Therefore every other SCC is a
singleton, and arbitrary-order existence reduces to one strong child.
The classical strong minimum is375 at order12, so all the negative
targets above lie inside the previously complete strong-core range.
This promotes inherited finite nonexistence to all-order nonexistence.
[Proof, exact marked-class series, and threshold implication](induction/PROOF.md).
The primary [Busch minimum theorem](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
was checked directly; it is an inherited external input, not a new bound.

## An explicit source grammar and its first new obstruction

Every H3 source is TT_a -> C3 -> TT_b. Its complete one-vertex extension
language now has a two-integer compiler. A and B encode the two named
flank cuts, with the right cut complemented and reversed. Let g be the
inherited one-apex recurrence

    g(0)=g(1)=1,
    g(2q+1)=g(q),
    g(2q)=g(q)+q+(q mod2), q>=1.

For a new vertex beating exactly one triangle vertex, its count is

    C(A,B)=12AB+6A+3B+3g(A)+3g(B)
                         +2(A mod2)+(B mod2)-1.

Beating two triangle vertices exchanges A,B. The uniform cases have a
separate formula retaining module reentry. The lower bound
C>=12AB+6A+3B+5 gives a terminating inverse for any target, with actual
binary cut words as certificates. This is an all-order generator and
inverse for a restricted source family, not a conjectured scalar update.

Combining the exact H1 and H3 inverses shows that their one-birth images
first miss **55** beyond the global exceptions7,21. Thus every H55
tournament has every vertex-deletion count at least5. The bound is sharp:
row masks `(54,124,56,113,96,64,5)` give H55 and deletion counts
`(9,25,15,19,29,29,5)`.
[Complete formulas, proof, inverses and witness](attachment/PROOF.md).

The independent audit caught and repaired a draft claim that uniform
module contact makes each path visit the module once. A six-vertex
example has H27 while three times its collapsed count is15. The exact
omitted term is6A(B+1), or its swapped version. The missing information
is the exterior return route, not just which direction the new vertex
faces the module. [Audit and correction](history/ATTACHMENT_AUDIT.md),
[MISTAKE-756](../../../01-canon/MISTAKES.md#mistake-756-2026-11-10-codex-decomposition16-odd-spectrum-unpublished-attachment-compiler-draft----uniform-module-contact-was-mistaken-for-one-contiguous-visit).

## Retaining several named ports

For r named roots retain the whole deletion cube h_S, one count for each
subset of roots removed. If gcd_S h_S=1, removing only the unmarked
singleton SCCs leaves at most r marked SCCs and at most

    r+floor(log_5(H(T)^3))

vertices. Their marked-class generating series is a finite sum of
polynomials divided by (1-z)^(k+1), k<=r. With two roots, a nonzero
determinant H(T)H(T-x-y)-H(T-x)H(T-y) forces the roots into the same SCC.
Zero determinant does not imply separation: a strong example has cube
(15,3,5,1). [Joint-core proof and loss ledger](JOINT_CORES.md).

This keeps simultaneous deletions and names. It does not replace the
full joint attachment observer, named endpoint tensors, or F49/F343
and81/27 controls needed by the earlier native construction languages.

## Recovering the decomposition of sixteen

The exact historical observable behind **1–5–4–5–1** has not been located
in the bounded repository and chat recovery. It is not the ordinary
four-bit rank row1–4–6–4–1, nor a divisor-rank row for a number with
sixteen divisors. These distinctions do not refute a different intended
grading. A clarification remains useful; no arbitrary regrouping is
substituted for the original construction.

The actual recovered bridge is the four-bit fold

    (p1,p2,q,r) -> (p1+p2,q,r),

giving the twelve divisors of p²qr and the ten proper-divisor gradings
**3–4–3** and **4–2–4**. Its symmetry is exact, but it forgets a named
bit. Two genuine fixed-Hamiltonian-path tournaments in one folded fibre
have H9 and H5. The sixteen blue order-five presentations also have six
distinct H values, so no five-bin observer of that carrier can determine
H by the bin alone. [Recovered sources and literal controls](history/RECOVERY.md).

The recovered six arithmetic roots1,23,25,63,65,67 and the family
8*3^m+6*2^m-5 with two full return levels remain useful prior constructions.
Their addresses, clocks and joint ports matter; a shared decomposition
size does not supply a target-preserving map.

## Reproduction and evidence scope

Run from repository root:

    python3 -B 05-knowledge/results/decomposition16_odd_spectrum_20261110/reproduce.py

Use `--rerun` for fresh normal/optimized executions of every checker.
The default validates saved receipts, pinned inputs, script hashes,
manifest and local links. Small mathematical universes are fully
enumerated; the earlier order-ten positive census and complete strong
inventory are explicit inherited dependencies, not silently rerun here.

The independent negative-edge audit recounts263 selected strong classes
and2,171 vertex deletions using literal directed words from a different
audited inventory. The source compiler has1,032 all-cut controls and
an independent implementation. Two- and three-root kernel checks use
separate counting methods. [Session inheritance, incoming connections
and next obligations](SESSION.md).
