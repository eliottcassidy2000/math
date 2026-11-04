# Targeted comparison with the two-cycle update at 2cfc7bba

## 1. A new positive rearming witness is independently verified

The new source provides the order-thirteen rows

    128,9,3,85,39,15,567,4478,639,1215,511,2047,3967.

They form a strong tournament T with H(T)=12875, sole exit0→7 and e_7(T−0)=2. The independent verifier checks the named arcs, strongness, Hamilton count, deletion endpoint count, the correct marked image Phi(T,0), and its two Hamiltonian cycles. The actual C2 child is strong, has order fifteen and H=38629.

Consequently 12875 now joins the previously verified 12877 as having a smaller-order digit-two ready realiser. Together with the already established weak-allowed preparation reduction, this supplies the stated weak-allowed +4 route for BOTH depth-two leaves. It does not supply a strong order-fifteen tournament with H=12875 and digit two, or the analogous strong order-fifteen 12877 preparation.

The distinction is checked concretely: putting a two-name transitive prefix before this thirteen-name witness yields an order-fifteen tournament with the same H and digit, but it is NOT strong. The strong C2 child above has the desired NEXT value 38629; it is not an H-neutral order-fifteen rearming of 12875.

## 2. Which bijection is actually available

The earlier corrected marked theorem is unchanged: Phi identifies strong positive-ready marked tournaments with strong sole-exit marked tournaments whose DELETED CORE is strong and whose Hamilton-cycle count equals the digit. Strongness of the cycle tournament alone does not imply that image condition. The previously exhibited four-name bad marking still has a nonstrong inverse.

For digits one and two and order n≥4, the separately proved canonical-choice repair gives exactly one admissible sole-exit mark. It therefore gives a genuine UNMARKED correspondence between digit-two ready classes and strong tournaments with exactly two Hamiltonian cycles AND a sole-exit vertex.

This explains why the new ready-class counts can agree with the numbers of two-cycle tournaments even when some cycle tournaments have two sole exits. Only their unique admissible exit is the canonical inverse mark; counting all sole-exit markings would overcount. The old broad marked-bijection sentence is not made true by the new data.

The regular five-name tournament with rows(6,12,24,17,3) is an independently checked boundary. It has exactly two Hamiltonian cycles, every outdegree is two, and their edge sets are disjoint. Thus it has no sole exit and is outside the ready correspondence. A separate five-name two-cycle tournament with an admissible sole exit is also retained; its inverse is strong and digit-two ready. This matches the source's distinction between two total two-cycle classes and one ready class at order five, without rerunning its full census.

## 3. Finite source statements and single-run status

The targeted incoming README reports exactly-two-cycle class counts 0,2,5,18,60,191 at orders 4–9. It reports that every such tournament at orders 6–9 has a sole exit and that each pair of cycles differs in three one-sided arcs. The order-nine file is explicitly named twohc_all_9.single_run.out and reports:

- 144 unique-cycle classes;
- 191 two-cycle classes;
-all 191 with three arcs in the first cycle absent from the second;
- 107 classes with one sole exit and 84 with two.

These are retained SOURCE-REPORTED finite results. This comparison does not rerun or independently certify that census. In particular our own normal/optimized spot checks do not convert the upstream single-run receipt into a two-mode or independently reproduced census.

The reported higher ready-class counts 589,1775,5257,15360 at orders 10–13 do not by themselves count ALL two-cycle tournaments at those orders. That identification additionally needs the unproved claim that every relevant two-cycle tournament has a sole exit, together with the ready-generator completeness inputs. No all-order Fibonacci/Lucas formula is established by the inspected source or by this note.

## 4. What the segment-swap observation means

There is an elementary conditional combinatorial fact: if two directed Hamiltonian cycles differ in exactly three arcs of either cycle, deleting those three arcs from the first leaves three nonempty directed path segments. Their shared internal arcs are also present in the second. The second cycle must reconnect the same three segments in the other cyclic order. Cutting the shared segment containing the marked vertex at that vertex expresses this as

    P S1 S2 R  ->  P S2 S1 R.

Thus three-edge difference implies a segment-swap representation. It does NOT prove that every tournament with exactly two Hamiltonian cycles has three-edge difference. The regular five-name example has difference five and is a concrete failure of an unqualified all-order version.

The additional short-segment restriction and the assertion that deleting either swapped segment leaves a unique-cycle tournament are separate observed properties of the tested ready classes. They are not consequences of the three-edge combinatorial fact alone. The inspected swapshape.py checks a proposed swap against the full second permutation and reports failures. glue.py derives segments using that pattern and counts cycles after deletion; its interpretation relies on the matching swap check. Neither script supplies an all-order structural or enumeration proof.

The source's finite gluing counts at orders 8–10 and its extension of the cycle-family conjecture are useful motivation. They must remain labeled finite evidence or conjecture. The canonical admissible-mark repair makes their relationship to ready classes precise but does not fill the missing all-order sole-exit, gluing, or enumeration arguments.

## 5. Reproduction and limits

Run python3 verify.py and python3 -O verify.py and compare outputs. The standard-library verifier checks only the new explicit12875 witness and its actual child/padded forms, the regular order-five exception, one admissible five-name example, and the earlier four-name marked counterexample. controls.json retains their physical names, rows, cycles and counts. All failures use explicit exceptions.

INCOMING.json pins the six targeted source files. The publication allowlist includes those exact snapshots, this comparison, the verifier and normal/optimized receipts. No broad census, unrelated source search, or installation was performed.
