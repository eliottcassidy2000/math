# Audit status and reproduction

## Provenance and scope

The baseline is pinned in README.md and RESULT.md. The inherited theorem paths, including THM-4097/4111/4118/4128/4177, are recorded in RESULT.md. No new canonical theorem IDs are reserved by this worker.

The parent-appointed independent auditor has accepted on proof review:

* the common-child overlap graph and radius n-4;
* the hull/gap corollary;
* exact two-ear seam reciprocity, nested-cut lower bound, and seam lattice;
* the fixed-core triangle-rotation divisor;
* the arbitrary-parent SCC strong-completion criterion;
* the r-vertex and q-arc marked-boundary divisors, including the zero-walk caveat;
* the unique-pair/card reconstruction theorem, without a strongness hypothesis.

These are all-order proof claims, not conclusions inferred from enumeration.

## Finite audit ledger

### First and second phases

`audit.py` passes in normal and -O modes with byte-identical output. It checks all bases of orders2..4 and every two-ear attachment pair, every labelled strong parent through order6 for the constructive clone path, the complete labelled unit atlases at orders5/6, all H25/H27 two-edge bridges, and the substitution seam controls.

`phase2.py`, `fullstar.py`, and `deck_leaf_audit.py` each pass normal/-O byte comparison. The parent-appointed auditor independently reconstructed the complete strong order-seven deletion-card graph and the order-six repair census, and independently verified the H147 closed-orbit leaf and its four minimum three-arc repairs. Its artifacts are in the parent audit packet's `independent/` directory; they are not represented here as imports into the producer code.

### Symmetry and zero-rigidity continuation

* `symmetry_edit.py`: all strong quotient classes at orders3..5,500480 cut presentations,111 direct literal-child DP controls; normal/-O identical.
* `zero_escape.cpp`, k=6:35 classes/210 deletion-card orbits/860160 profiles/27525120 represented cuts; full -O3/-O0 output comparison passes.
* Parent-appointed independent selected k=6 check: quotientcode16, deletion0, all131072 literal cuts without orbit compression; target1049305833,parent28273779,gcd54,unique zero114689; four order18 children directly recounted. The auditor's O2/UBSan streams match. This is a selected-row independent check, not a full independent k=6 census.
* `zero_escape.cpp`, k=7:353 classes/2471 deletion-card orbits/40484864 profiles/2591031296 represented cuts, terminal PASS. **Producer FINITE-EXACT; the full census is not independently replicated.**
* `direct_child_audit.cpp`: selected Paley7 quotientcode85298; its order20 parent has H11160377792361, original order21 child has H123156134945235. All212 literal-child counts recover every one of the211 constant/linear/quadratic coefficients. The count recurrence is forward-extension DP, distinct from the producer's backward endpoint convolution.
* Selected k=7 -O3/-O0 replay uses the same Paley7 quotient. It is an optimization control, not an independent algorithm.
* `zero_closure_controls.py`: actual and formal2048-cut fields, all layer means and all pair-deletion lower bounds checked; normal/-O identical. The parent-appointed auditor independently confirms the formal perturbation and its nonrealizability firewall.
* `card_rigidity_audit.py`: all labelled quotients through4 plus all strong quotient classes5/6;116 quotient presentations,1575 deletion cards,46083 triple candidates,138249 pair candidates; normal/-O identical.

## Main replay

Run from this directory, in order:

```sh
python3 -B audit.py > audit.out
python3 -B -O audit.py > audit_optimized.out
cmp audit.out audit_optimized.out
python3 -B phase2.py > phase2.out
python3 -B -O phase2.py > phase2_optimized.out
cmp phase2.out phase2_optimized.out
python3 -B fullstar.py > fullstar.out
python3 -B -O fullstar.py > fullstar_optimized.out
cmp fullstar.out fullstar_optimized.out
python3 -B deck_leaf_audit.py > deck_leaf_audit.out
python3 -B -O deck_leaf_audit.py > deck_leaf_audit_optimized.out
cmp deck_leaf_audit.out deck_leaf_audit_optimized.out
python3 -B symmetry_edit.py > symmetry_edit.out
python3 -B -O symmetry_edit.py > symmetry_edit_optimized.out
cmp symmetry_edit.out symmetry_edit_optimized.out
g++ -O3 -std=c++17 -Wall -Wextra -Werror zero_escape.cpp -o zero_escape
./zero_escape 6 quotients6.txt > zero_escape6.out
./zero_escape 7 quotients7.txt > zero_escape7.out
g++ -O0 -std=c++17 -Wall -Wextra -Werror zero_escape.cpp -o zero_escape_O0
./zero_escape_O0 6 quotients6.txt > zero_escape6_O0.out
cmp zero_escape6.out zero_escape6_O0.out
./zero_escape 7 maxquotient7.txt > zero_escape7_selected.out
./zero_escape_O0 7 maxquotient7.txt > zero_escape7_selected_O0.out
cmp zero_escape7_selected.out zero_escape7_selected_O0.out
g++ -O3 -std=c++17 -Wall -Wextra -Werror direct_child_audit.cpp -o direct_child_audit
./direct_child_audit 7 85298 > direct_child_audit7.out
python3 -B zero_closure_controls.py > zero_closure_controls.out
python3 -B -O zero_closure_controls.py > zero_closure_controls_optimized.out
cmp zero_closure_controls.out zero_closure_controls_optimized.out
python3 -B card_rigidity_audit.py > card_rigidity_audit.out
python3 -B -O card_rigidity_audit.py > card_rigidity_audit_optimized.out
cmp card_rigidity_audit.out card_rigidity_audit_optimized.out
```

The k=7 full run used roughly335MB peak array storage and about16 minutes in this environment. Storage for the entire research packet remains below2MB. Runtime is not a theorem or a platform guarantee.

`SHA256SUMS` binds the delivered source, notes, and completed outputs; compiled binaries and SHA256SUMS itself are excluded.
