# Independent audit: EDGE_ORDER packet

Verdict: PASS for the stated exact construction, arbitrary nonempty private-module formulas, and limiting counterexample to a positive universal forward-edge sector fraction. No eventual odd-spectrum occupancy, tail threshold, or uniform gluing inequality is proved.

## Frozen scope

Author confirmed the following exact immutable inputs via collaboration before completion:
- EDGE_ORDER.md: 29db7536acfa23c4bcd0648bf425c54052ccb4644d0908b8b97fbe4666ae7979
- verify_edge_order.py: 658d6c62acfd143609a5a6f28f5a1b846349fc1b30b28572ff29f156f2c5b7cc
- EDGE_NORMAL.json and EDGE_OPTIMIZED.json: cb59e77684866b92fc5b7f4973f6f76c7d5c11a52362f394fe0888e8d1fba990

Other scratch explorations are explicitly excluded. Author files were not modified.

## Proof audit

1. Every unordered vertex pair receives exactly one orientation. In a path having y before x, x cannot be first and y cannot be last. Their sole external in/out neighbors force the literal segment y,b,x. Outside that segment, each of P and Q is split into its ordered before/after paths; edges P -> Q make these independent, giving c_P c_Q.
2. If x is not first in a forward path, b immediately precedes it. Then y must be last because its only outneighbor has already occurred. No Q vertex can occur before b: it cannot reach b directly or return to P. Thus the prefix is a P-path, and the remaining P-path precedes one spanning Q-path after x, giving c_P h_Q.
3. Dually, if y is not last, x is first, and the count is h_P c_Q. These first two cases cannot intersect. If x is first and y last, the middle path is P,b,Q with one Hamilton path in each module, giving h_P h_Q. Nonempty modules are needed for this case and are explicitly assumed.
4. Addition gives the stated product (c_P+h_P)(c_Q+h_Q). Each factor is indeed the Hamilton-path count of the directed triangle blowup with that module and two singleton vertices: the two singleton vertices either occur consecutively between the split module paths, or the module lies between them, yielding c+h. No additive closure or occupancy consequence is established by this identity.
5. For transitive modules each induced subset has exactly one Hamilton path, so h_P=h_Q=1 and c_P=2^p,c_Q=2^q. The cycle condition forces y,b,x in cyclic order; all remaining vertices follow their unique transitive order. There is exactly one directed Hamilton cycle modulo rotation and the tournament is strong because it has a spanning directed cycle.
6. At p=q=m, f/H=(2^(m+1)+1)/(2^m+1)^2 tends to zero, while r=2^(2m)=(f-1)^2/4. This is an all-order consequence of the proved formulas, not an extrapolation from finite data. It refutes a positive universal lower fraction for paths placing x before y when x -> y, including within the C1 subclass. It says nothing about the distinct event that the edge x -> y itself appears consecutively in a path.

No mathematical correction was required. The proof candidate's status line remains unchanged pending root's publication decision.

## Reproducibility

The author's verifier was rerun normally and with python3 -O; the resulting files AUTHOR_RERUN.json and AUTHOR_RERUN_O.json match the respective frozen receipts byte for byte. Those checks cover 16 transitive pairs p,q=1..4 and all 121 pairs of labeled modules of orders 1..3.

The independent check.py uses separately encoded edge predicates and literal permutation enumeration, not the author's builder or subset DP. RESULT.json records its exact inputs and results: 121 arbitrary module pairs of orders 1..3, plus 9 transitive module pairs p,q=1..3 with unique-cycle checks. These are finite-exact implementation controls only; the preceding decomposition proves the all-order statements.

Reproduce from the workspace root:
python3 session202-threshold-audit/edge-independent/check.py > session202-threshold-audit/edge-independent/RESULT.json
