# Verification receipt

Run:

    python3 -B module_contexts.py > module_contexts.normal.json
    python3 -B -O module_contexts.py > module_contexts.optimized.json
    cmp module_contexts.normal.json module_contexts.optimized.json

Both executions exit0 and the outputs are byte-identical. `require` raises exceptions and remains active under `-O`.

- All1099 labelled tournaments of orders1..5: canonical-subset path-cover DP equals literal permutations/cuts.
- All1099 labelled contexts of quotient orders1..5, hole fixed at label0: first coefficient equals H(Q), two-visit positivity agrees with independent reachability, inactive vertices have no repeat coefficients.
- Seven distinct factor profiles through order4 times1099 contexts:7693 direct substituted-graph checks.
- Threshold signatures at B=1,3,7,15,31,127:23079 applicable large-input checks.
- Triangle decoder shells s=1..5; padding hostile orders1..8; inherited same-H masks40/76.
- All-depth closures at B127: K1/J/C3 has58 states and31 values; all seeds through order5/J/C3 has65 states and34 values.
- All31+34 positive expression certificates reconstructed and independently counted at vertex level; maximum witness order12.

The bounded contexts in the finite audit do not prove the all-order theorems; PROOF.md supplies the bijections, bounds, and finite-algebra completeness argument. The grammar closure itself is exhaustive over all expression depths because its sufficient state space is finite and the algorithm stops only at a fixed point.

`probe.py` and `inherited_THM3121.md` are research scratch, not required for reproduction. The standalone main script imports only the Python standard library.

Supplementary complete width-five grammar:

    python3 -B width5_probe.py > width5.normal.json
    python3 -B -O width5_probe.py > width5.optimized.json
    cmp width5.normal.json width5.optimized.json

- All1024 labelled order-five quotients tested for proper modules and literal isomorphism; prime-bank masks10,41,76, multiplicities120,120,24.
- Inherited Moon–Busch f(11)=225>127 justifies total-order<=10 pruning only for strong quotient outputs.
- Three prime operations evaluated2929 times each over the completed saturation rounds;75 direct small-graph profile controls.
- Closure86 states,49 positive values, with49 directly checked adjacency certificates (maximum order11).
- Outside width-five H31 witness at order6 checked by DP and literal permutations.
- Hash-stable tie breaking was made explicit when selecting positive witnesses; normal/-O output now byte-matches.

Tiny width controls in the main script check all64 order-four quotients against four J/C3 recipes and count the order-six H19 seed/operator hostile by literal permutations.

Independent audit: a separate agent imported no producer code for its core tests. It rederived literal small profiles, used a separate forward-length quotient-word engine, and independently reproduced the58/65-state closures and the86-state complete width-five closure, including all attained/missing lists, all three prime-five classes and masses, and literal H31. It approved the all-order context/algebra proofs and the strong-order pruning. Separate audit artifacts are held in the session's independent-audit lane rather than this payload directory.
