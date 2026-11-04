# Independent audit: evolving ports and persistent clone recurrences

2026-11-04. **PASS**, with no mathematical correction required.

## All-order proof inspection

The full named path-matrix birth formula is correct. Splitting a path at the newborn proves the scalar port updates; multiplying in the squarefree ring removes all quadratic newborn terms. The four endpoint locations in the individual joint update are correct, and the inner-defect formula uses the actual current in-neighbor set. This is a full-state transport, not a recurrence on only three scalar responses.

The clone formula is also correct. Intervening old paths must start in the old out-neighbor set and end in the disjoint old in-neighbor set, so each consumes at least two vertices. Permuting the disjoint middle paths gives a free action and proves (r−1)! divides A_r. The ordered nonempty clone-run assignment contributes exactly r! S(t,r), with transitive order forced inside each run.

The eventual-recurrence argument is valid on a tail. Distinct positive exponential bases remain independent after any finite threshold, and the largest base has coefficient exactly A_R, which cannot be cancelled by lower indices. Therefore pure ternary recurrence forces A_r=0 above3, A_2=3A_3, and A_1=3A_3−d/2. The divisibility makes A_3 even; Rédei parity makes A_1 odd. Hence any integer persistent digit must satisfy d≡2 modulo4, excluding +4. This is a necessary condition only, and does not establish +2 realizability.

The proof correctly allows finite transients and correctly excludes only lengthening one fixed transitive module in a fixed external context. No exclusion of changing contexts, modules, or nontransitive interiors is justified or claimed.

Minor clarity request sent to the author: nonempty W is sufficient for a strong full rooted attachment when the suffix is strong, not by itself. The author is clarifying that hypothesis. The independent induction packet gives the sharper exact terminal-SCC criterion.

## Independent numerical checks

`independent.py` imports none of the author's code. It enumerates permutations of the outside vertices and every admissible segment-length tuple, directly counting old path covers before, between, and after clone runs. It checks 285 rooted tournaments through order4 and 1,425 clone counts (t=1,...,5) by a separately written subset Hamilton-path DP, also checking parity and the factorial divisibility.

It independently verifies the mixed example coefficients [9,10,2] and its formula 2·3^t+4·2^t−5 through t=7.

It additionally verifies the author's new pure ternary control on rows [26,24,67,116,36,71,19], cloning physical6: coefficients [145,150,50,0] and counts 145,445,1345,4045,12145,36445 for t=1,...,6. Thus H_t=50·3^t−5 and H_(t+1)=3H_t+10. This positive example shows the spectral obstruction is selective, not absence of ternary modes.

Normal and optimized independent outputs byte-match. The initial author verifier replay also byte-matched normally and under -O: 33,864 joint-port birth cases, 1,098 full-matrix birth cases, 402 strong birth guards, and 1,140 clone response checks. Subsequent author additions of the pure ternary control are separately covered by the independent script and may update the author receipt.

No installations or forbidden-file reads were used.
