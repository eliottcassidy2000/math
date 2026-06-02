        # Message: oracle-2026-06-02-S559o: the A+B=C summand graph = the LRC pinch denominator; addition-shadow=pinch, multiplication-shadow=sieve, equal on the rationals (HYP-2066)

        **From:** oracle-2026-06-02-S?
        **To:** all
        **Sent:** 2026-06-02 18:56

        ---

        User asked me to search the repo's prior work on even/odd and on multiplication/addition, recalling an additive tournament construction: nodes A,B,C with arcs (A,C),(B,C) iff A+B=C. I found the object and recorded a bridge to the live LRC pinch frontier that nobody had written down.

THE CONSTRUCTION = THE SUMMAND GRAPH (oracle-2026-05-03, 07-reflections/summand-graph-fermat-zeckendorf.md): the digraph on N with a->n and b->n iff a+b=n (a!=b). Known: node n's incoming-pair-count is floor((n-1)/2), PARITY-SPLIT (odd n -> (n-1)/2; even n -> n/2-1, the midpoint pair n/2+n/2 excluded). Its Section 4 already proves 'the summand graph IS the tournament staircase' (node n <-> delta_{n-2}; depth-k multiplicity 2^{k-2} = the tournament binary tree; H-values on Zeckendorf/Fibonacci sums; the forbidden {7,21} recur as nodes). Sequel oracle-2026-05-05 adds Fibonacci/Lucas summand chains.

THE NEW BRIDGE (recorded, HYP-2066). opus-S557's pinch lemma: the loneliness radius is M=r/s with s=(v_a+v_b)/gcd = the REDUCED PAIR-SUM, attained at the pinch time t*=m/(v_a+v_b). oracle-S555o's pinch pigeonhole runs over 'the floor(n/2) pairs summing to n' -- which are EXACTLY the incoming pairs of summand-graph node n. So:
   summand-graph node C  =  the LRC pinch denominator s = v_a+v_b;
   its incoming pairs     =  the pinch-pairs;
   its pair-count floor((C-1)/2) = the number of pinch times at denominator C.
The user's additive tournament on N literally INDEXES the LRC pinch times. The parity split is the parity of the pinch: an EVEN denominator C=a+b loses the midpoint pinch (the apex pair C/2+C/2) -- the very apex / Z_14 zero-divisor that obstructs n=14=2*7 (opus-S559, HYP-2063).

ADDITION vs MULTIPLICATION = THE TWO WITNESS SHADOWS (codex-S365 + THM-369). codex-2026-05-31 (natural-operation-digraphs-and-product-sum-s365.md) collapses the two-input operation systems to simple digraphs: the ADDITIVE shadow (a->z,b->z iff a+b=z) = the TRANSITIVE tournament (x->z iff x<z); the MULTIPLICATIVE shadow (xy=z) = the DIVISIBILITY DAG (x->z iff x|z, sparse). In LRC: the additive shadow is the PINCH (t=m/(a+b)); the multiplicative shadow is the SIEVE (THM-369: t=1/q lonely iff no speed divisible by q). And oracle-S555o PROVED the rational pinch IS the sieve -- so the two shadows COINCIDE on the rationals, and the FINE pinch (q>n, S556o) is exactly where the additive (summand) structure has witnesses the multiplicative (divisibility) structure cannot see = the open core. Same +/x spine as the S548 hyperoperation tower (succ->+->x, prod c_i = exp(-sum log c_i); log welds them; rationals = where the weld is exact).

PARITY thread surfaced: pair-count parity = pinch-count parity (even C drops the apex pair); Burnside even cycle types fix 0 tournaments, all-odd fix 2^pairs (A000568 is an all-odd sum, x2 is the killer); doubled prime 2q = 'parity-completion of the prime additive basis' (S546o); 14=2*7 is our hinge; the even-fold measure reduction (HYP-2065) is the x2 acting on the summand structure.

This is a search/synthesis session (no new computation); the contribution is surfacing the prior work the user remembered and recording the unrecorded link summand-graph <-> pinch/sieve. New HYP-2066; reflection 07-reflections/the-summand-graph-node-is-the-lrc-pinch-denominator-addition-and-multiplication-shadows-s559o.md.

HANDOFF: (1) use the summand-graph incoming-pair lattice of FINE denominators q in (n,Cn] to organize the fine-regime pinch pigeonhole (the open S555o handoff -- the additive shadow is the right index set); (2) restate HYP-2063's apex zero-divisor purely as 'the even node's missing midpoint pair'; (3) Zeckendorf/Lucas summand chains (oracle-2026-05-05) as a SPARSE pinch family.

        ---

        *Reply by writing to `agents/oracle/inbox/` or run `python3 agents/processor.py --send --to oracle`*
