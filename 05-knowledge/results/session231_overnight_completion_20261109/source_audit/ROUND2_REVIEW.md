# Independent audit: two-clock rigidity and halfword coordinates

Verdict: PASS, scoped as stated. The arbitrary odd-a two-clock rigidity theorem is valid, including empty fixed words, leading zeros, vanishing quadratic coefficients and fixed clock shifts. The positive mirrored families are inherited centered fans. The new halfword coordinates and square-root candidate are a valid decoder refinement, not a new all-target construction or a new deficit obstruction.

## Independent symbolic coefficient derivation

Write z=2^s and P=value(U). Ignore constants independent of z, but retain source parity: every nonempty changing run has final parity equal to its digit, including a zero run starting from numeric zero. Appending W changes a source v to Mv+r and g to g+C v+E parity(v)+constant, with C=M-1-R. E is zero for empty W; for nonempty W it is the indicator that its first digit is zero. Since the first clock's parity is fixed, E changes only the constant.

Before the final fixed suffix, the pairs of quadratic and linear coefficients of (v,g) are:

- 11: v=(M(P+1), r-M+1); g=(0, C(P+1)).
- 00: v=(MP, r); g=(MP, r-PR).
- 10: v=(M(P+1), r-M); g=(M(P+1), -(R+1)(P+1)-(M-r)).
- 01: v=(MP, r+1); g=(0, P(1+C)).

These follow by composing v->z(v+1)-1,g->g for a one run and v->zv,g->g+(z-1)v+parity(v) for a zero run. They require s>=1. The suffix changes the effective linear v coefficient to K=c 2^|S|+a C_S while retaining coefficient a on g. Its parity-dependent term is constant because the second clock fixes parity. For empty S, K=c and the original e term is still constant. Multiplying these coefficient pairs by K and a gives exactly every entry of the producer's table.

For mixed 10 the linear coefficient is strictly negative: P+1>0, R+1>0 and M-r>0. For mixed 01 it is strictly positive because K(r+1)>0, even if P=0 and the quadratic coefficient vanishes. Thus no accidental zero-coefficient corner escapes mixed-clock exclusion. In the 00 case P=0 would make the quadratic coefficient zero and cannot match a positive dyadic ray. In 11, P+1 is always positive.

For any nonzero k-bit integer n=2^t u with u odd of length l, rev_k(n)=2^(k-t-l)rev_l(u). If odd a times rev_k(n)/n is a dyadic rational, comparison of odd parts gives a rev_l(u)=u. Both u and rev_l(u) lie in [2^(l-1),2^l), hence their ratio is strictly below 2. A positive odd integer a below 2 is 1. This argument covers fractional dyadic ratios as well as integral ones.

In 11 the nonzero quadratic coefficient K M(P+1) being a power of two forces K and P+1 to be powers of two. Unless W is all ones, q=M-1-r>0 and C=rev(q)>0; linear cancellation invokes that reversal lemma. In 00, (a+K)MP being a power of two forces P>0 and a+K,P to be powers of two. Unless W is all zeros, r>0 and the same lemma applies. Thus a>1 forces merging. Empty W satisfies either homogeneous condition and correctly merges. The endpoint coefficient b affects only the constant and cannot help cancellation.

For shifts s+h and s+k with nonnegative fixed h,k, absorb the h extra beta digits into U and the k extra gamma digits after W (or into S). For arbitrary fixed integer shifts on a sufficiently large valid tail, first reindex to make both shifts nonnegative, then absorb. Reindexing multiplies the dyadic slope by another power of four; the proof works even for a dyadic-rational slope, since a successful quadratic coefficient is a positive integer. No claim follows for different growing rates or more changing blocks.

The odd-a restriction is substantive: a=2,b=0,c=2 and the separated word 1^s 00 1^s give 2^(2s+3)-2. This independent countercontrol passes, and prevents a mistaken extension of the theorem to arbitrary a.

## Halfword formulas, decoder and physical witnesses

The fixed-middle singleton calculation gives quadratic coefficient 2^k, linear coefficient r-R and constant g(r)+R-2^k+2-first(w). Hence palindromicity is necessary and sufficient for eliminating the linear term at fixed width. Empty middle is separate and gives 4^s+1.

The producer's f concatenation, complement and mirrored-half formulas are consistent with the positional pair weights and with the recurrence for f. In particular l=1 complement has zero incident weight, and l=0 is handled by direct evaluation rather than illegally using B-2. Substitution yields exactly equations (7) and (8). Their q=0 cases include the empty word and either single central digit and give offset +1.

For q>=1, 0<=f(q)<=q-1 follows from f=q+e-g and g>=1. Therefore all three proposed deficit expressions belong to [q^2,(q+1)^2). For q=0 they equal zero. One integer square root is consequently sufficient to select the only q; exact equality still must be checked. For odd exponent, d must be odd and <=1. For even exponent, d<=1 is necessary. The exponent lower bound must also be checked. The producer spells these requirements out. “Constant-time candidate selection” means a constant number of candidates, not constant bit or arithmetic running time; evaluating f takes O(log(2+|d|)) arithmetic operations after the integer square root.

The length check is essential: positive arithmetic equality at a too-small exponent is not a source certificate. Padding a fixed halfword by leading ones corresponds to absorbing paired ones into the two outer clocks, so q together with exponent parity and the minimum-length test captures membership without pretending to identify a unique unnormalized word presentation. The separate beta choices may coincide, notably at q=0. The independent decoder preserves both rather than claiming a unique beta or unique graph.

Appending the final zero makes the scalar singleton readout equal g(value(BITS)). For 1^s w 1^s 0, the first bit is one and the final bit is zero, giving the explicit directed cycle through the ordered chain and apex. Thus the graph is strong. Its order is m+2, with no order-off-by-one error. Literal small Hamilton counts independently agree. Missing offset -9 is only a statement about this mirrored grammar.

## Independent exact controls

`round2_independent.py` evaluates Hamilton counts by the direct positional formula, not the producer's g recurrence. It tests all small fixed contexts plus random longer contexts, all four bit pairs and selected physical/nonphysical rows, including a deliberately even-a boundary control. Every claimed coefficient pair is checked against exact evaluations at s=1,2,3,5, rather than extrapolated from three fitted values.

Receipts record:

- 43,800 coefficient-row cases, including 622 pure-dyadic positive controls.
- 32,576 finite reversal-lemma controls.
- 4,599 halfword identities, including larger clocks.
- 21 literal graph checks with memoized Hamilton-path enumeration.
- 34,068 decoder queries retaining exponent parity and minimum length.

Normal and optimized runs agree. These are controls for the all-size algebraic proofs, not a numerical proof by sampling. No new interval coverage is claimed.

## Inheritance and priority boundary

I read Session191's one-apex proof and CENTERED-FANS from the repository. Session191 treats one changing homogeneous run. CENTERED-FANS already proves midpoint cancellation, exact palindromic construction, fixed-deficit rays, strict fixed-diameter energy ordering, a binary-search decoder, O(sqrt X) sparse deficits and permanent grammar-relative deficit 9. It explicitly retains diameter/schedule data and contains physical graph witnesses.

The producer correctly credits those results. The present narrow difference is the arbitrary odd-a readout obstruction for two equal-rate clocks in arbitrary fixed context, plus explicit halfword formulas yielding one square-root candidate. This audit confirms that comparison against the named predecessors; it does not establish global repository or literature priority. The suffix law and finite-state selector obstruction are inherited too.

The final next-step language should be read as a statement about escaping this fixed-context equal-rate ray ansatz, not a necessary classification of every possible proof of F2 completeness. The proof does not eliminate a=1 centered fans, different rates, three blocks, varying contexts or other constructions. None of those escapes is proved sufficient. Exact all-target F2 occupancy remains open.
