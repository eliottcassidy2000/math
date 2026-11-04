# Independent review of the transitive ordered-C3 interface

**PASS.** The frozen proof is `cde272a68fb8c752080871b34c8e86529ab5e7eb6a82e9bfe2a07a7340a3c61c`; verifier `ee696dfe4d46eb1bd3c497ddcc3a08af41fbb5232dc9a220759cf5ede5cd2b20`. No mathematical correction was needed.

## Exact interface and polynomial

The physical quotient rows `(4,13,8,17,7)` give `X->A->Y`, `X->Y`, both X and Y above F and below L, `F->L`, and both F and L above A. This verifies the transitive attachment, rather than the cyclic attachment of the positive theorem. The eight interface quantities in the proof follow by selecting nonempty subsets of the two transitive chains. In particular `sigma=PR`, `mu=tau=0`, and `kappa=P+PR`; omitting a constant empty segment would change the answer.

The displayed rational expression follows from the complete four-term local triangle kernel. The producer checks it after putting over the common denominator, without numeric fitting. Its truncation in the marked X/Y run variables is safe because all later factors have nonnegative exponents. Only powers through four of the return term can contribute after the leading `xy` to marked degrees at most three.

Independently, `audit.py` does not import the producer algebra. It constructs every bounded run word directly in the exterior quotient X,Y,F,L, with exact start/end colours; builds its endpoint, return and path-matrix series; and substitutes them into the original four-term triangle kernel. It then applies the actual marked triangle sequence `(0,1,2,2)`. The complete polynomial agrees at every exponent in the larger box `(3,3,8,8)`, and has exactly the displayed 23 terms, with no surviving higher P/R exponents. The rational support argument supplies the all-order bound beyond this box.

Since `[F]P^s=s!S(a,s)` and `[L]R^t=t!S(b,t)`, the extraction has the ordered-cover normalization claimed, with no extra s or t factor for a prescribed marked slot. The complete physical local-triangle and quotient-converse finite differences agree, including their signs and physical labels.

## All-order sign argument

Every positive coefficient has both exponents at most four. The only negative coefficient outside the upper exponent quadrant is `-96 P^5 R^3`; removing it can only increase the value for actual positive transitive sizes. The previously proved strict Stirling minors imply the required monotonicity of `p_n(s)/p_n(4)` as n increases: it decreases for s<4 and increases for s>4. Their products are nonnegative, so termwise comparison in the stated directions is valid.

At the base size six, the positive sum is `2080753736`, and the retained negative sum is `-2560896000`. Their sum `-480142264` is strictly negative. Since `p_a(4)p_b(4)>0` for a,b>=6, the quantitative bound proves strict negativity throughout the unbounded stated range, not just for the finite tested sizes. The omitted contribution at six is `-93312000`, giving the full value `-573454264`.

## Root families and scope

The physical map into W5 is correct, and the quotient is prime. Therefore every proper original module is contained in one child. Scalar equality restores the unmarked chain profiles even for crossing interval families inside a transitive child. The marked triangle children have no nontrivial proper modules, and their physical endpoint profiles are invariant under converse. Consequently selected A restores the same negative quartic; unselected A restores the whole-converse response on the exposed-anchor queries. These are the precise restoration hypotheses used, without assuming arbitrary marked-profile invariance in the chains.

The six compressed exterior blocks and marked outdegree two are both even, so no anchor parity correction is missing. The union of sixteen quartic corners and seven anchor queries has exactly 22 members, since `(a={0},c=empty)` is shared; the scalar query is already included. Finite root-family controls cover all 32 whole-child subsets only. The all-original-module assertion follows from the inherited structural and profile-restoration theorems, not from treating those 32 cases as a census of every interval family.

The result is a whole-core root-family theorem with the stated ordered-SCC padding. No arbitrary strong outer lift is claimed. The smaller positive coefficients are sign-change controls, not response collisions.

## Independent computation

Both independent modes are byte-identical. All 598 retained weighted physical counts are reproduced with a separately written predecessor recurrence, followed by inclusion-exclusion for transitive ordered surjections. All 144 smaller expanded-name counts are also reproduced with a full predecessor subset DP. All nine local/quotient quartics and 32 whole-child battery rows agree. The largest queried full tournament has 27 names, and `27!` fits in the unsigned 128-bit weighted arithmetic; smaller literal cases fit in their stated range. The standard compiler is used, with no installation or external arithmetic package.
