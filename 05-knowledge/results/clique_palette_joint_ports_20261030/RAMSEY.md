# Sharp growth and an exact capacity recurrence

**Status: PROVED within the specified substitution class, using P2's exact polytope certificate; table FINITE-EXACT.**

Let \(\mathcal C\) be the class generated from the empty tournament and a singleton by arbitrary nonuniform substitutions into a transitive pair J2, a directed triangle C3, and Paley7 P. Empty ports are allowed, so the class is hereditary. This is a declared class, not the class of all tournaments and not the earlier T5/T7/T9 carousel class.

## A polytope criterion for an exponent

For a finite support family \(\mathcal S\) containing the singletons, define
\[
A(a)=\max_{S\in\mathcal S}\sum_{i\in S}a_i,\qquad a_i\ge0.
\]
Suppose its indicator polytope is downclosed and has a presentation
\[
Q=\{x\ge0:b_j\cdot x\le1\text{ for all }j\},\qquad b_j\ge0.
\]
Then the anti-blocking polytope
\[
R=\{a\ge0:A(a)\le1\}
\]
is the coordinatewise downward closure of \(\operatorname{conv}(\{0\}\cup\{b_j\})\).

Here is a direct proof of the reverse inclusion. If a nonnegative point a is outside that downward-closed convex hull D, separate it from D. The separating functional can be made nonnegative: negative coordinates never help the support function of a downward-closed subset of the nonnegative orthant, while replacing them by zero can only increase its value on a. For this functional x≥0 let \(t=\max_j b_j\cdot x\). Boundedness of Q gives t>0 for nonzero x. Now x/t is in Q, whereas \(a\cdot(x/t)>1\), contradicting a∈R. The forward inclusion follows directly from the inequalities.

**R1 (template criterion).** If γ>1 and \(\sum_i b_{j,i}^{\gamma}\le1\) for every j, then
\[
\boxed{\sum_i a_i^\gamma\le A(a)^\gamma\qquad(a_i\ge0).}
\]
Normalize A(a)=1. The increasing convex function \(\sum a_i^\gamma\) attains its maximum on D at one of the displayed generators, where it is at most 1. The zero case is immediate. This gives a reusable finite-template test, not merely a balanced-power calculation.

## Exact Paley exponent, including unbalanced substitutions

Take \(\gamma=\log_3 7\). By P2, the positive normalized facet normals for P are
\[
e_i,\qquad \tfrac12\mathbf1_C\ (C\text{ cyclic}),\qquad \tfrac13\mathbf1_7.
\]
Their γ-power sums are respectively \(1,\,3/2^\gamma,\,1\). The middle quantity is strictly below 1: \(7^3>3^5\) gives γ>5/3, and \(2^5>3^3\) gives \(2^{5/3}>3\).

For J2, the normalized support polytope normals are e_i. For C3 they are e_i and \(\tfrac12\mathbf1_3\). The same inequality holds, strictly at a C3 half-triple.

**R2.** Every nonempty \(T\in\mathcal C\) satisfies
\[
\boxed{|T|\le\alpha(T)^{\log_3 7},\qquad
\alpha(T)\ge |T|^{\log_7 3}.}
\]
Induct on a substitution tree. For child transitive sizes a_i and orders n_i,
\[
|T|=\sum_i n_i\le\sum_i a_i^\gamma
 \le \left(\max_{S\text{ transitive in the template}}\sum_{i\in S}a_i\right)^\gamma
 =\alpha(T)^\gamma.
\]
The leaf and empty cases are immediate. This handles every unequal choice of children.

Strict convexity in R1 makes its nonzero equality cases precise. For P, equality requires either a single supported coordinate or all seven a_i equal and positive. For J2 or C3 only a single supported coordinate can give equality. A strict coordinate decrease or a convex mixture of distinct generators makes the power-sum inequality strict. Induction therefore shows:

**R3 (equality classification).** After deleting empty children and suppressing single-child nodes, equality in R2 occurs exactly for the singleton and the iterated uniform powers \(P^{[r]}\). They have order \(7^r\) and transitive size \(3^r\). All equality children at a P node have the same transitive size, hence the same depth by induction.

This is a sharp family-specific Erdős–Hajnal-type exponent. It is not a proof of the general Erdős–Hajnal conjecture or an unrestricted Ramsey improvement.

## All-budget recurrence and certified finite values

Let
\[
F(k)=\max\{|T|:T\in\mathcal C,\ \alpha(T)\le k\},\quad F(0)=0,\ F(1)=1.
\]
For k≥2,
\[
F(k)=\max\left(
F(k-1),\
\max_{\substack{Q\in\{J2,C3,P\},\,a_i\in\{0,\ldots,k-1\}\\
|\{i:a_i>0\}|\ge2,\ \max_{S\text{ transitive in }Q}\sum_{i\in S}a_i\le k}}
\sum_i F(a_i)\right).                                      \tag{R4}
\]

Necessity follows by reducing away single-child nodes and applying P1. Because every pair of quotient vertices is transitive, any nontrivial node has every positive child budget strictly smaller than k; the recurrence is well founded. Conversely, substitute an attaining child for each F(a_i). Its actual alpha may be smaller than a_i, which only improves feasibility. This proves both attainment and the recurrence by induction.

The exact outputs, indexed by k=0,…,12, are
\[
0,1,3,7,9,13,21,25,33,49,51,55,63.
\]
For every r≥0, R2 and the balanced construction give the closed subsequence
\[
\boxed{F(3^r)=7^r.}
\]

The producer uses vertex-transitivity of the legal-support systems to put a maximum budget at port 0. For J2 this is a symmetry of the support system, not an automorphism of the directed pair; it is legitimate because R4 depends only on that system. A separate checker enumerates all labelled budget vectors without this symmetry reduction through k=6, agrees on every value, and validates every stored construction through k=12. No numerical logarithm is used in the exponent proof.

The previous carousel capacity \(A006046(k)\) remains correct for its own palette. Paley7 already beats it at k=3, where the earlier capacity is 5 and this one is 7. That boundary motivates the present class; it does not invalidate the earlier restricted theorem. The full mixed class containing both Paley7 and all earlier carousel templates needs additional weighted template tests.

**Next decisive test.** Determine the anti-blocking power-sum maxima for those additional templates with exact facets. If each is at most one at γ=log_3 7, R2 extends unchanged. If one exceeds one, its maximizing weights specify the first nonuniform candidate to test, retaining actual child budgets.
