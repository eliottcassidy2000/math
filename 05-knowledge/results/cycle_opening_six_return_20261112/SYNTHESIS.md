# Full-count returns: the repaired step, the repeatable families, and the missing induction

**PROVED IN NOTES / FINITE-EXACT. Chan–Pak remains OPEN.** The results below were internally audited by independently working lanes. They are not a literature-priority claim or an external referee certification. The current task concerns Hamiltonian **paths**; H is their full count.

The principal advance is an actual all-size source return that restores defect one after a target-changing excursion. It keeps the same physical root and incoming port, and its proof controls all Hamiltonian paths. Two other constructions give indefinite neutral returns inside explicitly classified families. What remains missing is an unlimited source-return grammar at the prescribed affine targets, sufficient to cover the odd spectrum outside the known forbidden values 7 and 21.

## 1. The five-core source now survives one exact target-changing return

Let k≥1. There are strong tournaments on k+7 vertices with the same full count

\[
H_k=45\cdot2^k+7
\]

and the same physical indegree-one root, connected by a literal two-arc switch changing its defect from \(6\cdot2^k+1\) to 1. The parent is the substitution

\[
C[TT_2,TT_k,TT_2,1,1],\quad C=(22,4,8,3,14),\quad H(C)=9.
\]

Thus the five-core connection is an explicit quotient map with homogeneous modules, not an association of equal cardinalities. The switch is an involution at every k.

From its defect-one side, two standard births give the exact target

\[
3H_k+2=135\cdot2^k+23,
\]

but defect \(24\cdot2^k+5\). A specified replacement on the same k+9 vertices restores defect one at **that very target**. It retains the root, its unique incoming neighbor, every internal edge of TT_k, and all but six fixed arcs and three module boundaries. Its complete changed-arc list and retained-name map are in [MODULE_RETURN](compiler/MODULE_RETURN.md). An alternative [module-migration recipe](returns/MODULE_MIGRATION.md) uses 3k+3 changed arcs while retaining the root and changing its incoming neighbor. The two target families are explicitly isomorphic, with the substituted module mapped correctly.

The general theorem is stronger than transitive inflation. For any nonempty tournament module M, write u=H(M), v=P_M(2). The two stages have

| Stage | Full H | Root defect |
|---|---:|---:|
| High-defect five-core source | 45v+7u | 6v+u |
| Two-arc neutral repair | 45v+7u | u |
| Two standard births | 135v+23u | 24v+5u |
| Parent replacement | 135v+23u | u |

These are exact ordered-run identities; no uncounted sector remains. For TT_k, u=1 and v=2^k. For a nontransitive M, the offset is 2u rather than 2, which is the exact realization guard.

The next target is \(405\cdot2^k+71\). We do **not** have a uniform ready replacement there. The first two templates cannot simply be reused: their constant residues exclude numerical reentry modulo45 and135. At the first successor135·2^k+23, the full fixed-parent port slice is already ruled out at every k; the proved replacement escapes by changing that parent. Larger parent changes or a different interface are necessary. The [ancestry boundary](compiler/ANCESTRY_BOUNDARY.md) makes this structural: the raw excursion has two indegree-one vertices, while each repaired source has only one. The replacement leaves the strict indegree-one-port birth language, although it remains an ordinary adaptive birth from a parent with a higher-indegree port.

## 2. What is completely classified

Each row is a complete classification in its stated universe; the union is not a classification of the whole tournament metagraph.

| Universe | Complete result | Composition and boundary |
|---|---|---|
| Indegree-one extensions of one fixed named parent | Level sets of \(Q_a=E_a(2)-s_a\) are exactly the full-H fibres; each is a clique under two-arc predecessor switches. | Defect labels are \(s_a\). Old vertices/arcs and the physical newborn are retained. Parent rewrites are outside this slice. |
| One non-port substituted module in the complete ordinary parent bank through order8 | 225 coefficientwise neutral charts, 78 with strong parents. | Each certificate works for every nonempty substituted tournament. This counts charts, not nonisomorphic infinite families. |
| Prepared L/T constructor words | Exact three-coordinate state; uniform neutral A/B macros; complete-state diamonds ABA=BAA and ABB=BAB. | All depths are known. Counts have a permanent factor51, so this family cannot supply the missing odd-spectrum coverage. |
| Fixed-target macro addresses with prescribed prefix multiplicities a,b | Johnson graph \(J(a+b,a)\square K_2\) for the specified neutral bit switches; adjacent swaps give a distributive-lattice refinement. | For a fixed strong source, canonical birth peeling proves distinct addresses are nonisomorphic. Additional neutral edges outside these bit switches are not classified. |
| Unique-start funnel observer | Six coordinates are necessary and sufficient even for affine rational state models. | Every periodic empty/one-tail program misses unlimited \(R\mapsto3R+2\), even at period ends. This does not forbid larger or nonperiodic grammars. |
| Positive integral scaled monomial six-cycles | Prime-wise valuation hexagons; \(\prod_{p^c\parallel C}(1+3c(c+1))\) ordered states at reservoir C. | One fixed state; every other state has exact period six. This is an arithmetic classification, not a classification of all tournaments with those counts. |
| Marked block-resource tournaments | Autonomous hexagon map, exact six-return, fixed full H, and a strong wrapper preserving the same root/port. | Its whole uniform-module response is fixed; the wrapper's defect is zero, so it does not perform the positive-defect repair above. |

The incoming [Session293 conveyor](../session293_cycle_returns_20261113/opening/PROOF.md) supplies another complete all-order slice: neutral directed cycle-opening births with a moving entrance port. Its feedback obstruction and our parent-changing repair are integrated in [INCOMING](INCOMING.md).

## 3. A repeatable neutral family with a known metagraph

For an indegree-one root use \((x,b,d)=(H-H(T-r),\delta,H(T-r))\). The actual one-vertex constructor L and two-vertex constructor T have

\[
L(x,b,d)=(x+2b,d-b,x+d),
\quad T(x,b,d)=(2x+2b,d,2x+2b+d).
\]

Set \(q=2x+2b+d\), \(A=TTTLL\), \(B=TTLTL\). Both add eight vertices and yield full H=51q, with states

\[
A:q(24,9,27),\qquad B:q(22,12,29).
\]

A simultaneous four-arc switch exchanges these outputs, retaining their final root and its predecessor. Future macros multiply q by93 or97 according to the current type. At macro depth m the complete state list is

\[
q_0\,93^a97^b v_E,\qquad a+b=m-1,\quad E\in\{A,B\}.
\]

There are exactly2m continuation states, m distinct H-values, and \(\binom{m-1}{a}\) formal words per state. The full-H distribution assigns weight \(2\binom{m-1}{a}\) to \(51q_0 93^a97^{m-1-a}\). This compresses an exponential address enumeration into a linear-size exact state census.

The general unique-ancestry and automorphism mechanism is inherited from [ANCESTRY A1–A3](../adaptive_metagraph_reformulation_20261103/ANCESTRY.md); the present application classifies the selected macro family by full H and neutral rewrites. For a strong source, the newest vertex is canonically recognizable: it is the unique indegree-one vertex whose sole predecessor also has indegree one. Peeling births recovers every new vertex position, and the private four-arc gadgets then recover every address even from its unlabelled isomorphism class. Consequently the 2^m addresses are exactly 2^m nonisomorphic strong tournaments in this family, with the stated binomial full-H multiplicities. With a,b fixed, neutral exchanges of two prefix choices and the final type generate the Johnson-product graph above. It has \(2\binom{a+b}{a}\) vertices, degree ab+1, and diameter min(a,b)+1. Restricting prefix exchanges to adjacent positions gives the rectangle-ideal graph times K2, of diameter ab+1. See [FIBRE_GRAPH](returns/FIBRE_GRAPH.md) for proof and exact distance enumerators.

The broader [binary source language](returns/BINARY_LANGUAGE.md) also closes exactly: the same birth can mark its newborn (L) or retain its old root (K), with K(x,b,d)=(2x+2b,0,d). At c added vertices it gives2^(c−1) ordinary graphs and2^c rooted graphs for c≥1; the original L/T restriction has Fibonacci count F_(c+1), since T=KL costs two vertices. These selected counts are smaller than the already known unrestricted adaptive-descendant census. Their useful refinement is the exact three-coordinate full-H observer and the specified neutral-fibre laws, not a new general ancestry theorem.

At six macro stages, a=2,b=3 gives J(5,2), the line graph of K5, times K2. The Petersen graph is the complement of the J(5,2) factor. This is a precise five/six connection through addresses and exchanges. It is not an identification of these address states with the five tournament cores. The older square-sum completion Johnson graphs share the same fixed-inventory exchange mechanism, with different objects and different preserved predicates.

## 4. Six as a period, a lift, and an information requirement

The recovered quadratic orbit is

\[
-7/4\longmapsto5/4\longmapsto-1/4\longmapsto-7/4,
\quad f(x)=x^2-29/16.
\]

The repeated triple in the prompt omits5/4. This particular quadratic map has no rational six-cycle. Its interpolating matrix has cube−I and sixth power I, so the lifted vectors have order six while the projective points have order three.

A different rational map on pairs,

\[
L(a,b)=(b,b/a),
\]

has exact period six at every positive rational pair other than (1,1). Its full positive integral scaled recurrence is classified by finite valuation hexagons. The autonomous resource controller

\[
F_m(r,s)=(2m-s,r+s-m)
\]

on \(0\le r,s\le2m,\ m\le r+s\le3m\) turns that action into actual named tournament returns. Retaining complementary blocks keeps the full count constant. A fixed cyclic wrapper makes the tournament strong and retains the same indegree-one root and incoming port. Exact period six belongs to the named resource-and-cut state: forgetting the cut can identify phases, even when the labeled unmarked graph is unchanged. This is genuine iteration, but it preserves the entire external response and hence cannot reset a missing defect.

By contrast, six in [the counting compiler](compiler/SIX_STATE_MINIMALITY.md) is a minimal information dimension. Exact observation determinants552960 and811008 prove that neither linear nor affine rational models with fewer coordinates can reproduce all legal continuations. That theorem is not an order-six dynamical statement.

Finally \(K(a,b)=(b,(1+b)/a)\) has period five on positive rationals and a golden-ratio fixed point over the reals. Removing +1 changes its period to six. As an operation on odd full tournament counts, K fails immediately: an integral next value would be even. The comparison is an exact arithmetic bridge with an explicit realization obstruction, rather than a claimed cluster-algebra tournament construction.

The raw, Terras, and odd-only 3x−1 clocks give different periods for the familiar cycles at minima1,5,17: respectively (2,5,18), (1,3,11), and (1,2,7). The [dynamics proof and recovery ledger](dynamics/PROOF.md) retain these clocks, the corrected circulant dictionary, and the primary literature's quadratic-period hypotheses.

## 5. Near-squares and enumeration: the actual consequence

For block types A,B with h-values a,b and two-bin values alpha,beta, the strong six-return carrier has full count

\[
N_m=(ab)^{2m}+(\alpha\beta)^{2m}.
\]

This is an odd square plus an even square. Primality requires m to be a power of two and the two bases to be coprime; these conditions are not sufficient. For the H5/H17 blocks, m=1 gives51,387,449=29·1,771,981. The structural unit-square branch forces a=b=1 and hence transitive blocks, reducing to the Fermat-type form \(1+2^{2m(|A|+|B|)}\). It gives actual tournament realizations, not new primes or a proof of Landau's fourth problem.

The strict ancestry slice also has an exact size recurrence. If u_n counts ordinary strong tournaments with exactly one indegree-one vertex and t_n those with exactly two, then t_4=1 and t_n=u_(n−1)+2t_(n−1) for n≥5. The complete inherited bank through order8 gives u_8=2209 and t_8=236, hence **t_9=2681** without a new order-nine enumeration. Its generating-function decomposition retains the unknown one-root core inventory U(x) explicitly: S(x)=(1−x)/(1−2x)·(x^3+U(x)). This is a precise structural refinement of tournament enumeration, not a solution of that core inventory.

The exact new enumeration now includes an ordinary strong-tournament subfamily: from a fixed strong rooted source of order n, there are exactly 2^m nonisomorphic outputs at order n+8m. For that fixed source the size series is x^n/(1−2x^8). Automorphisms of a positive-depth output are exactly the rooted automorphisms of the source. For a rooted triangle all outputs are asymmetric, giving2^m(3+8m)! labeled members at positive depth. These are proved structured contributions to A000568 and its labeled family refinement, not new global coefficients. The continuation-state compression and Johnson fibres retain exact count multiplicities at every depth.

A135244 is an arithmetic sequence: the largest integer whose proper-divisor sum is the given target, or zero when no such integer exists. It is not a tournament enumeration. The inherited [extremal aliquot theorem](../uniform_return_metagraph_20261111/arithmetic/PROOF.md) gives A135244(N+1)=N^2 exactly when N is prime. Applied to N=N_m, this supplies a precise conditional bridge from the strong return carrier to the aliquot maximum. It does not decide primality or improve the certified arithmetic range. No larger global A135244 range or new A002496 terms are claimed here.

## 6. The remaining proof obligation

The most concrete next target is a uniform source at405·2^k+71 with compatible retained names and a law that returns to a reusable family. The latest bounded sparse-parent search has no defect-one or defect-two repair at the smallest target881 within its declared radius-three universe. A separate [all-defect radius-two census](returns/GROWING_SUCCESSOR.md) finds only defects0 and213, with213 already the raw excursion defect; its seven positive single-module lifts retain the same108·2^k−3 defect law. These are local obstructions, not global nonexistence theorems. A uniform repair may need more parent arcs, more than one changing module, larger defect, or a different constructor.

A proposed successful induction must preserve the actual target and expose a usable next interface, including the coordinates observed by its next operation. The five-core replacement proves that this can be done beyond the fixed-parent obstruction. The neutral diamond language proves that complete-state reentry can compose indefinitely. Combining those two strengths without the factor51 restriction remains open.
