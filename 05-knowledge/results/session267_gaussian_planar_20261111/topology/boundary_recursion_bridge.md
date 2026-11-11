# Session 267: a realizable boundary grammar, rather than a planar-topology analogy

**Status:** new conditional reduction architecture; elementary deductions proved below. The geometric compression, chart-recognition, and moment-transport conjectures are OPEN. No proof of JC(2), no degree bound for all Keller maps, and no newly verified classification of infinity is claimed.

## 1. What is inherited, and what is new

The working object is a complex polynomial Keller map

    F=(P,Q): A²_C -> A²_C, det JF ∈ C*,
    d(F)=[C(x,y):C(P,Q)].

Its complex dimension is two and its real dimension is four. The boundary curves in a compactification have complex dimension one. The planar drawing of a dual resolution tree is a combinatorial representation, not a reduction of C² to a flat real surface. The topology relevant here is the topology of curve complements and their covers, plus algebraic extension over the boundary.

Closest inherited proved mechanism: THM-2633, `01-canon/theorems/THM-2633-derangement-character-obstruction-and-d4-keller-exclusion.md`, gives finite affine survivors over every target divisor and normal generation of monodromy by sheet-fixing affine inertia. This is considerably stronger than an unmarked picture of a discriminant curve.

Closest corrected near miss: THM-1365, `01-canon/theorems/THM-1365-galois-reduction-jc-fixed-point-bridge.md`, has a scope pin: generic deck transformations do not automatically preserve the affine source chart. Thus an algebraic symmetry of a cover does not automatically yield a polynomial automorphism of A².

Useful sidecar: `05-knowledge/reference/CORE-PAPERS-PLANAR-JC.md` distinguishes the actual affine complement, actual line at infinity, and local links. It records the classical exclusion of generic degrees 2, 3 (Orevkov) and 4 (Domrina). These concern d(F), not max(deg P,deg Q).

Recursive inspiration: HYP-9033, `05-knowledge/hypotheses/HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid.md`, retains valuation, norm, index, and chart data in a fixed higher-dimensional tower. We use its research design, not its higher-dimensional claims as input to a planar theorem. HYP-9286, `05-knowledge/hypotheses/HYP-9286-local-system-dictionary-for-tournaments-unit-distances-collatz-and-keller-maps.md`, motivates keeping the permutation local system rather than merely its scalar statistics.

**New contribution here:** a precise distinction between three meanings of recursion, a realizability-preserving pumping criterion, an explicit size/degree consequence, and two non-vacuous local terminal mechanisms. This turns “no hidden high-degree anomaly” into a finite list of identifiable missing theorems.

## 2. Three kinds of recursive structure must not be confused

1. **Resolution recursion.** Successive point blowups produce a decorated boundary graph; Newton–Puiseux expansions describe branches. This is available for each individual curve/map, but its depth and integer decorations need not be uniformly bounded.
2. **Cover recursion.** A block system of the transitive monodromy action produces an intermediate function field and a factorization of finite covers/normalizations. It does not identify an intermediate affine surface as A².
3. **Polynomial recursion.** F=G∘H with G,H polynomial Keller self-maps of A². Here d(F)=d(G)d(H), and the factors are genuine members of the same problem.

Tournament substitution has an honest counterpart only at level 2: a block is a module for a permutation action. Level 3 requires a separate affine-source recognition theorem. Resolution-tree substitution is neither of these automatically.

Do not force an orientation on a dual graph and call it a tournament. The actual transferable structure is substitution with a preserved interface, an associative composition law, and an obstruction predicate. If a tournament comparison is used, its vertices should be genuine modules/blocks and the missing binary observable must be supplied first.

## 3. The full boundary passport

Let S_F be the nonproperness set and U=A²\S_F. Over U the map is a connected finite étale cover of degree d. Normalize A² in C(x,y), obtaining Xbar -> A². By Zariski's Main Theorem the original source embeds as an open subset A² -> Xbar. The inclusion of this particular open source is part of the data.

A usable boundary passport retains:

- A resolution of the target pair (P², closure(S_F) ∪ L_infinity), including the **actual** line at infinity.
- The source compactification and the rational-map resolution, with every exceptional divisor identified as included in or excluded from the affine chart.
- Divisorial valuations of P,Q, a primitive field generator, the Jacobian form, and relevant conductor/index ideals; chart transition functions and unit information.
- The transitive permutation representation π₁(U) -> Sym(Ω), not just conjugacy classes or cycle counts. Boundary meridians, their conjugating paths, and identifications at gluing interfaces must travel with it.
- For each generic affine boundary divisor, the subsets of sheets that survive in the affine source, together with ramification and residue-field degree of the other branches.
- Puiseux/contact data at singular points and infinity, and the actual specialization of surviving sheets there.

The passport is finite for a fixed F, but no uniform bound on its size is presently supplied. A finite list of unlabelled tree shapes does not control arbitrarily large valuations, Puiseux denominators, sheet sets, or gluing monodromy.

## 4. Inherited finite-survivor principle and a terminal certificate

### Proposition 4.1 (inherited mechanism, proof included for the bridge)

For every irreducible target divisor D, a dense open subset of D has a finite affine inverse branch. Consequently the inertia generator σ_D on Ω fixes at least one sheet.

Proof. If D=V(h), dominance implies h∘F is a nonconstant polynomial. It has a nonempty zero divisor in A². Since F is étale, it is quasi-finite, so an irreducible component of that zero divisor maps to a one-dimensional constructible subset of D, hence densely. After removing finitely many bad points of D, an affine inverse branch is unramified and supplies an inertia-fixed sheet. This is the THM-2633 mechanism, not a new claim. ∎

### Proposition 4.2 (normal-inertia terminal test)

Let G act faithfully and transitively on Ω. If σ fixes a point and the cyclic subgroup I=<σ> is normal in G, then σ=1. Thus a Keller passport for which every generic affine boundary inertia subgroup is normal in G has d=1.

Proof. I fixes some ω. For ω'=gω and i∈I,

    iω' = g(g^(-1)ig)ω = gω = ω'.

So I fixes every sheet and is trivial. Apply this to each boundary component. The inherited normal-generation/purity result says these inertia groups normally generate G, so G=1. Equivalently the finite normalization becomes unramified in codimension one, purity makes it finite étale over A², and the trivial étale fundamental group forces degree one. Finally a birational Keller self-map is an automorphism. ∎

This is a local group-theoretic certificate, independently testable from a passport. It is not a conjecture that all actual inertias are normal. A broader recursive program may use this as one kind of terminal leaf without claiming it covers every leaf.

**Hostile control.** In the degree-three action of S₃, a transposition fixes one point, and transpositions normally generate S₃; its order-two subgroup is not normal. Thus “each meridian fixes a sheet” does not prove JC. This abstract control does not purport to realize a planar Keller map.

**Another control.** A cyclic cover t↦t^m of C* has arbitrary degree and one simple repetitive topological pattern. Its inertia is regular and has no fixed sheet. It shows why finite branch shape alone has no degree content and why the affine-survivor mark matters.

## 5. Polynomial factorization is a valid recursion, but needs recognition

### Proposition 5.1 (minimal counterexamples are composition atoms)

If planar JC fails, choose a counterexample of least generic degree d>1. It has no polynomial Keller factorization F=G∘H with d(G)>1 and d(H)>1.

Proof. The field degree is multiplicative. Both factors would have smaller generic degree; by minimality they would be automorphisms, a contradiction. ∎

More generally, a decomposition into polynomial factors of a Keller map automatically has Keller factors: det JF=(det JG)∘H · det JH is a nonzero constant in C[x,y], so both factors on the right are units, and dominance of H makes det JG constant. This does **not** turn an intermediate field into polynomial factors.

### Open bridge A: affine-source recognition for recursive blocks

For a specified class of block decompositions extracted from the full boundary passport, prove that the intermediate cover admits an affine chart isomorphic to A² containing the image of the original affine source, and that both induced morphisms extend as polynomial maps between the designated A² charts.

This is intentionally not asserted for arbitrary block systems. A first meaningful target is a family of boundary subtrees with fixed contact type, primitive valuation lattice, and explicitly trivialized conductor/gluing obstructions. The output must construct two coordinate generators, not merely show rationality or contractibility of the intermediate surface.

**Why this is a substantive bridge:** it compares two different categories, finite covers and polynomial affine-space maps. A rational function field does not by itself identify the coordinate ring of the required chart; a normal affine rational surface is not automatically A².

## 6. A precise finite-state theorem that really controls high degree

The following is an abstract theorem about a proposed encoding, not an assertion that Keller maps already admit it.

### Definitions

A **bounded labelled grammar** for a class of Keller maps gives each map a finite rooted tree T with:

- at most b children per vertex, where b≥1 is fixed;
- a state in a fixed set Q of N elements at each vertex;
- a positive integral local multiplicity w_v≤W, W≥2;
- enough supplementary source-chart, valuation, and gluing data to reconstruct a valid member of the class;
- a proved estimate d(F)≤∏_v w_v.

The phrase “supplementary data” is not a license to hide unbounded degree. The local estimate must remain valid after retaining those data. Numerical decorations can be unbounded only if proved irrelevant to this estimate and to descent, for example because a specified polynomial automorphism absorbs them.

Say the grammar has **badness-preserving pumping descent** if whenever a nonautomorphic map has a repeated state on one root-to-leaf path, the repeated segment produces either:

(P) a polynomial Keller factorization with both factor degrees >1; or

(R) a new nonautomorphic polynomial Keller map F' with 1<d(F')<d(F).

Both alternatives must be realized algebraically, including the affine chart. A formally smaller decorated tree is not sufficient.

### Theorem 6.1 (conditional bounded obstruction theorem)

Assume every planar Keller map admits such a grammar and pumping descent. If JC(2) fails, it has a counterexample of generic degree at most

    B = W^(1+b+...+b^(N-1)).

For b=1 this is W^N. Therefore exclusion of generic degrees 2 through B implies JC(2).

Proof. Take a least-generic-degree counterexample. Alternative (P) is excluded by Proposition 5.1 and (R) by minimality. No root-to-leaf path repeats a state, so every path has at most N vertices. The arity bound gives at most 1+b+...+b^(N-1) vertices. Multiply the local multiplicity bounds. ∎

### What this does and does not say

It proves that high-degree exceptions cannot be irreducibly new **if** the stated geometric realization and bounded-label hypotheses are established. It does not bound degrees of all maps, or even degrees of all possible counterexamples: iterating one counterexample could produce unbounded generic degrees. It bounds the least obstruction and hence reduces exclusion to a bounded generic-degree problem.

Even “bounded generic degree” is not a finite coefficient search. There are continuous moduli, arbitrarily large coordinate degrees under automorphism changes, and possibly unbounded boundary presentations. The terminal exclusion must handle those families uniformly. A second finite-moduli or normal-form theorem would be needed for a literal finite census.

### Concrete new conjecture B: bounded essential passports

For planar Keller maps, after explicitly reversible polynomial source/target automorphism changes, the full boundary passports admit a bounded labelled grammar with badness-preserving pumping descent.

To keep B from becoming a tautological renaming of JC, a successful proof must provide Q,b,W and the local transformation rules before claiming coverage; validate the rules on broad classes of non-Keller generically finite maps as well; and exhibit exactly where the Keller identity supplies descent. If no finite local alphabet survives those controls, abandon this version and retain a well-founded, potentially infinite-state descent program instead.

**A useful split:** B1 is a bounded interface alphabet; B2 is bounded sheet-weight per essential vertex; B3 is algebraic realization of pumping; B4 is strict descent or polynomial factorization. None follows from the others.

## 7. A second route: a potential instead of a finite-state census

A finite-state theorem is only one strategy. Given a directed graph of allowable essential boundary transitions, seek a rational potential h on states and positive cost c(e) such that every transition of a surviving obstruction satisfies

    h(head(e)) - h(tail(e)) ≥ c(e) > 0.

A directed cycle is then impossible by summation. When the state set is finite, this implies bounded depth; with infinitely many weighted states, a bounded potential plus a uniform positive cost still does so. The boundedness and positive gap must be proved, not read off from sampled transitions.

This supplies a concrete place to connect a moment theorem: prove a **signed, realizability-sensitive boundary cost** that survives gluing. Plain counts of branches or discriminant squareclasses lose the source chart and generally do not furnish c(e).

### Open bridge C: Gaussian-to-boundary positivity, dimension pinned

Construct, for each allowable recursive interface in a planar Keller passport, an explicit polynomial observable f_e in exactly two Gaussian variables and an explicit reconstruction identity from its Gaussian moments to a boundary cost c(e). Prove:

1. The observables are honest polynomials in the asserted Gaussian model; poles, extra sheet variables, and conditioning have been eliminated lawfully.
2. The chosen GMC(2) theorem supplies the precise moment implication needed, with all test-polynomial quantifiers preserved.
3. The reconstruction is invariant under the allowed chart changes and sees the marked affine survivors.
4. Costs add under the actual algebraic gluing; either strict positivity or a classified zero-cost move results.
5. Zero-cost cycles are absorbed by polynomial automorphisms or yield a valid smaller polynomial Keller map.

This is a proposed interface, not a claimed implication GMC(2)⇒JC(2). It is stronger and more useful than asserting that a scalar moment is “like” a branch count, but its construction is currently missing. If the construction uses four or more Gaussian variables, it no longer provides the requested GMC(2) bridge without a separate dimension-preserving reduction.

## 8. Exact toy controls and the information they retain

### 8.1 Unbounded coordinate degree without any obstruction

    T_m(x,y)=(x+y^m,y),  m≥1.

The determinant is one, inverse is (u−v^m,v), and generic degree is one. Coordinate degree m is unbounded. Any “no very high degree” statement must specify generic degree or a normalized coordinate complexity. A grammar should absorb these maps as automorphism moves, not count their visible polynomial size as essential obstruction.

### 8.2 Repeated state with unbounded arithmetic label

A one-vertex graph labelled m can encode the covers z↦z^m. The unlabelled state count is one but degree is arbitrary. This refutes any attempt to deduce a degree bound from finitely many unweighted resolution types. B2 is indispensable.

### 8.3 Abstract finite-survivor monodromy is not enough

For every n≥3, the standard action of S_n is transitive, and transpositions fix n−2 sheets while normally generating S_n. Thus the inherited fixed-sheet normal-generation gate allows arbitrarily large groups. A high-degree exclusion needs a genuinely planar algebraic input; group theory alone has not supplied it.

### 8.4 Block systems are not affine charts

An imprimitive cover has an intermediate field by Galois theory. Its normalization is a well-defined finite algebraic object. None of these facts supplies two global coordinate functions generating the needed affine ring. The destructive shortcut is to replace “an intermediate normal surface” by “another A²”; bridge A exists precisely to pay this debt.

## 9. Working concept board and decisive next steps

- **Inertia / normal generation:** exact terminal certificate available. Next test: identify Puiseux/contact configurations forcing a boundary inertia subgroup to be normal in the full monodromy, not merely in a local subgroup.
- **Block substitution / affine recognition:** composition degree is exact. Next test: one explicitly parameterized family of intermediate normalizations; construct its A² coordinates or identify a conductor obstruction.
- **Boundary grammar / pumping:** Theorem 6.1 supplies an honest target. Next test: one repeated blowup/Puiseux motif, with all gluing and source-chart marks retained, and a literal polynomial transformation realizing its removal.
- **Moment cost / gluing:** no transport theorem yet. Next test: two-chart model with a nontrivial transition and a two-variable polynomial observable whose reconstructed cost is invariant under it.
- **Dimension profile:** planar complex geometry, higher-dimensional tower examples, and Gaussian variable count remain distinct. Next test: write the dimension of every source, target, and auxiliary space next to each proposed transform.

The first recommended deliverable is one algebraically realized pumping move, not an enumeration of thousands of unmarked trees. The second is a source-chart recognition theorem for one meaningful family. Those would turn the tournament/fractal analogy into actual reduction machinery.

## 10. Implication architecture, with every missing arrow visible

    GMC(2) in a precisely specified formulation
          + OPEN Gaussian-to-boundary reconstruction C
          + OPEN gluing/zero-cost realization
        -> boundary cost or pumping descent

    OPEN bounded-interface and bounded-weight statements B1,B2
          + realized pumping descent B3,B4
        -> bounded least-counterexample generic degree (Theorem 6.1)

    bounded least-counterexample generic degree
          + uniform exclusion of all resulting terminal families
        -> JC(2).

Separately:

    actual boundary passport with all affine inertia normal
        -> d=1 -> polynomial automorphism.

The latter arrow is already rigorous but is a terminal test, not a coverage theorem. The former chain is a new research program with named, inspectable obligations. It gives a precise interpretation of “there is no anomaly that first occurs only at unimaginably high degree”: every putative obstruction must either descend to bounded complexity or violate one of these explicit structural bridges. No claim about flat-surface classification can substitute for those bridges.

## 11. Concrete all-degree obstruction: one path rule generates every S_d

The finite-state issue can be made exact without a hypothetical Keller example.

### Proposition 11.1 (adjacent-transposition path family)

For every integer d≥3, put Ω_d={1,…,d}, G_d=S_d, and label the edges of the path

    1—2—3—…—d

by σ_i=(i,i+1), i=1,…,d−1. Then:

(a) each σ_i fixes exactly d−2 sheets;

(b) the σ_i generate G_d transitively;

(c) each <σ_i> is nonnormal in G_d, while its normal closure is all of G_d;

(d) deleting the numerical labels and global sheet identities leaves one repeated local edge rule, but d is unbounded.

Proof. (a) is immediate. Any transposition (i,j) is obtained from adjacent swaps by conjugation along the path; all transpositions generate S_d, proving (b). Conjugating σ_i by a permutation sending its unordered pair to a different pair gives a different transposition, outside its order-two subgroup. All transpositions are conjugate, so its normal closure is S_d. This proves (c). The construction extends by appending a single edge and sheet at each stage, proving (d). ∎

This is not a realization theorem for algebraic covers, much less Keller maps. It proves a specific negative result about the **information content** of an unlabelled recursive grammar: arbitrary sheet degree is compatible with the same adjacent-swap motif, finite-survivor condition, transitivity, and fixed-point normal generation. It does not prove that finite-state grammars are impossible. It proves that a valid degree-bounding grammar must include the labels/sidecars or possess an additional algebraically realized descent rule.

The local pattern does not become harmless merely because the global graph is a tree, a path, bounded valence, or recursively generated. In this control the obstruction to the normal-inertia terminal test remains at every edge and every degree.

## 12. Exact composition of sheet passports and marked affine survivors

Fix transitive groups H≤Sym(n) and K≤Sym(m), with m≥2. The imprimitive wreath product H wr K acts on pairs (i,j), 1≤i≤m, 1≤j≤n. Write an element in the convention

    τ(i,j)=(h(i), g_i(j)),   h∈K,  g_i∈H.

### Proposition 12.1 (local wreath ledger)

Its fixed-sheet count is exactly

    |Fix(τ)| = Σ_{i∈Fix(h)} |Fix(g_i)|.

More generally, if (i_0,…,i_(ℓ−1)) is an h-cycle, define the return permutation

    R = g_(i_(ℓ−1)) ∘ … ∘ g_(i_0).

An R-cycle of length k produces a τ-cycle of length ℓk. The resulting cycle lengths are independent of which block in the h-cycle is used as starting block, because its return permutation changes by conjugacy.

Proof. A fixed pair must have h(i)=i and g_i(j)=j. After ℓ applications, the block first returns and the inner sheet has moved by R. It returns to its original pair after exactly ℓk steps. ∎

Thus inertia data do compose exactly. They require the full return permutations and gluing labels; the separate unordered cycle-count histograms are not sufficient to multiply the g_i correctly.

For an **actual polynomial factorization** F=G∘H, the generic sheets of F are grouped by sheets of G. Around a chosen target divisor, retain A, the set of outer inverse branches that remain in the designated affine middle chart. For i∈A let B_i be the inner sheets that remain affine along that actual outer branch. Then the finite affine survivors of the composite are exactly

    Survivors(F) = disjoint union_{i∈A} ({i}×B_i).

In particular, |Survivors(F)|=Σ_{i∈A}|B_i|. This is simply the fibre identity

    F^(-1)(y) = disjoint union_{q∈G^(-1)(y)} H^(-1)(q)

with branch labels retained. It applies at the generic point under consideration after discarding special points where the local bookkeeping changes. It does not infer that every inertia-fixed sheet is affine: an unramified branch in the finite normalization may still lie outside the designated affine source.

For abstract covers the same rule defines a marked-sheet composition only when the chart/survivor markings and gluing are supplied. Arbitrary marks may fail algebraic realizability. Accordingly, the executable mark example below is expressly an information-loss control, not an invented Keller example.

## 13. An exact wreath theorem: the fixed-sheet gate sees only the outer action

For a finite permutation group A, define

    N_fix(A) = subgroup generated by {a∈A : Fix(a) is nonempty}.

This subgroup is normal because the generating set is conjugacy-invariant. The inherited Keller condition is N_fix(A)=A.

### Theorem 13.1 (fixed-sheet normal generation under full substitution)

For the full imprimitive wreath action above,

    N_fix(H wr K) = H^m ⋊ N_fix(K).

Equivalently the maximal quotient left by the fixed-sheet gate is

    (H wr K)/N_fix(H wr K) ≅ K/N_fix(K).

Proof. Any element of H supported in one block fixes all sheets outside that block, since m≥2. These elements generate the full base group H^m, so H^m⊆N_fix(H wr K). Every h∈K having a fixed block has a pure outer lift that fixes all n sheets in that block; these lifts generate a copy of N_fix(K). This gives one containment. Conversely, if τ has a fixed sheet, its outer action fixes the containing block. Hence the projection to K of every generator of N_fix(H wr K) lies in N_fix(K), giving the reverse containment. ∎

**Consequence 13.2.** Full substitution can hide arbitrarily complicated inner monodromy from the fixed-sheet normal-generation test. In particular H wr S₃ satisfies N_fix=the whole group for every finite transitive H, because transpositions generate S₃ and fix a block in its degree-three action. Taking H cyclic of arbitrarily large order already gives unbounded sheet degree with an inner regular action that would fail the Keller gate if it were realized as a separate Keller factor.

This is a strong warning with a positive use: a proof that recognizes the inner block as an actual polynomial Keller factor would immediately exclude this entire cyclic-inner family. The recognition theorem, not another scalar monodromy statistic, is the missing bridge.

**Exact degree-six contrast.** In their natural actions:

- C₂ wr S₃ has order 48 and N_fix of order 48, so passes the gate.
- S₃ wr C₂ has order 72 and N_fix of order 36, so fails the gate; the residual quotient is the outer C₂.

The order of substitution therefore matters even when total degree is the same. An unordered list of factors, branch multiplicities, or recursion depths cannot replace the ordered block tree.

The theorem is about the **full** wreath product. Actual composite monodromy can be a proper subgroup. For such a subgroup, the fixed-count and cycle-return formulas remain exact, but the equality for N_fix need not hold because the individual supported base generators may not be present.

## 14. Positive class beyond individually normal inertia: nilpotent block quotients

### Lemma 14.1 (nilpotent actions cannot be fixed-sheet normally generated)

If a nontrivial finite nilpotent group A acts faithfully and transitively on a set with more than one element, then N_fix(A) is proper in A.

Proof. Choose a point stabilizer J<A. Every element fixing some point lies in a conjugate of J, so N_fix(A) is the normal closure of J. Choose a maximal proper subgroup M containing J. Every maximal subgroup of a finite nilpotent group is normal, so the normal closure of J lies in M and is proper.

For completeness the normal-maximal fact follows from the normalizer condition. Every proper subgroup B of a finite nilpotent group has N_A(B)>B: induct on the nilpotence class using the nontrivial center Z(A). If Z(A) is not contained in B, a central element outside B enlarges its normalizer. If Z(A)⊆B, apply the inductive assertion to B/Z(A) in A/Z(A) and lift. Applying this to maximal M forces N_A(M)=A, so M is normal. The abelian base case is immediate. ∎

### Proposition 14.2 (Keller exclusion for nilpotent block quotients)

A nontrivial Keller monodromy action cannot be nilpotent. More strongly, any nontrivial transitive action induced by Keller monodromy on a system of blocks has nonnilpotent image.

Proof. The inherited theorem gives N_fix(G)=G. For a block system, a permutation fixing a sheet necessarily fixes the block containing that sheet. Therefore the induced image K is generated by fixed-block elements, so N_fix(K)=K. Lemma 14.1 excludes a nontrivial nilpotent K. Taking singleton blocks proves the first assertion. ∎

This requires no assertion that individual inertia subgroups are normal, so it is a useful additional terminal class. It also avoids all intermediate-chart recognition: a forbidden quotient of the actual monodromy is sufficient. It is an elementary corollary of the inherited normal-generation theorem, not a claim of a new literature-first exclusion.

Examples include any nontrivial p-group monodromy, cyclic regular endpoint, or recursion with a nontrivial nilpotent **outer block action**. Inner nilpotent blocks do not suffice, as C₂ wr S₃ demonstrates. For a general block subgroup the condition concerns the faithful image on blocks, not the possibly nonfaithful group before quotienting its kernel.

## 15. Executable finite controls and independent all-degree proofs

Files:

- `permutation_recursion_controls.py`
- `permutation_recursion_controls.json`

Reproduction:

    python3 topology/permutation_recursion_controls.py
    python3 -O topology/permutation_recursion_controls.py

All checks use explicit `require` with `RuntimeError`, never Python `assert`; optimization cannot suppress them. Both runs were executed and their output compared byte-for-byte.

The finite universe consists of adjacent-transposition groups in degrees 3 through 7, both indicated degree-six wreath actions, the order-eight dihedral group on square vertices, and one abstract marked-survivor control. The outputs include:

- Path groups of orders 6,24,120,720,5040, each with d−2 fixed sheets per generating edge.
- Wreath fixed-generation orders 48 and 36 as above; the fixed-count formula tested on every group element.
- Dihedral square action: total group order 8, fixed-generated subgroup order 4. This is a nonregular nilpotent action excluded by Proposition 14.2.
- A formal mark set {2} strictly smaller than the inertia-fixed set {2,3,4,5}. This checks that markings are not silently equated with fixed points; it claims no geometric realization.

The proofs in Sections 11–14 are independent of these finite tests and establish their asserted all-degree statements. No finite computation is used as evidence that every monodromy group is classified, that any abstract passport is realizable, or that the planar conjecture is solved.

## 16. A literal tournament realization: Paley substitution masks regular modules

The even-order groups in Sections 11–15 are arbitrary permutation controls, not tournament automorphism groups. Indeed every finite tournament automorphism group has odd order: an involution would exchange some pair of vertices and hence reverse its directed edge; Cauchy's theorem then excludes even group order. To make the tournament connection intrinsic rather than cosmetic, use the following odd-order construction.

### Theorem 16.1 (Paley–cyclic module masking)

Let q≥7 be a prime power with q≡3 mod 4. On F_q define the Paley tournament P_q by

    x -> y  iff  y−x is a nonzero square.

Let r≥3 be odd, and let T_r be the cyclic tournament on Z/rZ with

    i -> j  iff  j−i mod r ∈ {1,…,(r−1)/2}.

Let K_q be the affine subgroup

    K_q = {x -> a x+b : a∈(F_q*)², b∈F_q},

and let C_r act on T_r by translation. Then:

1. K_q acts transitively by automorphisms of P_q, has odd order q(q−1)/2, and satisfies N_fix(K_q)=K_q.
2. The regular action of C_r has N_fix(C_r)=1.
3. The lexicographic tournament P_q[T_r] admits the explicit transitive automorphism subgroup

       W = C_r wr K_q,
       |W| = r^q q(q−1)/2,
       N_fix(W)=W.

Thus literal recursive tournament substitution can conceal a regular, fixed-point-free inner action inside an action wholly generated by elements fixing vertices.

Proof. Since q≡3 mod 4, −1 is not a square, so exactly one of y−x and x−y is a square. Multiplication by a square preserves this distinction, as does translation. This proves the tournament and group-action assertions; transitivity follows from translations.

Every affine map with a≠1 fixes exactly one point, b/(1−a), while x↦a x fixes the origin. There exists a square a≠1 because (q−1)/2≥3. Put s_a(x)=a x and t_b(x)=x+b. The product of the two fixed-point elements t_b s_a t_(−b) and s_a^(-1) is translation by (1−a)b. As b ranges over F_q, these products give every translation. All square multipliers themselves fix the origin. Hence fixed-point elements generate all of K_q. In the regular cyclic action only the identity fixes a point, proving assertion 2.

For the lexicographic tournament the vertices are (x,i). Between distinct outer vertices x,y the orientation is inherited from P_q, independent of i,j; within one block {x}×Z/rZ it is inherited from T_r. Independent cyclic translations within each block and every outer K_q automorphism preserve all these edges. They generate the semidirect product C_r^q⋊K_q with the stated order and faithful transitive action. Theorem 13.1 gives N_fix(W)=W. ∎

The directed relation is intrinsic here: vertices, oriented pairs, and modules are all explicit. Each block is a genuine tournament module because every outside vertex has the same edge orientation toward all vertices of that block. The operation is literal lexicographic substitution. Nothing has been manufactured by arbitrarily orienting a resolution graph.

### Small exact control

For q=7 and r=3, the outer squares are {1,2,4}. The outer group has order 21 and is generated by x↦x+1 and x↦2x. The tournament P_7[C_3] has 21 vertices, and the explicit wreath subgroup has order

    3^7·21 = 45,927.

The script checks the 21-element outer group and its fixed-generated subgroup, the regular inner C₃ and its trivial fixed-generated subgroup, the orientation of every pair of distinct vertices in the 21-vertex lexicographic tournament, and edge preservation by each of nine wreath generators. The subgroup order is derived from the semidirect product, and its fixed-generation conclusion from Theorem 13.1; the script does not pretend to enumerate all 45,927 elements.

The normal and optimized runs were repeated after adding this control and their JSON outputs again agree byte-for-byte.

### Exact transfer ledger

- **Source:** lexicographically substituted tournaments with their actual automorphism action and module partition.
- **Target:** abstract transitive marked-cover permutation passports with a block partition.
- **Map:** retain the group action and block system; forget the tournament edges.
- **Preserved:** group composition, permutation degree, fixed-point counts, normal generation, and the wreath-module structure.
- **Lost:** the orientation relation; in the reverse direction, algebraic realizability, inertia placement, affine source charts, ramification, and the actual marked survivor subsets are absent altogether.
- **Conclusion transported:** recursive substitution can mask a regular inner factor from an outer fixed-point-normal-generation test. This is a proved obstruction to that particular classification strategy.
- **Conclusion not transported:** existence of a polynomial Keller map, or a procedure for removing its inner module as a polynomial factor.

A useful positive target is now sharply localized: prove affine-source recognition for a class of algebraic cover modules resembling these homogeneous substitution blocks. If a regular inner block genuinely becomes a nontrivial polynomial Keller factor, the inherited fixed-sheet gate excludes that factor. The tournament example explains why the module must first be recognized algebraically; it does not provide that recognition.
