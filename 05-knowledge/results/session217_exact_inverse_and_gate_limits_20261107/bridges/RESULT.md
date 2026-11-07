# Adaptive uniform-module gates reduce to multiplication

Status: **PROVED deduction from inherited run-support and substitution laws; FINITE-EXACT controls.** This is a structural reduction, not a construction proving Chan–Pak completeness, and not a literature-priority claim. No new theorem ID is reserved.

## Scope

A template is a fixed finite tournament Q with k distinct labeled holes. Substitute arbitrary nonempty tournaments U_1,...,U_k into those holes, with uniform external incidences, and leave the other vertices singleton. A menu is any fixed set of such templates, finite or infinite. Selection may depend arbitrarily on the entire input graphs, not merely their Hamiltonian-path counts. Every input occurs once in every template. Empty inputs, deletion of inputs, nonuniform incidences are outside this theorem.

Write h_i=H(U_i). A scalar gate F is implemented by the menu if, for every input tuple, at least one menu output has count F(h_1,...,h_k). There is no regularity or algorithmic restriction on the selector.

**Adaptive-menu theorem.** If a menu, even an infinite one, implements F, then

    F(h_1,...,h_k) belongs to {C_j product_i h_i : j is a safe template}

for every attainable positive scalar tuple. A template is safe precisely when every hole is a singleton strongly connected component of its quotient. Here C_j=H(Q_j). Conversely, a scalar function that pointwise selects these safe-template values is implementable by selecting the corresponding safe template.

In particular no menu can implement either

* F(x,y)=x+y+1, or more generally ax+by+c with a,b>0, on all nonempty tournament inputs when the output is required to be positive; or
* F(x)=ax+b with b nonzero, on every attainable x, whenever these outputs are positive.

The plus-one example avoids the irrelevant parity objection to x+y. The obstruction is not merely parity or a same-order seed obstruction.

## Proof

### Inherited one-hole test

For a rooted quotient (Q,v), let c_r count quotient skeletons with r occurrences of v. The inherited positive run expansion is

    H(Q[U at v]) = sum_(r>=1) c_r F_U(r),

where F_U(r) counts ordered spanning covers by r nonempty directed paths. In particular c_1=H(Q). Session216's matching-support theorem says c_2>0 iff the root out-to-in bipartite graph has an edge. This condition is equivalent to v lying in a nontrivial strong component: if the graph has no edge then N^-(v) => v => N^+(v), with every in-neighbor beating every out-neighbor, so v is singleton in the condensation. Conversely a triangle through v puts it in a nontrivial component. This establishes both directions without needing a new Hamiltonian-cycle theorem.

If v is singleton in the condensation, every Hamilton path uses the substituted module in exactly one run. If not, c_2>=1. Since F_TT_m(1)=1 and F_TT_m(2)=2^m-2, we have the exact two-point criterion

    H(Q[TT_2 at v]) = H(Q) iff v is a singleton strong component,

and in the other case H(Q[TT_m at v])>=H(Q)+2^m-2.

This fixed-template criterion is an immediate corollary of the Session216 theorem, explicitly not a new independent result.

### Padding that preserves arbitrary scalar inputs

For a nonempty U and m>=0, put U^(m)=U => TT_m. The SCC/ordinal product law gives H(U^(m))=H(U). Let

    Q_U(2)=sum_(S subset V(U)) H(U[S]) H(U[V(U) minus S]),

using H(empty)=1. An ordered bipartition of an ordinal join splits independently in its two factors, hence

    Q_(U^(m))(2) = 2^m Q_U(2),
    F_(U^(m))(2) = 2^m Q_U(2)-2H(U)
                      >= (2^(m+1)-2) H(U).

The last inequality uses the two partitions with an empty side. Thus this padding holds H fixed while forcing the second cover coordinate to infinity.

### Elimination of every unsafe template

Fix any input tuple U_i and pad every input simultaneously by TT_m. In a template with an unsafe hole i, retain only skeletons visiting hole i twice and every other hole once. Such skeletons exist: they are exactly the positive c_2 skeletons of the one-hole quotient at i, with the other holes regarded as ordinary vertices. The multivariate positive run expansion therefore gives

    H(Q_j[U_1^(m),...,U_k^(m)])
      >= c_2(i) F_(U_i^(m))(2) product_(l!=i) h_l
      >= (2^(m+1)-2) product_l h_l.

This tends to infinity. The desired scalar F(h_1,...,h_k) is fixed under padding. The displayed lower bound is **uniform over every quotient in the menu**: it uses only c_2(i)>=1. Therefore every unsafe template exceeds this target simultaneously for all sufficiently large m, even if the menu is infinite. At least one safe template must then realize it. For a safe template the SCC product law gives exactly C_j product h_i. This proves the claimed pointwise safe-template restriction. In particular product_i h_i divides F(h_1,...,h_k), and its quotient must itself be a tournament count. A finite menu additionally restricts this quotient to a fixed finite set. Conversely these safe templates do implement their indicated products.

### Why affine addition does not escape by adaptive selection

The attainable scalar spectrum is unbounded: ordinal sums of directed triangles have H=3^t. In one variable the divisibility consequence says x divides ax+b. For an affine map with integer a,b and b nonzero, take x=3^t>|b|; this is impossible. More generally, for real a,b, the ratios a+b/3^t must be integers for all t. They converge to a, so they are eventually constant integers, forcing b=0.

For F(x,y)=ax+by+c with a,b>0 fixed real constants, choose x=y=3^t. Every safe output is at least xy, whereas (a+b)x+c<x^2 for sufficiently large t. Thus no menu can implement the affine addition. In the parity-compatible special case x+y+1, x=y=3 already asks for 7, below the minimum safe output 9. This is an obstruction to the gate, not an assertion that its target values are all missing from the unrestricted spectrum. For example x=y=9 asks for the attainable value 19 (witness outgoing masks [40,1,3,6,15,30]) but the gate's safe outputs are at least81.

## What the reduction changes, and what it does not

The inherited fixed-root barrier excludes scalar-only arithmetic in a single template. The deduction above excludes an arbitrary *adaptive* bank, even an infinite one, when it can inspect full inputs and deliberately choose a different template after every padding. This is a stronger quantifier order than ruling out a fixed universal context.

A constructive arithmetic program must escape at least one explicit premise: retain additional endpoint/path-cover coordinates; use nonuniform external incidences; omit or duplicate input modules rather than using each exactly once; or restrict the permissible input family so ordinal padding is forbidden. The three-core binary construction is **not excluded**: its chosen cuts give nonuniform incidence to the apex, and it does not claim to transform arbitrary inputs from their scalar H alone. Existing O/Z/J carrier searches also retain extra state and are not subsumed by this no-go.

No positive construction reaching new scalar counts or improved occupancy bound was obtained in this lane. Cycle/path conversion and scalar semigroup ideas did not produce a verified escape, so they are not advertised as progress toward completeness.

## Reproduction and independent controls

Run `python3 verify.py > VERIFY.json` and `python3 -O verify.py > VERIFY_OPT.json`, then compare outputs. The script exhausts all 5405 rooted labeled tournaments of orders 1 through 5 and all 16215 literal TT blowups of sizes 1 through 3. It separately computes Hamiltonian paths by vertex-subset DP, the full one-hole run kernel by permutations and gap permissions, and singleton-SCC status by transitive closure. It checks the exact expansion, two-point criterion, and exponential lower bound. There are 4190 unsafe rooted contexts. Explicit exceptions, not removable assertions, enforce controls.

The arbitrary-input padding identity and multivariate argument are proved above; finite verification is not substituted for their quantifiers. The verifier additionally checks both identities for arbitrary-input ordinal padding on every labeled tournament of orders 1 through 4 and padding lengths 0 through 3 (300 checks).

The separate `verify_multi.py` constructs literal two-hole substitutions with simultaneous padding. Its 976 cases cover all quotients of orders 2 and 3, sixteen specified order-4 quotients, every pair of holes, singleton/C3 input pairs, and padding lengths 0 and 1. It checks both safe factorization and the uniform unsafe lower bound with the separate one-run baseline retained. Normal and optimized outputs match. This independently exercises the multivariate use of the run kernel; the arbitrary-quotient-size quantifier remains supported by the proof's integer lower bound c_2>=1.

## Inheritance ledger and concept board

* Closest proved mechanism: Session216 `../../session216_overnight_constructions_20261107/star/PROOF.md`, Theorem 1 (matching-number run support). The new deduction begins at arbitrary adaptive selection, not at c_2 positivity.
* Complete substitution kernel: `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md` (ordered run-cover expansion).
* Prior scalar failure: `01-canon/theorems/THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant.md` (C3 scalar failure); its historical open exact-kernel wording is superseded by THM-3121.
* Ordinal product law: `01-canon/theorems/THM-1860-c3-le-H-reduced-to-strongly-connected.md`, retained SCC section. No dependence on its explicitly superseded historical open triangle-bound discussion.
* Canonical hostile: a directed triangle rooted at any vertex; TT_1 and TT_2 both have scalar H=1 but outputs are 3 and 5.
* Corrected near miss: a selector can switch templates, so refuting each fixed gate separately does **not** suffice. Simultaneous padding makes all unsafe outputs exceed one fixed target at once.
* Least-used sidecar: ordinal padding preserves arbitrary scalar H while its second cover coordinate diverges.

Concept board: scalar multiplication (safe SCCs), additive arithmetic (blocked by multiplicative divisibility), adaptive templates (uniform lower bound defeats even unbounded menus), full cover profiles (padding separates them), and nonuniform endpoints (permitted escape, no new construction here).
