# Session 58: honest point-count and Frobenius bridges

Status: PROVED below from the explicitly named inherited counting theorem and standard rationality; FINITE-EXACT controls are separate. No F1 descent or Riemann Hypothesis result is claimed. Prepared 2026-10-25.

## 1. Inheritance and live concept board

Closest proved mechanism: Session 57, `prime_plane_frames.md`, sections 4–5: for odd prime p, a fixed affine F_p-plane A in the additive group of F_q, q=p^n, n>=2, has C_p q unordered unique-sum p-by-p frames, C_p=(p+1)p^(p-1)-p(p+1)/2, and C_p(q-p²) disjoint frames. Each gauge-equivalence class contains q frames, and precisely p² forbidden translations. Choosing an orientation is essential: translation is not an intrinsically defined action on unordered pairs.

Canonical hostile: p=2, where 2t=0 and the disjointness formula fails. Corrected near miss: replacing a varying F_p-vector dimension by field extension of a single variety without descent data. Least-used relevant sidecar: the permutation of deleted points under Frobenius.

Board: (1) additive gauge classes and forbidden translations; (2) fixed-base extension counts; (3) Frobenius cycles on finite deleted sets; (4) converse involution on tournament types; (5) polynomial absolute zeta versus actual descent; (6) simultaneous Frobenius lifts.

Read-only corpus search of canon, hypotheses, results, and reflections for F1/absolute zeta/Borger/Deitmar/Soulé/Guinand located no developed absolute-geometric carrier relevant to this count. A historical reflection, `07-reflections/the-fundamental-truths-in-every-language.md`, says that every tournament over F1 is transitive, referring to HYP-1539. This supplies no definition/base-change construction and is not inherited as a theorem. `THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md` is relevant for signed/fixed-point bookkeeping, not an F1 model. Its 'blueprints' paragraph concerns unrelated user attachments, not Lorscheid blueprints.

## 2. Fixed-base frame realization (PROVED)

Fix k=F_Q, Q=p^d, d>=2, and A=C+U contained in k, where U has F_p-dimension two. For each of the C_p unordered gauge classes choose one ordered representative (B_i,D_i) with elements in k. Set

S_i = { (d-b)/2 : b in B_i, d in D_i },
X_i = A¹_k minus S_i.

The unique-sum property implies |S_i|=p². Indeed, d-b=d′-b′ implies b+d′=b′+d, and uniqueness of representations forces b=b′ and d=d′. Each point is k-rational. The map t -> {B_i+t,D_i-t} identifies X_i(F_{Q^m}) with the disjoint frames in that gauge class over F_{Q^m}. For completeness, the full free-gauge argument is as follows. A translation stabilizing both ordered factors with parameter u!=0 makes u a nonzero element of both B_i-B_i and D_i-D_i; this contradicts unique sums, since those difference sets intersect only at zero. If two parameters t and s give the same unordered pair through swapping, then B_i+t=D_i-s, so D_i is a translate of B_i. Their difference sets would then coincide and contain a nonzero difference (p>1), again contradicting uniqueness. Therefore t=s in the nonswapped case and the swapped case is impossible. Equivalently, translations and factor-swap act freely on ordered frames, with classes of size 2|k| and unordered classes of size |k|. The same argument works over every extension field. Every extension-field frame belongs to a base-field shape class because after subtracting one element from each factor, both factors lie in the fixed U. Hence the disjoint union X of the X_i counts all the frames:

#X(F_{Q^m}) = C_p(Q^m-p²),
Z(X,t) = ((1-t)^(p²)/(1-Q t))^(C_p).

In fact the deleted configurations can be made identical. Normalize every representative so that 0 belongs to B_i. Then B_i is contained in U and D_i is contained in A, so D_i-B_i is contained in A; its cardinality p² forces D_i-B_i=A. Therefore every S_i equals A/2, and X is noncanonically a disjoint union of C_p copies of the single complement A¹_k minus A/2. This stronger statement retains the representative-choice caveat.

This is a chosen parameter-space realization, not a claimed canonical moduli scheme. Different gauge choices change charts by affine transformations; no intrinsic unordered additive-group action is asserted. Sum-label, square, or ambient-graph restrictions are NOT included unless independently proved invariant and imposed geometrically.

Alternative repair over the prime field: freeze the additive ambient space V=F_p^d and its fixed plane. Realize every gauge class as A^d_{F_p} minus p² F_p-rational points. Over F_{p^m} it has p^(dm)-p² points and zeta (1-t)^(p²)/(1-p^d t). This is a dimension-d variety, not a dimension-one variety over F_p. The identical displayed rational functions do not identify their bases or geometry.

## 3. The stronger varying-dimension obstruction (PROVED, using zeta rationality)

For an odd prime p and positive integer C, no finite-type variety Y/F_p satisfies

#Y(F_{p^n})=C(p^n-p²) for every n>=2.

Proof: write a=#Y(F_p)-C(p-p²). The defining exponential series gives

Z(Y,t)=((1-t)^(p²)/(1-pt))^C exp(a t).

Zeta rationality makes exp(a t) rational. If a nonzero rational function R has R'/R=a constant, comparison at infinity gives R'/R=O(1/t), so a=0. Thus #Y(F_p)=C(p-p²)<0, a contradiction. This also excludes finite-type schemes by reduction and stratification. More generally, changing finitely many early extension counts cannot repair an eventually impossible exponential-polynomial count: the correction is exp(P(t)), rational only for constant P, and here P(0)=0.

This is stronger than observing the missing n=1 case. Even declining to specify that case cannot repair the proposed geometric interpretation.

Concrete Frobenius-stable repair: remove the roots of x^(p²)-x from A¹/F_p. These consist of p rational closed points and (p²-p)/2 degree-two closed points. Therefore

#Y(F_{p^n})=p^n-p^gcd(n,2),
Z(Y,t)=(1-t)^p (1-t²)^((p²-p)/2)/(1-pt).

The desired p^n-p² appears only for even n. For p=3, n=3 the correct count is 24, not 18. Base extension to F_{p²} splits every puncture and recovers the fixed-base formula.

## 4. Converse as Frobenius: a genuine but deliberately limited bridge (PROVED)

For any finite set Omega with involution c, let N=|Omega| and f=|Fix(c)|. Define

Z_c(t)=exp(sum_{m>=1} |Fix(c^m)| t^m/m).

Since fixed counts alternate between f and N,

Z_c(t)=(1-t)^(-f)(1-t²)^(-(N-f)/2).

For every prime power Q this is the Hasse–Weil zeta of the explicit finite étale F_Q-scheme consisting of f copies of Spec F_Q and (N-f)/2 copies of Spec F_{Q²}. Its geometric Frobenius permutation is c up to isomorphism. Conversely these Frobenius orbit lengths recover N and f, and Burnside gives |Omega/<c>|=(N+f)/2.

Take Omega to be tournament isomorphism types of a fixed order and c to be converse. This gives a literal arithmetic realization of the complement-orbit census. For labelled tournaments with at least two vertices, converse has no fixed point; for unlabelled types it can have many. They must never be interchanged.

Exact enumeration through n=5 gives (N,f,number of 2-cycles): (1,1,0), (1,1,0), (2,2,0), (4,2,1), (12,8,2). At n=4 the zeta is (1-t)^(-2)(1-t²)^(-1). This bridge is a standard finite-permutation/finite-étale construction applied to the census, not a novelty claim about the general construction.

Common mechanism with section 3: a finite Frobenius permutation with a_j cycles of length j has finite-set zeta product_j (1-t^j)^(-a_j). Removing it from affine d-space reverses these factors and adds (1-Q^d t)^(-1). The missing sidecar is exactly the cycle index, not a vague geometric analogy.

## 5. What an absolute-zeta interpretation does and does not provide

For a fixed formal counting polynomial N(T)=C(T-r), its polynomial absolute-zeta transform is s^(Cr)/(s-1)^C. The normalization follows directly by replacing t by Q^(-s) in ((1-t)^r/(1-Qt))^C and multiplying by (Q-1)^(C(1-r)) before the formal limit Q->1. For the frame chart r=p², p stays fixed in this formal calculation; it cannot simultaneously vary as the residue characteristic.

N(1)=C(1-r)<0 rules out interpreting N(1) as the cardinality of an ordinary finite set. It does not rule out signed Euler characteristics, virtual classes, open varieties, or all F1 frameworks. Conversely, polynomial counts alone do not supply an F1 model. The fixed-characteristic chart supplies neither a single integral family with uniform all-prime counts nor descent data.

A particularly cheap lambda-structure probe: lift the p=3 grid to the nine integer points {0,1,2}² in A²_Z. The usual coordinatewise squaring map sends (-1,0), outside the grid, to (1,0), inside. Therefore this naive Frobenius lift does not restrict to the complement. This refutes only that candidate lift, not existence of all possible lambda structures. The obstruction is backward stability of the deleted locus.

## 6. Source–target / loss ledger

- Frame -> punctured affine chart: t maps to the translated pair. Preserves uniqueness and disjointness with fixed A; loses canonical orientation and moduli functor. Sidecar: chosen representatives, characteristic, base field, and square restrictions. Test: enumerate forbidden differences and ensure cardinality p².
- Varying dimension -> one F_p-variety: proposed count sequence, no valid map. Rationality plus n=1 negativity refutes it. Repair: freeze dimension, or enlarge the base and restrict extension degrees.
- Frobenius finite set -> zeta: cycle-index product. Preserves all extension fixed counts; loses positions and incidence of punctures. Sidecar: embedded configuration. Test: same number of punctures with different cycle lengths.
- Tournament types -> finite étale scheme: converse cycles become degree-one/two closed points. Preserves exactly (N,f), extension traces, complement-orbit count; loses arcs, Hamiltonian paths, automorphism groups, metagraph adjacency, and all non-converse structure. Test: identity versus transposition on a two-element set has the same N but different zeta.
- Counting polynomial -> absolute rational function: formal normalized limit. Preserves coefficients/virtual weights; destroys characteristic dependence and does not manufacture descent, intersection positivity, or a Riemann-zeta spectral realization.

## 7. Hypotheses/tests and outcomes

H1: A single F_p-variety realizes the entire Session-57 sequence for n>=2. REFUTED by section 3; at p=3,C=30 rationality would force -180 rational points.

H2: A fixed-base extension-compatible realization exists. PROVED in section 2. Positive control: Q=9 gives 30(9^m-9), m>=1; base count zero is permitted. Hostile: the F_p model with F_{p²} punctures has parity-dependent counts instead.

H3: Complement counts admit a genuine Frobenius realization. PROVED in section 4 and FINITE-EXACT through order five. Hostile: using labelled fixed-point counts in the unlabelled formula fails already at n=2.

H4: The obvious power-map lambda structure survives grid deletion. REFUTED by the explicit square-map witness. Strongest survivor: the closed points remain a useful deletion arrangement but an all-prime lift needs new construction.

H5: Negative N(1) disproves F1 geometry wholesale. REFUTED as an inference: the conclusion only applies to a literal finite-set cardinality. An actual descent claim remains OPEN/unconstructed.

Reproduction: `python3 verify_bridges.py` creates `verification.json`, using only the Python standard library. All verification gates raise explicit exceptions and remain active under Python -O. It checks p=3,5,7 and extension degrees 1–8; exhausts all labelled tournament masks and all vertex relabellings through order five, then directly iterates converse; checks the lambda hostile. Algebraic proofs above, not finite tests, establish the unrestricted claims.

## 8. Primary-source pins, checked 2026-10-25

- Deitmar–Koyama–Kurokawa, *Counting and zeta functions over F1*, section 1, https://arxiv.org/html/1306.3681v2 . Imported: local zeta definition, uniform all-prime Soulé condition, coefficient product for polynomial absolute zeta. Their torsion extension is broader than naive polynomial specialization. No theorem here establishes frame descent.
- Borger, *Lambda-rings and the field with one element*, introduction p.2, https://arxiv.org/pdf/0906.3146 . Imported: for schemes flat over Z, lambda structure is expressed by commuting prime-indexed Frobenius lifts. One finite-field Frobenius is not this datum. The RH discussion is motivation, not a proof.
- Deitmar, *F1-schemes and toric varieties*, abstract/version record, https://arxiv.org/abs/math/0608179 . Scope: monoid-scheme model and toric relationship. No universal exclusion theorem for other F1 models is imported.
- Lorscheid, *Blueprints – towards absolute arithmetic?*, https://arxiv.org/abs/1204.3129 . Checked bibliographic scope only; no technical blueprint theorem is used.
- Connes–Consani, *Geometry of the arithmetic site*, https://arxiv.org/abs/1502.05580 . Their abstract identifies arithmetic-site geometric points and Frobenius correspondences with an adelic construction related to the completed Riemann zeta. This is a specific geometry, not a claim that every counting polynomial proves RH. No positivity theorem is imported here.
- Stacks Project, section 58.2, https://stacks.math.columbia.edu/tag/04JI . Imported: étale schemes over a field correspond to continuous Galois sets. Section 4 also proves its instance directly with finite fields.
- Dwork, *On the rationality of the zeta function of an algebraic variety*, Amer. J. Math. 82 (1960), 631–648; primary journal/translation record https://www.mathnet.ru/eng/mat201 . Imported classical theorem: zeta rationality. Bibliographic record checked; no new audit of Dwork's proof claimed.

The RH boundary is exact: these constructions give only finite permutation eigenvalues and affine Tate factors. They supply no identified operator whose full spectrum is the nontrivial Riemann zeros, no trace identity for the completed Riemann zeta, and no positivity mechanism forcing real part 1/2. The user's geometric direction is a legitimate research program, but these particular counts test its prerequisites rather than solve RH.
