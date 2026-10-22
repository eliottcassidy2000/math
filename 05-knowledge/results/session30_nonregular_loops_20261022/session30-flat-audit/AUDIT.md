# Independent audit: flat almostregular formal universality

## Verdict and scope

**PASS**. The harmless boundary correction identified below has been applied. The claimed coefficient 65/8 is justified for the complete finite second-order predictor. This is not a theorem about the exact Hamiltonian count or an infinite sum of cores. No missing invariant enters the displayed second-order coefficient under the stated hypotheses.

The producer theorem was independently checked by direct asymptotic derivation, new exact symbolic arithmetic, a fresh evaluation of the audited integer core tensors, and numerical controls not importing producer code. The core tensors themselves are an explicitly reused, separately audited dependency; this audit does not claim a new derivation of those tensors.

## Analytic checks

Put B=S²+(n−1)I. Its entries are bounded and its diagonal is zero. Thus ||S||²=||−S²||=O(n), and all rational propagators are analytic on any fixed complex disk for all sufficiently large n, uniformly over the matrix class. Their resolvents have uniformly bounded operator norms. Enlarging the disk before applying Cauchy's estimates justifies every fixed Euler derivative of the errors. There is no appeal to differentiating an unproved pointwise asymptotic.

The diagonal fourth-power estimate is (S⁴)ii=(n−1)²+O(n), and the off-diagonal fourth-power estimate is O(n). The third power has all entries O(n). These estimates give the stated diagonal expansions through n⁻², off-diagonal G,H errors O(n⁻³), and O errors O(n⁻²). Using only an operator estimate for the cubic term would have been insufficient; the proof correctly uses the entrywise S³ bound.

The theta off-diagonal sum is fixed by the exact identity sum(B)=n(n−2). Figure-eight, diagonal theta, off-diagonal theta, terminal, and lollipop contributions independently reproduce

A₁(z)=−(z²+7z⁴/2)/n+(12z⁴−8z⁶/3)/n²+O(n⁻³).

The dumbbell contribution is O(n⁻³), not generally estimated sufficiently sharply by replacing its vectors with unrelated arbitrary vectors. The producer correctly separates their constant part.

The vacuum uses tr(S⁶)=−n⁴+O(n³), and all traces starting with S¹⁰ give O(n⁻⁴). This yields the stated expansion of X₀.

## Finite second-weight contractions

After summing terminal labels first, a weight-two core has E=v+2 scalar factors and v remaining indices. Error insertions cost a factor n⁻¹ under the producer assumptions. Nonloop identity factors force equality and also cost n⁻¹. Collision removal is valid with error O(n⁻³).

For an odd residual edge, restore unrestricted indices after recording that collision error. Conditional on the other indices, the contraction is fᵀSg, with bounded f,g. Thus it is O(n³ᐟ²), rather than O(n²), and contributes O(n⁻⁵ᐟ²). The restriction to distinct labels should not be left inside a purported separable bilinear form; the preceding collision removal is the justification. Empty residual signed graphs survive exactly when all d parities are even.

Independent exact evaluation of all 26 audited tensors gives precisely four nonzero polynomials: 12z⁶, 14z⁶, 3z⁴, and 4z⁴. Hence A₂=(7z⁴+26z⁶)/n²+O(n⁻⁵ᐟ²). Symbolic conversion independently gives κ₂/e=1−5/(2n)+65/(8n²)+O(n⁻⁵ᐟ²).

## Perturbations and triangle packing

For skew E with E1=0 and bounded row absolute sums, the same bound holds for columns. Consequently SE, ES, E² have uniformly bounded entries, d stays fixed, and Ed is uniformly bounded. The class is closed under the claimed perturbations.

A maximal vertex-disjoint triangle packing leaves a transitive tournament. Its row sums directly give ||U||²≥(m²−1)/3, whereas principal compression gives ||U||≤||S||. Therefore m=O(√n), and reversing independently any subset gives exponentially many distinct labeled matrices. No unlabeled/isomorphism-count conclusion follows.

Boundary correction requested: for a deleted DRT, write ||S||²≤n+1 rather than equality without qualification. At n=2, deleting one vertex from C₃ gives norm squared 1, not 3. The claimed leftover bound m≤√(3n+4) remains valid. Root applied this correction before publication; the corrected document hash is recorded in MANIFEST.json.

## Independent controls

`audit.py` constructs deleted Paley tournaments and makes TWO independently ordered rounds of maximal vertex-disjoint triangle reversals. The rounds may overlap in vertices. This tests bounded per-vertex perturbations beyond a single disjoint packing. Orders are 10,22,46,82,130,198; tested z values are 1, 0.7+0.4i, and 1.2i. Score vectors are exactly preserved, perturbation row absolute sums are at most 8, and the saved controls have max|Sd|≤9 and max|B|≤18.

The 18 finite scalar checks test the n³-scaled A₁ remainder. They are consistency controls, not proof of a uniform estimate. They all passed. Full receipts are in `receipt.json`; no large-order exact Hamiltonian or full A₂ enumeration is implied.

## Stronger energy-condition lemma

The proof actually permits replacing max|Sd|=O(1) by the weaker condition

**||Sd||₂²=O(n).**

Keep almostregularity and bounded off-diagonal S². Write w=Sd. Then

||G1||₂=O(||w||₂/n²)=O(n⁻³ᐟ²),
1ᵀG1=z²/n+O(||w||₂²/n⁴)=z²/n+O(n⁻³).

For D=c1+r, c=−2z²/n and ||r||₂=O(n⁻³ᐟ²), Cauchy–Schwarz gives |rᵀG1|=O(n⁻³). Thus the lollipop still has remainder O(n⁻³); the dumbbell remains O(n⁻³), with its cross term O(n⁻⁴). The path correction is O(||w||₂²/n⁴)=O(n⁻³). All vacuum and other A₁ estimates are unchanged.

For each terminal kernel, its row-sum remainder begins with a bounded rational function of Q applied to t²w. Therefore its Euclidean norm, and hence every entry, is O(n⁻³ᐟ²). This is weaker than the producer's O(n⁻²) terminal error, but still sufficient: replacing a leading n⁻¹ terminal factor by this error changes a weight-two contraction from O(n⁻²) to O(n⁻⁵ᐟ²). There are finitely many cores and factors. The residual-edge bound is already O(n⁻⁵ᐟ²), so the total A₂ error order is unchanged. The identical coefficient 65/8 follows under this energy hypothesis.

The energy condition is also preserved by the same bounded-per-vertex perturbations, since ||Ed||₂≤L√n. This lemma is a sufficient-condition extension, not a claim that every almostregular tournament satisfies the energy bound. Uniformity constants depend on the entrywise B bound, the energy constant, and the fixed z disk.
