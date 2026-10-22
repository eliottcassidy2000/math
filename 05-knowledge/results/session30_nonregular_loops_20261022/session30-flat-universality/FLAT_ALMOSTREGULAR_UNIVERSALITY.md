# Formal second-order universality for flat almostregular tournaments

## Scope and theorem

This is a theorem about the **finite-resummed first- and second-order predictors** defined by the complete audited linear-forest core expansion. It is not an asymptotic theorem for the exact Hamiltonian path count. In particular it does not supply the missing uniform remainder bound at z=1.

Let S=S_n be a skew tournament sign matrix of even order n, and put d=S1. Suppose, with constants independent of n,

1. d_i belongs to {−1,+1};
2. |(S²)_ij|≤C for i≠j;
3. |(Sd)_i|≤C.

Then the formal predictors obey

- κ₁/e = 1−5/(2n)+20/(3n²)+O(n⁻³);
- κ₂/e = 1−5/(2n)+65/(8n²)+O(n⁻⁵ᐟ²);
- y₁ = 3/2−13/(3n)+O(n⁻²);
- y₂ = 3/2−139/(24n)+O(n⁻³ᐟ²).

Here y_j=n(1−κ_j/ρ), ρ=det(nI−S)/|det(inI−S)|, and the predictors are exactly those in `session30-decorated-engine/predict.py`. No coefficient has been fitted. Constants implicit in all bounds depend only on C and the fixed compact z-domain, not on the tournament.

The often separately imposed condition ||S||op=O(√n) is redundant: B=S²+(n−1)I has zero diagonal and bounded off-diagonal entries, so ||S||op²=||−S²||op≤n−1+||B||op=O(n).

## Propagator estimates, with controlled errors

Write t=z/n, Q=tS, K=S², and B=K+(n−1)I. Every entry of B is bounded, ||B||op=O(n), and (B²)_ii=O(n). In particular

- (K²)_ii=(n−1)²+O(n);
- (K²)_ij=O(n) for i≠j;
- (S³)_ij=O(n), since S³=−(n−1)S+SB;
- ||K||op=O(n), and ||Q||op=O(n⁻¹ᐟ²).

These estimates justify the Neumann expansions uniformly on every fixed closed complex disk |z|≤R, once n is sufficiently large. Higher tails are bounded in operator norm. Entrywise estimates below are sharper where needed.

Define

G=−Q²(I+Q²)⁻¹, H=Q²(I−Q²)⁻¹, O=Q(I−Q²)⁻¹, D_i=(H−G)_ii.

The diagonal expansions are

G_ii=z²/n+(z⁴−z²)/n²+O(n⁻³),
H_ii=−z²/n+(z²+z⁴)/n²+O(n⁻³).

For example, G=−t²K+t⁴K²+O_op(n⁻³), and the displayed K² diagonal estimate proves the stated diagonal remainder. For off-diagonal entries,

G_ij=−t²B_ij+O(n⁻³), H_ij=t²B_ij+O(n⁻³), O_ij=tS_ij+O(n⁻²).

The G,H errors use (K²)_ij=O(n) and the operator tail O(n⁻³). For O, the cubic term has entries O(n⁻²), and its remaining tail starts with t⁵S⁵, whose operator norm is O(n⁻⁵ᐟ²), which is also O(n⁻²). Thus the required errors are uniform entrywise.

The source hypotheses give the additional bounds

G1=O_∞(n⁻²), 1ᵀG1=z²/n+O(n⁻³), D_i=−2z²/n+O(n⁻²).

Indeed K1=Sd has bounded entries. In G1=−t²Sd+t⁴K(I+t²K)⁻¹Sd, the second term has Euclidean norm O(n⁻⁵ᐟ²), so the claimed entry bound follows. Also 1ᵀK1=−||d||²=−n and 1ᵀK²1=||Sd||²=O(n); the scalar resolvent remainder is bounded using this last quadratic form, giving 1ᵀG1=z²/n+O(n⁻³).

## Complete A₁ through order n⁻²

Use the audited exact expressions

X₁=Σ_i[2G_ii H_ii−(3/2)G_ii²]
   +Σ_ij G_ij[H_ij²−O_ij²−(2/3)G_ij²]−DᵀGD,
T=−z²(d/n)ᵀ(I+Q²)⁻¹(d/n),
L=−2DᵀG1,
A₁=X₁+T+L.

The first, figure-eight sum is

−7z⁴/(2n)+(7z⁴−3z⁶)/n²+O(n⁻³).

The diagonal part of the theta sum contributes z⁶/(3n²)+O(n⁻³). On its off-diagonal part the cubic G,H terms total O(n⁻⁴), and

−G_ij O_ij² = t⁴B_ij+O(n⁻⁵).

This is where a possible nonuniversal invariant could have entered. It does not: almostregularity gives the exact identity

Σ_{i≠j}B_ij=1ᵀK1+n(n−1)=n(n−2).

Consequently the off-diagonal theta term is z⁴/n²+O(n⁻³).

The dumbbell is O(n⁻³). To see this without an overly weak operator estimate, write D=c1+r with c=−2z²/n and ||r||∞=O(n⁻²). The main quadratic term is c²1ᵀG1=O(n⁻³), the cross term is bounded with G1=O∞(n⁻²), and rᵀGr=O(n⁻⁴).

Hence

X₁=−7z⁴/(2n)+(8z⁴−8z⁶/3)/n²+O(n⁻³).

For T, the leading term is −z²/n. The correction is bounded by a constant times |t|²||Sd||²/n²=O(n⁻³), using skew-symmetry and the resolvent operator bound. Thus T=−z²/n+O(n⁻³).

Finally L=4z⁴/n²+O(n⁻³): use D=c1+r, the scalar sum 1ᵀG1 above, and Σ_i r_i(G1)_i=O(n⁻³).

We obtain the same full subleading coefficient as in deleted doubly regular tournaments:

A₁(z)=−(z²+7z⁴/2)/n+(12z⁴−8z⁶/3)/n²+O(n⁻³).

All these statements hold uniformly on a slightly larger disk than the one desired, so Cauchy's estimate also gives the same bounds after any fixed number of Euler derivatives z∂z.

## A₂ leading contraction remains universal

Use the audited 26 weight-two core tensors in `session30-decorated-engine/symbolic-kernels.json`, with four independent edge propagators R₊, R₋, Q₊₊, Q₊₋. Their required expansions are

R_±=n⁻¹(±zS−z²I)+O_entry(n⁻²),
Q₊₊=n⁻¹z²I+O_entry(n⁻²),
Q₊₋=−n⁻¹z²I+O_entry(n⁻²).

These follow from the K and S³ bounds above, because all four kernels are rational functions of Q analytic near zero. Terminal sums must be performed before size estimates. They satisfy

R_±1=±zd/n+O_∞(n⁻²), Q_st1=O_∞(n⁻²).

For the first identity, the remainder begins at t²Sd; terms of order three and higher are bounded using ||Sd||₂=O(√n), ||S||op=O(√n), and the resolvent bound. The same argument applies to Q_st, which starts at degree two. Column sums have the corresponding reversed leading sign.

After summing ell leaves, a weight-two core has v=V−ell free vertex indices and E=v+2 propagator factors. Replacing one factor by an entrywise O(n⁻²) error changes the contraction by O(n⁻³). A leading identity propagator between distinct core vertices forces their equality and also costs O(n⁻³); only identity loops can survive at order n⁻².

The remaining nonloop factors are S_ij/n, with bounded d vertex colors. Remove collisions of distinct free indices at cost O(n⁻³). Even edge multiplicities then disappear since S_ij²=1. If a simple signed edge (a,b) remains, fixing all other indices leaves a sum S_ab f(a)g(b) with bounded f,g. Its absolute value is at most ||S||op||f||₂||g||₂=O(n³ᐟ²), gaining n⁻¹ᐟ² relative to naive n² counting. This bounds that entire contraction by O(n⁻⁵ᐟ²). If no signed edge remains, any odd d vertex color annihilates the leading term since Σ_i d_i=0.

The same exact finite tensor calculation as in `session30-deleted-drt/derive_a2_leading.py` therefore applies without alteration. Its four surviving cores give:

- bouquet of three loops: 12z⁶;
- a loop at each of two vertices and two parallel joining edges: 14z⁶;
- four parallel edges: 3z⁴;
- one loop with two terminal leaves: 4z⁴.

Therefore A₂(z)=(7z⁴+26z⁶)/n²+O(n⁻⁵ᐟ²). The tensor enumeration is finite and exact; the class extension lies in the preceding remainder proof, not in a fitted or sampled extrapolation.

## Vacuum and falling-factorial conversion

The determinant vacuum X₀=logρ(z) has only traces of powers S^(4j+2). Its first term is exactly z²(1−1/n). The next term uses

tr S⁶=−n(n−1)³−3(n−1)tr B²+tr B³=−n⁴+O(n³).

Here trB²=O(n²) and |trB³|≤||B||op trB²=O(n³). Terms starting at power S¹⁰ total O(n⁻⁴). Hence

X₀(z)=z²(1−1/n)+z⁶/(3n²)+O(n⁻³),
ρ/e=1−1/n+5/(6n²)+O(n⁻³).

The audited falling-factorial conversion uses D=z∂z, c_j=D^jX₀ and exponential Bell polynomials b_j:

p₁=(b₂−b₁)/2,
p₂=(3b₄−2b₃−3b₂+2b₁)/24,
shift=[D²A₁+(2c₁−1)DA₁]/2,
Y₁=A₁+p₁/n,
Y₂=A₂+A₁²/2+(A₁p₁+shift)/n+p₂/n²,
κ₁=ρ(1+Y₁), κ₂=ρ(1+Y₁+Y₂).

At z=1, the now-proved universal expansions give

Y₁=−3/(2n)+13/(3n²)+O(n⁻³),
Y₂=35/(24n²)+O(n⁻⁵ᐟ²).

Multiplication by the displayed ρ/e proves the theorem.

## A genuinely larger, stable class

Let S satisfy the theorem, and let S'=S+E be another tournament sign matrix with E1=0 and every absolute row sum at most L, independent of n. Since E is skew, the same holds for columns and ||E||op≤L. Then d'=d. Each entry of SE and ES has absolute value ≤L; each entry of E² has absolute value ≤L². Therefore off-diagonal (S')² remains uniformly bounded. Moreover

S'd=Sd+Ed, ||Ed||∞≤L.

Thus S' satisfies the same theorem. This permits O(n) total reversed arcs, not just a fixed finite number. For example, reverse any vertex-disjoint collection of directed triangles. Each reversal preserves scores, and the perturbation has absolute row sums at most four regardless of how many triangles are reversed. More generally, a bounded per-vertex number of score-preserving reversals suffices.

Starting from a deleted doubly regular tournament gives a concrete infinite source of examples. The result does not assert that arbitrary almostregular tournaments satisfy these flatness assumptions; almostregularity alone does not control off-diagonal S² or Sd.

### There are linearly many available score-preserving flips

This is not merely a finite-sample phenomenon. Take any maximal vertex-disjoint packing of directed triangles in S, and let m be the number of unused vertices. The induced remainder has no directed triangle and hence is a transitive tournament. For the skew sign matrix U of a transitive tournament on m vertices, the row sums are m−1,m−3,…,1−m, so

||U||op² ≥ ||U1||₂²/m = (m²−1)/3.

Principal compression gives ||U||op≤||S||op=O(√n); therefore m=O(√n). The packing has (n−O(√n))/3 triangles. Reversing any subset preserves all scores and has perturbation row sums at most four. The resulting family has at least 2^((n−O(√n))/3) distinct **labeled** matrices, all with the same uniform formal constants. This is not a lower bound on their number of isomorphism classes.

For a deleted doubly regular tournament, ||S||op²≤n+1, so m≤√(3n+4). Reversing the entire packing flips at least n−√(3n+4) arcs. Some reversals may still produce deleted doubly regular tournaments; non-DRT membership must be checked separately and is not asserted for every subset.

## Checks and remaining obligation

`check_controls.py` constructs deleted Paley matrices and reverses a greedy vertex-disjoint family of directed triangles. It verifies score preservation, bounded perturbation degree, and reports the defect from the deleted-DRT matrix identity. Nonzero defect proves these controls are not deleted doubly regular tournaments. The n=6 control retains zero defect; all saved controls from n=10 onward have nonzero defect. At n≤10 it also compares the independent scalar A₁ expression against all five first-order kernel contractions and evaluates the complete second-order predictor. Larger-order controls check A₁ only. `verify_formal.py` independently checks the exact symbolic A₁ and falling-factorial coefficient arithmetic. No numerical check is used to justify a uniform theorem.

The remaining exact-H obligation is unchanged: bound the aggregate of all omitted cores and falling-factorial terms at z=1 uniformly in this class. Neither these entrywise contraction bounds for finitely many cores nor the finite controls establish that infinite remainder bound.
