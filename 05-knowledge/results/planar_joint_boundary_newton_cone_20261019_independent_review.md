# Independent review: joint-boundary Newton completion

2026-10-19. **ACCEPTED ANALYTICALLY**, including the convergence-conditional polynomial-limit corollary. The reviewed scope is the stated Newton cone and the specified weight topology; no general planar Jacobian assertion is inferred.

## 1. The completion is a well-defined complete topological algebra

The requirement that only finitely many supported monomials have weight 2j-i at most any fixed bound is stronger than a mere lower weight bound. It supplies both a minimum weight and finite truncations. Addition is immediate. For a product of elements with lower bounds a,b, every contributing pair of total weight at most M lies in the finite truncations at M-b and M-a. This proves that multiplication is defined, stays in J, and has the claimed filtration estimate.

A weight-Cauchy sequence stabilizes coefficientwise below every fixed cutoff; that stabilized truncation agrees with one actual element of the sequence and hence is finite. The resulting series is therefore in J and is the limit. Finite weight truncations show polynomial density. Convergent sequences have eventually uniform lower weight bounds; together with the multiplication estimate this gives joint multiplication continuity, not only continuity of multiplication by a fixed element.

For fixed t exponent j all terms have weight at most 2j, so their x polynomial is finite. The coefficient inclusion J into K[x][[t]] is injective. Weight at least N forces t exponent at least ceil(N/2), proving source continuity. Under boundary substitution x=-1/s,t=-s², weight is exactly the resulting s exponent. Local finiteness therefore proves convergence, finite negative support and multiplicativity of epsilon:J to K((s)).

The excluded series sum over n>=0 of x^(2n)t^n illustrates the necessary restriction: its source t-adic series exists, but all terms have weight zero and their boundary evaluations cannot be summed formally. Thus source convergence alone cannot replace convergence in J. The proof does not make that substitution.

## 2. Exponentials and boundary transport

The x derivative increases weight by one and the t derivative decreases it by two; both are continuous on J. Since the nonconstant Hamiltonian part has minimum weight at least two, each of the two Poisson-bracket terms increases input weight by at least one. This proves topological nilpotence of the derivation. The exponential, inverse and group law follow from the finite number of relevant operator powers below each cutoff. The same estimates justify multiplication and Poisson-bracket identities for the exponential, even though the filtration is not ideal-adic.

The independent source completion also supports the exponential because H-c is divisible by t². Continuity and matching polynomial partial sums identify the J and source operations; finite compositions inherit this compatibility. Thus source-polynomial coordinate images really are the same polynomial elements of J before applying epsilon.

The boundary velocity v_R=s⁴R(-1/s,-s²) has order at least two by the monomial condition. The displayed coordinate identities prove epsilon delta=(v_R d/ds) epsilon on polynomials. Polynomial density and continuity extend this to J. The boundary derivation increases s-order by at least one, so its exponential converges and intertwines term by term. No closed-kernel or q-adic shortcut is required.

## 3. Polynomial rigidity and closure corollary

A polynomial output has determinant one and fixes t=0 pointwise. The cited Gwoździewicz fixed-line theorem applies after algebraic closure. The automorphism preserves the line ideal, so T=at; determinant one and pointwise line fixing force T=t and X=x+f(t), f(0)=0. The retained boundary identity forces f(-s²)(2-s f(-s²))=0. The second factor has constant term two, hence f=0. The resulting polynomial identity descends to the original characteristic-zero field.

The closure corollary is also valid. Coordinate convergence in J passes determinants by continuous differentiation and multiplication, fixed-line restrictions by the continuous source embedding, and the q-boundary identity by continuous multiplication and epsilon. The resulting polynomial limit therefore meets exactly the same fixed-line theorem and shear argument. This deduction needs neither a limiting inverse nor convergence of maps on arbitrary inputs. It makes no claim that arbitrary infinite words converge or that weaker source convergence suffices.

## 4. Edge cases and checker inspection

The examples R=x and R=x² genuinely extend the prior K[xt,t] coefficient family. For R=x², v_R=s² and its formal boundary flow s/(1-lambda s) are correct. For R=x³, the Hamiltonian weight drops to one and the induced boundary velocity is -s; failure of weight increase marks failure of this particular convergence argument, not a counterexample to rigidity. Dropping q allows the stated genuine polynomial shear from H=t².

The standalone checker source was inspected. Its checks remain active under Python optimization, its finite partial-word computations are labeled as such, and the sparse-support tests corroborate multiplication, derivative truncation and evaluation interfaces. The checker does not purport to prove the infinite support or closure assertions; those follow from the analytical arguments above. The universal proof required no correction; an illustrative common-coordinate formula was corrected as recorded below.

## 5. Final proof and independent reproduction receipt

The accepted final proof SHA-256 is `1e35ad3181038db63e17992c149a877e2e0ce51112767b04117d4aa17ba01dcb`. The final changes add review and reproduction references and clarify the x³ edge: source t-adic exponentiation still exists, while the weight-completion and induced-boundary order-growth arguments fail. This clarification is correct and does not change the theorem. The standalone checker was independently run normally and with `python -O`; both completed successfully on 726 checks, and their output matched each other and the staged canonical JSON byte for byte. The reproduced JSON SHA-256 is `e560c82b13b2718f645eeaeacc536766e76d8d9f068edaf0a397b8d84d08aedd`. The final review is complete with no unresolved mathematical blocker.

## 6. Corrected illustrative common-coordinate formula

A denominator omission in the initial illustrative R=x² formula was identified after the first review and corrected before publication. Direct substitution of x=s/tau and t=tau gives H=s²(s²+tau)/tau, not s²(s²+tau). Consequently delta(s)=-s⁴/tau and delta(tau)=-4s³-2s tau. These corrected expressions were independently verified. The previous review had missed the omission. It affected only the example explaining why the older tau-order completion is insufficient; the universal weight-completion proof and its closure corollary are unchanged. Four dedicated symbolic checks were added, and the corrected checker was independently replayed in both modes with the 726-check count and final hashes recorded above. The corrected final text is accepted.
