# The cubic excellent-Frobenius congruence improves sharply to modulo27

Status: PROVED deduction from the explicitly inherited modular-polynomial framework. This resolves Open Problem2 of the inherited root note under that framework; it is not an independent reproof of its modularity and Fricke degree lemma. No external novelty claim is made.

## Statement and dependency boundary

Use the inherited notation F=G², G(t(q))=θ(q), s²=Q=(1+t)(1−27t), G(0)=s(0)=1, and the excellent lift t^σ(t(q))=t(q³). Then

F(t) ≡ (1+4t+48t²) F(t^σ) (mod27).

The modulus is sharp: the coefficient of t³ in F/F(t^σ)−(1+4t+48t²) is760−4=756=27·28, not divisible by81.

Inputs inherited from sources/inherited-root-note.md are the modular definition of T_p, its exact Fricke eigenspaces and pole-degree bounds (LemmaT.1), the eta/theta parameter identities, integrality of the formal series, and the exact root identity T_p(G/G^σ)=0. The mod-p² vanishing lemma, TheoremB, and the proof of TheoremE are not needed for the deduction below.

## Exact polynomial T3

At p=3, χ(3)=−1. The inherited degree bounds force:

T_3(x)=27+A x²+B s x³+(C+D t)x⁴.

In particular the x coefficient vanishes: its odd Fricke component would have polynomial degree at most−1. The x² coefficient is constant, the x³ coefficient is s times a constant, and the x⁴ coefficient has degree at most1. At q=0, θ=θ(q³)=s=1 and t=0, so

T_3(x)=(1−x)(3+x)³=27−18x²−8x³−x⁴.

Thus A=−18,B=−8,C=−1. To find D, let z=θ(q³)/θ(q), and e_j be the elementary symmetric functions in the three values θ(ζ^k q^(1/3))/θ(q), k=0,1,2. Newton identities give, through q¹,

z=1−2q+O(q²),
e1=3−6q+O(q²),
e2=3−36q+O(q²),
e3=1−22q+O(q²).

The leading x⁴ coefficient is −z e3=−1+24q+O(q²). Since t=q+O(q²), D=24. Therefore exactly

T_3(x)=27−18x²−8s x³+(24t−1)x⁴.

This is an identity of modular-function coefficients, established by the inherited finite degree spaces and the displayed exact coefficients, not by truncation alone. verify.py independently forms the Newton power sums from θ, constructs t and s from eta products, and checks all five polynomial coefficients through q^80.

## Solving its unit root modulo27

Let R=G/G(t^σ) and X=R/s; these are integral unit series with constant1. The exact root equation is

27−18R²−8sR³+(24t−1)R⁴=0.

Divide by the unit R² and reduce modulo27. Since R=sX, it becomes

Q[(1−24t)X²+8X]+18 ≡0 (mod27).                         (*)

Modulo3, Q is a unit and (*) becomes X(X−1)=0. Since X is a unit, X≡1 modulo3. Write X=1+3a. Modulo9, substitution in (*) gives

3(1+t)(a+t)≡0 (mod9),

hence a≡−t modulo3 and X≡1+6t modulo9. Therefore write X=1+6t+9b. Substitution in (*) now gives the exact reduced residual

9[(1+t)b−t−t³]≡0 (mod27).

Consequently b≡t(1+t²)/(1+t) modulo3, where1+t is an invertible formal series. The verifier independently expands (*) as a polynomial in b,t and checks this residual coefficient by coefficient.

Finally Q≡1+t modulo27 and

R²=QX²
 ≡(1+t)(1+12t+9t²+18b)
 ≡(1+t)(1+12t+9t²)+18t(1+t²)
 ≡1+4t+21t²
 ≡1+4t+48t² (mod27).

Since R²=F/F(t^σ), multiplication by the integral unit F(t^σ) proves the claimed congruence. The difference between the penultimate polynomial expansion and1+4t+21t² is exactly27t+27t³.

## Precision and scope controls

- The lift remains t^σ=t(q(t)³), not t³. It agrees with t³ modulo9 but cannot be silently replaced modulo27.
- This theorem does not prove C_3²G≡G+9θ_t²G modulo27, which still needs its own Cartier/lift-error argument.
- It also does not prove the OPEN higher-depth tower. Those packets are unchanged.
- No finite coefficient count is used as a substitute for an infinite modular identity: the exact T3 identity uses finite-dimensional coefficient bounds, and the remaining congruence is formal algebra.

## Search and verification record

A bounded search of the existing math repository and session50–53 publication trees found the original root note's Open Problem2 and finite observations modulo27, but no later proof of this particular improvement. This is a repository-status check, not a literature novelty certification.

The exact checks use only Python's standard library and survive python -O. They independently check the full T3 coefficient identities through q^80, the resulting excellent-Frobenius congruence through q^80, and the symbolic residual. replay.py checks recursive exact file membership and all hashes before reproducing both fixtures. hostile.py rejects a corrupted proof, a missing fixture and an undeclared nested file.
