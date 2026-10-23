# Session37 independent arithmetic audit

Verdict: **ACCEPT**, subject to keeping the explicit domain and information-loss fences below. The initial draft and the later regular-two-sink/minimal-order additions have been independently proved and tested. No substantive mathematical error was found. The classical transcendence results are cited inputs, not newly proved or computationally certified. No Git state was changed.

Audited producer deliverable: `session37-arithmetic/ARITHMETIC_CUT_COMPLETION.md` and its exact-control script. Audit artifacts are outside the repository in `session37-arithmetic-audit/`.

## 1. Core field construction

For an integral skew tournament matrix S, K=-S²=SᵀS is an integral symmetric positive-semidefinite matrix. If μ is a nonzero eigenvalue, its monic irreducible polynomial over Q divides the characteristic polynomial of K. Every conjugate is consequently a nonnegative real eigenvalue. Zero cannot be a conjugate of nonzero μ, since an irreducible polynomial with root zero is x. Thus μ is a totally positive algebraic integer. Every embedding of F=Q(μ) is real. At every embedding, -μ is strictly negative, so adjoining sqrt(-μ) gives a quadratic extension E/F with no real embeddings. This is exactly a CM field.

The exclusion μ=0 is essential: Q(0,sqrt(0))=Q is not a CM field. Odd-order tournaments have one zero mode; for n=1 this is the only mode.

For s=i sqrt(μ), with positive real sqrt(μ), put t=(1+s)/(1-s). Then

- t+t⁻¹=2(1-μ)/(1+μ),
- s=(t-1)/(t+1),
- c=13-6(t+t⁻¹)=(1+25μ)/(1+μ),
- μ=(c-1)/(25-c).

The denominators are nonzero, t lies in the open upper semicircle, and 1<c<25. These identities prove Q(t)=E and Q(c)=Q(μ)=F, not merely a coarse arithmetic-type correspondence. Complex conjugation interchanges t and t⁻¹ and has fixed field F. Every conjugate of t has modulus one. **This alone does not imply torsion**, since t need not be an algebraic integer. Rational c need not be an integer, e.g. DRT19 gives c=119/5, although rational μ must be integral.

## 2. Full matrices versus spectral quotients

The rational orthogonal Cayley matrix U=(I+S)(I-S)⁻¹ exists, since a real skew matrix has no eigenvalue 1. Its inverse is its transpose. Direct functional calculus gives

C=13I-6(U+U⁻¹)=(I+25K)(I+K)⁻¹.

Two inverse maps are important:

S=(U-I)(U+I)⁻¹,
K=(C-I)(25I-C)⁻¹.

Hence full labelled U retains the full labelled tournament S. Full C retains the full Gram matrix K. K forgets at least the distinction S versus -S, whereas taking a spectrum also discards eigenvectors and labels. Statements of Hamiltonian-path blindness belong to the spectral multiset, determinant functions and scalar log/cut observers. The audit does not assert that full U loses this information, or that the T3/C3 collision has identical full K matrices.

The independent controls check both inverse maps on every labelled tournament through n=5.

## 3. Torsion, normalized angles and the first non-torsion order

Fix θ=2 arctan sqrt(μ) in (0,π), so t=e^(iθ). If θ/π is rational, t is a root of unity. If t has infinite order, θ/π cannot be rational. If it were algebraic irrational, the value exp(iπ·θ/π)=t of (-1)^(θ/π) would be transcendental by Gelfond–Schneider, a contradiction. Thus the normalized angle is rational or transcendental, never algebraic irrational. The branch log(-1)=iπ is explicit. Every other argument branch changes θ/π by an even integer, and reciprocal branches negate it, so the distinction is branch-independent.

For rational positive μ, root-of-unity t makes t+t⁻¹ a rational algebraic integer in (-2,2), hence one of -1,0,1. Solving gives μ=3,1,1/3; conversely these do give roots of unity. For rational tournament energies only μ=1,3 survive integrality. This classification must not be extended to irrational μ. The transitive T4 supplies the hostile: μ=3±2sqrt(2), but t is a primitive eighth root of unity and E=Q(ζ8).

The order-five minimum is valid by proof, not merely census:

- n=1 has no nonzero mode.
- n=2 has energy 1, with phase order 4.
- n=3 has energy 3, with phase order 3.
- At n=4, the skew characteristic polynomial is x⁴+6x²+p², where p=s12s34-s13s24+s14s23 is an odd sum of three signs. Thus p is ±1 or ±3. The energies are 3±2sqrt(2) or the repeated value 3; the phases have order 8 or 3.
- At n=5, take C3 followed by two ordered sinks: C3→a→b, with every old vertex beating both new vertices. The upper signs are all +1 except s02=-1, mask 1021 in lexicographic upper-edge enumeration. Its characteristic polynomial is x(x²+3)(x²+7), so K has energies 0,3,3,7,7. The μ=7 phase has infinite order and transcendental θ/π.

The finite census finds 640 labelled order-five tournaments with this 3/7 spectrum and 384 with energies 5±2sqrt(5), each doubled, plus zero. The latter phases are fifth roots of unity. These counts are finite-exact; they are not needed for the minimal-order proof.

## 4. Regular two-sink spectral insertion

The stronger operation also passes. Let S be an odd-order regular tournament sign matrix, S1=0, and form

S' = [[S,1,1],[-1ᵀ,0,1],[-1ᵀ,-1,0]].

The subspace W=1⊥ of the old coordinates is invariant and the old action survives unchanged. Since S is skew and S1=0, no old/sink mixing occurs on W. The complementary invariant space has basis (1_n,e_a,e_b), in which S' acts by

[[0,1,1],[-n,0,1],[-n,-1,0]].

This block has characteristic polynomial x(x²+2n+1). The orthonormal version replaces n and 1 by sqrt(n) in the corresponding off-diagonal entries, as in the producer's statement. Since the old regular S has exactly one zero mode,

χ_(S')(x)=χ_S(x)(x²+2n+1).

The inserted energy is 2n+1, with multiplicity two; its field is Q(sqrt(-(2n+1))), its charge is 25-12/(n+1), and for n≥3 its normalized angle is transcendental. n=1 is a genuine exception: the output is transitive T3 and the new energy 3 is torsion.

Regularity is necessary. Appending the same two ordered sinks to nonregular transitive T3 gives transitive T5, whose energies 5±2sqrt(5) do not include 7. The output of the regular-base operation is not generally regular, so the statement does not license automatic iteration. Independent exact controls verify the characteristic-polynomial insertion for regular cyclic bases n=1,3,5,7,9,11, including non-DRT bases.

## 5. Doubly regular specialization and parity

For an actual tournament with K=nI-J, the nonzero spectrum is n with multiplicity n-1, and the zero-mode multiplicity is one. The associated field is Q(sqrt(-n)) and charge is 25-24/(n+1). At n=3 the phase has order 3; every admissible DRT n≥7 gives non-torsion phases and transcendental normalized angles. There is no existence claim at every odd n. At n=1 there are no nonzero modes, so substituting μ=n into a mode statement is not licensed.

The rank parity assertion is correct: every even principal sign-skew matrix reduces modulo 2 to J-I. In even dimension this matrix is invertible over F2. Therefore even-order S is invertible; odd-order S has a nonsingular (n-1)-principal minor and skewness forces its determinant to vanish. Thus nullity(S)=nullity(K)=n mod 2. This is linear-algebraic parity, not a claim about physical spin.

## 6. Formal occupancy functions and logarithms

The identities

F(z)=Σ tr(Λ^k K)z^k=det(I+zK),
B(z)=Σ tr(Sym^k K)z^k=det(I-zK)⁻¹

are valid exterior/symmetric power identities. F is a polynomial; B is initially a formal series and analytically converges for |z|<1/ρ(K). Its coefficients and those of F are nonnegative integers. The weights are μ, not thermal Boltzmann factors e^(-βμ). The terminology describes formal one-particle occupancy; no quantum-statistical physical model has been established.

For n≥2, ρ(K)≤(n-1)². Thus rational 0<z<(n-1)⁻² guarantees rational F>1 and B>1, with log B>log F>0. All these positive nonzero logarithms are transcendental by Hermite–Lindemann: otherwise exponentiating a nonzero algebraic logarithm would contradict its rational value. The same proof works with positive algebraic z in the convergence domain.

The positive even-sector difference must be log B-log F, not the reverse. Precisely,

2O=log(FB),
2E=log(B/F)=-log det(I-z²K²).

The arguments FB and B/F are rational greater than one. Multiplying a nonzero transcendental logarithm by the nonzero rational 1/2 preserves transcendence. This does not imply algebraic independence among O,E,log F,log B.

The fermionic logarithm itself remains positive and transcendental for every positive rational z, without the Bose radius restriction. Its Taylor series about zero, and the common even/odd trace split, do require their stated convergence domain. Keeping the shared safe radius is sufficient, not sharp.

The Bose radius cannot be dropped just because its rational formula exists. For the n=2 tournament K=I, at z=3 the algebraic expression B=(1-3)^(-2)=1/4 is positive but log B is negative and the occupancy series diverges. At n=1 or z=0 all stated logarithms vanish; the transcendence conclusions do not apply.

## 7. Rational cuts and completion

For rational A>1, r=(A-1)/(A+1) lies in (0,1). The rational partial sums

L_N=2Σ_(j=0)^(N-1) r^(2j+1)/(2j+1),
U_N=L_N+2r^(2N+1)/[(2N+1)(1-r²)]

satisfy L_N<L_(N+1)<log A<U_(N+1)<U_N and their width tends to zero. The lower bound is the positive atanh series truncation; the upper bound follows by replacing every tail denominator by 2N+1, with strictness because later denominators are larger. Direct subtraction proves U_N decreases. Thus {q in Q: q<L_N for some N} is exactly the lower Dedekind cut of log A.

These are error-certified rational enclosures, not a finite transcendence certificate. The external transcendence theorem does essential work. Likewise a transcendental cut has no supremum in the real algebraic numbers, establishing their failure of Dedekind completeness. Every order-field automorphism of R fixes all rational cuts and is the identity; the logarithm/Cayley constructions are not such automorphisms.

For transitive T4 at z=1/18, exact controls recover

F=(433/324)²,
B=(324/217)²,
O=log(433/217),
E=log(104976/93961).

Both last arguments exceed one. No decimal approximation is used to support their arithmetic type.

## 8. Recovered prior statements and novelty fence

The audit directly read the old theorem files in `math/01-canon/theorems/`:

- THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density.md proves that scalar comparison x_i<x_j loses external arithmetic labels and limit information, including finite-prefix continuation hostiles. Its precise scope is the scalar-order observer, not every intrinsically arithmetic tournament construction.
- THM-4095-exact-arithmetic-field-transport-gapped-pair-margins-and-order-tournament-blindness.md proves Q(F_S(t)-β)=Q(t) for the piecewise-affine LRC observer and separates fixed-observer density from optimized gapped spectra. It is a model of field-preserving maps, not a proof of this spectral construction.
- The Session35 elliptic charge note already used the nonintegral rational trace obstruction and the Gelfond–Schneider normalized-angle argument. These are inherited mechanisms, not new transcendence theorems.
- LEM-004 provides the skew orientation/formal arctanh dictionary. Its abbreviation CM in sign-function tables means completely multiplicative, unrelated to CM fields. Its quasirandom/DRT limit qualifications are not imported into the present all-tournament finite matrix identities.

The lawful new-to-session synthesis is the explicit integral-tournament → totally positive energy → CM pair → field-preserving charge → rational occupancy values → transcendental cut chain, its sharp n=5 non-torsion threshold, and the regular two-sink spectral insertion. Literature-wide novelty is not established.

The transitive T3 and directed cycle C3 both have χ_S=x(x²+3), K spectrum {0,3,3}, and identical F,B and all spectral scalar completions, while their directed Hamiltonian-path counts are 1 and 3. This proves the claimed spectral blindness. It does not conflict with THM-4088/4095, and it does not transfer any unresolved Hamiltonian, LRC, zeta, VOA or physical statement.

## 9. Classical source check

Author-hosted sources were freshly opened and the relevant statements checked:

- J. S. Milne, Algebraic Number Theory, version 3.08, section on CM fields: https://www.jmilne.org/math/CourseNotes/ANT.pdf . The definition matches the field construction above.
- Richard Elman, Lectures on Abstract Algebra, Corollary 73.8 and Theorem 74.29/Corollary 74.30: https://www.math.ucla.edu/~rse/algebra_book.pdf . These respectively supply the nonzero algebraic exponential obstruction and the algebraic-base/irrational-algebraic-exponent obstruction used above.

No deeper claim is attributed to those sources. The displayed all-order matrix, field, extension and cut proofs are independent elementary arguments using these cited inputs only at their explicitly marked transcendence steps.

## 10. Independent finite verification

`independent_exact_checks.py` does not import producer code. Matrix multiplication, determinant, inverse and Faddeev–LeVerrier characteristic coefficients use Python integers/Fractions. SymPy is used only for exact polynomial factoring, rational root isolation and polynomial equality in the extension controls. Explicit exception checks remain active under `python -O`.

The run enumerates all 1,099 labelled tournaments n=1 through 5, tests all full-matrix Cayley and charge forward/inverse maps, exact partition inequalities, nested rational cuts, even-product determinant identity, skew rank parity, and positivity of all conjugates in the finite spectrum. Additional controls cover transitive Cayley powers through n=12, irrational-energy torsion, rational torsion candidates, Paley3/7/11/19, the T3/C3 Hamiltonian collision, n=1 and outside-radius hostiles, and the regular-two-sink insertion.

Normal and optimized output must match byte for byte. The JSON result and SHA256 manifest document the final executed revision. None of these finite checks purports to test transcendence numerically or replace an all-order proof.

## 11. All-eight extension refinement and unbounded spectral ambiguity

The later refinement is **ACCEPTED**. Let the three quotient blocks be T,a,b, and choose any of their eight tournament orientations. All edges between a singleton and T have the selected uniform direction. For a regular n-vertex base, W=1⊥ still survives. With signs α,β,γ in {±1}, the residual matrix in (1_n,e_a,e_b) is

[[0,α,β],[-nα,0,γ],[-nβ,-γ,0]],

whose characteristic polynomial is x(x²+2n+1), independent of the signs. Thus all eight extensions have the same entire S spectrum and K spectrum, not just the inserted energy.

Define pc_T(2) as unordered covers of all old vertices by two nonempty directed paths, including singleton paths. Equivalently, writing H(U)=1 for a singleton induced tournament,

pc_T(2) = (1/2) Σ_(∅≠A⊊V(T)) H(T[A]) H(T[V(T)\A]).

This definition counts actual directed paths, not merely partitions of the vertices. Then the six transitive quotient orientations have H_ext=H(T). A path cannot return to an earlier quotient block, so the T vertices form one run in the unique quotient order.

For either cyclic quotient orientation, e.g. T→a→b→T, a spanning path has either one or two T-runs. One run can occur in the three orders T,a,b; a,b,T; b,T,a, contributing 3H(T). Two runs must have the form P,a,b,Q, contributing one for each ordered two-path cover, or 2pc_T(2). Therefore

H_cyclic=3H(T)+2pc_T(2).

The path formulas themselves need no regularity; regularity is needed for the common spectral polynomial. For the directed C3 base, H=3 and pc2=3, so the two kinds of order-five extensions have Hamiltonian counts 3 and 15 despite identical full spectral data, including the added μ=7 and c=22.

An elementary strengthening follows from the same counting. Cut every Hamiltonian path of T at one of its n-1 gaps. Each unordered two-path cover permits at most two concatenations, so

(n-1)H(T) ≤ 2pc_T(2).

Consequently

H_cyclic ≥ (n+2)H_transitive.

For every odd n, a regular cyclic tournament exists, with explicit forward Hamiltonian path 0,1,...,n-1. The ambiguity ratio among these equal-spectrum pairs is therefore unbounded as their order N=n+2 increases; it is at least N. In particular, no fixed universal multiplicative approximation factor from the spectral observer alone is possible. For a symmetric multiplicative estimate satisfying H/C≤estimate≤CH on both members, necessarily C≥sqrt(N). This conclusion concerns unrestricted tournaments; the constructed extensions are not generally regular and cannot be represented as a regular-only impossibility theorem.

Every cyclic extension is strongly connected for nonempty T: any old vertex reaches any other via a,b, and the singleton blocks communicate with the old block around the quotient cycle. Every transitive-quotient extension is reducible. The pair therefore does **not** establish an unbounded ambiguity ratio inside a fixed strong class or a fixed score-sequence class. All eight share skew-S/energy-K/phase/charge spectra and determinant/log observables, not full Cayley matrices, full charge matrices, or necessarily score sequences. They are not asserted to have the same 0/1 adjacency-A spectrum: A=(J-I+S)/2, and the extended tournaments are not regular, so the skew spectrum does not by itself determine the adjacency spectrum. Thus the ratio obstruction does not contradict a result using the separate adjacency characteristic polynomial and permanent.

The independent script checks all eight orientations for regular cyclic bases n=1,3,5,7, using a separately written subset-DP path counter and permutation cross-checks at small orders. It verifies the common characteristic polynomial, exact path formula, strong/reducible distinction and cut-count ratio bound in every case. The general proof above, not these bounded controls, supplies the unbounded ratio.

## 12. Final source-revision checks and freeze

The final producer note's P163 arithmetic specialization also passes an independent check. Trial division through 12 verifies primality; all 162 nonzero quadratic-character autocorrelations equal -1. This gives K=163I-J, t=(-81+i sqrt(163))/82 and c=1019/41 exactly. The receipt is `paley163_final_check.json`. No class-number-one, Heegner, modular, VOA or RH assertion is needed for that calculation.

The raw angle θ is already transcendental for every nonzero algebraic mode: if θ were algebraic, iθ would be nonzero algebraic while e^(iθ)=t is algebraic, contradicting Hermite–Lindemann. The sharp five-vertex threshold concerns **normalized θ/π**, not θ itself. This distinction is correctly present in the final producer note.

Final arithmetic additions are frozen after the all-eight/strongness/ratio checks. The strong/reducible and score-class caveats were passed to the producer and are present in its final text. Further theorem additions require a new audit rather than silently inheriting this acceptance.

Final prose correction: the displayed even logarithm is written as log(B/F)=log(det(I-z²K²)^(-1)) when invoking arguments greater than one. The equivalent form -log det(I-z²K²) is also correct, but its un-inverted determinant argument is less than one. This is a wording correction only, not a failed substantive claim.
