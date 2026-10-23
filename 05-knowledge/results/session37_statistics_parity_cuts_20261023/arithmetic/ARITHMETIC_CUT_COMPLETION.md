# Tournament spectral fields, Cayley charges, and exact logarithmic cuts

Status: PROVED elementary matrix/field/series results; the explicitly identified Hermite–Lindemann and Gelfond–Schneider inputs are CITED. FINITE-EXACT controls are recorded separately. These are operational connections and new-to-session syntheses, not claims of new classical transcendence theorems or literature-wide novelty.

## 1. Inheritance, portfolio, and precision

The inherited source is `session35-publication/clean-repo`, read at commit `989b11f2ad12`, with a clean working tree and cached `origin/main` equal to HEAD. No repository or Git state was changed. The read-only startup packet ran. Targeted message inspection was blocked by absent `.machine-id`; no registration or state mutation was attempted.

Closest proved mechanism: session35's exact reciprocal quotient `C(t)=13−6(t+t⁻¹)`, together with the actual elliptic charge orbit and its CM field. Canonical hostile: THM-4088 already proves that bare order tournaments cannot recognize rational, algebraic-irrational, or transcendental limits. Corrected near miss: real and transcendental are not opposite predicates, and invariant measure does not prove equidistribution of a chosen doubling orbit. Least-used sidecars: the exact generated-field equality in THM-4095, LEM-004's odd formal logarithm, and the retained order/analytic completion rather than another orientation statistic.

Anchor: an exact tournament-to-arithmetic map. Niche: odd/even logarithms of exterior/symmetric power traces. Wildcard: the inherited Coulomb coordinate applied to genuine tournament spectral modes.

Concept board: (1) skew orientation and parity; (2) positive spectral energy; (3) CM conjugate pair and real quotient; (4) finite rational partition data; (5) analytic completion and Dedekind cuts; (6) spectral information loss. The source object and every map below are explicit. No intrinsic tournament is manufactured from tied scalar values.

## 2. The rational operator and its parity

Let T be a labelled tournament on n vertices, with sign matrix S having S_ij=1 when i→j, S_ji=−S_ij, and zero diagonal. Define

    K = −S² = SᵀS.

K is a positive semidefinite symmetric integer matrix. Its trace is n(n−1). The eigenvalues of S are conjugate pairs ±iλ_j with λ_j>0, plus zeros; each nonzero energy μ_j=λ_j² occurs twice in K, counting multiplicity.

For a tournament, nullity(S)=nullity(K)=n mod 2. Indeed an even-order sign-skew principal matrix reduces modulo 2 to J−I. For even m this matrix is invertible over F₂: if (J−I)v=0, all coordinates equal their sum, and summing gives their sum=m times itself=0. Thus every even-order principal determinant is odd and nonzero. For even n this proves full rank; for odd n a principal (n−1)-minor is nonsingular while skewness forces determinant zero. The n=1 case is K=0. This is a genuine parity effect, not a physical spin assertion.

## 3. A CM field attached to every nonzero mode

Theorem A. If μ>0 is an eigenvalue of K, then μ is a totally positive algebraic integer. The field F=Q(μ) is totally real, and

    E = F(sqrt(−μ))

is a CM field, namely a totally imaginary quadratic extension of F.

Proof. The minimal polynomial of μ divides the integer characteristic polynomial of K, whose roots are all real and nonnegative by symmetry and positivity. Therefore every conjugate of μ is nonnegative. None is zero: otherwise the irreducible minimal polynomial of the nonzero μ would have zero as a root. Hence every conjugate is strictly positive. At every real embedding of F, −μ is negative and cannot be a square. Adjoining its square root has degree two, and every embedding of the resulting field is nonreal. This proves the claim from the definition of a CM field, without a classification theorem. QED.

Choose λ=√μ>0 and set

    t = (1+iλ)/(1−iλ),
    c = C(t) = 13−6(t+t⁻¹) = (1+25μ)/(1+μ).

Then |t|=1, t≠±1, and 1<c<25. More precisely,

    t+t⁻¹ = 2(1−μ)/(1+μ),
    μ = (c−1)/(25−c),
    sqrt(−μ) = (t−1)/(t+1).

Consequently

    Q(c)=Q(μ)=F,        Q(t)=E,
    complex conjugation t↔t⁻¹ fixes exactly F in E.

This is the precise real/complex conjugate-pair duality. Arithmetic type is a separate coordinate: c is rational iff μ is rational; a rational eigenvalue μ of the integer matrix K is necessarily an integer. Irrational spectral charges are algebraic irrational, never transcendental.

There is a matrix-level map, not only an eigenvalue construction:

    U=(I+S)(I−S)⁻¹ ∈ O_n(Q),
    C_U=13I−6(U+U⁻¹)=(I+25K)(I+K)⁻¹.

The inverses exist because S has imaginary spectrum and K is positive semidefinite. Transposing proves Uᵀ=U⁻¹. The displayed charge identity follows by adding the two Cayley fractions. C_U is rational symmetric, with the above c-values as its spectrum. The full matrices have exact inverse reconstructions:

    S=(U−I)(U+I)⁻¹,
    K=(C_U−I)(25I−C_U)⁻¹.

Thus full U retains the labelled tournament, and full C_U retains K. Only after passing to spectra or scalar determinants does the later spectral-blindness assertion apply. When n is odd, the unique zero mode gives t=1 and c=1. Every nonzero mode gives the open elliptic interval (1,25).

Source→target: tournament S → K → individual field-decorated mode μ → reciprocal pair {t,t⁻¹} → real charge c. The last step preserves Q(μ) exactly but forgets the conjugate sheet; λ>0 selects one sheet. K and its spectrum also lose much more of the original tournament, quantified below. These are not positive-real mirror pairs: the session35 formulas X(t),Y(t) on t>1 do not apply to this unit-circle branch.

## 4. Rational or transcendental normalized angle, with an exact torsion boundary

Write t=e^(iθ), where θ=2 arctan√μ lies in (0,π). If t is a root of unity, θ/π is rational. Otherwise θ/π is transcendental.

Proof of the second assertion. An algebraic irrational β=θ/π would make the value exp(iπβ)=t of (−1)^β transcendental by Gelfond–Schneider, contradicting t algebraic. Infinite order rules out rational β. There is no algebraic-irrational middle possibility. This does not infer transcendence by multiplying two transcendental numbers. In fact the raw angle θ is transcendental for every nonzero mode by Hermite–Lindemann, since e^(iθ)=t is algebraic and θ≠0. The sharp torsion and five-vertex distinctions concern the normalized angle θ/π, not θ itself.

For rational μ>0, t has finite order iff

    μ ∈ {1/3,1,3}.

Indeed a root of unity makes t+t⁻¹ a rational algebraic integer, hence an integer in (−2,2), necessarily −1,0,1. Solving 2(1−μ)/(1+μ) in that set gives precisely the three values, and the corresponding t are explicit roots of unity. Since a rational tournament energy is an integer, its finite-order exceptions are exactly μ=1 and μ=3. Irrational μ can also produce roots of unity; it would be false to extend this rational-μ criterion indiscriminately.

Doubly regular specialization. For any tournament satisfying

    S1=0,        −S²=nI−J,

all nonzero energies are n, with multiplicity n−1. Standard doubly regular tournaments have this identity; the elementary results here require the displayed identity rather than an existence claim for every n. For n≥3 its modes have

    E=Q(sqrt(−n)),       c_n=25−24/(n+1).

The n=3 mode is torsion; for every admissible n≥7 the angle divided by π is transcendental, since its rational energy n is not 1 or 3. For Paley7 specifically,

    t=(-3+i√7)/4,        c_7=22,        E=Q(√−7).

A compact D=−163 recovery is equally explicit. The Paley tournament P₁₆₃ is valid because 163 is prime and 163≡3 mod4; its quadratic-character correlation verifies K=163I−J. The present map gives

    μ=163, E=Q(√−163), t=(−81+i√163)/82, c=1019/41.

The exact companion checks all 162 nonzero correlation shifts, without enumerating tournaments of that size. The old D=−163 motif is thus reached through a specified spectral field. Class-number-one/Heegner status is not used, and confers no extra modular, VOA or RH consequence here; the same construction applies to every admissible doubly regular order.

Thus an exact rational charge can carry a provably transcendental normalized angle. This extends the inherited session35 charge-38/5 example through a natural tournament operator, not through numerical coincidence.

Hostile and complementary control: for the transitive tournament on four vertices,

    spec(K)={3+2√2,3+2√2,3−2√2,3−2√2}.

Its positive-energy charges are 13±6√2, and its CM field is Q(√2,i)=Q(ζ₈). Here the parameters t are eighth roots of unity, so their normalized angles are rational despite irrational charges. This supplies both directions of the rational-charge/angle non-equivalence.

## 4A. Two ordered sinks force a new transcendental normalized angle

Theorem B. Let T be any regular tournament of odd order n≥3. Form T⁺ by appending two ordered sinks a,b: every old vertex points to both, and a→b. Then all old nonconstant skew modes survive, and the new spectral block contributes K-energies

    0, 2n+1, 2n+1.

The new positive mode has the imaginary quadratic CM field and rational charge

    Q(√−(2n+1)),        c=25−12/(n+1),

and its angle divided by π is transcendental.

Proof. Regularity means S1=0, so W=1⊥ is invariant and the constant cross-edges couple trivially to W. In the orthonormal basis (1/√n,a,b) of the remaining three-dimensional subspace, the skew matrix is

    [[0, √n, √n], [−√n, 0, 1], [−√n, −1, 0]].

Its characteristic polynomial is x(x²+2n+1), either by expansion or because the sum of the three squared upper entries is 2n+1. Hence the positive energy is 2n+1. The field and charge follow from Theorem A. Since 2n+1≥7 is a rational integer distinct from 1 and 3, the torsion criterion above and Gelfond–Schneider prove the angle assertion. QED.

This is an explicit orientation-preserving extension operation that creates a new arithmetic mode. It does not claim to recover the whole extended tournament from the old spectrum; regularity and the chosen constant-incidence interface are load-bearing sidecars. The output need not be regular, so this is not an automatically iterable closure operation.

## 4B. Sharp first occurrence at five vertices

Theorem C. Five is the smallest tournament order for which a nonzero Cayley spectral angle divided by π can be transcendental.

For n=1 only the zero mode occurs. For n=2 the positive energy is 1 and the Cayley modes are fourth roots. For n=3 the characteristic polynomial of S is x(x²+3), so the nontrivial modes are third roots. For n=4,

    det(xI−S)=x⁴+6x²+pf(S)²,
    pf(S)=S₁₂S₃₄−S₁₃S₂₄+S₁₄S₂₃ ∈ {−3,−1,1,3}.

Thus the energies are either 3 repeated four times, or 3±2√2 each twice. Their Cayley modes are respectively third roots or eighth roots of unity. These analytic formulas prove the lower bound independently of the enumeration.

At n=5 use T⁺ from Theorem B with T=C₃. In zero-based labels its sign matrix is

    [[0,1,−1,1,1], [−1,0,1,1,1], [1,−1,0,1,1],
     [−1,−1,−1,0,1], [−1,−1,−1,−1,0]].

It has K-characteristic polynomial x(x−3)²(x−7)². The μ=7 mode has charge 22 and transcendental normalized angle; tr(U)=−3/2. In the supplied script's lexicographic upper-edge bit convention this is mask 1021. The exact all-labelled n=5 census is:

- 640 tournaments: K-polynomial x(x−3)²(x−7)², with an infinite-order Cayley mode.
- 384 tournaments: K-polynomial x(x²−10x+5)², with energies 5±2√5; their Cayley modes are fifth roots of unity.

These counts are FINITE-EXACT, not a count of isomorphism classes. The minimal-order claim has the displayed analytic proof and explicit witness. In particular rational energy 7 can produce a transcendental normalized angle, while irrational energy 3±2√2 produces a rational normalized angle: the naive arithmetic-type↔torsion pairing fails in both directions.

## 5. Exterior/symmetric powers: the lawful fermionic/bosonic terminology

On V=Qⁿ with the endomorphism K, define the formal fugacity partition functions

    F_T(z)=Σ_(k=0)^n tr(Λ^k K) z^k = det(I+zK),
    B_T(z)=Σ_(k≥0) tr(Sym^k K) z^k = det(I−zK)⁻¹.

These identities follow by diagonalizing K and respectively choosing each eigenvalue at most once or with unrestricted repetition. The exterior/symmetric occupancy interpretation is the exact algebraic reason to use the terms fermionic/bosonic. This is a formal one-particle weight model, not the thermal formula det(I+exp(−βH)), and not a physical assertion about a tournament.

All coefficients of F_T and B_T are nonnegative integers. For F_T, this follows from elementary symmetric functions of nonnegative eigenvalues and integrality of the determinant. For B_T, the eigenvalue expansion gives nonnegative real coefficients and the inverse of its integer polynomial denominator gives integer coefficients recursively.

The exact reciprocal law is

    F_T(z) B_T(−z)=1.

For n≥2 put R=(n−1)². The norm estimate ||S||₂≤sqrt(||S||₁||S||∞)=n−1 gives 0≤K≤RI. At any rational fugacity

    0<z<1/R,

both F_T(z) and B_T(z) are rational and greater than one. The logarithms

    L_F(z)=log F_T(z),         L_B(z)=log B_T(z)

are positive real transcendental numbers: a nonzero algebraic logarithm would have transcendental exponential by Hermite–Lindemann, whereas its exponential here is rational. The same conclusion holds for all positive algebraic fugacities in this domain, replacing rational partition values by algebraic ones. Rationality is retained below because it makes every cut endpoint exact rational.

## 6. Even and odd trace sectors survive analytic completion

Set a_m=tr(K^m), positive integers for every m≥1 when n≥2. Absolute convergence in |z|<1/R gives

    L_F(z)=Σ_(m≥1) (−1)^(m+1) a_m z^m/m,
    L_B(z)=Σ_(m≥1) a_m z^m/m,
    L_F(z)=−L_B(−z).

Therefore

    O_T(z)=(L_B+L_F)/2 = Σ_(m odd) a_m z^m/m,
    E_T(z)=(L_B−L_F)/2 = Σ_(m even) a_m z^m/m.

This is the promised even/odd duality under z↔−z. It is degree parity in the logarithmic expansion, distinct from vertex-order parity and from exchanging the tournament with its converse.

For positive z in the stated domain, O_T and E_T are strictly positive. Moreover,

    2O_T=log[F_T(z)B_T(z)],
    2E_T=log[B_T(z)/F_T(z)]=log(det(I−z²K²)⁻¹).

Both log arguments are rational greater than one. Thus all four L_F,L_B,O_T,E_T are real transcendental, with no unproved claim about algebraic independence. The sum and difference do not cancel because every nonzero eigenvalue contributes a positive odd/even series. The excluded n=1 or z=0 has all logarithms zero, so that boundary is essential.

For the transitive T₄,

    F_T(z)=(1+6z+z²)²,
    B_T(z)=(1−6z+z²)⁻².

At z=1/18, a fully explicit three-layer chain is

    integer sign matrix
      → algebraic irrational energy 3±2√2
      → rational F=(433/324)², B=(324/217)²
      → real transcendental logarithms.

In particular O_T=log(433/217) and E_T=log(104976/93961), both positive real transcendental. The field extension at the spectral step and the analytic operation at the logarithmic step are exactly where the arithmetic distinctions change.

## 7. Explicit Dedekind cuts, not decimal evidence

There are two complementary exact constructions.

Trace construction. Since a_m≤nR^m, write q=zR<1 and let s_N=Σ_(m=1)^N a_m z^m/m. Then

    s_N < L_B(z) ≤ s_N + n q^(N+1)/[(N+1)(1−q)].

The upper error tends to zero. Restricting the sums to even or odd m gives corresponding rational cut bounds for E_T and O_T. The analogous absolute remainder bound works for L_F. This uses the actual trace coefficients, not merely their order.

Universal logarithm construction. For any rational A>1 set r=(A−1)/(A+1), so 0<r<1, and define

    ℓ_N=2 Σ_(j=0)^(N−1) r^(2j+1)/(2j+1),
    u_N=ℓ_N+2r^(2N+1)/[(2N+1)(1−r²)].

The identity log A=2 arctanh r gives

    ℓ_N < ℓ_(N+1) < log A < u_(N+1) < u_N,
    u_N−ℓ_N → 0.

The upper bound follows by replacing all tail denominators by 2N+1; strictness follows because later denominators are larger. The upper endpoints decrease by direct subtraction. The lower Dedekind cut is exactly

    D_(log A)={q∈Q: there exists N with q<ℓ_N}.

Thus a finite combinatorial object supplies a computable real transcendental cut, with a certified rational error budget. One may use A=F,B,FB,B/F; division by two gives O_T,E_T. No floating-point transcendence inference occurs.

All real numbers admit rational cuts, including rational numbers, algebraic irrationals, and transcendental reals. A rational a has Q\D_a with a least element a; an irrational does not. This infinite-cut distinction does not make a finite prefix an arithmetic certificate. THM-4088 already proves that every finite increasing rational prefix admits all three limit types. The present positive replacement is to retain a defining operator, a convergent series, and a tail bound, then apply a proved arithmetic theorem.

The real algebraic field contains √2 and every algebraic spectral energy above, but it is not Dedekind complete: the rational cut for log(433/324) has no supremum in that field, because such a supremum would equal its real limit and therefore be transcendental. Its completion is the real numbers. This is a precise completion distinction, rather than the false dichotomy “real versus transcendental.”

There is also a strict no-go for a literal field duality: every ordered-field automorphism of R is the identity. It fixes Q, preserves all comparisons with rationals, and hence fixes each Dedekind cut. Therefore no order-and-field-preserving involution can exchange rational and irrational or algebraic and transcendental reals. The logarithm/exponential and Cayley maps used here are different typed operations; they are not such field automorphisms.

## 8. What the spectrum forgets

The transitive tournament T₃ and directed three-cycle C₃ both have

    det(xI−S)=x(x²+3),
    spec(K)={0,3,3},
    F(z)=(1+3z)²,        B(z)=(1−3z)⁻².

Nevertheless their numbers of directed Hamiltonian paths are respectively 1 and 3. Consequently the entire partition pair, all traces, all logarithmic cut functions, and all spectral CM/charge modes fail to determine the Hamiltonian-path count. The same full spectral observer even misses the presence of a directed 3-cycle in this smallest witness.

S↦K already identifies converse tournaments, since S and −S have the same square. Scalar spectral data lose eigenvectors and vertex labels in addition. No inequality or reconstruction for the Hamiltonian count is inferred from the arithmetic completion. A marked response, incidence data, or another target-bearing sidecar would be needed.

The boundary is useful: this is a valid bridge between tournament spectra, arithmetic fields and analytic completions, not a replacement for the original tournament.

## 8A. All eight block orientations: identical arithmetic, different paths

Theorem D. Fix a regular tournament T of odd order n≥3. Add singleton blocks a,b and orient every pair of blocks uniformly, allowing all eight quotient tournaments on {T,a,b}. Every extension has exactly the same skew spectrum, positive-energy spectrum, CM fields, Cayley phase multiset, charge multiset, and partition/logarithmic functions. Nevertheless:

- For each of the six transitive quotient choices, H(extension)=H(T).
- For each of the two cyclic quotient choices,

      H(extension)=3H(T)+2pc_T(2),

  where pc_T(2) counts unordered spanning covers of T by two nonempty vertex-disjoint directed paths, with singleton paths allowed.

Spectral proof. Write the three upper block signs as α,β,γ∈{±1}. The old W=1⊥ remains invariant. The residual orthonormal block is

    [[0, α√n, β√n], [−α√n,0,γ], [−β√n,−γ,0]],

whose characteristic polynomial is x(x²+2n+1) independently of all three signs. Thus for every choice

    det(xI−S_extension)=det(xI−S_T)(x²+2n+1).

This proves spectral equality, not equality of the eight matrices U, C_U, or K themselves. These extensions are skew/energy-cospectral: the shared spectra belong to S and K=−S², not necessarily to the 0/1 adjacency matrix A=(J−I+S)/2. Adjacency characteristic and permanental polynomials can distinguish them. In particular the five-vertex H=3 versus H=15 pair does not contradict corrected THM-506 recovery from adjacency characteristic plus permanental data through order eight; those are different, richer observers.

Path proof. In a transitive quotient, an original-vertex run cannot be interrupted by a singleton and then resume: the singleton's uniform orientation with T would need to point both ways. Hence all original vertices are contiguous, and their unique block order leaves exactly H(T) paths. For a cyclic quotient label it T→a→b→T. There are three one-run words T,a,b; a,b,T; b,T,a, each contributing H(T). The sole two-run word is P,a,b,Q, where P,Q partition the original vertices into two directed paths. Each unordered cover gives two ordered choices (P,Q), so it contributes 2pc_T(2). There are no three-run words because a→T is forbidden. This independently reproduces the inherited THM-3121 cycle-word kernel, with its factorial convention explicit.

For T=C₃, H(T)=3 and pc_T(2)=3. The eight extensions therefore have path counts 3 (six choices) and 15 (two choices), while all share

    spec(K)={0,3,3,7,7},
    F(z)=(1+3z)²(1+7z)²,
    B(z)=(1−3z)⁻²(1−7z)⁻².

The μ=7 CM field, charge 22, and transcendental normalized angle are identical in both path-count classes. Thus the arithmetic lens has a universal, operation-level information-loss certificate, not merely an accidental small cospectral pair.

The ambiguity is unbounded, not uniformly small. Cutting every Hamiltonian path at each of its n−1 gaps gives (n−1)H(T) incidences with unordered two-path covers. Every cover supports at most two concatenations, so

    2pc_T(2) ≥ (n−1)H(T),
    H(cyclic extension)/H(transitive extension) ≥ n+2.

Regular carousel tournaments exist for every odd n, and have a Hamiltonian path by their consecutive labels. Hence these skew/energy-cospectral pairs have unbounded Hamiltonian-count ratio as their order N=n+2 grows. Equivalently any single positive estimate depending only on this spectral data has symmetric multiplicative error at least √N on one member of the pair. This concerns the stated spectral observer, not algorithms retaining the full labelled matrix U. The inequality was identified in the independent arithmetic audit and follows from the displayed incidence proof. The cyclic quotient extension is strongly connected for every nonempty base, since routes may pass through T→a→b→T; the transitive quotient extension is reducible. Therefore this does not yet prove an unbounded ambiguity ratio inside a fixed strong-connectivity class or a fixed score sequence. Those targets would need new witnesses.


## 9. Reproduction and provenance

Run `python3 -B session37-arithmetic/checks.py > session37-arithmetic/checks.json` from the scratch root. The script enumerates every labelled tournament through n=5, checks rank parity and exact rational partition inequalities, checks the matrix Cayley/charge identity through n=4, checks symbolic logarithmic derivatives through degree 12 on seven declared controls, checks the smallest spectral/Hamiltonian collision, and verifies Paley Gram identities at primes 3,7,11,19,23. A symbolic three-dimensional block verifies the all-order extension; all eight block orientations are checked for regular carousel bases of orders 3,5,7, including the exact two-path-cover Hamiltonian formula. The P₁₆₃ case is checked by 162 character-correlation sums, not a large tournament census. It also constructs exact rational log-cut intervals for the explicit T₄ example. No finite census is the proof of the general theorems.

Inherited repository sources, relative to `session35-publication/clean-repo`:

- `01-canon/theorems/THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density.md`: order-only blindness and finite-prefix hostiles.
- `01-canon/theorems/THM-4095-exact-arithmetic-field-transport-gapped-pair-margins-and-order-tournament-blindness.md`: exact generated-field transport as a model for the present Q(c)=Q(μ) statement; not a dependency of its proof.
- `01-canon/theorems/LEM-004-tournaments-are-odd-functions.md`: skew orientation, odd arctanh formal logarithm, and DRT spectral scope. Its “CM” abbreviation in multiplicative sign-function tables means completely multiplicative, not CM field; these must not be conflated.
- `01-canon/theorems/THM-4666-class-numbers-are-ranking-excesses-of-kronecker-tournaments-and-local-optimality-of-the-natural-paley-order-excludes-a-real-zero.md`: recovered P₁₆₃/class-number provenance only; its class-number and L-function conclusions are not dependencies of the present spectral field map.
- `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`: inherited cyclic block-word/path-cover kernel; the exact two-block specialization is reproved above, with unordered covers and the factor 2 distinguished.
- `01-canon/theorems/THM-416-cm-density-quantum-totient-ladder.md`: recovered CM-field setting and definition only; no density-quantum classification or lattice conclusion is imported.
- `05-knowledge/results/session35_dimensions_dualities_20261023/lie/COULOMB_MIRROR_DUALITY.md` and `COULOMB_POWERING_ADDENDUM.md`: reciprocal quotient and branch arithmetic. The present Cayley spectral map lies in their elliptic, not positive-real-mirror, branch.
- `05-knowledge/results/session35_dimensions_dualities_20261023/charge-dynamics/ELLIPTIC_CHARGE_DYNAMICS.md`: normalized-angle transcendence mechanism and its no-doubling-genericity guardrail.
- `00-navigation/META-PATTERNS.md`: “Type every analogy and every implication,” “Controlled forgetting and unlabeled quotients require a sidecar,” and “Search the statement before the method” are applied; no new meta-pattern is promoted from one lane.

Classical arithmetic inputs were checked in Richard Elman's author-hosted `Lectures on Abstract Algebra`, https://www.math.ucla.edu/~rse/algebra_book.pdf : Theorem 73.7 (Lindemann–Weierstrass), Corollary 74.25 and Example 74.26(2) (nonzero logarithms of algebraic numbers), and Theorem 74.29/Corollary 74.30 (Gelfond–Schneider). These are CITED inputs; the new work is the explicit operator/field/cut application.

No Riemann Hypothesis, Collatz, general LRC, new VOA realization, physical-statistics classification, or generic equidistribution claim is made.
