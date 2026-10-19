# Pythagorean level counts: an exact ternary digit law and a natural boundary

19 October 2026. **PROVED**, using the **CITED** Pólya–Carlson theorem; exact finite controls accompany the proofs. Independently [reviewed](review-root/REVIEW.md) for the stated scopes, including the count-tower inversion checked by the coordinating reviewer. No historical-priority claim is made. The constant occurring below is proved irrational, but its transcendence remains open.

## 1. Result and inherited mechanism

Let c_n count the words of length n in the Berggren ternary tree whose primitive Pythagorean triple (a,b,c) satisfies c²>4ab. Equivalently, these are the primitive shapes giving positive arithmetic-progression square-sum triangle labels. Set c_0=0. The existing [ternary-depth theorem](../pythagorean_ternary_depth_limit_20261019/TERNARY_DEPTH_LIMIT_AND_FINITE_MEMORY.md) proved

    c_n / 3^n → ρ = 0.06039052312486240241155650348537631073311780098403…,

with ρ irrational and |c_n−ρ3^n|≤2. The [positivity-language theorem](../dynamics_labeling_entry_20261019/BERGGREN_POSITIVITY_NONREGULAR.md) had already excluded every exact finite-state recognizer.

This continuation supplies a stronger, exact arithmetic obstruction:

1. If a_1,a_2,… are the spatial ternary digits of the lower positivity wall, and b_n is their length-n integer prefix, then

       c_n = 2b_n + a_(n+1),
       c_n − 3c_(n−1) = a_(n+1) − a_n ∈ {−2,−1,0,1,2}.

2. The discrepancy improves to the strict bound

       |c_n − ρ3^n| < 2/3.

3. Writing A(z)=Σ_(j≥1) a_j z^j and C(z)=Σ_(n≥0)c_n z^n,

       C(z) = (1−z)A(z) / [z(1−3z)].

   The germ C extends meromorphically throughout |z|<1, with exactly one pole, at z=1/3. The unit circle is its natural boundary.

4. Consequently C is transcendental over C(z) and is not D-finite. The sequence c_n satisfies no nonzero linear recurrence with polynomial coefficients, even one valid only eventually. The same conclusion holds for every arithmetic-progression subsequence c_(mn+r), with m≥1 and r≥0.

The exact generating-function identity, rather than an extrapolation from the counts, carries the proof. A finite-state obstruction alone would not establish the generating-function conclusion: nonregular languages can have rational counting series.

### Inheritance and connection contract

- Closest proved mechanism: the two ambiguous cylinders and irrational ternary-depth limit in the two linked Pythagorean notes.
- Recovered underused mechanism: the Pólya–Carlson step in [the September natural-boundary/Collatz note](../procgen_continuous_20260925_natural_boundary_mahler_harmonic.md), §1.2. Its basin object is different. No Collatz coverage claim is transported.
- Canonical hostile: integer counts with rational density, c_n=(3^n−1)/2, have bounded defects but a rational generating function. Irrationality is load-bearing in the decimation argument.
- Corrected near miss: the finite interval c_n/3^n≈ρ cannot itself prove a natural boundary. The integer defect series, its exact evaluation, and the published dichotomy are necessary.
- Preserved predicate: the exact Pythagorean positivity count at each depth.
- Map: pair the root probes 2 and 3 by the slope involution; encode their position in the wall cylinder as the next spatial ternary digit.
- Lost information: one isolated level count forgets which words are accepted. In contrast, the complete count tower in this fixed two-probe family recovers the wall itself, by Corollary 1.2 below; that recovery depends on the family and orientation sidecars.
- Sidecar: orientation of each Möbius prefix and the actual lower wall, not an unlabelled branching count.
- Cheapest decisive test: direct Euclid-pair enumeration through depth ten versus the boundary-only count compiler.

## 2. Definitions and the two-root identity

Use the outer-first convention and put

    D=(1,∞),  L(t)=2−1/t,  M(t)=2+1/t,  R(t)=t+2.

For u=u_1…u_n let f_u=f_(u_1)∘…∘f_(u_n). Reversal is a bijection on words of each length, so this convention leaves c_n unchanged. The actual finite tree has root slope 2. Its positivity region is

    (1,μ) ∪ (λ,∞),
    μ=2−√3+√6−√2,
    λ=2+√3+√6+√2.

The rational Möbius involution S(t)=(t+1)/(t−1) satisfies

    S∘L=R∘S,  S∘R=L∘S,  S∘M=M∘S,
    S(λ)=μ,  S(2)=3.

Thus, on swapping L and R in every word, the condition f_u(2)>λ becomes f_(swap u)(3)<μ. Therefore

    c_n = #{|u|=n : f_u(2)<μ} + #{|u|=n : f_u(3)<μ}.        (1)

The probe 3 is an auxiliary rational slope, not a second primitive root node. Equation (1) is valid because of the explicit conjugacy and the word bijection.

## 3. The next digit is exactly the boundary correction

At depth n, order the 3^n branch-cylinder intervals f_u(D) from left to right. Let b_n be the zero-based index of the unique cylinder containing μ. Its spatial ternary digits satisfy

    b_0=0,  b_(n+1)=3b_n+a_(n+1),  a_j∈{0,1,2}.

Spatial order matters: an increasing prefix orders its children L,M,R, while a decreasing prefix orders them R,M,L. The M branch reverses orientation.

Every cylinder strictly left of μ contributes both probes in (1), and every cylinder strictly right contributes neither. The unique wall cylinder contributes the number of its two internal subdivision endpoints lying left of μ. Those endpoints are precisely f_u(2) and f_u(3), in one order or the other. Hence that number is its spatial child index a_(n+1).

### Theorem 1. Exact digit-count law. PROVED.

For every n≥0,

    c_n = 2b_n+a_(n+1).                                     (2)

For every n≥1,

    δ_n := c_n−3c_(n−1)=a_(n+1)−a_n.                         (3)

Proof. The preceding cylinder count proves (2). Substituting b_n=3b_(n−1)+a_n into (2) gives (3). At n=0 the wall lies in the leftmost first-level cylinder, so a_1=0=c_0. No sample can equal μ, since samples are rational and μ has degree four. ∎

The first counts are

    0,0,0,2,5,15,44,132,396,1188,3566,10698,32094,… .

The defect values 0, ±1 and ±2 all occur in the exact checked prefix. A long apparently stable recurrence is a trap: δ_n=0 for every 11≤n≤56, but δ_57=−1. Constant-digit runs explain such plateaux without making them permanent.

### Corollary 1.1. Improved discrepancy. PROVED.

Let B=Σ_(j≥1)a_j/3^j=ρ/2. Then, for every n≥0,

    |c_n−ρ3^n|<2/3.                                        (4)

Proof. Write 3^nB=b_n+u_n, where

    u_n = Σ_(j≥1) a_(n+j)/3^j.

Its first digit is a_(n+1), so a_(n+1)/3 < u_n < (a_(n+1)+1)/3. Both inequalities are strict: a terminating or eventually all-2 expansion would make B rational. By (2), c_n−ρ3^n=a_(n+1)−2u_n. For a_(n+1)=0,1,2 its possible intervals are respectively (−2/3,0), (−1/3,1/3), (0,2/3). ∎

### Corollary 1.2. The complete count tower recovers its wall. PROVED.

More generally, for any irrational α>1, define

    c_n(α)=#{|u|=n:f_u(2)<α}+#{|u|=n:f_u(3)<α}.

The same proof gives c_n(α)=2b_n(α)+a_(n+1)(α). Consequently

    a_(n+1)(α)=c_n(α)−2Σ_(j=0)^(n−1)c_j(α).                (4a)

Thus counts through depth n recover the first n+1 spatial ternary digits exactly. The entire count sequence recovers the infinite branch itinerary, using the two-state orientation decoder. The shrinking-cylinder lemma then recovers α uniquely. In this precisely specified family, the complete count tower is therefore injective in its irrational wall. This does not assert that arbitrary languages can be reconstructed from their counting sequences.

Proof of (4a): use a_1=c_0 and telescope (3). Every further assertion follows from the inherited orientation decoder and shrinking-cylinder property. ∎

## 4. Exact generating function and analytic obstruction

Set

    A(z)=Σ_(j≥1)a_j z^j,
    C(z)=Σ_(n≥0)c_n z^n,
    D(z)=(1−3z)C(z).

The prefix recurrence gives Σ_(n≥0)b_n z^n=A(z)/(1−3z). Summing (2) gives the formal identity

    C(z)=2A(z)/(1−3z)+A(z)/z
        =(1−z)A(z)/[z(1−3z)],                            (5)
    D(z)=(1−z)A(z)/z.                                   (6)

The apparent singularity at zero is removable. The coefficients of D are δ_0=c_0=0 and the bounded integers δ_n in (3). The formulas hold analytically for |z|<1/3, and the right sides give the continuation thereafter.

### Why the digit series is not rational

The inherited wall theorem proves B irrational. For completeness, the mechanism is: rational B has eventually periodic base-three digits; the two-state orientation rule then makes the branch itinerary eventually periodic. Shrinking cylinders identify equal infinite itineraries, so some nonidentity rational Möbius word fixes an inverse iterate of μ. This forces μ to have degree at most two. But μ has minimal polynomial

    X^4−8X^3+2X²+8X+1

and degree four, a contradiction. The shrinking-cylinder and degree arguments are proved in the linked inheritance notes.

If a rational power series has coefficients in a finite set, its coefficient sequence is eventually periodic: beyond its polynomial numerator, its fixed-order constant-coefficient recurrence determines the next coefficient from a tuple of finitely many previous coefficients. There are finitely many such tuples. Repetition therefore forces a periodic tail. This argument works even if the putative recurrence coefficients are complex.

Thus A cannot be rational, because an eventually periodic digit sequence evaluates to a rational number at z=1/3. Its radius is exactly one: bounded digits give radius at least one, while radius greater than one would force the integer coefficients eventually to vanish.

### Theorem 2. Natural boundary and non-D-finiteness. PROVED with cited input.

The circle |z|=1 is a natural boundary of A and D, and of the meromorphic continuation of C. Inside that disk, C has a unique pole, a simple pole at z=1/3 with residue −ρ/3. In particular C is transcendental over C(z) and not D-finite.

Proof. The Pólya–Carlson theorem says that an integer-coefficient power series of radius one is rational or has the unit circle as a natural boundary. Apply it to A, which is not rational.

Multiplication by the nonzero rational factors in (5)–(6) transfers the obstruction at every boundary point except possibly z=1. A continuation through a neighborhood of z=1 would also continue through nearby unit-circle points different from 1, which is impossible. Thus the entire circle remains a natural boundary.

A is analytic in |z|<1, and at z=1/3 the numerator in (5) gives

    lim_(z→1/3) (1−3z)C(z)=2A(1/3)=ρ>0.

This proves the simple pole and its residue, and shows that there are no other poles in the open unit disk. A D-finite function satisfies a linear differential equation with polynomial coefficients; it can be continued along paths avoiding the finitely many zeros of the leading coefficient. Such a function cannot have an entire circle of singularities. Every algebraic function is D-finite, so algebraicity is also excluded. ∎

The equivalence between D-finite series and polynomial-coefficient linear recurrences follows directly by coefficient extraction in the differential equation, and conversely by summing the recurrence. Hence there is no such nonzero recurrence for c_n, including an eventually valid one. This does not rule out nonlinear recurrences, arithmetic algorithms, or efficient boundary computation.

**Published input and primary pin.** F. Carlson, “Über Potenzreihen mit ganzzahligen Koeffizienten,” Math. Z. 9 (1921), 1–13. The exact theorem used here is restated as Theorem 2.2 in J. P. Bell and J. C. Lagarias, “3x+1 inverse orbit generating functions almost always have natural boundaries,” Acta Arith. 170.2 (2015), 101–120, DOI https://doi.org/10.4064/aa170-2-1; [primary article PDF](https://www.impan.pl/shop/publication/transaction/download/product/82190), pp.107–108. An [English translation of Carlson's original paper](https://people.math.ethz.ch/~airibar/Polya_Carlson.pdf) is included as an appendix to Aitor Iribar López's exposition. These sources were checked 19 October 2026. The original theorem is cited, not claimed as new here.

## 5. No arithmetic-progression subsequence repairs the recurrence

The following general mechanism is useful beyond this tree.

### Lemma 3. Irrational exponential density with bounded discrepancy. PROVED with cited input.

Let B≥2 be an integer and let x_n be integers satisfying x_n=βB^n+O(1), where β is irrational. Then

    X(z)=Σ x_n z^n

has meromorphic continuation to |z|<1, a simple pole at z=1/B, and the unit circle as a natural boundary. It is not D-finite.

Proof. Put e_n=x_n−βB^n, which is bounded. Then

    X(z)=β/(1−Bz)+Σe_nz^n.

Thus the claimed continuation and pole are established. The regularized series E(z)=(1−Bz)X(z) has bounded integer coefficients x_n−Bx_(n−1), apart from its initial term, and satisfies E(1/B)=β. If E were rational, its finite-valued integer coefficients would be eventually periodic, and its value at 1/B would be rational. Contradiction. Its radius is exactly one, as an integer series of larger radius would be a polynomial. Pólya–Carlson gives the natural boundary. Division by 1−Bz preserves it, and the same finite-singularity argument excludes D-finiteness. ∎

Apply this with

    x_n=c_(mn+r),  B=3^m,  β=3^rρ.

Equation (4) gives the bounded discrepancy and β is irrational for every m≥1,r≥0. Thus every such subsequence has the same analytic obstruction, with its pole moved to z=3^(−m). This rules out attempts to repair an exact polynomial-coefficient recurrence merely by selecting even depths, odd depths, or any other fixed depth progression.

## 6. Exact computation and scope

The boundary compiler follows one irrational wall, retaining its integer Möbius matrix, orientation and ternary prefix. Equation (2) requires only the next child digit. Producing all counts through depth N uses O(N) arithmetic steps on integers with O(N) bits; this is an arithmetic-operation statement, not a unit-cost bit-complexity claim. The supplied implementation stores the returned lists and repeatedly constructs small rational samples for clarity.

Run:

    python3 -B verify_ternary_natural_boundary.py
    python3 -B -O verify_ternary_natural_boundary.py

The two modes produce identical JSON. The declared controls comprise:

- The exact boundary count and digit-difference formulas through depth 2048.
- An independently coded chronological Euclid-pair tree through depth ten: 88,573 total nodes, checking the direct integer predicate c²>4ab.
- Telescoping D(1/3) partial sums through depth 128.
- All fixed residue subsequences for strides 1 through 8 within the compiled prefix.
- The rational-density hostile and the false stabilization plateau at depths 11–56.

Finite tests check normalization, word order, arithmetic and identities; they do not numerically prove a natural boundary or non-D-finiteness. Those conclusions come from the exact digit law, the inherited irrationality proof and the cited theorem.

The result makes no claim about transcendence of ρ, irrationality of odd zeta values, Euclid–Mullin or split-sum completeness, the global Collatz problem, or the hypotenuse-cutoff density (which remains the different value one third). A function being transcendental does not imply that its value at a particular rational argument is transcendental.
