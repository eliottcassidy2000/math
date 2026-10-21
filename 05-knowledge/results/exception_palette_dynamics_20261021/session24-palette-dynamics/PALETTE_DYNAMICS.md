# Twisted palettes, asymptotic rarity, and labelled signed exceptions

Session 24, 2026-10-21. **PROVED relative to the independently audited structural classification and cited classical dependencies; FINITE-EXACT controls independently replayed.** No claim of external literature priority. All colors are named, may be unused, and are preserved individually. They are neither proper/acyclic colors nor colors identified up to a color permutation.

## 1. Inheritance and scope

The independently audited Session 24 structural theorem classifies the tournaments of order at least two with no signed-regular binary coloring as

    X[C3^[k]],  X in {C3,R5,P7}, k>=0,

where R5=Z5{1,2}, P7=Z7{1,2,4}, and C3^[0]=K1. Its proof file is `../session24-structure/STRUCTURE.md`; the approved pre-promotion SHA-256 is recorded in SOURCES.json. Its Sections 8 and 9 already give the arbitrary twisted C3-palette transfer and unlabelled generating function. These are inherited, not new claims here. The full-family statements below use the independently audited classification; the same statements about the explicitly displayed family do not require completeness.

The Session 23 `signed_tournament_palette_rigidity_20261021/session23b-stabilizers-review/REVIEW.md` gives the ordinary P7 named-color polynomial c^7−7c^3+6c and signed polynomial c^7−7c^4−7c^3+7c^2+6c by exact subgroup inclusion–exclusion. Those seed polynomials and the universal signed distinguishing bound three are inherited. This note extends to nontrivial converse involutions on the block palette, then derives all-depth labelled enumeration and quantitative asymptotics.

Live board: free orbit palettes; involution fixed points; modular substitution; sparse order support; probability versus orbit fraction. Cheap hostiles: f=0, q=0,1,2, the untwisted f=q case, and C3[R5] versus R5[C3]. No Hamiltonian-path multiplicity, metagraph equivalence, Collatz-convergence, or LRC consequence is inferred.

## 2. Uniform transfer with the missing involution coordinate

Let B be a nonempty strong self-converse tournament. For c named vertex colors, let q be the number of ordinary-regular Aut(B)-orbits of colorings. Converse induces a well-defined involution alpha on these orbits; let f count its fixed points. For B=K1 use q=f=c and alpha=identity, without pretending there is a faithful negative identity on its one vertex.

For X in {C3,R5,P7}, the copies of B are the canonical maximal strong modules of X[B] (the singleton convention covers K1), since X is prime. Thus the full ordinary action is the wreath action, and every anti-map has a quotient anti-map and local anti-isomorphisms. A global coloring is ordinary-regular iff every local coloring is ordinary-regular and its palette word is regular for Aut(X). Quotienting the local ordinary actions identifies the remaining converse action with a quotient anti-map combined with alpha in every coordinate.

Writing Q,F for the resulting ordinary-regular orbit count and its converse-fixed count, respectively, gives:

    C3: Q=(q^3−q)/3,             F=f(q−1).
    R5: Q=(q^5−q)/5,             F=f(q^2−1).
    P7: Q=(q^7−7q^3+6q)/21,     F=f q(q^2−1)/3.

For q=0 these formulas are still valid because f=0. More generally every possible abstract palette involution has 0<=f<=q and f≡q (mod 2). The formulas depend only on this conjugacy type, not on a particular labeling of palette elements.

Proof for C3/R5. The ordinary quotient group is cyclic of prime order p=3 or 5, so only constant words are nonregular, giving Q=(q^p−q)/p. A specified reflection has one fixed position and (p−1)/2 exchanged pairs; combining with alpha gives f q^((p−1)/2) fixed words. Exactly f are constant. A regular word has at most one twisted anti-stabilizer, since two distinct ones have a nonidentity ordinary product. There are p reflections, and an ordinary regular orbit has p representatives. Consequently F=f(q^((p−1)/2)−1).

Proof for P7. Its ordinary group is the affine group x↦a x+b, a in {1,2,4}, b in F7, of order 21. There are seven order-three subgroups H_t, centered at t. Each has cycle type (1)(3)(3), so fixes q^3 words. Two distinct H_t generate a group containing a nonzero translation (compare maps x↦2x−t), hence a transitive group, and their common fixed words are exactly the q constant words. A nontrivial stabilizer has order divisible by 3 or 7; order-seven invariance is already constant. Thus there are q^7−7q^3+6q regular words, proving Q.

Fix tau:x↦−x. Twisted tau-fixed words number f q^3. Their ordinary stabilizer is normalized by tau: the palette involution commutes with position permutations and cancels in the conjugation. If such a word is fixed by H_t with t≠0, it is also fixed by H_−t, a distinct subgroup, hence constant and already fixed by H_0. Therefore nonregular twisted-fixed words are precisely those fixed by H_0. The two nonzero H_0-orbits are exchanged by negation; choosing one color and its alpha-image for these orbits, and an alpha-fixed color at 0, gives exactly f q such words. Thus tau fixes f(q^3−q) regular words. All seven anti-involutions are conjugate reflections, and an ordinary-regular word fixed by any anti-map has that map as its unique anti-stabilizer and that map is an involution. Double-counting incidences and dividing by 21 yields F=7f(q^3−q)/21. This also explains why the fourteen anti-elements of order six contribute nothing on the regular locus.

The number of signed-regular orbits is (Q−F)/2. The actual number of signed-regular vertex colorings is

    |Aut(X[B])| (Q−F),    |Aut(X[B])|=|Aut(X)| |Aut(B)|^|X|.

These factors must not be replaced by the signed group order without also halving the orbit count. Substitution with f=q=c recovers the inherited seed polynomials, providing a normalization control.

## 3. Triangle tower dynamics and two distinct notions of rarity

Let B_k=C3^[k], n_k=3^k, and let q_k,f_k refer to c named colors. Then

    q_0=f_0=c,
    q_(k+1)=(q_k^3−q_k)/3,
    f_(k+1)=f_k(q_k−1),
    a_k=|Aut(B_k)|=3^((3^k−1)/2).

The pair (2,2) is fixed: binary palettes never escape the exception. The scalar polynomial P(q)=(q^3−q)/3 has derivative P'(2)=11/3. This is a repelling fixed point in the real scalar extension, not a statement about arbitrary integer iterations or Collatz dynamics.

For integer c>=3, q_1>=8 and q_(k+1)>q_k. Moreover

    f_(k+1)/q_(k+1)=(f_k/q_k) 3/(q_k+1).

In particular f_k/q_k <= (3/(c+1))^k→0. Almost all ordinary-regular ORBITS are therefore signed-free. Since all ordinary-regular orbits have the same cardinality, the same relative statement holds among ordinary-regular colorings. It does NOT say most of the c^n arbitrary vertex colorings are regular.

Indeed let R_k=a_k q_k and p_k=R_k/c^(3^k). Direct substitution gives

    R_(k+1)=R_k^3 (1−q_k^−2),
    log p_k /3^k = S_k := sum_(j=0)^(k−1) 3^(−j−1) log(1−q_j^−2).

For c>=2 define the negative constant

    lambda_c = sum_(j>=0) 3^(−j−1) log(1−q_j^−2) < 0.

Monotonicity q_j>=q_k for j>=k and −log(1−x)<=x/(1−x) show the rigorous tail bound

    0 <= S_k−lambda_c <= −log(1−q_k^−2)/(2*3^k)
                           <= 1/[2*3^k*(q_k^2−1)].

Thus for c>=3,

    1 <= p_k/exp(lambda_c 3^k) <= exp(1/[2(q_k^2−1)]) → 1.

Ordinary regularity is exponentially rare in vertex count, despite asymptotically overwhelming signed-freeness conditional on ordinary regularity. For k>=1 the signed-regular probability is p_k(1−f_k/q_k), hence also asymptotic to exp(lambda_c 3^k). The case k=0 must not use a formal signed identity as a faithful vertex action.

As an optional orbit-growth normalization, beta_c=log c−(log 3)/2+lambda_c is positive for c>=3: lambda_c>=−1/[2(c^2−1)] and the resulting lower bound is positive already at c=3 and increasing. Hence

    q_k ~ sqrt(3) exp(beta_c 3^k).

For c=2 the exact formula is p_k=(2/sqrt(3))(sqrt(3)/2)^(3^k), while all signed-regular counts are zero. It is a useful hostile: the asymptotic multiplicative prefactor one above depends on q_k→infinity.

All-depth results for X[B_k] now follow by applying the transfer once to (q_k,f_k). More strongly, all three exceptional towers have the SAME exponential rate per vertex for every c>=3. If d=|X| and n=d3^k, their ordinary-regular probability is p_k^d h_X(q_k), where

    h_C3(q)=1−q^−2,
    h_R5(q)=1−q^−4,
    h_P7(q)=1−7q^−4+6q^−6.

Each h_X tends to one. Also F/Q tends to zero by the transfer formulas and f_k<=q_k. Therefore both ordinary-regular and signed-regular probabilities for X[B_k] are asymptotic to exp(lambda_c n), independently of the exceptional seed X. This equality of exponential rates is a statement about these three families, not a universal tournament coloring law.

The dependence on f cannot be dropped: abstract palettes with the same q but different f have identical ordinary counts and different signed counts.

## 4. Exact labelled enumeration and a non-D-finite EGF

For d=|X| in {3,5,7}, a_X in {3,5,21}, and n=d3^k,

    |Aut(X[B_k])| = a_X a_k^d = a_X 3^((n−d)/2),
    number of labelled copies on [n] = n!/[a_X 3^((n−d)/2)].

Canonical strong-module preservation is essential here; otherwise a displayed wreath subgroup would give only a lower bound for the full automorphism order. The three order families are disjoint. Together with the independently audited classification these formulas count ALL labelled signed binary exceptions, not merely examples.

Let ell_n be this labelled count, zero off the displayed support, and let L(z)=sum ell_n z^n/n!. Then

    L(z)=z^3/3+z^5/5+z^7/21+L(z^3/3).

Indeed bottom triangle substitution sends T of order m to T[C3] and multiplies its automorphism order by 3^m; the coefficient of z^(3m) is correspondingly (ell_m/m!)/3^m. This coefficient substitution is not the ordinary unlabelled equation E(z)=z^3+z^5+z^7+E(z^3), inherited from the structural note.

Both E and L are non-D-finite over any characteristic-zero field containing their coefficients. Their coefficients have infinitely many nonzeros and arbitrarily long zero runs before later nonzeros. For example between 7*3^k and 9*3^k there are no supported orders. If a coefficient sequence satisfied a polynomial-coefficient recurrence of finite order r, after all integer roots of its nonzero leading coefficient, any run of r zeros would force all later entries to vanish. This contradicts the support. Order r=0 would itself force eventual vanishing. D-finite power series have P-recursive coefficient sequences, giving the claim. This proves neither a natural boundary nor any Hamiltonian-path-count gap theorem.

## 5. Exact checks and loss ledger

Reproduce from this directory:

    python3 verify_palette_dynamics.py --output receipt.json
    python3 -O verify_palette_dynamics.py --output receipt-O.json
    cmp receipt.json receipt-O.json

The independent checker constructs every full seed automorphism/anti-automorphism group by testing all vertex permutations. It enumerates all q-color words for q=0..6 and every involution type (q,f); ordinary regularity is checked against every nonidentity ordinary element, and fixed ordinary orbits by direct canonical orbit representatives. A second reflection-incidence count is compared. For q<=4 every regular word is tested against every anti-map, checking uniqueness and raw signed counts. All checks use explicit exceptions, so Python optimization does not remove gates. Exact integer/rational checks cover tower recursions for c=2..6 and depths 0..4, and labelled recursions for each seed through depth five. These finite tests corroborate but do not replace the proofs.

The source-to-target map is local coloring→ordinary orbit palette→outer regular palette orbit. It preserves existence of ordinary/signed stabilizers, but loses orbit sizes and converse behavior unless a_B and f are retained. The passage from an isomorphism class to its labelled count loses canonical decomposition unless the full automorphism order is proved. The passage from regular-orbit fractions to all-coloring probability loses the normalization R_k/c^n; the binary fixed point is the cheapest decisive hostile.
