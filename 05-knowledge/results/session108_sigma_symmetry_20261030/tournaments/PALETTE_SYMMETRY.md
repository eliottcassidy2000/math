# Palette symmetry is a joint quotient, not division by 168

Session 108, 2026-10-29. **PROVED extension with exact finite checks.** The named-color Burnside formula and color-and-edge refresh construction are inherited from Session 107 `colored/COLORED_RECURRENCE.md` and its credited concurrent count-polynomial lane; the rational Weyl decomposition is inherited from Session 105. We add the explicit finite-palette-group quotient, source stabilizers, and its observable modes. No new prime-distribution theorem, faster enumeration algorithm, literature-priority claim, or open-problem resolution is asserted.

## 1. Three actions that must not be conflated

Let the eight feasible answer colors be the dual space C=(F2^3)*, with a response at squareclass port v given by z(v). For the seven nonzero squareclass ports, GL3(2) permutes ports and acts contragrediently on answers. Adjoining the sigma(3)=4 port adds v=0, a coordinate identically zero; every code-coordinate permutation fixes it. The order of GL3(2) is (8−1)(8−2)(8−4)=168.

1. **Named ports/colors:** no group quotient is authorized by the arithmetic observation itself. Distinct named answers remain distinct.
2. **Linear code-coordinate symmetry:** port maps v→Av induce z→z A^−1. They fix the zero answer and preserve all linear relations. Marking the port associated with 19 means preserving its nonzero functional on C, restricting to an order-24 subgroup. It does not mean marking one answer color. The two stabilizer descriptions must not be silently substituted.
3. **Affine palette gauge:** maps z→Az+t act on the eight answer colors and preserve the feasible color *set*. These are legitimate symmetries of the abstract uniform colored-tournament model, but not automatically physical symmetries of named prime roots, their incidences, or the prime carrier. In particular a translation changes the zero answer. Affine translations of *port labels* are a different operation: they do not preserve the homogeneous code without adjoining the constant row (the first-order affine code).

The fixed-base limiting distribution from Session 107 can be pushed forward along any declared quotient; invariance of the limiting uniform law does not construct a permutation of arithmetic primes realizing that gauge. Finite prime prefixes need not have it.

## 2. Joint Burnside theorem

For a finite nonempty palette C of size q and any finite permutation group G acting on it, let a_n^G count colored n-vertex tournaments modulo both vertex relabeling and one global palette element. Set a_0^G=1 and a_−1^G=0. For an all-odd partition λ of n, let

    e(λ)=sum_i (λ_i−1)/2 + sum_(i<j) gcd(λ_i,λ_j),
    z_λ=product_m m^(c_m) c_m!,
    F_g(l)=#{c in C : g^l(c)=c}.

Then, in every degree,

    a_n^G = (1/|G|) sum_(g in G) sum_(λ all odd)
                  (2^e(λ)/z_λ) product_i F_g(λ_i).       (1)

**Proof.** Apply Burnside to S_n×G on labeled colored tournaments. An even cycle of a vertex permutation swaps an antipodal pair under a power, so cannot fix an oriented tournament, independently of palette behavior. For odd cycles the number of fixed tournaments is 2^e, by the inherited pair-orbit count. Along a vertex cycle of length l, a combined fixed coloring propagates by g or g^−1 depending on action convention, returning consistently exactly when g^l fixes the initial color. Thus that cycle contributes F_g(l), independently of the others. Summing vertex permutations by cycle type gives (1). QED.

This includes g≠identity terms and permits nonconstant colors around a vertex cycle. Neither a_n(q)/|G| nor replacing q by the number of color orbits is correct. Example: GL has two singleton-color orbits but a_2^GL=5, not the named-two-color count 4; and a_1(8)/168 is not even integral.

## 3. Attachment, deletion, and refresh survive the quotient

On rational vector spaces with joint-orbit bases, define U by all q·2^n attachments and D by all n deleted cards, retaining multiplicities. Relabeling and palette permutations biject choices, so these are well defined. The same local proof as Session 107 gives

    D_(n+1) U_n = 2 U_(n−1) D_n + q 2^n I,
    E_n=U_n/(q 2^n),                   [D,E]=I.           (2)

Thus the inherited rational harmonic decomposition applies verbatim:

    h_j^G=a_j^G−a_(j−1)^G,
    V_n^G=direct-sum_(j=0)^n E^(n−j) ker D_j,
    spectrum(K_n=E_(n−1)D_n/n): 1−j/n, multiplicity h_j^G. (3)

Equivalently these are the G-invariant harmonics of the named-color model, or its rational coinvariant quotient. The identification is not an integral lattice assertion. Both the color and incident edges are refreshed; retaining old colors defines a different chain.

Uniform labeled colored tournaments push forward to a reversible stationary measure. If Stab(T) is the stabilizer in S_n×G of a representative, then

    pi_n([T]) = n! |G| / (q^n 2^binom(n,2) |Stab(T)|).   (4)

This uses the **combined** stabilizer, which may contain a vertex permutation compensated by a nontrivial palette element. Replacing it by a color-preserving automorphism group is wrong. Reversibility also follows directly by lumping the reversible uniform labeled refresh chain. Irreducibility and self-loops survive quotienting.

## 4. First modes, sparsity, and exact gaps

If G has r orbits on C, then a_1^G=r and h_1^G=r−1. The invariant centered orbit counts are explicit backward eigenfunctions with eigenvalue 1−1/n.

* **Named eight colors:** r=8, h_1=7, gap 1/n.
* **GL3(2):** orbits {0} and C\{0}; h_1=1 and gap 1/n. An eigenfunction is N_0−n/8. Linear symmetry does not erase the zero-versus-nonzero contrast.
* **Marked nonzero port:** the order-24 stabilizer has orbits {0}, the three nonzero answers annihilating that port, and the four answers evaluating to one. Hence h_1=2 and gap 1/n. The port answer bit alone does not distinguish zero from the other three annihilators.
* **AGL3(2) on answer colors:** r=1, so h_1=0. At n=2 the unique tournament edge orders its two endpoints, and two-transitivity gives precisely two pair-color orbits (equal/different). Thus a_2=2, h_2=1, and the gap is exactly 2/n for n≥2. At n=1 its state space is a singleton, so there is no nontrivial relaxation mode.

For the affine quotient an explicit first surviving eigenfunction is

    F(T)=sum_(i<j) 1[color_i=color_j] − binom(n,2)/8.

A refreshed vertex meets any fixed pair with probability 2/n, after which equality has probability 1/8. Consequently K_n^T F=(1−2/n)F. The function is nonzero on quotient states for n≥2, proving the mode exists directly.

For a sparse character description put χ_v(z)=(−1)^(z(v)). The seven nonconstant Walsh characters span the named singleton contrasts. GL permutes them transitively; their sum is 8·1[z=0]−1, giving the surviving zero contrast. A marked-port stabilizer has two orbits on nonzero character indices, giving two invariant contrasts. Palette translations multiply χ_v by χ_v(t), killing every nonconstant degree-one invariant. Thus a quotient can remove observables and accelerate its own relaxation, without speeding the original named process or improving an arithmetic mixing bound.

Affine invariants retain higher joint information: AGL3(2) is transitive on ordered triples of distinct colors, but four distinct colors split according as z1+z2+z3+z4=0 or not. In characteristic two this sum condition is translation invariant; in the second case the three differences from z1 are a basis. These are the two orbits. Forgetting named coordinates therefore does not mean independent per-vertex zero storage; relational memory survives.

## 5. Source-free counts need palette stabilizers

A universal source is unique. Its color can be moved to one representative c of its G-orbit; thereafter the remaining tournament is quotiented only by G_c, not all of G. Therefore for n≥1

    source_bearing_n^G = sum_(c representatives of C/G) a_(n−1)^(G_c),
    b_n^G = a_n^G − sum_c a_(n−1)^(G_c).               (5)

This is a bijective argument, not a division heuristic. It recovers q a_(n−1) for named colors. For GL there are two terms: a_(n−1)^GL plus a_(n−1)^(stabilizer of nonzero answer). The latter fixes an answer, not a port. For AGL, transitivity reduces the answer to a_(n−1)^GL. Already at n=2, AGL has two source-bearing classes, whereas a_1^AGL=1.

The harmonic dimension is still a_n^G−a_(n−1)^G and is generally not b_n^G. If G fixes a distinguished color c0, adding a source of that color gives a bijection from previous-level G-classes; classes without such a source positively count h_n^G. If no color is fixed, this argument is unavailable. Source-deletion δ remains well-defined and satisfies δU−Uδ=qI, since all q newly sourced attachments contribute even when some share a quotient class. Its kernel contains signed differences of source lifts and is not the source-free coordinate subspace. Source stripping retains correlated colors along the prefix under a single global group element; the named identity A=B/(1−qt) cannot be reused for quotient counts.

## 6. Exact checks and limits

| Palette group | order | a_0,...,a_4 | h_1,h_2,h_3 |
|---|---:|---|---|
| named | 1 | 1,8,64,688,11008 | 7,56,624 |
| GL | 168 | 1,2,5,24,188 | 1,3,19 |
| marked-port stabilizer | 24 | 1,3,11,70,756 | 2,8,59 |
| affine palette | 1344 | 1,1,2,8,48 | 0,1,6 |

`verify_palette.py` builds GL and AGL explicitly as permutations of eight vectors and builds the marked-port subgroup by preserving its bit functional. It evaluates (1) through n=8 with rational arithmetic. Independently, it enumerates all labeled colored tournaments and direct S_n×G orbits through n=3, constructs U,D, checks (2) through input degree 2, and checks exact spectra, stationary measures, detailed balance, and (5) through degree 3. Named eight-color values are computed by the formula, not independently enumerated in this checker. Finite tests do not replace the all-degree proofs.

    python3 tournaments/verify_palette.py > tournaments/verification.json
    python3 -O tournaments/verify_palette.py > tournaments/verification.optimized.json
    cmp tournaments/verification.json tournaments/verification.optimized.json

Requires SymPy (tested 1.14.0); guards use explicit exceptions. Higher-order formula evaluations are not a faster class generator or a new scalar recurrence. No publication was performed by this lane.

## 7. Relations become objects: exact binary-subspace classification

This reframing was supplied by the root integration lane and is proved and independently checked here. It identifies what the gauge quotient remembers, not just how many orbits it has.

Represent an ordered n-tuple of colors in F2^r as the columns of Z∈F2^(r×n). The palette action GL_r is left multiplication; the vertex action permutes columns. Then

    GL_r\{color tuples} ↔ {W≤F2^n : dim W≤r},
    [Z] ↦ rowspace(Z).                                  (6)

**Completeness proof.** Left multiplication preserves rowspace. Conversely any two r-row matrices with the same rowspace reduce by invertible row operations to the same reduced row-echelon matrix (including their zero rows); composing the row operations gives an element of GL_r relating them. Every subspace of dimension ≤r is represented by choosing a basis and padding with zero rows. Column permutations carry W to its coordinate-permuted subspace. Thus joint orbits are exactly **tournaments decorated by binary subspaces, modulo simultaneous vertex-coordinate permutation**. Equivalently retain the relation code W^⊥={u:Zu=0}; it records exactly which selected colors XOR to zero. Rank alone is not this object.

For affine palette transformations Z→AZ+t·1_n, the complete invariant is instead

    W_hat = span(1_n, rows Z),
    1_n∈W_hat,       dim W_hat≤r+1.                     (7)

For n≥1, normalize the first color to zero by subtracting it from every column. The normalized rows all vanish in coordinate 1, and their span equals W_hat intersect {x:x_1=0}. Hence equal augmented rowspaces give equal normalized rowspaces, and (6) supplies a linear transformation between the normalized tuples. Undoing the translations proves affine orbit completeness. Conversely an affine transformation preserves W_hat. Every subspace containing 1_n and of dimension ≤r+1 has a complement to span(1_n) of dimension ≤r; its basis supplies Z. At n=0 there is just the empty tuple and the zero space, handled separately. The dual affine relation code is W_hat^⊥: precisely the zero-sum color relations supported on an even number of vertices. This augmentation is on the **color matrix**; it does not claim that translating the original named port labels preserves their homogeneous answer code.

For the stabilizer of a fixed nonzero port functional λ on F2^r, the complete invariant is

    (W,w),      W=rowspace(Z),       w=λZ∈W,
    dim W≤r,    and w=0 implies dim W≤r−1.              (8)

To prove completeness, choose coordinates λ=(1,0,...,0). Then Z has first row w, and the permitted row transformations are

    [ 1  0 ]
    [ b  B ],       B∈GL_(r−1), b∈F2^(r−1).

If w=0, arbitrary invertible operations on the remaining rows classify their rowspace, and at most r−1 independent rows exist. If w≠0, work modulo span(w): the remaining r−1 rows span W/span(w). Equality of this quotient rowspace permits invertible B transforming one tuple of quotient rows to the other, by the same row-equivalence argument as (6), now with values in F2^n/span(w). Differences of the actual rows are multiples of w; the vector b removes them. This proves completeness and, by choosing bases, existence for all pairs satisfying (8). In particular preserving the observed port-answer word w does not mean preserving a palette point.

Deletion of a vertex acts on these objects by the corresponding coordinate projection (and on w by the same projection); it commutes with the classification. Attachment must retain the allowed extension subspace, not merely its dimension. For the affine object, the all-ones vector projects to the shorter all-ones vector. Thus the quotient admits exact operations on relational objects, while forgetting W or W_hat down to rank generally destroys information needed for extension.

### The small numbers, with their types kept separate

For r=3, the unique two-vertex tournament has no automorphisms and F2^2 has five subspaces. Hence a_2^GL=5, and since a_1^GL=2, h_2^GL=3.

At n=3 there are two tournament types. The transitive triple has trivial automorphism group; its decorations are the 16 subspaces of F2^3. The directed triangle has automorphism group C3 acting by cyclic coordinate rotation; those 16 subspaces have eight orbits. Therefore

    a_3^GL=16+8=24,          h_3^GL=24−5=19.

The affine augmented-space condition gives five transitive decorations and three cyclic ones, hence a_3^AGL=8. For the marked-port pairs (8), it gives 50 transitive decorations and 20 cyclic ones, hence a_3^marked=70. These are structurally derived, independently enumerable counts.

The integers **5, 3, 19** thus occur respectively as a degree-two class count and degree-two/degree-three harmonic dimensions in this specific GL_3 quotient. The last integer is not thereby the prime label 19, an order modulo 19, a divisor sum, or a shared arithmetic mechanism. The coincidence of numerical values does not identify these differently typed objects.

`verify_relations.py` uses no color-matrix orbit enumeration and no Burnside formula. It enumerates XOR-closed subsets of F2^n for n≤3, filters the augmented and marked conditions, and directly computes cyclic-coordinate orbits at n=3. Its normal and optimized receipts agree byte for byte:

    python3 tournaments/verify_relations.py > tournaments/relations.verification.json
    python3 -O tournaments/verify_relations.py > tournaments/relations.verification.optimized.json
    cmp tournaments/relations.verification.json tournaments/relations.verification.optimized.json
