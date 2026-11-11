# Physical Fourier filters for rooted length responses

**PROVED:** the length-response product, exact universal-filter criterion,
dyadic obstruction, and parity-vanishing order. **FINITE-EXACT:** an explicit
eight-vertex universal mod-eight block, with two independent counts of every
induced Hamiltonian count. **OPEN:** all dyadic filter orders and their use in
an all-target Chan–Pak controller. No external novelty is claimed.

## 1. Why retain path length?

The spin/Fourier theme becomes a precise question: can an ACTUAL tournament
operation erase selected Fourier channels of a source, without signed weights,
virtual projectors, or a division that has no combinatorial realization?

The inherited source state keeps the total U of paths starting at a selected
vertex, weighted by the Hamiltonian count of their complement. THM-5129's
head correction, at its stated pre-audit status, identifies the parity split.
We do not need its proof: the constructions below start from literal paths.
The earlier all-size inputs are Rédei parity, the ordinal product, and
THM-002 (`01-canon/theorems/THM-002-ocf.md`). The two-colour cover product is
related to THM-3166 (`01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md`),
but the variable here records VERTEX LENGTH in one slot, not number of colours.

The corrected near miss is extending the bare-apex formula through its missing
endpoint head. The hostile here is believing approximate Fourier mixing is
exact erasure. The least-used coordinate is the full length-residue vector.

## 2. Two positive polynomials and their actual operation

For a rooted tournament (B,y), define

    P_(B,y)(z) = sum_(directed paths R starting at y)
                         H(B - V(R)) z^|R|.

The path R is nonempty; H(empty)=1. For an unmarked tournament C define

    K_C(z) = sum_(S subset V(C)) H(C[S]) H(C-S) z^|S|.      (1)

Both have nonnegative integer coefficients. P(1)=U_B(y), K_C(1)=Q_C(2),
and K_C(0)=K_C's leading coefficient=H(C). This K is neither the usual
path-cover polynomial nor an adjacency characteristic polynomial.

**Theorem 1 (physical length convolution).** Append C ordinally after B,
keeping the old root y. Then

    P_(B=>C,y)(z) = P_(B,y)(z) K_C(z),
    K_(A=>C)(z) = K_A(z) K_C(z).                          (2)

*Proof.* A path beginning at y meets B first, in one nonempty path R_B, and
then C in one possibly empty path R_C. For its C-vertex set S there are
H(C[S]) choices for R_C. The complement is the ordinal join of the two
complements, so its Hamiltonian count is their product. Every choice glues
back uniquely. The same argument applied to both ordered path slots proves
the second identity. Empty subsets are necessary for both equations. ∎

In the ring Z[z]/(z^m-1), write p_r for the coefficient sum at lengths
congruent to r. Ordinal attachment is cyclic convolution by K_C's residue
vector. It is a literal positive integer operation, not an abstract Fourier
projection. Fourier characters merely diagonalize this operation.

## 3. Exact filtering is possible only at dyadic moduli

Call C a **universal m-filter** if the m entries of the length-residue vector
of P_(B=>C,y) are equal for EVERY rooted source (B,y). We mean exact equality
of finite integer counts. Put S_m(z)=1+z+...+z^(m-1), with m>=2.

**Theorem 2 (exact criterion and obstruction).** The following are equivalent:

1. C is a universal m-filter.
2. K_C(z) is divisible by S_m(z) in Z[z].
3. For every nontrivial m-th root of unity zeta, K_C(zeta)=0.

Such a C can exist only if m is a power of two and |C|>=m-1.

*Proof.* Divisibility by the monic S_m is equivalent to the Fourier zeros.
Modulo z^m-1 one has S_m(z) f(z)=f(1)S_m(z), so divisibility implies a constant
residue vector, and multiplication by any P preserves constancy. Conversely
use the actual initial source consisting of its single marked vertex: P=z.
If zK_C has constant residue vector, so does K_C. Polynomial division then
shows S_m divides K_C.

Rédei parity applies to every induced subtournament, including H(empty)=1.
Consequently each summand of (1) is odd, and

    K_C(z) == (1+z)^n (mod 2), n=|C|.                    (3)

Thus S_m modulo two must divide (1+z)^n. Write m=2^a u with u odd. In
characteristic two, if u>1, z^u-1 has a nontrivial root different from 1
in the algebraic closure, and this root remains a root of S_m; it cannot
be a root of (1+z)^n. Hence u=1. Then S_m=(z+1)^(m-1) modulo two, forcing
n>=m-1. Reduction loses no leading coefficient since S_m is monic. ∎

This rules out universal 3-, 5-, 6-, 10-, and all nondyadic length filters
at EVERY order. It is not a prohibition on a particular source already
balanced modulo an odd integer. The universal source quantifier does work.
The result also does not prohibit nonordinal operations or extra path states.

For a transitive tail TT_t, K=(1+z)^t. It kills only the character -1 among
roots of unity. Thus no amount of ordinary tail padding is a universal
m-filter for m>2. Approximate mixing is not a substitute for exact erasure.

## 4. Concrete filters with all initial cases supplied

A singleton has K=1+z and gives the universal 2-filter. A directed triangle
C3, with rows (2,4,1), has

    H=3, K=3(1+z+z²+z³)=3(1+z)(1+z²).                  (4)

It is a universal 4-filter on three vertices, meeting the lower bound.
Its four output residue counts are exactly 3 U_B(y).

The following explicit tournament F8 has eight vertices. Rows are bitmasks:
bit j in row i means i beats j.

    (234,236,161,228,111,192,132,16).                   (5)

Every unordered pair is oriented exactly once. Its exact polynomial is

    K_F8(z)=231(1+z) S_8(z)
           =231(1+z)²(1+z²)(1+z⁴),
    coefficients (231,462,462,462,462,462,462,462,231).
                                                               (6)

Therefore F8 is a universal 8-filter. It has H=231 and Q(2)=3696, and the
eight residue counts after attachment are each exactly 462 U_B(y).
Equation (6) is a finite explicit identity, not a pattern conjecture:
`filters.py` supplies both subset path-DP and independent permutation-path
enumeration for ALL 256 induced Hamiltonian counts entering (1). The
polynomial multiplication/factorization is exact integer arithmetic.
We claim neither smallest possible order nor unique filter. The theorem
only puts its order between the lower bound seven and our upper bound eight.

Discovery provenance: deterministic Random(263092), 1,000 proposed labelled
order-eight tournaments, stopping at zero-based proposal795, yielded (5).
This was an exploratory positive search, not a preregistered density test.
The retained constructor and certificate no longer need the search.

A preliminary distinct order-eight block (42,120,99,68,13,216,17,95) has

    K=3(1+z)²(1+z⁴)(115z²-2z+115).

It kills primitive eighth-root modes but NOT fourth-root modes. It is a useful
hostile: killing the highest Fourier level alone is not a universal eight-
filter. Its composition with C3 does give an eleven-vertex eight-filter;
(5) is the more economical complete solution found subsequently.

## 5. Exact source service and the paid offset

**Corollary 3.** Given ANY rooted source (B,y) with D=H(B), U=U_B(y), append
k copies of F8 ordinally, k>=1, retaining y. This literal operation produces

    D'=231^k D, U'=3696^k U,
    p'_r = (3696^k/8) U  for every r modulo8.            (7)

Its endpoint count E'=0, because an old vertex cannot be the last vertex
of a spanning directed path when a nonempty block follows it. Consequently
both its bare-apex and positive-tail closures have the unified law

    H_closed = 231^k D + 2^t 3696^k U, t>=0,             (8)

where t=0 means the bare apex. There is no fictitious zero-length head:
the head correction vanishes for this actual root by its endpoint geometry.
The response remains uniform after EVERY further unmarked ordinal attachment,
since cyclic convolution preserves a constant vector.

For the real initial singleton D=U=1 and k=1, (7) gives eight channels462;
the bare and first-tail counts are3927 and7623. With C3 instead, the analogous
numbers are D'=3,U'=12, four channels3, and closures15,27.

These are paid-offset services: the block multiplies D by231 (or by3 for C3).
It does not preserve unit offset, prescribe arbitrary target increments, or
turn a residue projection into a new vertex-selection gadget. The usefulness
is an exact reusable finite source class for constructions that consume
length residues. A new consumer attaching physical meaning to those channels
is still required for an all-target Chan–Pak argument.

## 6. The order of parity erasure records odd-cycle geometry

For an odd-cycle packing F in C, write k=|F| and u=n-sum_(gamma in F)|gamma|.
Colour each cycle entirely into one of the two path slots. Applying OCF to
the two induced parts of (1) and summing over their partition gives

    K_C(z) = sum_F 2^k (1+z)^u prod_(gamma in F)(1+z^|gamma|). (9)

Every nontrivial cycle here has odd length>=3; unused vertices are singletons.
No bijection between Hamiltonian paths and OCF packings is needed: apply
the inherited count identity separately to each induced part.

**Theorem 4 (exact parity-jet order).** Let

    kappa(C)=min_F [n-sum_(gamma in F)(|gamma|-1)].

Then the multiplicity of z=-1 as a root of K_C is exactly kappa(C). Moreover,
writing K_C=(1+z)^kappa L_C, one has the positive integer certificate

    L_C(-1) = sum_(F attaining the minimum) 2^|F| prod_gamma |gamma| >0. (10)

*Proof.* For odd ell, (1+z^ell)/(1+z) evaluates to ell at -1. Each term of
(9) has exact vanishing order u+k and a positive leading quotient. The
smallest orders cannot cancel, giving (10). The empty packing makes the
minimum nonempty and finite. ∎

Thus kappa is the smallest number of components in a spanning permutation
cover by odd directed cycles and singleton fixed points. It is additive
under ordinal composition by (2). In particular every nonempty C erases
parity; a transitive n-block has kappa=n, while C3 has kappa=1. The mod-eight
block above has kappa=2 and exact leading quotient924. It therefore erases
more Fourier characters than the three-vertex transitive block, while having
a LOWER parity-jet order than that block. Fourier depth and derivative depth
are distinct resources. This is the precise replacement for a universal
small-dimension threshold in the snippets.

## 7. Connection and remaining obligations

Map: rooted paths with complementary H -> positive length polynomial ->
finite residue vector -> Fourier diagonalization. Preserved: literal ordinal
responses, exact residues, base cost, endpoint head. Lost by Fourier/scalar
projection: path and vertex identities. Sidecar: literal block, root address,
integer residue vector, and source offset. Cheapest decisive tests: singleton
initial source; primitive-eighth-only hostile; exact B8 induced-path certificate.

A promising next question is whether universal filters exist for every2^a,
and whether they can be built without inflating source offset prohibitively.
The lower bound |C|>=2^a-1 means no fixed-size block can erase arbitrary
Fourier depth. A different question is to realize a chosen residue channel
as a nonnegative physical selector rather than equalizing all channels.
Neither is settled by an abstract idempotent or by the present filters.
