# Prime-order symmetry folds the path-colour polynomial

Let T be a tournament with an automorphism g of prime order p. Necessarily p is odd: a nonidentity involution would reverse an arc. Let F be the subtournament on the fixed vertices, and let O₁,…,Oₘ be the other orbits, each of size p. Define a directed weighted quotient D on these orbits: the weight dᵢⱼ of i→j is the number of vertices of Oⱼ beaten by any one vertex of Oᵢ, for i≠j. This is independent of the choice, and dᵢⱼ+dⱼᵢ=p.

A path in D has weight the product of its arc weights, and a singleton path has weight one. Let W_D(u)=Σⱼ wⱼuʲ, where wⱼ sums products of path weights over all unordered spanning covers by j directed paths. Empty graphs contribute 1. Write Q_T(t)=Σₖ cₖ(T)(t)ₖ, where cₖ counts unordered directed spanning path covers and (t)ₖ denotes a falling factorial, as in THM-3166.

**Fold law (as polynomials over the field with p elements):**

    Q_T(t) = Q_F(t) W_D(t^p−t) mod p.

Since m≥1, Q_T has the factor t^p−t modulo p. At each residue r, its derivative is −Q_F(r)w₁. For a fixed-point-free action the common multiplicity of all p residues is exactly the least j for which wⱼ is nonzero modulo p. This multiplicity can jump when the number of orbits changes even though the path-colour values at the p residues remain zero.

## Proof

Act by (g) on the unordered covers. Contributions of size-(p)
cover orbits vanish modulo (p). A fixed cover has some paths
individually fixed by (g), and other paths in (p)-cycles. A fixed
directed path has all vertices fixed because its ordered vertices are
rigid. Conversely no path in a (p)-cycle can contain a fixed vertex:
all its translates would contain that vertex. For the same reason a
path in a (p)-cycle cannot use two vertices from the same vertex orbit:
two translates would overlap. Its translates therefore use every vertex
of each vertex orbit it meets, exactly once.

It follows that the fixed cover consists of an independent cover of
(F) by (r) paths and a cover of the quotient by (j) directed paths,
each lifted to a (p)-cycle of paths. To lift a particular quotient path,
fix once and for all the first vertex in its first orbit. Its representative
then has (d_{ij}) choices for each next phase. There is no closing
condition because the object is a path. This constructs exactly one lift
per phase choice, and all lifts have exactly one representative with that
initial vertex. The total is thus the product of the weights specified
above. The corresponding cover of (T) has (r+pj) paths.

Finally (t)_(r+pj) ≡ (t)_r(t^p−t)^j as polynomials modulo (p)
by taking (j) consecutive complete residue blocks. Summing proves the
law. The derivative follows by the chain rule and W_D(0)=0 for m>0.

## New consequences from the skew quotient

Modulo (p), reversing a quotient path on (h) vertices multiplies its
weight by ((-1)^{h-1}). In a cover having an even-size path, reverse the
even-size path with the smallest vertex (under any fixed order of quotient
vertices). This pairs the covers with opposite weights. Consequently
(w_j≡0) unless (j≡ m mod 2): only covers entirely by odd-size
paths contribute. The top coefficient (w_m=1) comes from the singleton
cover. In particular, for any fixed-free action the fold polynomial is
nonzero, its roots have common multiplicity congruent to (m mod 2), and
for even (m) its first derivative vanishes at every residue. The
multiplicity is always at least two for even (m).

If all pairs of moved orbits have uniform incidence, every quotient
arc weight is 0 modulo (p). Then (W_D(u)≡ u^m). The converse is
false in general: already for (m=2), the sole nontrivial coefficient is
(w_1=d_{12}+d_{21}=p), independent of incidence. Thus an exact symmetry
quotient can miss all information about its first nonuniform interface.
A count-only quotient cannot be assumed to preserve an extension.

If (F) is empty and (m=3), put (a=d_{12},b=d_{23},c=d_{31}).
Summing the six directed orders gives
    W_D(u) ≡ u³ + 2(ab+bc+ca)u mod p.
For (p=3), exactly two nonzero incidences among (a,b,c)
force multiplicity one; with three nonzero incidences the multiplicity
is three if their signs all agree and one otherwise. With zero or one
nonzero incidence it is three. The signed information here is about
orientation of whole cross-orbit incidence, not just the order of the
automorphism group.

## A coarse count barrier for a large symmetry support

The majority orientation on the (p)-orbits and fixed singleton
vertices is a tournament. Take one of its directed Hamilton paths.
Each moved orbit is itself strong: its (p)-cycle is vertex-transitive,
whereas the strong components of a tournament are invariant in their
linear order. For every moved orbit let (h_i) be its internal
Hamiltonian path count. Every vertex is the first vertex of exactly
(h_i/p) internal paths, by transitivity. For each consecutive pair of
moved orbits in the quotient Hamilton path, there are at least
((p+1)/2) admissible choices of the next start for every choice of the
preceding endpoint. Fixed vertices have uniform incidences to any moved
orbit. Concatenation gives at least

    (∏ᵢ hᵢ) ((p+1)/(2p))^(m−r).

actual spanning paths, where (r) is the number of nonempty runs of
moved orbits in that quotient path (with no moved orbits interpret the
bound as 1). Since every (hᵢ ≥ p), it follows that for (m≥1),
(H(T) ≥ p((p+1)/2)^(m−1)). Thus a fixed target has only logarithmically
many moved (p)-orbits, independently of the total number of neutral
singleton vertices. For example, an automorphism of prime order
(p>H(T)) is already ruled out by the free path action; this stronger
inequality charges additional moved orbits individually. One can replace
(p) in each (h_i) bound by any valid strong-tournament floor at order
(p) (see THM-4741 for Busch's cited floor). The count barrier uses
majority paths and does not require the weighted polynomial.

## Source trail and proof boundary

The generating function is the repository's inherited falling-factorial
transform THM-3166. LEM-003 gives freeness on *one* spanning path;
THM-4747 counts rotation-fixed *cycles*. This fold counts all invariant
**covers**, where the action can permute distinct directed paths. The
phase weights are the same kind of interface data that survive the
marked symmetry isolation of Session 306, but do not give a way to
realize all odd (H): at (t=1) the congruence only sees (H=0mod p),
consistent with the free action. In particular, this theorem supplies
obstructions and an exact quotient observer, not a proof of Chan–Pak.

## Prime targets cannot borrow any symmetry

Suppose H(T)=p for a prime p≥5. If the automorphism group were nontrivial,
freeness on Hamiltonian paths would force its order to be p. Its generator
has a moved orbit O of size p. The restriction to O is a circulant on
the additive group of order p, with (p−1)/2 permitted increments.
For each permitted increment, the constant-increment orders starting
at each of the p vertices give p distinct Hamiltonian paths. Hence
H(T[O])≥p(p−1)/2>p. Ordinary tournament path insertion extends each of
these paths to a Hamiltonian path of T retaining its restriction to O;
distinct orders on O cannot coalesce. Contradiction.

In particular, constructing witnesses for arbitrary large primes in
Chan–Pak necessarily uses asymmetric tournaments; growing a permanent
global group cannot handle even this sparse subproblem. This does not
prohibit symmetric temporary carriers or a construction that breaks
symmetry at its final step. For p=3 there are symmetric examples:
a cyclic triangle with arbitrary transitive source/sink padding.

More generally each prime p≥5 in the order of the automorphism group forces
p to divide H and H≥p(p−1)/2, by Cauchy's theorem and the same orbit proof.
If the support of the mod-p quotient weights has several connected
components, its cover polynomial is a product across those components.
Reversal cancellation forces each odd-order component to contribute at
least one path and each even-order component at least two. Thus the
multiplicity at every field residue is at least the number of odd
components plus twice the number of even components. Precisely, for an
arbitrary fixed set F it is the multiplicity for Q_F at that residue plus
the vanishing order of W_D at zero; the two polynomials are monic and
t^p−t has a simple root at every residue.

## Prime-power symmetry: exact independent length tiers

Now let g have order p^a, with F its fixed set. For e=1,…,a let D_e
be the same phase-weighted quotient on the orbits of exact size p^e
(using an empty quotient when no such orbit occurs). The exact law is

    Q_T(t) ≡ Q_F(t) ∏_{e=1}^a W_{D_e}((t^p−t)^(p^(e−1))) mod p.

Indeed only fixed covers survive reduction modulo p: all other cover
orbits have size a positive power of p. If a path in a fixed cover has
orbit size p^e, then g^(p^e) fixes all its vertices by directed-order
rigidity. No vertex can have shorter orbit, since then two distinct
translates of the path would share it. No two vertices of the path
can be in the same vertex orbit, by the same overlap argument. Thus
each translate orbit lies inside a single length tier, and the
quotient path choices in each tier lift independently by the same
phase argument. A quotient cover with j_e paths contributes p^e j_e
paths, so the falling factorial breaks into exactly
Σ_e p^(e−1)j_e residue blocks. This proves the formula.

The exact multiplicity at any residue r is consequently

    ord_r(Q_F) + Σ_e p^(e−1) ord_0(W_{D_e}).

This makes a size change precise: a visible symmetry orbit of length
p^e is not just e independent prime actions. Its entire tier has a
Frobenius-scaled weight p^(e−1), and coupling between tiers disappears
in this particular characteristic-p observer. It may still matter
to the integer Hamilton count.
