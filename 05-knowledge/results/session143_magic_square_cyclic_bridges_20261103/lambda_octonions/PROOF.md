# Signed cubes, fixed subalgebras, and a square-magic half-turn lift

**Session 143; 2026-11-03. Status: PROVED elementary identities and explicitly
bounded exact checks.** The intended meaning of the user's “lambda cube” is
**UNRESOLVED** in the inspected repository evidence. The conditional type-theory
comparison below does not resolve that provenance question. The main algebraic
mechanisms are classical, not claims of new octonion theory. No LRC, magic-square
existence, tournament-spectrum, or sporadic-group conjecture is settled here.

## 1. Inheritance and connection board

The closest proved mechanism is the repository's signed XOR construction and
Hamming-code gauge kernel: `arithmetic_braids2_20260917_squarefree_symmetry.md`
§5 and `arithmetic_braids2_20260917_fano_code.md` §§1–2. Its exact sidecars are the
basis labels, product-output XOR, identity, diagonal squares and sign cochain.
The Paley/Fano note §S13 supplies an explicit relabelling from Z/7 to this table.
The canonical hostile is an independently oriented Fano plane: only 16 of its
128 orientations are the inherited alternative octonion tables. The corrected
near miss is “eight vertices imply octonions”; no multiplication follows from
that count. The least-used sidecar is the eight-character subgroup fixing the
multiplication table while changing signs of its named basis vectors.

The live board has five objects: Boolean inclusion cube; binary squareclass
space; signed octonion product; the nine-cycle; and the perimeter-root space of
a 3×3 magic-square problem. We retain exact maps and state their losses below.

## 2. The missing octonion coordinate is a sign cochain

Let V=F₂³ and label its elements by integers 0,…,7, with the first bit least
significant. Let e₀=1. Use the Cayley–Dickson convention

    (a,b)(c,d) = (ac − d b*, a* d + cb).

Its multiplication is

    e_x e_y = α(x,y)e_(x+y),  α(x,y)=(-1)^f(x,y),
    f(x,y)= Σ_(i≤j) x_i y_j
             + x₂x₃y₁ + x₁x₃y₂ + x₁x₂y₃.                 (1)

Coordinates and exponents in (1) are in F₂. Recursive doubling gives (1);
`verify.py` independently compares the 64 signed products from recursion with
this polynomial. The same classical cochain construction appears in
[Albuquerque–Majid, §4, Corollary 4.3 and Proposition 4.4](https://arxiv.org/html/math/9802116v1).
The formula above is checked from the declared recursion rather than inferred
from a displayed formula or an unrecorded basis convention. The Fano/cube picture
and the distinction between an ordinary cocycle and the octonion twist are
also described by [Baez, §2.1](https://math.ucr.edu/home/baez/octonions/node4.html).

**Proposition 1 (exact obstruction).** For all a,b,c∈V,

    (e_a e_b)e_c = (-1)^det(a,b,c) e_a(e_b e_c).           (2)

**Proof.** The exponent of the ratio is

    δf(a,b,c)=f(a,b)+f(a+b,c)+f(b,c)+f(a,b+c).

The bilinear summand of f cancels. The three cubic summands produce the six
monomials a_i b_j c_k with {i,j,k}={1,2,3}, namely the determinant in
characteristic two. This symbolic polynomial identity proves (2). The checker
has a second path that computes all 512 signed products recursively. ∎

The additive algebra associator is **not** this sign. It is

    [e_a,e_b,e_c]=(e_a e_b)e_c−e_a(e_b e_c),

which is zero at determinant zero and ±2e_(a+b+c) at determinant one. In this
convention

    (e₁e₂)e₄=+e₇,   e₁(e₂e₄)=−e₇,   [e₁,e₂,e₄]=2e₇.    (3)

There are 168 ordered independent triples, or 28 unordered noncollinear
imaginary triples. The remaining 344 ordered triples have zero associator.
Among the 35 triples of distinct imaginary labels, seven are Fano lines.
These counts concern basis triples only, not arbitrary octonion vectors.

A basis-sign change e_x↦(-1)^g(x)e_x changes f by

    f'(x,y)=f(x,y)+g(x)+g(y)+g(x+y),  g(0)=0.

The added terms cancel from δf. Thus no such gauge can remove (2). The gauges
that fix the table are exactly the eight linear characters g(x)=u·x; this is
the inherited 128/8=16 orientation count, not a new classification claim.

## 3. A real cube correspondence: subalgebras, not basis elements

For a linear subspace W≤V, set A_W=span_R{e_x:x∈W}.

**Proposition 2.** A_W is a unital subalgebra, and it is associative iff
`dim_F2 W≤2`. Its real dimension is 2^dim W. Inclusion and intersection of W's
are preserved, and algebraic generation of two A_W's corresponds to W+W'.

**Proof.** Product labels add and signs never vanish, so multiplication closes
exactly on addition-closed support. Every product of supported basis elements
can generate its sum label; hence the assertion about joins. If dim W≤2, every
triple determinant is zero, so (2) and trilinearity imply associativity. If W
has dimension three, choose an independent triple and apply (2). ∎

The four levels are R; seven labelled copies of C; seven labelled copies of H;
and O. For a two-dimensional W, choose independent a,b, put i=e_a, j=e_b and
k=e_a e_b: the associative table has i²=j²=k²=−1 and ij=−ji=k.

Fix any independent labels b₁,b₂,b₃. Then

    S ⊆ {1,2,3}  ↦  A_(span{b_i:i∈S})                 (4)

is an inclusion-lattice isomorphism onto eight chosen coordinate subalgebras,
with dimensions 1,2,4,8 by rank. There are 28 such chosen Boolean sublattices:
the 7·6·4 ordered bases, divided by their 3! reorderings. A translated cube
face not containing zero need not be a unital subalgebra. The full lattice of
binary-supported subalgebras has 16 vertices, not eight.

**Conditional λ-cube interpretation.** If the intended object is Barendregt's
cube, the three optional dependency rules index subsets S; its arrows represent
inclusions of systems, as in [Barendregt, *Introduction to generalized type
systems*, §2](https://homepages.inf.ed.ac.uk/wadler/papers/barendregt/pure-type-systems.pdf).
Map those corners by (4), after explicitly choosing a correspondence between
the three rules and b₁,b₂,b₃. It preserves inclusion, meet and join **as the
labelled Boolean poset**. It does not map terms, typing judgements, proofs,
normalization, or expressive power. Rule union is idempotent, unlike XOR:
`S∪S=S`, while `x+x=0`. The basis-element map alone cannot preserve that
operation. This is a transparent structural comparison, not a type-theory
transport theorem.

## 4. A stronger squareclass/Galois correspondence, and its exact loss

Let K=Q(√d₁,√d₂,√d₃), for independent rational squareclasses d_i. Write
r_x=∏√d_i^x_i. These eight elements form a Q-basis and

    r_a r_b = q(a,b) r_(a+b),  q(a,b)=∏ d_i^(a_i b_i).     (5)

The function q is symmetric and satisfies

    q(a,b)q(a+b,c)=q(b,c)q(a,b+c).                         (6)

For each bit, the integer exponent on both sides is a_i b_i+a_i c_i+b_i c_i
−2a_i b_i c_i. This proves (6) without computation.

Let O_Q=span_Q{e_x}, the rational octonion algebra inside real O. The Galois
sign action σ_u(r_x)=(-1)^(u·x)r_x transports under the Q-linear isomorphism
K→O_Q, r_x↦e_x, to the explicit rational octonion automorphism

    D_u(e_x)=(-1)^(u·x)e_x.                              (7)

Indeed characters multiply, so D_u(e_a e_b)=D_u(e_a)D_u(e_b), and D_uD_v=D_(u+v).
Extending D_u R-linearly gives an automorphism of real O as well. Thus this is
an actual common action, not a group-order coincidence. For each subgroup H≤V
the fixed supports in K and O_Q are H^⊥. Consequently the sixteen intermediate
fields and sixteen rational fixed octonion subalgebras have matching inclusion
lattices and Q-dimensions. Their real octonion spans have the same R-dimensions.

But the linear map is not an algebra map: (6) has trivial associator, whereas
(2) does not. For d=(2,3,5), even r₁²=2 differs from e₁²=−1. No Q-algebra
isomorphism K→O_Q can exist, independently of this chosen map, because
associativity is preserved by algebra isomorphisms. The same obstruction
separates K⊗_Q R from real O. Passing only
to binary labels loses precisely coefficients q/α, including this obstruction.
A multiquadratic Galois cube therefore supports a genuine character-action and
subobject-lattice bridge, not an identification of algebras.

## 5. The literal nine-cycle cannot be the eight-label cube

The sequence

    4,8,3,7,2,6,10,5,9,4

is translation by +4 in Z/9, using representatives {2,…,10}; it visits every
residue once. Every group homomorphism C₉→V is trivial: the image h obeys
9h=0, but 9h=h in V. More generally, if a C₉ action is by permutations of
eight labels, each cycle length divides nine and is at most eight, hence is
one or three. An equivariant image of the transitive nine-cycle has at most
three points. A faithful nine-cycle-to-eight-basis-label transport is therefore
impossible. This does **not** rule out order-nine continuous linear actions on
an eight-dimensional space; that is a different target and needs its own map.

## 6. Positive synthesis: an octonionic lift of a square-magic half-turn

Index grid cells row-major by 0,…,8. The center is cell 4. Index the eight
perimeter roots clockwise by

    cell(a): 0,1,2,5,8,7,6,3,  for a=0,…,7.

Thus opposite cells have labels a and a⊕4. Encode their roots as
x=Σ x_a e_a, retain the center root t separately, and let

    J(x,t)=(e₄x,t).                                      (8)

Left multiplication by e₄ is the signed permutation

    e_a ↦ e_(a⊕4) times signs (+,+,+,+,−,−,−,−).           (9)

Its square is −I on the eight perimeter coordinates, either directly from
(9), or by alternativity e₄(e₄x)=(e₄²)x=−x. Consequently J has order four on
the nine-root space, with J²(x,t)=(−x,t). It preserves Σx_a²+t².

Let sq send the nine roots to their coordinate squares and let ρ be 180-degree
grid rotation. Then

    sq ∘ J = ρ ∘ sq.                                    (10)

Each output squared perimeter coordinate equals the opposite input squared
coordinate; the central square is fixed. Rotation permutes the three rows,
three columns and two diagonals. Hence (10) preserves the **actual eight-line
magic-equality predicate**, the magic sum, total square norm, rationality,
integrality, nonzero entries, and distinctness of the squared entries. This is
an explicit C₄→C₂ lift of a familiar symmetry of the magic-square equations.

**Failure boundary and non-novelty.** J is a linear norm-isometry, not an algebra
automorphism: it sends octonion 1 to e₄, not 1. It can change root signs, so a
positive-root convention is not preserved before taking squares. It merely
rotates the same nine squared values and cannot produce a new multiset or an
existence proof. Every paired-coordinate involution admits some orthogonal
signed lift with square −I; octonions give this specified lift, not a unique
magic-square mechanism.

**Norm/opposite-pair hostile.** The nine distinct positive integer roots

    13 23 35
    47 65 79
    85 89 91

have every opposite squared pair summing to 2·65²=8450, and their perimeter
square norm is 8·65²=33800. The center row, center column and two diagonals
all sum to 3·65²=12675. But the outer row sums are 1923 and 23427, and the
outer column sums are 9603 and 15747. Thus even the correct octonion norm
plus all opposite-pair equations is insufficient for magic. This statement
concerns the displayed arrangement, not an uncomputed census of rearrangements.

Nor do all eight XOR translations preserve the magic equations. In this
labelling only 0 and 4 do. As a single explicit hostile, the positive real
roots of the ordinary Lo Shu square have square array

    8 1 6       XOR-1 sends it to       1 8 7
    3 5 7                              4 5 6
    4 9 2                              3 2 9.

The new row sums are 16,15,14, not constant. This uses real square roots and
is **not** a rational/integer magic-square-of-squares witness. The checker
exhibits a failure for each of the other six nonzero translations. No F₂³
magic-preserving action, octonion automorphism action, or G₂ action is claimed.

## 7. Repository recovery and corrections

The recorded bounded local scan examined 29,025 locally present canon/results
`.md` and `.out` files below 150,000 bytes each (203,356,314 bytes total). It found
no exact “lambda cube”, “λ cube”, or “Barendregt” phrase. This is a bounded
negative result, not a claim about every file, commit or conversation.

There are distinct nearby meanings:

* `lucas_mirror_gap_spaces_20261028/LUCAS_MIRRORS.md` calls its toggle graph the
  Lucas cube Λ_n. At n=5 it has 11 vertices and 15 edges, not eight vertices.
* `sixth_clock_branches_20261004.py` uses the test-message string “lambda cube”
  for a scalar algebraic identity in its `sextic()` function, not for a cube
  object. This hit alone cannot resolve the user's intended meaning.
* Inherited pair-coverage/eigenvalue uses of lambda are not type-system cubes.

The docstring of `lrc_octonion_associator_fano_macmini_S74.py` says the octonion
associator is ±1 off lines. Equation (3) corrects that normalization for the
ordinary additive associator; ±1 belongs to the multiplicative ratio. Its
sector-cumulant experiment is a different observable. The companion stored
output reports doubling invariance on only 2/20 inner triples, so it does not
supply an octonion-to-LRC predicate transport. Historical claims in snapshots
are provenance, not newly endorsed assertions.

## 8. Reproduction and stopping boundary

Run `python3 verify.py` and `python3 -O verify.py`: both outputs must equal
`certificate.json`. `python3 verify.py --mutate-sign` must fail at the independent
Cayley–Dickson comparison. The checker uses stdlib integer arithmetic and
explicit exceptions, not optimization-removable assertions. It tests the
consequence objects: product signs, associators, subalgebra supports, character
actions, field cocycles, the given nine-cycle, and every magic-line permutation.

The recovered proposition is now precise: shared binary labels can carry an
inclusion lattice and an explicit diagonal character action; their signed
multiplication remembers the determinant associator. The one square-magic
half-turn also has a concrete root-space lift. The unresolved λ-cube referent
and any requested stronger transport require additional definitions/evidence;
repeated dimensions or spectral coincidences are not substitutes.
