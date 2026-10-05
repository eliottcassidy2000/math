# Small structures that can be reused recursively

2026-10-04. **PROVED:** the elementary identities and scoped compression
mechanisms below, with proofs in the three linked packages. **FINITE-EXACT:**
their explicitly bounded controls. **OPEN:** completion for every positive
Collatz source. This note synthesizes reusable constructions; it makes no
literature-priority claim and no claim that numerical resemblance identifies
two mathematical operations.

## Inheritance and the live board

The closest proved mechanisms are the
[word-doubling trace law](quadratic_escape_rank_atlas_20261004.md#6-the-first-and-third-quadratics-act-on-repeated-word-presentations),
[ordered affine carry decoder](geometry_collatz_drift_carriers_20261004.md#3-same-axis-means-one-primitive-pattern-different-axes-retain-order),
and [155 mod 2048 common-future dependency](collatz_checkpoint_reroute_20261004.md).
The tournament side inherits
[THM-291, Mode B](../../01-canon/theorems/THM-291-mode-b-multilinear-recursion.md)
and [THM-447, skew-Sylvester doubling](../../01-canon/theorems/THM-447-skew-sylvester-doubling.md).
The current incoming paid-controller work keeps the original source rank
through anchor changes; that is relevant evidence for keeping a proof's
boundary information through recursion.

Canonical hostiles are valuation words 12 and 21, whose identical slopes
hide different carries, and two strongly connected four-vertex presentations
whose same named edge flip produces different classes. The corrected near
miss is treating a square law for an observable as a square law for the whole
object. The least-used useful sidecars are the marked affine anchor, the
tournament presentation gauge, and the supplied terminal root certificate.

The live board has six concepts: **ordered pair; diagonal/mixed split;
composition; quotient; retained boundary data; terminal certificate.**
Anchor: usable Collatz proof components. Niche: a small tournament memory
alphabet. Wildcard: the commuting/anticommuting versions of the same
four-product expansion.

## 1. One square, three different uses

For J(z) = z + z^(-1), arrange the four products as a square:

| | z | z^(-1) |
|---|---:|---:|
| z | z^2 | 1 |
| z^(-1) | 1 | z^(-2) |

The diagonal sum is J(z^2); the mixed sum is 2. Therefore

    J(z)^2 = J(z^2) + 2.

This explains exactly why word doubling appears as lambda -> lambda^2 in
one coordinate and J -> J^2 - 2 in another. The two mixed positions are the
source of the correction. For unequal inputs the same construction gives

    J(a)J(b) = J(ab) + J(a/b).

The relative coordinate replaces the constant correction; for positive a,b
it equals 2 precisely when a=b. On the unit circle this is the elementary
cosine product identity: the same four products split into sum and difference
angles. A coefficient-2 correction and an operation that adds two vertices
still have different types.

The general two-by-two trace identity is

    (tr M)^2 = tr(M^2) + 2 det(M).

It follows directly from expanding a two-by-two matrix. The representation
split is three-dimensional symmetric square plus one-dimensional exterior
square, not a global decomposition into a matrix square and two trivial
representations. The distinction matters as soon as factors fail to commute.

There is also a precise connection to the earlier three quadratic portraits.
Retain the determinant together with the trace: matrix squaring induces

    (t,d) -> (t^2-2d,d^2).

At d=0, 1/2, and 1 the first coordinate is respectively t^2, t^2-1, and
t^2-2. The determinant fibers d=0 and d=1 are invariant. The middle fiber
is not: d=1/2 becomes 1/4, so the next correction is 1/2. For the rational
matrix M=[[0,-1/2],[1,0]], the actual traces are 0,-1,1/2, whereas repeated
t->t^2-1 would give 0,-1,0. This is a minimal explicit test of the missing
coordinate. Keeping (t,d) gives a valid recursive summary for every rational
two-by-two matrix; it still does not recover the matrix's marked carry.

The user's earlier base-phi identities fit this same cell exactly. Let
phi=(1+sqrt(5))/2 and C=[[1,1],[1,0]]. Its eigenvalues are phi and -1/phi,
so det(C)=-1. Direct multiplication gives

    C^4=[[5,3],[3,2]],   trace=7,  determinant=1,
    C^5=[[8,5],[5,3]],   trace=11, determinant=-1.

Their eigenvalue sums give phi^4+phi^(-4)=7 and
phi^5-phi^(-5)=11. These are precisely phi^8+1=7phi^4 and
phi^10=1+11phi^5. Squaring the same two small carriers now gives

    trace(C^8)=7^2-2=47,
    trace(C^10)=11^2+2=123.

The two mixed products each equal the determinant, so their sign explains
the minus-two/plus-two split. This transports an exact algebraic mechanism
between the phi examples and the word-matrix calculation. Primality of 7
or 11 is unnecessary for the identity, and the transport supplies no integer
Collatz legality condition by itself.

The same four-product expansion has three reusable outcomes:

1. **Reciprocal scalars:** the mixed products equal 1 and add to 2.
2. **Skew-Sylvester matrices:** the two mixed products cancel.
3. **Affine maps with different anchors:** the two orders retain a translation
   difference; taking a character discards it.

For the second case, set

    A = [[1,1],[1,-1]],  B = [[0,1],[-1,0]].
    A^2 = 2I, B^2 = -I, AB + BA = 0.
    D(M) = A tensor M + B tensor I.

Expanding four terms proves the inherited tournament identity

    D(M)^2 = I_2 tensor (2M^2-I).

If i*mu is an eigenvalue of the skew tournament matrix M, its squared real
magnitude z=mu^2 evolves by z -> 2z+1, so z+1 doubles. Iteration gives

    D^d(M)^2 = I_(2^d) tensor [2^d M^2-(2^d-1)I].

This is an exact bridge back to the earlier 2x+1 theme, including 63 from
z=0 after six iterations. It does not depend on interpreting 63's factors
as a convergence certificate. Here the recursion is matrix block algebra.

For the third case, if f(x)=a*x+b and g(x)=c*x+d, then

    f(g(x))-g(f(x)) = (a-1)d-(c-1)b.

The common slope ac alone cannot recover this difference. These letters f,g
are explicitly defined affine maps in this paragraph, not an interpretation
of the earlier incompletely parenthesized operation quartet.

Proofs, equality boundaries, hostiles, and exact controls:
[four-channel algebra](reciprocal_four_channel_kernel_20261004.md) and
[tournament recursion](tournament_recursive_four_state_20261004.md).

## 2. Name the coordinate that is changing

| Construction | Size coordinate | Other exact update |
|---|---|---|
| Source -> T -> sink | n -> n+2 | Hamiltonian-path count H stays H |
| Two-copy ordered join T -> T | n -> 2n | H -> H^2 |
| Self-substitution T[T] | n -> n^2 | No general H -> H^2 law |
| Skew-Sylvester D(T) | n -> 2n | z -> 2z+1 for squared spectral magnitude |
| Repeated valuation word ww | word length r -> 2r | lambda -> lambda^2, J -> J^2-2 |

Thus the user's square/add-two connection has an exact local mechanism, and
several distinct global realizations. In the ordered join, every Hamiltonian
path uses the first block completely before the second, giving the product
law. Source/sink enclosure forces only one new endpoint on each side and
preserves the interior choice of path. The operations need not have the same
update on every observable.

The Mode B audit also repairs a printed boundary count in THM-291: each arm
has n-1 free tiles, with one apex, totaling 2n-1. The original path bijection
and restriction identity are unchanged; the correction lineage is in
[MISTAKES.md](../../01-canon/MISTAKES.md).

## 3. A four-state description versus a four-channel processor

A scalar character loses carry. The full four-dimensional tensor operator
does not: for a positive marked carrier M=[[P,B],[0,Q]], entries of M tensor M
include P^2, PB, and Q^2, which recover P,B,Q. Further,

    R(MN)=R(M)R(N),   R(M)=M tensor M.

This gives a genuine compositional four-channel processor. Repeating the
same word squares the same four-by-four operator; its dimension stays four.
Its integer entries grow, so this is a fixed number of coordinates, not a
fixed number of bits. Taking only its trace crosses the information-loss
boundary. For positive valuation words, an alternative small coordinate
set is the typed pair (J,alpha), where alpha=B/(Q-P); this also reconstructs
the ordered word and its arithmetic source guard.

Four-vertex tournaments give a different useful alphabet: the transitive
class, a cycle plus sink, a source plus cycle, and the strongly connected
class. Their H values are 1,3,3,5. An ordered join of four-vertex blocks
faithfully encodes a word over this alphabet: the intrinsic ordered strongly
connected components give cuts after each four vertices. Repetition can be
stored as a shared subtree. This is an exact graph encoding of an ordered
word, not a theorem that Hamiltonicity implies Collatz convergence.

The presentation sidecar is essential for predicting arbitrary tile flips.
With the fixed path 0->1->2->3, let a,b,c reverse the remaining edges
(0,2),(1,3),(0,3). Presentations ab and c are both strongly connected.
Flipping a gives b, which is source plus cycle, and ac, which stays strongly
connected. A single isomorphism-class label has discarded the necessary
edge address.

More generally, in the fixed-path cube of any order n, any quotient that
preserves class and all named flips must distinguish every presentation.
Indeed, translating two distinct bit vectors u,v by u sends the first to
the unique transitive presentation 0 and the second to nonzero u XOR v.
The two classes differ. A chosen four-state halfsection does support its
two specified flips; the full family of flips demands the full gauge.

The reusable rule is elementary but strong. If q is a summary and every
allowed operation F has a well-defined operation Fbar satisfying

    q(F(x)) = Fbar(q(x)),

then induction transports every composition to the summary. If two states
have the same summary but different next summaries, either add the missing
coordinate or restrict the allowed operations. This identifies precisely
when recursive reuse is justified.

## 4. A completed Collatz family built from a reusable component

The earlier dependency is

    H(n)=(729n+669)/1024,  n=155 mod 2048.

It is a smaller common-future dependency, not itself an actual forward
Collatz valuation word. Write

    v=(1,2,1,1,1,2),   W=(3,2,1,1,1,2).

For the affine maps of these words,

    F_v(n)=4H(n)+1,       F_W F_v = F_v H.

The second identity is the recursive transport rule. The exact repeat
condition for m applications is

    v_2(295n-669) >= 10m+1.

The repeat count, this guard, and the terminal source c=H^m(n) are the
boundary data. If c has a supplied first-hit word (b,tail), the source has
the actual first-hit word

    v W^(m-1) (b+2) tail.

The repeat component adds 6m odd edges and 10m halvings. A typed certificate
verifier checks the closed endpoint and terminal proof without materializing
the repeated prefix. A separate expansion reconstructs the inherited root
certificate and checks the literal route. Fuel without the terminal proof
is only a prefix certificate.

There are completed instances at every repeat depth m>=1. Choose b=2s with

    4^s = 2302 * 295^(-1) mod 3^(6m+1),
    c=(2^b-1)/3 > 669/295,
    n=[669+1024^m(295c-669)/729^m]/295.

The target is 1 mod 3; powers of 4 run through that entire principal-unit
group, so suitable s exist and form an arithmetic progression. The congruence
makes 729^m divide 295c-669, and 1024=729 mod 295 proves n is an integer.
The source is positive odd, n>c>1, and its binary repeat guard holds. Since
c goes directly to 1 with valuation b, this constructs a completed source
with word v W^(m-1) (b+2) for every m. The note gives the first-hit audit and
the elementary principal-unit proof. In fact these sources have exactly m
available repetitions: b>=4 gives v_2(295c-669)=1, hence initial fuel 10m+1.

The first two chosen phases have b=888 and b=786750. Literal integer routes
are checked for those two, including the 786750-bit second source. Higher
depth plans have exact symbolic and modular readers; their compact
descriptions must not be confused with small expanded integers.

See [the dependency kernel proof and program](collatz_recursive_dependency_kernel_20261004.md).
This is an all-depth completed family. It does not assert that an arbitrary
supplied integer has this form, or that a terminal certificate can always be
found by a terminating search.

### A concurrent connection: powers of three enter every repeat depth

Incoming [HYP-9175, powers-of-three shadows](../hypotheses/HYP-9175-powers-of-three-shadow-only-cycle-points-1-or-3-mod-8-resisting-share-constant.md)
separates its proved modular group fact from conjectures about coverage. The
proved fact supplies a useful new test here. The kernel anchor 669/295 is
3 mod 8. For k>=3, the powers of 3 modulo 2^k are exactly the units that are
1 or 3 mod 8: the order is 2^(k-2), as the elementary valuation formula
v_2(3^(2^j)-1)=j+2 shows for j>=1.

Consequently, for every m>=1 there is exactly one exponent class

    e = e_m mod 2^(10m-1),
    3^e = 669 * 295^(-1) mod 2^(10m+1),

whose positive powers-of-three sources admit at least m kernel repetitions.
The first class is the inherited e=483 mod 512; higher classes refine it.
The program lifts these exponent digits using modular powers without
expanding the enormous source integers. This proves arbitrary repeat-depth
applicability inside a previously hostile source family. It does not prove
that the resulting terminal children have root certificates. The constructed
completed family above and these power-of-three guard families have separate
quantifiers; their intersection has not been established.

## 5. What the connection preserves, and the next decisive test

| Source -> target | Preserved predicate | Lost information / retained repair | Cheapest hostile |
|---|---|---|---|
| Reciprocal pair -> trace | Exact square identity | Carry; retain alpha or full operator | Words 12 versus 21 |
| Marked word matrix -> full tensor operator | Composition and exact marked carrier | No carrier loss on the positive type; bit sizes grow | Recover P,B,Q from entries |
| Tournament presentation -> isomorphism class | Unlabeled tournament properties | Named edge; retain gauge or a closed halfsection | ab,c followed by flip a |
| Four-symbol word -> ordered graph join | Symbol order and decoding | H alone loses order; retain graph | Exchange the two diamonds |
| Repeated dependency -> proof node | Common future and exact prefix legality | Terminal home status; retain terminal first-hit proof | Correct fuel with wrong or absent terminal |

The next useful target is a library of these typed repeat components and a
selector that seeks a proved terminal family while preserving exact source
identity. Every proposed new component should supply an affine identity,
its integer legality cylinder, a rank or bounded repeat fuel, and a terminal
interface. Compose two different components only after testing the interface
guard on the actual output; matching scalar slopes is insufficient.

The current advance is reliable reuse across arbitrarily many repetitions.
The remaining Collatz problem is coverage by completed components, not the
algebraic ability to repeat a component once its hypotheses are satisfied.

## 6. Audit and reproduction

The three packages received independent proof, type, and hostile-example
audits. The parent session separately reran each program with Python's normal
and optimized interpreters and compared both outputs with its saved output;
all three matched. The new determinant-coordinate identity also received an
independent rational-matrix check, and the first power-of-three guard phase
was independently enumerated over all 512 exponent classes.

Run the three scripts linked by the package notes with `python -X utf8 -B`,
then with `python -X utf8 -B -O`. The LF-normalized output SHA256 values are:

| Package | SHA256 |
|---|---|
| reciprocal_four_channel_kernel_20261004 | f41e95fb52dda817bbf21d07fc158480f69b1df668f0e088c616340fc22e2a66 |
| tournament_recursive_four_state_20261004 | 19c998989f099302967c57c7fcb98e2945b35a4856b192eb71ee40b7994f9bc9 |
| collatz_recursive_dependency_kernel_20261004 | e980329e5c49fd7ec7b51f13fd744891dd7b14b279aaa906e85d3688868a5773 |

The exact finite universes and rejected controls are specified in each note.
These are mathematical proofs with reproducible integer/rational experiments,
not a Lean formalization claim.
