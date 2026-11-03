# All-prime boundary divisors and shared-slot geometry

2026-11-03. **PROVED, scoped**, except explicitly labelled finite controls
and open extensions. Same-agent proof review; independent computational
paths in `audit.py`. No external novelty or independent-agent audit claim.
Local labels AP1--AP7 belong to this file, not the global THM namespace.

## Inherited definitions and mechanisms

For a nonempty tournament T on n vertices, P_j(T) counts unordered covers
by j nonempty directed paths; C_j=j!P_j counts ordered covers. Empty paths
are permitted in the r labelled slots counted by

    F_T(r)=sum_j C_j(T) binom(r,j).

Write H(T)=C_1 and c(T)=gcd_j C_j. The fixed-divisor identity

    c(T)=gcd_{r=0,...,n} F_T(r)=gcd_{r>=0} F_T(r)

and its preservation by uniform module substitution are inherited from
G3--G4 in
`../coverage_machine_induction_20261103/gluing/NOTE.md`. The binomial
coefficient/sample transformations are inverse integral triangular
matrices. In particular c(T)|H(T), c(T)|n!, and c(T) is odd.

The identities F_(A ordinal-sum B)=F_A F_B and
F_T(-r)=(-1)^n F_T(r) are inherited from
`../session159_constructions_20261103/cycle_profiles/NOTE.md`.
That note derives reciprocity both from odd-cycle packings and from
reversing vertex permutations. These are integer power-basis polynomials.

The fixed divisor c(T) is not the potentially larger gcd of repeated
cyclic-wrapper Hamilton counts: the latter uses (r-1)!F_T(r). We never
substitute one for the other. Equal F does not imply equal named-port
responses; the two five-vertex seeds in TM1--TM3 of
`../digit_memory_all_bases_20261103/TOURNAMENT.md` remain hostile controls.

## AP1. A simultaneous prime budget for arbitrary old uniform modules

Let an old tournament contain M pairwise disjoint nonempty vertex sets
T_1,...,T_M, each a module relative to all old vertices outside it. The
remaining old vertices and the orientations between modules are arbitrary.
Add q vertices, with arbitrary edges to every old and new vertex, retaining
every old edge. Put c_i=c(T_i), t=min(2q,M), k=M-t, and

    D_t = gcd_{U subset [M], |U|=k} product_{i in U} c_i.       (AP1)

The empty product is 1. Then D_t divides the new raw Hamilton-path count.
For every prime p simultaneously,

    v_p(D_t) = sum of the k smallest members of (v_p(c_i))_i.

**Proof.** Fix the sequence of maximal old-module runs, individual other
old vertices, and new vertices of a Hamilton path. At most 2q module
occurrences touch a new vertex: each new vertex has at most one predecessor
and one successor. Hence at most 2q distinct modules are touched. For an
untouched module with j occurrences, every ordered j-path cover fills its
runs legally, independently of other untouched modules. Its external edges
are uniform; its internal factor is exactly C_j(T_i), divisible by c_i.
The choices in touched modules and other old vertices give another integer
factor, possibly zero. There are at least k untouched modules, so every
such contribution is divisible by a product from (AP1), and hence by D_t.
Sum the contributions. For the valuation formula, minimizing a sum over
k-subsets selects the k smallest valuations. QED.

This generalizes the common-content-3 case of TM2, whose all-size sharpness
construction is already proved. We do **not** claim that every heterogeneous
layout attains D_t, nor that a single subset minimizes the valuations for
all primes. D_t is exactly the gcd of the available products, rather than
an assertion that one product divides the entire sum.

**Controller consequence.** If M old modules all have c_i divisible by 3,
then any old-edge-preserving attachment producing 3H+2 or 3H+4 must have
q>=ceil(M/2). This is a restriction on this attachment menu. Native
complements changing old edges can escape its hypotheses. It is not an
impossibility theorem for a return grammar with such complements.

## AP2. Exterior powers, K4 and ideals without unique factorization

Let A=diag(c_1,...,c_M). Its kth determinantal divisor, the gcd of its
k-by-k minors, equals D_t. Equivalently it is the gcd of the matrix entries
of the genuine linear map exterior^k A. If the Smith factors of A are
s_1|...|s_M, then D_t=s_1...s_k.

This follows directly: the nonzero minors of a diagonal matrix are the
products indexed by k-subsets, and the same products are the diagonal
entries in the exterior basis. No determinant sign or real spectral
approximation is used.

For M=4,q=1,k=2 there are exactly six coordinates e_i wedge e_j, identified
with K4's six edges. This is an exact index/linear-map correspondence.
It does not turn an arbitrary four-module source into a four-vertex
tournament or recover its forgotten attachments.

The actual modules used in `check.py` have contents

    (3,15,63,55), with Smith factors (1,3,15,3465).

The six pair products are 45,189,165,945,825,3465; their gcd is 3. These
modules are respectively C3, the regular cyclic tournament of order5,
and the quadratic-residue tournaments of orders7 and11. Literal rows and
all profiles are saved in `receipt.json`. Only the order7 and11 cases use
quadratic-residue orientation; the order5 rule is specified separately.

There is also an ideal version. Give internal module edges weights in a
commutative ring R, inter-module and other old edges weight1, and arbitrary
new-edge weights in R. Let I_i be the ideal generated by all internal
weighted C_j(T_i). The same grouping proves that the weighted Hamilton
sum belongs to

    J_t = sum_{|U|=M-t} product_{i in U} I_i.                  (AP2)

No UFD assumption is needed. The ideal is
Fit_t(direct_sum_i R/I_i): a presentation has one row per module, with
generators of I_i supported in that row, so its (M-t)-minors are precisely
the indicated products (up to sign and zero). This is the standard Fitting
ideal construction and direct-sum rule, not a new general algebra theorem;
see [Stacks, Lemma 15.8.4](https://stacks.math.columbia.edu/tag/07ZA).

Over a Dedekind domain, prime-ideal valuations give the same sorted formula.
Over an arbitrary ring, retain the ideal itself; lists of rational-prime
valuations can lose information. The repository's norm19 unit-direction
example distinguishes different prime ideals above the same rational19
(`../cubic_adelic_20261016e.md`), illustrating why the finer labels matter.

## AP3. The arithmetic phase profile and a sharp nine-vertex emergence

For each prime p retain the function

    phi_(T,p)(x)=v_p(F_T(x)), x in Z_p,

with v_p(0)=infinity. Then

    min_x phi_(T,p)(x) = v_p(c(T)),
    phi_(A ordinal-sum B,p)(x)=phi_(A,p)(x)+phi_(B,p)(x).

To prove the first equality, binom(x,j) is p-adically integral on Z_p
(approximate x by nonnegative integers). Thus all F_T(x) have valuation
at least min_j v_p(C_j). Conversely finite differences express each C_j
integrally in F_T(0),...,F_T(j); one of these finitely many samples attains
the minimum. This is the finite-polynomial case of the usual Mahler norm
identity: max_x |F_T(x)|_p=max_j |C_j|_p. For the general classical theorem
see [Conrad, Mahler expansions, Theorem4.1](https://kconrad.math.uconn.edu/blurbs/gradnumthy/mahlerexpansions.pdf).

Taking only the minimum forgets where it is attained. Therefore

    c(A)c(B) divides c(A ordinal-sum B),

and strict divisibility can occur when the low-valuation phases disagree.
Independent slot counts x,y would instead have gcd_{x,y} F_A(x)F_B(y)
equal to c(A)c(B). The constraint here is equality of the slot counts.

**Explicit emergence.** Take outgoing bitmasks, vertices numbered from0,

    A=(8,1,3,6),             F_A(x)=x^4+4x^2;
    B=(24,1,3,6,14),         F_B(x)=x^5+8x^3+2x.

Both have c=1: (F_A(1),F_A(2))=(5,32) and
(F_B(1),F_B(2))=(11,100). Modulo5 the first polynomial vanishes at
0,1,4, and the second at0,2,3. Consequently their product vanishes on
all five residues. Its fixed divisor is **exactly5**: it divides H=55
and9!, ruling out11, and the finite-difference calculation confirms5.

This is the smallest total order of two nonempty tournaments for which a
prime divides the joined fixed divisor but divides neither input divisor.
Here is a finite-case proof of minimality, independently exhaustively
checked by permutations and by subset covers through order5:

* The only possible primes at total order<=8 are2,3,5,7, since c|n!.
* The prime2 cannot divide any c, since H is odd.
* Modulo3, reciprocity and F(0)=0 imply 3|c iff 3|H: the only nonzero
  residues are1,-1. Thus3 cannot first arise in a product.
* For5 to arise, one factor must have 5|H and therefore order>=4.
  If the other factor has order<=2, its F is x or x^2, nonzero on all
  nonzero residues; it cannot supply a missing zero. The only remaining
  splits at total<=8 are4+3,4+4,5+3 (and reversals).
  Every tournament through order4 has F(2) nonzero mod5. At order5 the
  only profiles with H divisible by5 have (H,F(2))=(5,64),(15,108),
  (15,120); the last already has c=15 and is excluded by the hypothesis.
  Thus the remaining products have F(2) nonzero mod5.
* For7 to arise, one factor has 7|H and hence order>=6: the exact H set
  through order5 is {1,3,5,9,11,13,15}. The other factor has order<=2
  and cannot supply a nonzero-residue zero.

The enumeration is a complete finite proof component (1,099 labelled
tournaments), not an extrapolation about larger orders. The minimality
claim concerns prime **support emergence**, not every possible increase
of a previously present prime's valuation.

For odd p, reciprocity pairs x and -x, so the first-level phase profile is
indexed by (p-1)/2 nonzero antipodal pairs plus0. Prime3 has one such pair;
prime5 has two and admits the complementary coverage above. At higher
prime powers the zero set alone is inadequate; retain valuations and
residue lifts. A table mod p^e is exact because F has integer coefficients.

## AP4. Exact marked counting after two exceptional vertices

Take M nonempty tournament blocks in ordinal order. Let the first a blocks
be side A and the rest side L, with1<=a<M. Add c with A->L->c->A. Add z
with c->z; z beats only the named vertex0 of each old block and loses to
its other vertices. These are literal graphs, not weighted averages.

Delete c,z from a Hamilton path. Its remaining pieces are three ordered
segments, each visiting any block at most once, in increasing block order.
Each block is therefore specified by a nonempty subset of the three
segments: seven choices. Enumerate both orders of c,z and these7^M masks.
Reject skeletons whose c incidences are illegal or whose z->c edge occurs.

For each block let C_j be its ordered j-cover count and S_j,E_j the counts
whose specified run starts/ends at its root. The root-selector weights are

|Incidence at z|Internal factor|
|---|---|
|Neither adjacent run belongs to this block|C_j|
|Only the incoming run belongs to this block|C_j-E_j|
|The outgoing run belongs to this block|S_j|

If both runs belong to the same block they are distinct runs. Specifying
the root at the outgoing start already excludes it from the incoming end,
so the last row is still correct. Multiply the block factors and sum.

This is a bijection: a skeleton plus internal ordered covers gives one
Hamilton path, and cutting that path recovers both. Thus it is exact for
all M, at most2*7^M skeletons after computing module profiles. It preserves
each actual block and named root. It is not a polynomial-time algorithm
for arbitrary tournaments, and additional arbitrary ports require richer
joint endpoint profiles.

The module profiles can themselves be computed by subset DP. If h_0(S)
counts paths on S starting at0 and P_j(U) unordered covers of U, then

    S_j=(j-1)! sum_{S containing0} h_0(S) P_(j-1)(V\S),

and the analogous reversed-edge formula gives E_j. Only j<=3 is needed.

**FINITE-EXACT.** For the3,5,7,11 modules of AP2 with split a=1, this gives
the28-vertex raw count 1,133,904,592,664,865, using2,739 positive skeletons.
Its valuations at3,5,7,11 are4,1,2,0. The universal AP1 divisor is3;
the additional factors in this graph are not claimed universally. Thirty-
three smaller literal endpoint-DP counts agree, including all eight
three-bit seed words. A separate permutation/cut audit checks every module
through order5, and suffix recursion checks four further reader graphs.

## AP5. All prime views require a common source and a height condition

For nonzero rational q, the finitely supported list (v_p(q))_p and its sign
reconstruct q by prime factorization. Multiplication is coordinatewise
addition. For addition at a **single** prime one also needs the units:

    v_p(a+b)>=min(v_p(a),v_p(b)),

with equality when the valuations differ, but possible extra cancellation
when they agree. The complete global valuation list still determines the
numbers: it permits reconstruction, addition and refactorization. The
point is that addition is not a coordinatewise rule on valuations alone.

Prime powers form rooted p-ary residue trees; every infinite coherent
branch is a p-adic integer. Across primes, CRT combines every finite set
of residue specifications. Ordinary integers lie diagonally in their
product, which contains more points than that diagonal.

**Height criterion.** Let M_j|M_(j+1) tend to infinity and eventually be
divisible by every positive integer; let r_j mod M_j be coherent. For a
fixed H>=0, the following are equivalent:

1. There is an integer n with |n|<=H in every residue class.
2. Every class has some representative in[-H,H].
3. All sufficiently late classes have such a representative.

For M_j>2H the representative is unique. Compatibility then forces all
later representatives to be the same, and coherence supplies the earlier
ones. This proves the equivalence. Hence bounded balanced representatives
characterize the ordinary-integer diagonal. Without a uniform bound,
finite CRT satisfiability does not do so.

**Hostile coherent point.** Take the2-adic coordinate0 and every odd-prime
coordinate1. It is idempotent in the product. Every finite CRT projection
has an integer representative, but a single integer satisfying all of them
would be divisible by every2^k, hence0, while also1 mod3: impossible.
This is not a counterexample to the Hasse principle for x^2=x: global
roots0,1 exist. It is an obstruction to gluing these *selected* local roots
into one prescribed rational/integer source.

For a nonzero rational source the real and p-adic norms satisfy the
classical product formula |q|_infinity product_p |q|_p=1; see
[Conrad, Ostrowski's theorem, Theorem2](https://kconrad.math.uconn.edu/blurbs/gradnumthy/ostrowskiQ.pdf).
That equation alone is not a diagonal-membership certificate. For example,
real coordinate1,2-adic coordinate3,and all other p-adic coordinates1
have norm product1, but cannot be the images of one rational: the real
coordinate forces that rational to be1. Units and compatibility survive
after every absolute value is forgotten.

## AP6. Divisibility tricks in every base are local coordinate operations

The all-base theorem DM1 and minimal-observer theorem DM2 in
`../digit_memory_all_bases_20261103/PROOFS.md` are inherited, not re-proved
as new. For base B>=2 and modulus m, write m=s*u, with s supported on
primes dividing B and gcd(u,B)=1. If L is least with s|B^L, divisibility
by s uses the last L digits; divisibility by u uses invertible weights
B^i mod u. Block sums and alternating block sums occur when B^ell=1
or -1 modulo u. For n=Ba+d,

    n=0 mod u iff a+B^(-1)d=0 mod u.

This separates a finite valuation tail from a periodic unit action. Digital
roots use only modulus B-1; decimal doubling's6-cycle mod9 is one unit
orbit, not a universal numerical geometry. At a prime p, appending a
digit is the affine map x->Bx+d on its residue tree; its valuation obeys
the cancellation rule in AP5. Base is a choice of observation/update
coordinates; the integer and its prime valuations are not base-dependent.

The tournament connection is exact at the level of integer functions:
F_T is an integer polynomial, so the same residue rules evaluate it in
every modulus, while its binomial coefficients retain actual path-cover
counts. It is **not** an assertion that divisibility is intrinsically a
tournament relation. Divisibility has ties/incomparabilities and forms an
exponent lattice (gcd=min, lcm=max); forcing a tournament orientation would
discard that structure.

For a concrete spectrum application, AP1 excludes a whole family of
insufficient return menus, and AP3 identifies phases to retain to detect
new obstructions. Neither proves closure under actual-source returns.
For Collatz, the inherited all-odd family n=q*2^k-1 reaches q*3^k-1 after
k shortcut steps, with ordinary clock2k. With q=1 it approaches -1
2-adically while its real size is unbounded. A proof based only on finite
local views cannot replace the missing source-specific height/descent
estimate with compactness. Convergence remains OPEN.

## AP7. The exact ideal of one-birth responses, and its coverage boundary

This section uses, rather than reclaims, **THM-4104**,
`../../../01-canon/theorems/THM-4104-selected-order-eleven-strong-ear-solid-interval.md`,
Lemma1.1. Its exact cut quadratic is also the observation mechanism of
`../mrb_square_probe_spectrum_20261103/README.md`; that packet already
recovers one pair probe for free and compiles a square-sum probe bank.

For a labelled old tournament retain endpoint counts A_u (ending at u),
B_v (starting at v), and K_uv, the number of ordered two-path covers with
first path ending at u and second starting at v. Let x_v=1 when the new
vertex beats v. Cutting a new Hamilton path at the newborn gives

    H_new(x)=H + sum_v L_v x_v - sum_{u<v}(K_uv+K_vu)x_u x_v,
    L_v=B_v-A_v+sum_u K_uv.                                 (AP7)

It follows that the ideal in Z generated by **all** one-birth counts is

    (gcd(H, (L_v)_v, (K_uv+K_vu)_(u<v))).

Indeed values are integer combinations of coefficients, and finite
differences at the empty, singleton and pair cuts recover all coefficients
integrally. The inherited constant-cut identity H_new(all)=H_new(empty)=H
recovers one missing quadratic coefficient from the others. Thus, for
n>=2, n+binom(n,2) scalar values including the baseline suffice even for
this exact ideal, not merely for a rational response span. This is a
direct corollary of the older integral reconstruction, not a replacement
for its stronger finite-horizon result.

This ideal retains the actual source; it can be much smaller than the
uniform-slot divisor. For the AP3 nine-vertex source H=55, c=5 but its
one-birth response ideal is (1). A literal cut mask3 gives1111, escaping5.
The512 cuts give354 distinct counts, minimum55 and maximum2513. Neither111
nor167 nor169 occurs. These are **FINITE-EXACT** source-specific failures,
not global forbidden Hamilton counts.

In particular, gcd1 does not imply that any one response is1, that every
odd number occurs, or that the two required ternary targets occur. Ideals
allow signed combinations of responses; a physical continuation does not.
Every one-birth count is at least H: each old Hamilton path admits an
insertion position for the new vertex (endpoints or a transition from a
vertex beating it to a vertex it beats); delete the new vertex to recover
the old path, so chosen insertions for distinct paths are distinct.

As a sanity family, the ordinal sum of M copies of C3 has one-birth response
gcd exactly3^max(M-2,0). AP1 gives the lower divisor. For M>=2 the actual
six-vertex source C3 ordinal-sum C3 with new vertex beating local names0,3
(mask9) has H=41. Put this seven-vertex extension in the ordinal sum with
the remaining M-2 triangles. Compare its count41*3^(M-2) with the all-one-
side birth count3^M: their gcd is3^(M-2). For M=1, counts3 and5 suffice.
The construction is not generally strong; stronger connected boundary-
sharpness families are already in the inherited attachment work.

The useful next problem is to retain both this ideal and the **attainable
response set with its next ports**. No divisibility invariant by itself
can supply the missing compatible return.
