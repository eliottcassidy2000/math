# Black deficits, strong centres, and the first Fibonacci centre sector

Status: **PROVED** for the factorization, grammar, capped-word family, and
black-support degree statements below. **FINITE-EXACT** for the complete
order-three through order-nine centre census. **OPEN**: the assertion that
the four-vertex strong tournament is the only strong black-deficit-two
centre at every order; all-order black connectivity and terminal-basin
confluence are not consequences of these counting formulas.

This is the deficit-centre lane of the 20 October continuation. Its closest
proved mechanisms are the SCC centre law in
[THM-2454, pure-blue alphabet and centre law](../../../01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md),
the completed pure-blue/rigid alphabet in
[the current pure-blue packet](../pure_blue_atom_classification_20261020/PURE_BLUE_COUNTS_AND_BRANCH_GEOMETRY.md),
and Theorems 4--6 of
[continuation width and symmetry](../continuation_width_symmetry_20261020/README.md).
The historical open classification sentences in THM-2454 are superseded by
the current packet, but its proved fibre factorization remains the input.
The canonical hostile is `C5 -> C5`: each component is pure-blue, whereas
the stack has `(h,b)=(9,3)` and six black half-edges. The corrected near miss
is treating a power of a centre series as an additive black deficit. The
least-used sidecar is the subset/path convolution `J(U)` below.

Live board: SCC palindromes; black half-edge mass; additive symmetry slack;
subset/path convolution; actual black-support neighbours; complete boundary
path data. The map from a tournament to `(h,b)` preserves black incidence
mass and loses the destinations of those incidences. Hamiltonian words and
their actual complement partners remain necessary for connectivity.

## 1. Three parameters that must not share one name

Write

    h(T) = H(T)/|Aut(T)|,
    b(T) = number of blue tilings in the ordinary class,
    beta(T) = h(T)-b(T).

For self-converse T, `h,b` are positive odd integers, so beta is a
nonnegative even integer. It counts black **half-edges**, not black lines.
For a self-converse node, ordinary and converse-merged fibres coincide.
For a non-self-converse merged node the total fibre is doubled; this note's
centre classification is always self-converse.

Independently, let `f(T)` be minimum feedback arcs, `s(T)` the number of
vertices moved by the whole automorphism group, and

    delta(T) = 3 f(T)-s(T).

This is the symmetry slack of the continuation-width packet. The capped
feedback-core deficiency `2k+1-h_heavy` in
[the strata packet](../defect_transitions_20261020/strata.md) is a third
quantity. In particular the four-vertex strong tournament R4 has

    (h,b,beta)=(5,3,2),       (f,s,delta)=(1,0,3).

The notation `D_2(z)` below means the series of **strong self-converse
centres with beta=2**. No inherited definition identifies it with a transfer
determinant, the capped-core deficiency, or a symmetry-slack series.

## 2. Exact SCC factorization and the first deficit layers

Every self-converse tournament has the unique SCC form

    T = A -> Z -> A^op,                                      (1)

where A is an arbitrary ordered SCC word, and Z is either one self-converse
strong SCC or empty. Set `q=h(A)`, `c=h(Z)`, `d=b(Z)`, with `c=d=1` for an
empty centre. Hamiltonian paths and automorphisms factor over the ordered
SCCs. A blue tiling chooses an arbitrary tiling on each left component,
its prescribed mirror on the right, and a blue central tiling. Thus

    h(T)=q^2 c,       b(T)=q d,
    beta(T)=q(qc-d)=q[(q-1)c+(c-d)].                          (2)

This is an all-order identity, not a counting ansatz. Since q is odd,

    q>1  implies  q>=3 and beta>=q(q-1)c>=6.                 (3)

Consequently every self-converse tournament of black deficit two or four
has a rigid half A. The rigid strong alphabet is exactly `{1,C3}`; put

    A(z)=1/(1-z-z^3).

For odd q, let `W_q(z)` count **ordinary tournament classes** A with
`h(A)=q`. Let `C_(c,d)(z)` count allowed centres with masses `(c,d)`, including
the empty centre in `C_(1,1)`. Let `D_beta(z)` count strong self-converse
centres of positive black deficit beta. Then the self-converse beta-stratum
series satisfies the exact formula

    B_beta(z) = sum_(q,c,d: q^2 c-q d=beta)
                    W_q(z^2) C_(c,d)(z).                    (4)

The q=1 term is `A(z^2) D_beta(z)` for beta>0. For each fixed positive beta,
only finitely many **q>1** correction triples q,c,d are possible: q divides
beta, `q(q-1)<=beta`, and `c<=beta/[q(q-1)]`, with `1<=d<=c`. There is no
such c,d bound in the q=1 term. Thus (4) reduces a deficit layer to its
strong-centre series and finitely many outer-mass corrections.

The first corrections are:

| beta | B_beta(z) beyond A(z^2) D_beta(z) |
|---:|---|
| 2,4,8,10,12,14,16 | 0 |
| 6 | W_3(z^2)(1+z+z^3) |
| 18 | W_3(z^2) z^5 |
| 20 | W_5(z^2)(1+z+z^3) |
| 24 | W_3(z^2) C_(3,1)(z) |
| 30 | 0 |
| 36 | W_3(z^2) C_(5,3)(z) |
| 42 | W_3(z^2) C_(5,1)(z)+W_7(z^2)(1+z+z^3) |

Here `C_(1,1)=1+z+z^3` and `C_(3,3)=z^5` follow from the proved pure-blue
classification. The potential beta30 correction would require a pure-blue
centre of mass five, which does not exist. `C_(3,1)` has not been declared
empty by a finite check. If the global conjecture `D_2=z^4` is proved, it
would imply `C_(3,1)=0` and `C_(5,3)=z^4`.

For beta0, (4) recovers exactly

    B_0(z)=(1+z+z^3+z^5)/(1-z^2-z^6).

For beta2 it gives the useful all-order equivalence

    B_2(z)=D_2(z)/(1-z^2-z^6).                               (5)

Thus `D_2=z^4` would classify all deficit-two nodes as rigid palindromic
padding around R4. Equation (5) is proved; that global identification of
`D_2` is presently open outside the sector proved next.

## 3. A complete infinite strong-centre sector

The symmetry theorem classifies every strong tournament of delta3 as

    C(U): x -> U -> y -> x,                                  (6)

where U is a word in singleton and directed-triangle SCCs, and `|U|>=2`.
Suppose U has r singleton and t triangle letters. Define

    J(U)=sum_(S subseteq V(U)) H(U[S]) H(U[V(U)\S]),

with the empty tournament assigned one Hamiltonian path.

**Proposition 1.** Throughout (6),

    |Aut(C(U))|=3^t,
    h(C(U))=1+2^(r+2t),
    f(C(U))=t+1,  s(C(U))=3t,  delta(C(U))=3.                (7)

**Proof.** Among maximum-score vertices, x is the unique vertex dominating
all the others: its score is |U|, and another vertex has that score only
when it is a source vertex of U. Similarly y is the unique minimum-score
vertex dominated by all the other minimum-score vertices. These descriptions
work for every `|U|>=2`, so every automorphism fixes x,y and restricts to an
arbitrary automorphism of U.

Hamiltonian paths with x preceding y consecutively in cyclic block order
are `x P_U y`, `y x P_U`, and `P_U y x`. Every remaining path has the unique
form `P_S y x P_(U\S)` for a nonempty proper S. Therefore

    H(C(U))=3H(U)+J(U)-2H(U)=H(U)+J(U).                      (8)

Both H and J multiply over ordered sums. Directly `H(1)=1,J(1)=2` and
`H(C3)=3,J(C3)=12`. Thus `H(U)=3^t`, `J(U)=2^r 12^t`, proving the mass
formula. An optimal order for U, placed between x and y, uses its t
feedback arcs and yx. Conversely each triangle of U requires an internal
feedback arc and the cycle `x,u,y` requires at least one arc outside U.
The automorphism statement gives the moved support. This proves (7).

**Proposition 2.** C(U) is self-converse exactly when the atom word U is
palindromic. For such U,

    b(C(U))=1                         if |U| is odd;
    b(C(U))=1+2^(r/2+t)               if |U| is even.         (9)

**Proof.** The intrinsic descriptions above force every anti-automorphism
to swap x and y and reverse the ordered SCC word U. This proves the
self-converse criterion. The paths `x P_U y` contribute exactly one blue
tiling, because a palindromic rigid word U has h=b=1. A blue path with
adjacent y,x must put those vertices at the middle two positions, since
its positional reversal swaps them. Such paths are impossible when |U| is
odd. The two unsplit end placements are also impossible for `|U|>=2`.

For even U write `U=W -> W^op`. Let `a=|Aut(W)|`. Its anti-involutions
swapping the two halves number a: an isomorphism between the halves
determines the inverse action on the other half. Fix one such involution.
A blue split path chooses one vertex from each of its paired orbits for
the segment before y, and the other segment is forced by reversal. If S
is the subset chosen from W, the number of choices of the first segment
is `H(W[S])H(W\S)`. Thus this involution contributes J(W) paths. Each
blue path has its unique positional anti-involution, so the normalized
number of split paths is

    a J(W)/a^2=J(W)/a=2^(#1(W)+2#C3(W)).

Adding the one unsplit blue tiling proves (9).

The case `U=1` is deliberately excluded. It gives C3, whose automorphism
group moves x,y and is three times larger than the formula would predict.
The case `U=C3` is allowed and gives the nonregular five-vertex centre of
mass `(5,1)`, not the regular C5.

## 4. Fibonacci polynomials determine the sector coefficients and degrees

For a half-word W define the effective weight

    ell=#1(W)+2#C3(W),
    p_ell(z)=sum_(W of effective weight ell) z^|W|
            =sum_(j=0..floor(ell/2)) binom(ell-j,j) z^(ell+j).

It follows either by the last letter or by choosing j triangle positions
that

    p_0=1, p_1=z, p_ell=z p_(ell-1)+z^3 p_(ell-2),
    sum_(ell>=0) p_ell(z) u^ell=1/(1-zu-z^3u^2).            (10)

This is a literal Fibonacci-word correspondence. At z=1 the number of
words is `F_(ell+1)`; vertex size remains a separate grading. No equality
between a tournament and an arithmetic object is inferred from equal counts.

Every self-converse centre in the entire delta3 sector belongs uniquely
to one of the following families:

| U | Range | black deficit | centre polynomial at fixed ell |
|---|---|---:|---|
| W W^op | ell>=1 | 2^ell(2^ell-1) | z^2 p_ell(z^2) |
| W 1 W^op | ell>=1 | 2^(2ell+1) | z^3 p_ell(z^2) |
| W C3 W^op | ell>=0 | 2^(2ell+2) | z^5 p_ell(z^2) |

For offset a=2,3,5, the coefficient of `z^(a+2ell+2j)` is exactly
`binom(ell-j,j)`; all other coefficients vanish. The polynomial has degree

    a+2ell+2 floor(ell/2)                                   (11)

and minimum exponent a+2ell. This determines every coefficient and degree
for the sector, rather than only its first samples.

The first nonzero black-deficit centre polynomials in this sector are

| beta | polynomial |
|---:|---|
| 2 | z^4 |
| 4 | z^5 |
| 8 | z^5 |
| 12 | z^6+z^8 |
| 16 | z^7 |
| 32 | z^7+z^9 |
| 56 | z^8+2z^10 |
| 64 | z^9+z^11 |
| 128 | z^9+2z^11 |
| 240 | z^10+3z^12+z^14 |

In particular **R4 is the unique beta2 centre in every order within the
entire delta3 sector**. The sum over all beta in this sector is

    z^2[(1+z+z^3)/(1-z^2-z^6)-1-z].                         (12)

The powers have a precise origin: q=`2^ell` counts the split blue choices,
and the black deficit is `q^2-q` for the even inner word. The odd inner
words have one blue tiling and deficit a pure power of two. This is a
square/difference structure arising from paired paths, not an assumed
additive deficit law.

There are genuine low-deficit centres outside delta3. The Paley tournament
on seven vertices has `(h,b,beta,delta)=(9,3,6,14)`. Thus a small beta does
not justify restricting to the first few symmetry-slack layers.

## 5. Squares, cubes, and fourth powers in the SCC word grammar

There are two valid power expansions, with different meanings.

First let `G_m(z)` count ordinary **strong** classes of mass m. In the
Dirichlet convolution algebra in the mass index, the ordered-word law is

    W=(1-G)^(-1),        G_1=z+z^3,  W_1=A(z).

For an odd prime p this yields, with each formula an identity of series,

    W_p     =A^2 G_p,
    W_(p^2) =A^2 G_(p^2)+A^3 G_p^2,
    W_(p^3) =A^2 G_(p^3)+2A^3 G_p G_(p^2)+A^4 G_p^3.       (13)

The power `A^(j+1)` counts the j+1 rigid gaps around j nonrigid strong
letters. These formulas combine with (4); they do not identify beta with
the number of nonrigid letters. The series G_m are actual strong-class
inventories, not universally known polynomials merely because the first
few coefficients were computed.

Second, symmetry slack **is** additive over SCCs. If S_j(z) counts strong
classes of slack j, the continuation-width theorem gives

    S_0=z+z^3,
    S_3=z^2[A(z)-1-z],
    S_4=z^5,
    S_1=S_2=S_5=0.

The full slack generating function is

    F(z,u)=1/[1-z-z^3-sum_(j>0) S_j(z)u^j].                 (14)

Consequently its first genuinely composite layers satisfy

    [u^6]F=A^2 S_6+A^3 S_3^2,
    [u^7]F=A^2 S_7+2A^3 S_3 S_4,
    [u^8]F=A^2 S_8+A^3 S_4^2,
    [u^9]F=A^2 S_9+2A^3 S_3 S_6+A^4 S_3^3.                (15)

Unknown strong-centre terms remain displayed. In particular, the
phase-locked two-triangle hostile at slack6 cannot be silently absorbed
into independent triangle letters. Formula (15) is the exact way to
separate new strong centres from powers of previously classified ones.
The companion [arithmetic lane](arithmetic.md) gives a different exact
Fibonacci codec for the restricted alphabet `{1,C3,R4}`; its variable
counts vertex size and its R4 insertion count is additive, unlike beta.

As a direct hostile, the ordered words with r identical R4 letters have

    h=5^r,
    b=5^floor(r/2) 3^(r mod 2),
    beta=2,20,110,600 for r=1,2,3,4.                        (16)

Thus `D_2(z)^r` counts r selected centres in that restricted language. It
does not count deficit 2r. The outer-pair squaring in (2), and the ordinary
ordered-word powers in (13)--(15), are different operations.

## 6. What the centre classification gives black normalization

**Proposition 3.** Every self-converse node with beta2 is a leaf of the
simple black-support graph: it has one distinct black neighbour, joined
by exactly two parallel black lines, and has no black loop.

**Proof.** The sigma-twin theorem
[THM-4731, self-lines, twins and depth](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md)
pairs each black line with a distinct black line having the same merged
endpoints. A loop and its distinct twin would use four half-edges at the
node. The only possibility with two black half-edges is a pair of
non-loop parallel lines to a single other node.

These are black-support leaves, **not** leaves of the full metagraph.
They are mixed and have positive odd blue degree. The proposition makes
`D_2` relevant to fragile branches of black normalization, but not to
terminal black components without path-partner information.

For a rigidly padded mixed centre `A -> Z -> A^op`, the fibre factorization
is particularly sharp: all choices are central path orbits. Removing the
padding preserves `(h,b)` but changes the complement operation by changing
the Hamiltonian path's length and exterior arcs. It is therefore not a
black-path contraction theorem.

There is already a concrete boundary hostile in
[the defect-transitions blue lane, Section 3](../defect_transitions_20261020/blue.md):
`1 -> R4 -> 1` has beta2, but a blue Hamiltonian path with fixed
`(first,second,last)=(0,1,5)` admits no black path with that same boundary.
Its other path orbits do include black paths. Thus the centre count alone
cannot justify a boundary-preserving blue bypass. At least the admissible
boundary path set or actual complement-partner map is required.

The rigorous reduction achieved here is:

1. deficit2 classification reduces exactly to strong centres, by (5);
2. within the complete slack3 sector those centres reduce uniquely to R4;
3. every such node is a black-support twig, by Proposition 3;
4. attaching that twig to a common terminal basin still requires an actual
   black-normalization path, and mixed blue bridge removal still requires
   a bypass preserving the needed boundary data.

No proof of universal black connectivity, of a global `D_2=z^4`, or of
the absence of non-leaf full-metagraph bridges is claimed here.

## 7. Exact computation and reproduction

Run from the repository root:

```text
python 05-knowledge/results/scc_basins_20261020d/centers.py
python -O 05-knowledge/results/scc_basins_20261020d/centers.py
```

The script uses explicit raising checks. It records SHA256 hashes for the
inherited complete ordinary-class tables at orders3--9. It independently
enumerates Hamiltonian words for **every self-converse class of beta<=20**
in those tables: the set of normalized adjacency words computes h, and
literal path-reversal tests compute b. Score-preserving permutation search
independently checks their automorphism groups; subset DP computes f.

The strong-centre inventory at small deficit is:

| order | beta2 | beta4 | beta6 | beta8 | beta10 | beta12 | beta16 | beta20 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 4 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 5 | 0 | 1 | 1 | 2 | 1 | 0 | 0 | 0 |
| 6 | 0 | 0 | 1 | 0 | 0 | 3 | 0 | 0 |
| 7 | 0 | 0 | 1 | 1 | 0 | 0 | 2 | 1 |
| 8 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| 9 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |

Thus `D_2(z)=z^4+O(z^10)` is **FINITE-EXACT**, independently of the
all-order delta3 sector proof. The script also checks all57 quiet words
with `2<=|U|<=9`, including nonpalindromic hostile controls; exact counts,
automorphisms, feedback minima, moved supports, and strong connectivity all
match (7)--(9). Seven explicit larger SCC stacks verify the nonlinear
deficit law, including `C5 R4 C5` at order14 and four R4 blocks at order16.
The Fibonacci coefficient identities are independently checked through
effective weight12. Raw canonical centre words and all receipts are in
[centers.json](centers.json).
