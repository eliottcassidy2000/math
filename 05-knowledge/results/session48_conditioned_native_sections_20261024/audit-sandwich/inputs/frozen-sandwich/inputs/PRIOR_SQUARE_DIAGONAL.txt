# Symmetry slack three: split-word paths, Fibonacci spectra, and blue diagonals

Status: PROVED relative to the inherited exact slack-three classification;
FINITE-EXACT literal controls are in check.py and its receipts. No claim to
classify every strong tournament of a given Hamilton-path orbit count is made.
This is a research connection packet, not a scarce theorem-ID reservation.

## Inheritance and scope

1. `05-knowledge/results/continuation_width_symmetry_20261020/README.md`,
   Theorem 6: for symmetry slack delta=3f-s equal to three, the unique
   non-neutral strong component is T(W)=C3[1,W,1], with W an ordered word
   in K1,C3 of total vertex size at least two. Its endpoints x,y are
   intrinsically fixed, x -> W -> y -> x, and Aut(T)=C3^t, where t is the
   number of triangle letters. Arbitrary neutral words may precede/follow.
2. `01-canon/theorems/THM-284-linear-coefficient-formula.md`: the proved
   single-backedge raw path count H=1+2^(skip-1). The same special case is
   recorded historically in `THM-115-h21-permanent-gap.md`. Neither the
   original one-backedge mechanism nor that special formula is new here.
3. `05-knowledge/results/deficiency_three_20261020/centers.md`: the exact
   center rule h(A -> C -> A^op)=h(A)^2 h(C),
   b(A -> C -> A^op)=h(A)b(C). Neutral words have h=1.
4. `05-knowledge/results/fixed_degree_metagraph_rationality_20261021/`
   proves fixed-h rationality abstractly with finite but partly unknown
   strong-atom polynomials. Here we compute one exact symmetry-slack slice
   of those polynomials. We do not replace them by this slice.

Here H is raw Hamilton-path count, h=H/|Aut| is the number of path orbits,
and b is the number of those orbits fixed by converse plus path reversal,
when the tournament is self-converse. Ordinary types are not merged with
converse unless explicitly stated. Put p(z)=1-z-z^3 and A(z)=1/p(z).
The spectral exponent e below is log2(h-1), NOT the center deficit h-b.

Live concept board: the inherited one-backedge subset split; independent
triangle rotations; exact symmetry-slack normal forms; finite strong-atom
polynomials at fixed h; the square-minus-diagonal involution; literal
metagraph incidence ports versus their unresolved target wiring.

## 1. An exact path-orbit alphabet

Write a=#K1 letters, t=#C3 letters in W, so |W|=a+3t>=2.

**Theorem 1.**

    H(T(W)) = 3^t + 2^a 12^t,
    h(T(W)) = 1 + 2^(a+2t).

More strongly, the path-orbit set is the disjoint union of one direct
orbit and the product of local alphabets

    {0,1} for each K1, and {0,1,2,3} for each C3.

**Proof.** A Hamilton path avoiding y -> x must start at x and end at y,
and its middle is a Hamilton path in W. Every Hamilton path in an ordered
join visits its factors in order. The unique path orbit of each K1/C3
factor therefore gives one direct path orbit and 3^t raw direct paths.
Every other Hamilton path is uniquely

    P_L, y, x, P_R,

where L is a subset of W, R is its complement, and P_L,P_R are Hamilton
paths in the induced subtournaments (the empty path has count one).
Every arc required at either end of y,x has the right direction. Thus

    H(T(W))=H(W)+sum_(L subset W) H(W[L])H(W[W\L]).

For an ordered join, both induced path counts factor letter by letter,
so the subset sum factors too. A singleton contributes 2. A triangle
contributes 3+3+3+3=12: the four terms are for |L|=0,1,2,3.
Under the C3 rotation action, for each fixed j=|L| there is exactly one
orbit of the local ordered pair of paths. For j=0,3 this is rotation on
the three Hamilton paths of the whole triangle; for j=1,2 it is rotation
on the three choices of the singleton side, with the other side's path
forced. Independent rotations of all triangle factors give exactly the
product of the local alphabets. This uses the FULL actual automorphism
group from the inherited classification, not a selected subgroup. QED.

The all-singleton slice t=0 recovers THM-284's one-backedge theorem.
A triangle is not merely three singleton choices: its normalized split
alphabet has four letters rather than eight. The resulting exponent is
number of singleton letters plus twice the number of triangle letters.

## 2. Explicit generating functions at fixed spectrum

The exponent is e=a+2t>=2. By the uniqueness part of the inherited
classification, two different ordered words W give different strong types.
With z recording vertices and u recording e, the strong ordinary series is

    S(z,u)=z^2 [1/(1-uz-u^2 z^3)-1-uz].                 (1)

The subtraction removes W empty and W=K1, for which the inherited
fixed-endpoint/aut-group conclusion fails. In particular W=K1 gives C3,
with actual h=1, NOT h=3; it has symmetry slack zero.

At a fixed e>=2 the strong polynomial is explicitly

    P_e(z)=sum_(t=0)^floor(e/2) binom(e-t,t) z^(e+t+2). (2)

Indeed a=e-2t, so there are a+t=e-t ordered letters, of which t are
triangles. The binomial coefficient is NOT binom(e,t): at e=2 the answer
is P_2=z^4+z^5, not z^4+2z^5. Summing over vertex order gives
P_e(1)=F_(e+1), with F_0=0,F_1=1. More precisely the polynomials satisfy
P_e=z P_(e-1)+z^3 P_(e-2) for e>=4, starting with

    P_2=z^4+z^5,
    P_3=z^5+2z^6,
    P_4=z^6+3z^7+z^8,
    P_5=z^7+4z^8+3z^9.

The extreme strong orders at fixed e are e+2 and e+floor(e/2)+2,
both attained. No assertion is made that these bounds hold outside
symmetry slack three.

Neutral outer SCC words preserve h, so the full ordinary slack-three
fixed-h=1+2^e series is

    F_(delta=3,e)(z)=P_e(z)/p(z)^2.                    (3)

Its reduced denominator over Q is exactly p(z)^2: p is an irreducible
cubic (rational-root test), and at its positive real root rho, P_e(rho)>0.
The bivariate full series is S(z,u)/p(z)^2. At u=1 this is exactly the
inherited slack-three OGF z^2(A(z)-1-z) A(z)^2, a required consistency check.

### An exact finite sampler, with its sampling measure stated

The fixed-e strong types are in bijection with compositions of e into 1
and 2: replace each 1 by K1 and each 2 by C3. Uniformly among these strong
isomorphism types, choose the first K1 with probability F_e/F_(e+1),
or the first C3 with probability F_(e-1)/F_(e+1), then recurse on the
remaining exponent (the remainder zero has one empty continuation).
This is an exact sampler of strong TYPES at fixed e, not labelled
tournaments, path orbits, or converse-merged types. Conditional on both
e and vertex order n, put t=n-e-2 and choose the t triangle positions
uniformly among e-t letter positions. No assumption about typical
unrestricted tournaments is involved.

### Balanced exponent cuts recover the older Fibonacci-product sequence

The sequence 1,2,6,15,40,104,273 and its identity as consecutive Fibonacci
products/sums of squares are inherited, not new: see
`05-knowledge/results/metagraph_spectrum_torus_session16_20261020.md`,
Section 5, including its matching interpretation and its warning about
unproved move-graph identifications. The following map gives a concrete
new realization inside the slack-three normal form.

Let R_r(z) count neutral words of exponent r by their vertex order, with
R_0=1, R_1=z, and R_r=z R_(r-1)+z^3 R_(r-2). Thus R_r(1)=F_(r+1).
For r>=1 choose two words U,V of exponent r and form

    W=U reverse(V),       T(W)=C3[1,W,1].

The exponent-midpoint cut is unique because every letter has strictly
positive exponent weight (1 for K1, 2 for C3). Therefore this map is
injective into actual strong isomorphism types; it produces exactly the
words admitting that midpoint cut. Converse swaps U,V, so these types
form a literal class-level square of size F_(r+1)^2, with diagonal
F_(r+1) and non-self-converse pairs counted by binom(F_(r+1),2).
This class-level square is DIFFERENT from Section 4's path-level square,
whose coordinate set at an even self-converse core has size 2^r.

A word of total exponent 2r either admits the midpoint cut or has a
unique triangle letter straddling it. In the latter case it is uniquely

    W=L C3 reverse(R),       exponent(L)=exponent(R)=r-1.

Thus the full strong fixed-exponent polynomial has the exact partition

    P_(2r)(z)=z^2 R_r(z)^2 + z^5 R_(r-1)(z)^2.          (11)

At z=1 this realizes F_(2r+1)=F_(r+1)^2+F_r^2, not merely a numerical
coincidence. In each square, converse swaps the two coordinates; its
diagonal gives respectively the even-order and odd-order blue channels
of Section 3. The midpoint condition is indispensable: at r=2,
W=(1,3,1) has exponent four but no exponent-two cut. There are four
balanced-cut strong types, whereas all exponent-four strong types number
five. It is the central-C3 straddler, not an exception to the formula.

Adjoin the ACTUAL EMPTY TOURNAMENT as one separate unit, not the invalid
slack-three expression C3[1,empty,1]. The number of this unit plus all
balanced-cut cores of half exponent 1 through r is

    1+sum_(e=1)^r F_(e+1)^2 = F_(r+1) F_(r+2),
    r=0,1,2,...: 1,2,6,15,40,104,273,... .             (12)

There is an explicit bijection with pairs of monomer/dimer tilings of
lengths r and r+1, the standard path-matching model. Given such a pair
(A,B), inspect B's first tile. If it is a monomer, remove it: the two
remaining tilings have equal length r, so convert their monomers/dimers
to K1/C3 words U,V and return the balanced-cut core. If it is a dimer,
remove it and swap the two tilings, obtaining lengths r-1,r; recurse.
At r=0 return the separate unit. This map is invertible: prepend a
monomer at the stopping level and undo each earlier step by swapping
and prepending a dimer. The output's half exponent identifies the
stopping level, so there is no missing history choice.

Connection boundary: this is an object bijection and an explicitly
tracked class-level converse quotient. It does not assert that a natural
matching move graph, a torus face-flip graph, or a metagraph-neighbor
relation is carried to another. The cumulative tiling-pair map is not
vertex-order or Hamilton-spectrum preserving: its stopping layer varies. The inherited 273-state edge-count
hostile remains in force. At the class level the retained coordinates
are the total exponent, the unique cut/straddling letter, and both words;
retaining only the Fibonacci number would erase the map.

## 3. Converse and actual blue path orbits

**Theorem 2.** T(W) is self-converse iff W is a palindrome as an ordered
K1/C3 word. When it is, its blue count is

    b(T(W))=1                         if |T(W)| is odd,
    b(T(W))=1+2^(e/2)                 if |T(W)| is even. (4)

The even case automatically has even e. For a full slack-three type
L -> T(W) -> R, self-converse means W is palindromic and R=L^op as an
ordered neutral word. In that case the same h and b apply.

**Proof.** The fixed intrinsic endpoints swap under converse, and W is
sent to its reversed/opposed SCC word. Both possible atoms are
self-converse, so the criterion is exactly a palindrome. On the product
alphabet of Theorem 1, the induced converse/path-reversal action is

    (j_1,...,j_l) -> (s_l-j_l,...,s_1-j_1),            (5)

where s_i is 1 or 3. This action is on FULL path orbits, so there is no
assumption that one selected anti-automorphism enumerates all blue paths.
The direct orbit is fixed. For split orbits, every mirrored letter pair
may choose its first j freely, forcing the second. A central letter has
odd size and cannot satisfy j=s-j. Therefore a palindrome with an odd
number of letters has no fixed split orbit; its vertex order is odd.
A palindrome of even length has a/2 singleton pairs and t/2 triangle
pairs and hence 2^(a/2)4^(t/2)=2^(e/2) fixed split orbits. The center
factorization with h(L)=1 proves the padding claim. QED.

Define E_e and O_e as the strong self-converse polynomials in even and
odd vertex order respectively. Their complete bivariate series are

    E(z,u)=z^2 [1/(1-u^2 z^2-u^4 z^6)-1],
    O(z,u)=z^2 [(uz+u^2 z^3)/(1-u^2 z^2-u^4 z^6)-uz]. (6)

Their fixed-e sum Q_e=E_e+O_e is finite. Mirrored neutral outer words give

    J_(delta=3,e)(z)=Q_e(z)/p(z^2),                    (7)
    M_(delta=3,e)(z)=[P_e(z)/p(z)^2+Q_e(z)/p(z^2)]/2. (8)

The first counts self-converse ordinary classes, the second converse-
merged classes. These formulae are exact, not a claim that every fixed-h
class belongs to delta=3. The blue-even and blue-odd channels in (6)
retain different b-values even when e agrees.

Examples:
- e=2, h=5: two strong types, W=(1,1) at order 4 with b=3, and W=(3)
  at order 5 with b=1. The latter has raw H=15 but Aut order 3; it is
  not cyclic R5, whose raw H=15 has Aut order 5 and h=b=3.
- e=3, h=9: W=(1,1,1) at order 5 has b=1. W=(1,3),(3,1) at order 6
  form a non-self-converse pair; their raw H=27 and Aut order is 3.

## 4. Square minus diagonal as a literal quotient map

For an even-order self-converse core write W=V V^op (word reversal,
with each atom replaced by its converse). Let X be the split alphabet
product of V, so |X|=N=2^(a(V)+2t(V)) and h=1+N^2.

Write a split state as (ell,r), where ell uses V's natural order and r
uses its mirrored partner's natural order. Replace r by its reversed
coordinatewise complement, called r*. Then the full split-state set is
literally X x X with coordinates (ell,r*), and action (5) becomes

    (ell,r*) -> (r*,ell).

Thus the N blue split states are the diagonal, the unique direct orbit
adds one more blue state, and nonfixed reflection pairs are in bijection
with unordered two-element subsets of X. The fundamental paired-black
local incidence degree is therefore

    (h-b)/2 = (N^2-N)/2 = binom(N,2).                 (9)

This is an explicit equivariant set map, not merely agreement between
a square-minus-diagonal numerical expression and a degree formula.
For odd-order self-converse cores there is no fixed split state, so the
same degree convention gives (h-b)/2=2^(e-1). A non-self-converse pair
has fundamental degree h=1+2^e, as in the inherited convention.

These bijections enumerate local path-orbit incidences. They do NOT
identify the complement-neighbor tournament, prove different neighbors,
exclude loops, determine edge multiplicities, or establish connectivity.
Reconstructing that wiring requires the literal source tournament and
Hamilton-path representative, followed by the nonpath-complement map.

**FINITE-EXACT wiring hostile.** The independent review's literal port
receipts compare W=(1,1,1,1) (core order 6) and W=(3,3) (core order 8).
Both have h=17,b=5,N=4 and six local black incidences. Nevertheless their
actual converse-merged nonpath-complement targets number respectively
five, with multiplicities [1,1,1,1,2], and six, with multiplicities
[1,1,1,1,1,1]. Thus even the smallest nontrivial square-minus-diagonal
case beyond N=2 does not force distinct neighbors. See the accompanying
`session23-slack-spectrum-review/normal.json`, `literal_port_receipts`.

**Stronger FINITE-EXACT same-order hostile, credited to the independent
review.** W=(1,3,3,1) and W=(3,1,1,3) have identical
(n,H,|Aut|,h,b)=(10,585,9,65,9), hence 28 local black incidences each.
Each has 28 distinct converse-merged complement destinations, all with
multiplicity one, but their two destination sets are DISJOINT. Thus even
vertex order, raw H, full automorphism order, h, b, and the multiplicity
profile together do not determine the neighboring classes. The exact
receipt is `session23-slack-spectrum-review/normal.json`, under `rows`
with `word=[1,3,3,1]` and `word=[3,1,1,3]`; each row's
`literal_port_receipts` retains the paths and `target_merged_canonical_key`
values. The optimized receipt independently repeats the same controls.

### Exact fundamental-degree support inside the slack-three layer

Let D_k^(3)(z) count converse-merged nodes of symmetry slack three and
fundamental paired-black degree k. The preceding formulas give the COMPLETE
support and exact series within this layer:

    k=2^(e-1), e>=2:          D_k^(3)=O_e(z)/p(z^2),
    k=binom(2^r,2), r>=1:     D_k^(3)=E_(2r)(z)/p(z^2),
    k=1+2^e, e>=2:            D_k^(3)=[P_e/p(z)^2-Q_e/p(z^2)]/2,
    all other k:             D_k^(3)=0.                         (10)

The three degree sets are disjoint: the first consists of powers of two
at least two; the second is 1 or a power of two times the odd number
2^r-1>1; the third consists of odd integers at least five. Each displayed
series is nonzero. Palindromic W with an odd central atom realizes every
first-family e; W consisting of 2r singletons realizes every second-family
r; asymmetric neutral outer padding realizes every third-family e.
For example

    D_1^(3)=z^4/p(z^2),
    D_2^(3)=z^5/p(z^2),
    D_4^(3)=z^5/p(z^2),
    D_6^(3)=(z^6+z^8)/p(z^2).

The degree-one formula recovers the inherited self-converse deficit-two
family. It does not include non-self-converse degree-one nodes, which have
h=1 and slack zero. Nothing here asserts that these are ALL nodes of
degrees two, four, or six across other slack layers. This gives fully
explicit infinite families inside Session22's general fixed-degree
rationality theorem, rather than replacing its unknown strong polynomials.

The reduced denominators in the first two lines of (10) are exactly
q(z)=p(z^2), and in the third line exactly p(z)^2 q(z). Hence the minimal
eventual constant-coefficient recurrence orders over Q are respectively
6 and 12, not merely upper bounds. This is a slice-specific statement.

For completeness, q is irreducible over Q. Negating it gives
f=z^6+z^2-1. Modulo 2, f=(z^3+z+1)^2 with the cubic irreducible, so
any proper monic integer factor would have degree three. Modulo 3,
f=(z^2+1)(z^4-z^2-1), both factors irreducible, which excludes degree
three. The quartic has no linear root; in a hypothetical product of monic
quadratics the constants must be 1 and 2, the cubic coefficient makes
the linear coefficients a,-a, and the linear coefficient forces a=0.
Their quadratic coefficient would then be 1+2=0, not -1. Gauss's lemma
finishes the irreducibility proof. Thus p and q are coprime. Every
nonzero E_e,O_e,Q_e is positive at sqrt(rho), where p(rho)=0; P_e(rho)>0.
These positive values prohibit cancellation of q or of either copy of
p in the respective expressions in (10).

## 5. Connection contract and hostile boundaries

Source: the inherited delta=3 normal form and its actual full Aut action.
Target: explicit fixed-h strong polynomials, rational padded counts,
converse channels, and local paired-black incidences.
Map: split paths at the unique endpoint backedge, quotient each atom by
its actual independent rotation, and retain the perletter occupation j.
Preserved: the complete Hamilton-path orbit at a fixed W; converse action;
h,b and the fundamental local degree; exact ordinary class enumeration
when W and the outer SCC words are also retained.
Destroyed if only e is retained: word order, vertex order, raw H, Aut,
self-converse status, b, and all external wiring. If W is retained,
raw H and Aut are recovered but neighbor data still need literal paths.
Needed sidecar: ordered W, outer SCC words, and a literal path/port for
any claim about metagraph targets.

Cheap hostile controls:
- W empty and W=K1 are not slack-three cores; their endpoints are not
  fixed as required and the normalization formula does not apply.
- At e=2, C3[1,C3,1] and cyclic R5 have the same raw H=15 but different
  automorphism orders and hence different h and b.
- W=(1,3) and (3,1) have equal e but are different, non-self-converse
  types, exchanged by converse.
- At symmetry slack six, the inherited synchronized pair of triangles
  has diagonal C3 rather than independent C3 x C3. The local product
  quotient proof must NOT be applied there.
- The fixed-e coefficient is binom(e-t,t), not binom(e,t).
- R5 -> R5^op is not neutral padding: it has h=9,b=3. Only words in
  K1,C3 are harmless here.

## 6. Reproduction and next directions

Run from this packet directory:

    python3 check.py --receipt receipt.json
    python3 -O check.py --receipt optimized-receipt.json

Both normal and optimized runs pass 20,745 explicit checks over 85 words
and 90,223 literal Hamilton paths. The checker imports no repository producer code. It constructs every W
with 2<=|W|<=10, enumerates all its literal Hamilton paths by DFS, checks
the raw total independently by subset DP, computes its FULL automorphism
and anti-isomorphism sets by adjacency-compatible backtracking, computes
feedback by exact subset DP, and verifies delta=3, the complete local
orbit alphabet, converse action, and b. Fixed-e coefficients and symbolic
rational identities are separately checked through e=14. Assertions are
not used, so Python optimization does not disable checks. The separate
midpoint partition is checked through half exponent seven, and the
old tiling-pair bijection through r=10, including all 12,816 target
objects at that final level.

New questions, not theorems claimed here:
1. Carry actual path representatives for the X x X coordinates into the
   complement target. Which unordered pairs share a neighbor, and which
   form loops? The square-minus-diagonal structure alone cannot answer.
2. Compare P_e with the unrestricted strong fixed-h atom polynomials.
   The computed slack-three slice is an explicit lower bound, not an
   asserted equality. Other symmetry slack layers may add atoms.
3. At slack six, keep the synchronized diagonal C3 action as the missing
   sidecar and ask for a finite transfer matrix on local split states.
   Independent letters have failed, but a finite-state quotient may survive.
