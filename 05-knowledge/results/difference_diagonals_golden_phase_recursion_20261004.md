# Difference families, golden phase trees, and the information on a diagonal

2026-10-04 (America/Denver).

**PROVED:** the all-denominator golden periodic-lift theorem, its binary and
ternary cycle trees, signed root-splicing families, and the elementary pair
and incidence identities specified below. **FINITE-EXACT:** the explicitly
bounded censuses in the companion notes. **CITED:** the classical zeta value
and the stated physical convention. **OPEN:** rational golden passports for
every integer input, Collatz convergence, and any physical derivation of the proposed coupling.
No claim of historical novelty is made for the classical ingredients.

## Inheritance and research board

The anchor is the user's family that retains its route home. The niche is
the difference diagonal and what it forgets; the wildcard is the connection
between primitive lattice density, prime quotients, and the proposed lambda.
The broader LRC(14) target remains OPEN; this session supplies no new LRC bound.

The closest proved mechanisms are the source-preserving return compiler and
reversible carry record in [the previous synthesis](golden_primes_carries_route_compiler_20261003.md),
and the exact signed lattice in [golden carriers](collatz_golden_carriers_20261004.md).
The canonical hostile is the rational cycle 1/13 sharing a denominator ideal
with the integer cycle -5. The corrected near miss is cancelling a primitive
coefficient pair as though it were a unit of the golden residue ring. The
least-used coordinate is the numerator's entire phase orbit, with the upper
boundary tag retained separately. The live board is:

1. a difference and its position along the diagonal;
2. a golden denominator ideal and numerator phase;
3. an ordered route with a certified root;
4. a quotient and its reconstruction data;
5. the diagonal, two triangular sheets, and typed flags;
6. primitive-point density and quadratic norm.

The companion notes supply proofs, hostile examples, reproduction commands,
and complete finite universes:

- [Golden periodic lifts and denominator trees](denominator_fibre_recursion_20261004.md).
- [Independent arithmetic realization filter](denominator_arithmetic_filter_20261004.md).
- [Difference-family operations and the clarified f/g laws](difference_family_operations_20261004.md).
- [Pair tiles, medial maps, and recoverable incidence](medial_pair_tiles_20261004.md).
- [Primitive visibility, golden norms, quadratic secants, and lambda](difference_visibility_norms_20261004.md).

During this session, incoming commit952b93f75 independently supplied
[an L-shaped fundamental-domain proof and the tower above76](difference_families_20261004.md),
plus [a rational-anchor return compiler](rational_anchor_returns_20261004.md).
The overlap is acknowledged explicitly: the phase correspondence and
difference-pair arithmetic have parallel proofs, not separate novelty claims.
The additional connections below use the numerator selector, the independent
arithmetic filter, and the compatible residue ring.

## 1. The ideal equality leads to a stronger, reconstructive invariant

The earlier prime connection was the exact equality

    (Phi_5(phi))=(phi^5-1)=(phi-4)=(11,phi-4)

in O=Z[phi]. It identifies the same quotient in the repunit, shift-clock,
and denominator constructions. The next step is to retain a numerator phase
in addition to that ideal. Set beta=phi^-1 and use the upper golden map

    B(x)=phi*x-1_(phi*x>1),       0<=x<=1.

For x=(a+b*phi)/q its coefficient pair changes by

    (a,b) -> (b-d*q,a+b),        d=1_(phi*x>1).

Modulo q the digit disappears, leaving the invertible matrix

    M=[[0,1],[1,1]].

**All-denominator theorem.** Every nonzero class of Q(phi)/O has exactly one
periodic representative under B. Its exact period equals the period of its
coefficient phase under M. The zero class has precisely three periodic
representatives: 0, and the cycle beta <-> 1.

This is stronger than an invariant: on periodic nonzero states the quotient
has an exact inverse. The proof combines a contracting conjugate coordinate
with a seven-vector overlap check. Periodic conjugates lie in [-phi,1]; for
x>beta the sharper lower bound is -beta. Any two periodic points in the
same class differ by one of 0, +/-1, +/-beta, +/-beta^2. The branch bound
leaves only the stated zero-class collisions. A bounded lattice orbit gives
existence. [Full proof and constructive inverse](denominator_fibre_recursion_20261004.md).

For q>1, the exact-denominator-q periodic points therefore number

    J_2(q)=q^2 product_(p|q)(1-p^-2).

There are q^2+2 periodic points whose denominator divides q, because the zero
phase has two extra lifts. For a transient state the phase still forgets its
present lift and digit; that loss has only been repaired on the periodic set.

The ideal-11 hostile now splits exactly. Evaluation at phi=4 modulo11 gives

    Theta(-5)=(-1+7phi)/11   -> phase5,
    Theta(1/13)=(1+4phi)/11  -> phase6.

Multiplication by4 has two nonzero five-cycles: the quadratic residues
{1,3,4,5,9} and the nonresidues. The two examples occupy different cycles.
Both still have the same denominator ideal. Thus the phase orbit resolves
the specific ambiguity that the ideal could not resolve, without turning
the ideal into a general integer-orbit certificate.

## 2. Two distinct recursive directions, each with a precise preserved target

At exact denominator 2^a every phase cycle has length 3*2^(a-1), and there
are 2^(a-1) cycles. At denominator 3^a every cycle has length 8*3^(a-1), and
there are 3^(a-1) cycles. Here a>=1. Reduction from p^(a+1) to p^a gives
exactly p child cycles over every parent, each p times as long.

There is an explicit selector. At parent denominator q=p^a, choose an
anchor v and its period L. Write M^L=I+qU and c=Uv modulo p. For p=2,3,
U is invertible, hence c is nonzero. The p^2 lifts v+qw return after L steps
by w -> w+c. The p children are exactly the translation lines, distinguished
by det(c,w) modulo p. This supplies a recursive address rather than only a
count of descendants.

The incoming independently audited tower continues the user's76 directly:
at q=4*19^k, k>=1, every primitive phase has period18*19^(k-1), and there
are240*19^(k-1) cycles. Each parent has19 children. Its stronger-than-order
identity is M^18-I=76*M^9: the leading matrix modulo19 is invertible even
on primitive nonunit phases. The same determinant selector applies with
p=19 and U=(M^L-I)/q. This extends the branch-address construction beyond
the two inert-prime towers, without cancelling a nonunit.

The first binary split is concrete:

    parent phi/2
        -> child phi/4
        -> child -1+3phi/4.

Both children have period6. Their parent operation is "double, reduce
modulo O, then take the unique periodic lift." They correspond to rational
arithmetic cycles through 1/29 and 5/7, respectively. Neither is an integer
cycle. This is the minimal test that prevents the golden tree from being
mistaken for a tree of new integer basins.

The second direction really does preserve an integer route. For any finite
positive valuation word w, write its odd-step affine map as

    U_w(n)=(3^p*n+S)/2^A,   p=len(w), A=sum(w), C=2^A+3S.

For each root r in {1,-1,-5,-17}, choose sufficiently large a0 satisfying

    2^(A+a0)*r=C modulo 3^(p+1).

With T=2*3^p, every member

    n_k=(2^(A+a0+Tk)*r-C)/3^(p+1), k>=0,

has the actual head w, then the valuation a0+Tk, then arrives at r. The
large starting exponent places every preterminal odd vertex outside the
known cycles; the final halving chain may enter a cycle before reaching r. The
recurrence is affine with multiplier 2^T; its golden coordinate contracts
by beta^T toward a shared finite-prefix boundary while retaining its
denominator ideal and phase orbit. For w=(1), concrete sources are
227 -> ... -> 1, -1821 -> ... -> -1, -285 -> ... -> -5, and
-3869 -> ... -> -17, each with exactly two accelerated odd steps.

| Direction | Operation | Preserved predicate | Necessary extra information |
|---|---|---|---|
| Denominator tree | q -> pq, then periodic lift | Golden periodicity and parent covering | Integer-realization test |
| Certified-root family | Increase the final valuation by 2*3^p | Actual integer route to the chosen root | Ordered head, root, first-hit convention |

Every finite head occurs in all four certified families. Consequently no
fixed finite head can alone decide which of these four root labels applies.
This does not exclude methods that retain the source integer or use an
unbounded certificate. Shared-prefix recursion and shared-root recursion
intersect, but preserve different predicates.

The four denominators {1,2,11,76} are therefore meaningful arithmetic labels,
not four objects canonically identified by their count with the four
isomorphism classes of tournaments on four vertices. At q=11 there are
already 13 golden cycles; at q=76 there are 240, all of period18. The signed
integer test selects a much smaller subset.

For the user's 105=3*5*7 the CRT decomposition gives a particularly clean
split. Modulo5 the four nonzero eigenline phases have period4, and the other
20 have period20. Combine these with the 8 phases of period8 modulo3 and
48 phases of period16 modulo7:

| Kind of phase modulo105 | Points | Period | Cycles |
|---|---:|---:|---:|
| Primitive but not a ring unit | 1536 | 16 | 96 |
| Ring unit | 7680 | 80 | 96 |

The common clock80 hides two equally numerous families of cycles. This
decomposition is proved by CRT and independently reproduced by the exact
golden lattice. All192 cycles fail the integer-realization test.

The new independent arithmetic filter enumerates every golden cycle for
q=1,...,40,64,76,81,105: 1065 cycles containing40730 periodic points. For a
word with s ones and e zeros, compose the ordinary branches to
(3^s*n+K)/2^e; its unique rational fixed point is K/(2^e-3^s). Exact replay
tests integrality and every parity edge. Exactly five integer cycles remain,
through0,-1,1,-5,-17. This is a complete result within the stated denominator
universe, not an exclusion at other denominators. The smallest denominator
with a rational-only cycle is q3, whose arithmetic root is1/11.

**Reuse the failed integer cycles.** Of the1060 rational-only cycles,48 are
negative. Every one supplies an expanding chart for the incoming rational
return compiler. There is a direct canonical choice: let r0<0 be the cycle
member of smallest absolute value. At a nonempty odd prefix,
r_i=M_i*r0+C_i<=r0 with C_i>0, so M_i>1. Thus this anchor automatically
satisfies every strict prefix-growth condition. The general cyclic-minimum
rotation argument is inherited from the incoming compiler; this canonical
negative representative makes its hypothesis transparent.

The filter now instantiates all48 charts at their least admissible repeat
count and exit budget. It saves each exact source cylinder and endpoint,
and independently replays three sources per cylinder:144 successful first
descents, always compared with the original source. All331 cyclic odd phases
of the51 negative cycles also pass the exact rotation test. These are
operational compiler instances; their mutual overlap and added coverage
relative to existing banks have not been classified.

One exact bridge is

    root -19/11, odd word(1,1,2), raw word1010100,
    Theta=(-4+20phi)/29, golden denominator29.

That cycle fails the integer-cycle gate, yet the guarded repeated chart
gives the integer first descent

    999 ->1499 ->2249 ->1687 ->2531 ->3797 ->89.

The general incoming compiler proves an all-height family after imposing
its exit budget; it does not identify the rational cycle with these integer
routes. This is the useful transfer: change the target predicate from
"is an integer cycle" to "supplies a guarded integer return chart."

## 3. Difference families have commutative arithmetic with exact memory

For every signed integer d define

    D_d={(a,b) in N_0^2 : b-a=d}.

This implements the requested integer-as-family viewpoint. Addition of
pairs descends to addition of differences. Componentwise multiplication
does not: (7,11)*(1,2) has difference15, while (8,12)*(1,2) has difference16.
The crossed product does descend:

    (a,b) tensor (c,e)=(ae+bc,ac+be),
    difference of the product=(b-a)(e-c).

It is associative, commutative, and distributive, with identity (0,1).
After quotienting by equal differences it realizes the ring Z. Before
quotienting, keep a translation register t=min(a,b). The raw pair is

    (max(-d,0),max(d,0))+(t,t).

The complete, lossless operations in these coordinates are

    (d,t)+(e,u)=(d+e,t+u+(|d|+|e|-|d+e|)/2),
    (d,t) tensor (e,u)=(de,|d|u+|e|t+2tu).

Thus the diagonal register plays the role that the carry polynomial played
for base-phi normalization: the quotient gives ordinary arithmetic, while
the register preserves the discarded representative. An operation history
still requires its own DAG; the final representative is not a factor tree.

The supplied f clarification has an unambiguous triangular component. With
T(n)=n(n+1)/2, A+B=N-2, A<B and d=B-A,

    T(N-1)=T(A)+T(B)+(A+1)(B+1),
    (A+1)(B+1)=(N^2-d^2)/4.

The allowed differences satisfy 0<d<=N-2 and d=N modulo2. At N=10,
A=3,B=5 gives 45=6+15+24; the equal-pair boundary A=B=4 gives
45=10+10+25. Moving A up and B down transfers d-1 units from the triangles
to the rectangle. This exact law is retained without silently choosing a
parenthesization for the unresolved string `(A+B)f1fAfB`.

Read literally over ordinary integers or reals, the g condition u=ABu
forces u=0 unless AB=1. It becomes a nontrivial universal relation in the
quotient Z/(AB-1): AB acts as the identity there. For AB=100 this is modulo99;
the factor10 itself acts as the identity modulo9. These are different from
the multiplication table modulo10. In the golden setting, the same
operation-minus-identity construction is O/(phi^L-1).

A different valid type is an invariant rational family r*c^Z under scaling
by c=AB. It has no distinguished first member. A forward family r*c^N has
a root but scaling drops that root. Thus retaining a route home requires
the base point and height; scale invariance by itself has erased them.

The periodic golden lifts also carry a genuine commutative addition:
take [x+y] in Q(phi)/O, then its canonical periodic lift. Choose0 for the
zero class. This yields the O-module (Q/Z)^2; multiplication by phi is B.
Ordinary multiplication of arbitrary classes does not descend, as
1/2~3/2 but their products with1/2 differ by1/2. Keep the beta/1 boundary
tag if the negative integer cycle through -1 must remain distinguishable.

There is a stronger structural boundary: **every biadditive product on
H=Q(phi)/O is zero**. Indeed H is both torsion and divisible. If n*a=0,
choose c with n*c=b. Then a*b=a*(n*c)=(n*a)*c=0. A novel distributive
multiplication cannot repair this particular quotient without changing
the retained object or its addition.

The useful alternative stores a compatible residue tower:

    v_a in O/(p^a),   v_(a+1)=v_a modulo p^a.

Levelwise addition and multiplication preserve compatibility, so complete
towers form a ring. Each phase point has p^2 lifts per level. For p2 or3,
over a primitive parent cycle those split into p child cycles, each with p
visits above a chosen parent phase: one branch coordinate and one relative
phase coordinate. The determinant selector stores the former; the latter
is needed for addition. Modulo2, 1 and phi lie in the same cycle, but
1+1=0 and 1+phi=phi^2 lie in different cycles. Multiplication alone does
descend to cycle orbits, because shifting either factor shifts the product;
it is the full ring structure that requires the relative phases. This
provides a precise choice of a recursively stored object supporting both
commutative operations. It still needs the separate arithmetic route test.

## 4. The dotted diagonal repairs a real graph-symmetry ambiguity

For an n-element set, the pair plane has

    n^2=n+2 binom(n,2).

The diagonal is the fixed set of transposition. The two strict triangles
are its two ordered sheets. A tournament selects one sheet at each ordinary
pair; optional loops occupy the diagonal separately. Values such as a+b or
ab are labels, not addresses: 1*4=2*2 already merges a strict pair with a
diagonal point.

On the six unordered pairs of a four-element set, adjacency by intersection
gives the octahedron J(4,2). It has48 automorphisms, while only24 come from
permuting the original four elements. Edge-complementation supplies the
extra symmetry. Add the four diagonal pairs and connect tiles when their
supports intersect. Each ordinary pair is now identified by its two
diagonal neighbors; the automorphism group becomes exactly S_4. The dotted
diagonal restores the identity of the original endpoints.

At five vertices there are10 ordinary pairs. Each belongs to two of the
five four-element stars, meets6 other pairs, and is disjoint from3. Hence
the 10*10 ordered pair-of-pair addresses split into

    100=10 equal +60 intersecting +30 disjoint.

The disjointness graph is the Petersen graph. This is a precise way that
3,4,5 and10 occur in one object. A separate address labelling is needed to
turn it into the numerical mod10 table; no such identification follows
from the count100. The original tournament wedge retains its own exact
split15=6 free arcs+4 forced Hamiltonian-path arcs+5 optional loops.

For a polyhedral map G, its medial map has vertices indexed by E(G), and
two types of faces indexed by V(G) and F(G). The underlying medial map is
shared by G and its dual. Retaining the two face types reconstructs the
original; forgetting them may add dualities. A typed flag graph goes one
level further: its three adjacency colors recover vertices, edges, and
faces as three kinds of two-color components. The intrinsic identity
r_0*r_2=r_2*r_0 describes the flag squares. This is a precise unification
of incidence types, with the types still available for reconstruction.

Vertex-, edge-, and face-transitivity can therefore be compared as vertex
actions on associated structures, but use the same specified map-symmetry
group. The full abstract medial-graph automorphism group can be larger, and
Euclidean symmetry can be smaller. An irregular tetrahedron has the same
abstract incidence but may have no nonidentity Euclidean symmetry. The
symbols F(G) and E(G) in the quoted diagram remain undefined here; the
argument uses explicitly defined medial and flag constructions instead.

## 5. The prime density and the three-color register meet at exact local counts

Along the difference-d diagonal, gcd(a,a+d)=gcd(a,d). For d>0 the coprime
density is phi_E(d)/d, where phi_E is Euler's totient. Across the whole
positive square it is instead

    lim_(N->infinity) #{1<=a,b<=N:gcd(a,b)=1}/N^2
       = product_p(1-p^-2)=6/pi^2.

The same local factors count exact-denominator golden phases. If q runs
through primorials, J_2(q)/q^2 tends to6/pi^2. This has a stated sampling
rule: there is no such universal limit along every q->infinity (prime
moduli give a fraction tending to1).

The user's two gap-four examples reveal a second filter:

    7+11phi=sqrt(5)*phi^5, norm5, coefficient gcd1;
    8+12phi=4phi^4,        norm16, coefficient gcd4.

The same difference has different content and ideal. In general
N(a+(a+d)phi)=a^2-ad-d^2; along a fixed family this is a quadratic prime
filter, not a prime-producing theorem. Inert, split, and ramified primes
give different unit counts. This is why primitive coefficient pairs cannot
always be cancelled as ring units.

The earlier red/black/blue Zeckendorf construction with marked extra copies
of1 retains its original ordered grammar. Its arithmetic register maps by
an invertible integral change to O. Modulo2, O/(2)=F_4: the three nonzero
states form one M-cycle; zero is a fourth state. Under the new periodic-lift
theorem, those three nonzero phases lift to the golden cycle

    phi/2 -> beta/2 -> 1/2 -> phi/2.

Thus the phase quotient connects the three-color arithmetic register to
the root of the binary cycle tree. It does not recover the extra-unit
markers; those remain explicit. The three children in the ternary tower
have a different source, the three translation lines in F_3^2.

The inherited Fourier reader also retains its precise map. For
chi_(r,s)(a,b)=exp(2*pi*i*(ra+sb)/m) and D_d(v)=Mv+d(1,0),

    chi_(r,s)(D_d v)=exp(2*pi*i*r*d/m) chi_(s,r+s)(v).

This transports complete characters and their phase; in real coordinates
it rotates sine and cosine together. The previous mod3 guard hostile still
shows why forgetting phase can forget a legal operation. At modulus2 the
additive characters are already real signs, so that particular loss does
not arise. The shared register is exact; identifying all triples or all
three-colorings would discard the input grammar.

## 6. Quadratic pair dynamics identifies the diagonal's differential role

For F_c(x)=x^2+c, let s=x+y and d=y-x. Then

    d_next=s*d,     s_next=(s^2+d^2)/2+2c.

The difference update cancels c but needs s, whose update retains c. The
diagonal d=0 carries the transverse derivative2x; off it the secant is x+y.
For a nontrivial periodic orbit the product of adjacent secants is1, by
telescoping the nonzero differences. The derivative product need not be1.

The inherited -29/16 example has cycle -7/4,5/4,-1/4, secant product1 and
derivative product35/8. At parameter -7/4 the parabolic three-cycle has both
products1 and a real cubic equation of discriminant49. After z=x+1/2,
its points are 2cos(k*pi/7), k=1,3,5, and the restricted map is angle
tripling. The full scoped proof and exact controls are in the
[visibility and secant note](difference_visibility_norms_20261004.md), with
[THM-4146's order-six lift](../../01-canon/theorems/THM-4146-rational-three-cycle-order-six-lift-horizontal-divisor-fibre-firewall.md).
The cubic field of discriminant49 and the golden quadratic field of
discriminant5 remain distinct; no Collatz conjugacy has been supplied.

The rational anchor just recovered supplies a literal connection to the
user's second parameter: its raw period7 and four halvings give
N(phi^7-1)/16=-29/16. Its denominator ideal is exactly
(phi^7-1)=(Phi_7(phi))=(29,phi-24), another prime-ideal identification.
This exact identity survives, but a proposed general
parameter N(phi^(p+A)-1)/2^A forgets the ordered carry. Words1113 and1122
both give -121/64 while their rational cycles, through -65/17 and -73/17,
are disjoint. The identity is not a reconstruction of either dynamics;
the [companion probe](difference_visibility_norms_20261004.md) records the
precise surviving statement and what the collision does and does not rule out.

The literal current lambda, (5phi^4)^(-1/6), is about0.5548555133. The older
1/[5(2phi^4)^(1/6)] is about0.1292805634. They are different expressions.
Both fit the exact scheme C^2*x^12-7C*x^6+1=0 with different C, using
phi^4+phi^-4=7. The old value's proximity to the conventional tree-level
Higgs coupling about0.13 is numerical. The companion records the primary
physical convention and the previous rejected Fourier identification;
no lattice or density argument here predicts a Higgs mass.

## 7. Focused continuation after these results

The promising next object is a **typed phase-and-route state**: a denominator
ideal, numerator phase, upper-boundary tag, arithmetic source, and a pointer
to a certified ordered suffix. Each component repairs a demonstrated loss.
An implementation can intern common suffixes while keeping the source and
edge legality checks local. The existing carry register and route compiler
already provide those local certificates.

The first finite phase/route intersection and the negative-anchor transfer
are now implemented in the independent filter above. Three decisive
continuation tasks follow:

1. Seek a structural description of the arithmetic fixed-point divisibility
   test along one complete denominator tower, using the saved smallest
   successful and failed examples. Increasing a census alone does not give
   this all-height rule.
2. Attach the determinant branch address to this filter and determine which
   child operations preserve integer realization or a certified suffix.
   The first binary split proves that unrestricted lifting does neither.
3. Rank the48 rational-only negative anchors using the inherited compiler,
   measuring genuinely new certified inputs and retaining the original source threshold.
   A larger inventory or shared-prefix count is not a coverage proof.

For a global Collatz conclusion, two separate obligations remain: establish
the required eventual periodicity/rational golden coordinate for every input,
and exclude every other integer-realized cycle in the resulting unbounded
denominator space. The present theorems give exact finite models and recursive
families; they do not discharge those obligations.

## Transfer audit

| Source -> target | Preserved | Lost or required extra coordinate | Cheapest decisive test |
|---|---|---|---|
| Raw pair -> difference | Crossed arithmetic | Diagonal translation t | (7,11) versus (8,12) |
| Golden point -> ideal | Denominator ideal along B | Numerator phase orbit | -5 versus1/13 |
| Nonzero periodic point -> phase | Point and period, reconstructively | Zero-class tag; transient lift outside domain | 0,beta,1 |
| Phase cycle -> higher denominator | Golden dynamics | Integer realization | Two q4 rational children |
| Certified splice -> source integer | Actual route and root | Ordered head and first arrival | Cycle-padding hostile |
| Pair tiles -> overlap graph | Intersection relation | Endpoints at the n4 exception | Octahedron48 versus S4 order24 |
| Typed medial map -> untyped map | Edge addresses | Vertex-face distinction | Self-dual tetrahedron |
| Fixed-gap count -> square density | Local gcd filters | Sampling measure | d1 density1 versus6/pi^2 |

These are the connections that survived the session's hostile tests. The
recursive structure is now explicit at the level of phases and certified
splices, with their distinct domains recorded rather than inferred from
matching small numbers.
