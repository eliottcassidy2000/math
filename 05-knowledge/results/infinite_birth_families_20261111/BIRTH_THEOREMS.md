# Complete families for one-vertex additions

**PROVED below**, with **FINITE-EXACT** coefficient tables and an explicitly
**CITED** order cutoff. No literature-priority claim or canon promotion.

## 1. Objects, inheritance, and the exact quotient

A birth is an isomorphism class of a pair (U,v), with v distinguished and
T=U-v its source. For one fixed representative T, births are the orbits of
subsets O of V(T) under Aut(T), where v points to precisely O. We keep
ordinary isomorphism, not converse identification. Two roots in the same
unpointed child count separately exactly when they lie in different vertex
orbits. Thus a birth count is neither an unpointed child count nor H(U).

The intrinsic relation is the tournament's arc relation, without ties.
Vertices are the actual tournament vertices; O is a named attachment side.
The quotient preserves rooted isomorphism and loses vertex names. Literal
adjacency and the marked root restore construction. The fixed-order
blue/black complement-line metagraph is a different relation.

Inherited mechanisms:

- Unique ordered strong-component decomposition, H multiplicativity, and
  finite strong alphabets at fixed H: [Session24, PROOF](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md).
- Canonical transitive-module families and their full group actions:
  [Session67](../session67_metagraph_enumeration_information_20261026/README.md).
  A fixed finite core bank misses almost all tournament classes.
- Davis/Burnside odd-cycle enumeration, already used in
  [the connected-diagram compiler](../sequence_gf_20261110/burnside/README.md).
- The complete attachment-orbit census and quadratic path response:
  [the preceding aspect packet](../metagraph_aspect_sequences_20261111/README.md).
- The common exterior divisor is already J1 of
  [the prior joint-deletion proof](../decomposition16_odd_spectrum_20261110/induction/PROOF.md).
  The coprime strong-child reduction and source-padding compatibility are
  already P1 of [the 55/105 source proof](../chanpak_55_105_golden_20261110/targets/PROOF.md).
  Section 5 below gives their marked multiplicity formula, not a new claim
  to the underlying coprime reduction.
- Neutral births are exactly the k+1 SCC cuts by the repaired
  [THM-1900, insertion-response-calculus](../../../01-canon/theorems/THM-1900-the-insertion-response-calculus.md).
  Accordingly C_(h,h)=0 and B_(h,h)=D_h.

The new formulation here separates the interval an addition absorbs and
compiles fixed-H transition families. It does not reattribute the inherited
SCC decomposition, Burnside formula, path insertion, or rational fixed-H
class series as new results.

## 2. B1: every birth has exactly one active interval

Write T=S1 -> ... -> Sk for its nonempty strong components. Let i be the
first component meeting O, and j the last component meeting its complement.
If i<=j, the strong component containing v is exactly

    U[{v} union Si union ... union Sj].

All earlier and later components remain unchanged. If i>j, the new vertex
is a singleton strong component inserted at the unique cut i=j+1. This
includes the empty source as its separate singleton birth.

Proof. For i<j, v reaches all of Si, which reaches every later component
through Sj; every vertex of Sj reaches a vertex pointing into v. These
paths make the entire displayed interval strong. For i=j, the old component
is strong and has both attachment directions, giving the same conclusion.
All earlier components point uniformly into v and all later ones receive
arcs uniformly from v, so none can join this component. If i>j, no component
is mixed and the orientations form the asserted cut.

Put alpha(S)=number of Aut(S)-orbits on vertex subsets. Automorphisms fix
each position in the ordered SCC list, even when components are isomorphic,
so Aut(T) is the product of their automorphism groups. Consequently:

    singleton births:               k+1;
    interval [i,i] births:           alpha(Si)-2;
    interval [i,j], i<j:             (alpha(Si)-1)(alpha(Sj)-1)
                                    product_(i<l<j) alpha(Sl).

The endpoint factors exclude the empty or full subset as appropriate.
The middle components have arbitrary subsets. These disjoint cases sum to
product_i alpha(Si), the complete number of attachment orbits. In particular
the statement is about actual orbit counts, not independent choices after
discarding an unknown group action.

## 3. B2: all-size interval-width series

Let A(x)=sum a_n x^n count ordinary tournaments including the empty one,
S(x) count nonempty strong tournaments (including the singleton),
Q(x)=sum_T alpha(T)x^|T|, and R(x)=sum_strong S alpha(S)x^|S|.
The usual ordered-word identities give

    A=1/(1-S),            Q=1/(1-R).

Here Q_n also counts tournaments on n+1 vertices with one distinguished
vertex, up to rooted isomorphism. If z records the number of old strong
components absorbed, the complete birth series is

    B(x,z)=A(x)^2 [1+z(R-2S)+z^2(R-S)^2/(1-zR)].       (B2)

Source order is the x grading; the added vertex is not counted in x.
The two external SCC words give A^2. The three terms classify zero, one,
or at least two absorbed components. Setting z=1 gives Q exactly.

For an odd partition lambda of n, let l(lambda) be its number of parts,
z_lambda=product_p p^m_p m_p!, and
e(lambda)=sum_i(p_i-1)/2+sum_(i<j) gcd(p_i,p_j). Then

    A_n=sum_lambda 2^e(lambda)/z_lambda,
    Q_n=sum_lambda 2^(e(lambda)+l(lambda))/z_lambda.     (B3)

The second formula is ordinary Burnside on a tournament with a subset:
a permutation fixing the tournament has only odd cycles, and its invariant
subsets are unions of cycles. It therefore has 2^l subset choices. Our
integer implementation multiplies by n! before division; no rounding is used.

Equations B2-B3 evaluate the entire interval-width triangle through any
finite order without generating the tournament classes. The saved triangle
reaches 64. This is an extra suite of exactly computable sequences, not a
claim that this implementation beats the repository's faster A000568 engines.

## 4. B3: the fixed-path-count family theorem

Let A_h(x) count ordinary tournaments of path count h, with the empty
tournament allowed only for h=1. Let C_(h,g)(x) count pairs (U,v) with
U nontrivial and strong, H(U)=g and H(U-v)=h, weighted by x^(|U|-1).
Finally let B_(h,g)(x) count all births with source count h and child count g.
Define

    D_d(x)=sum_(ab=d) A_a(x) A_b(x).

Then, at every order,

    B_(h,g)=1[h=g] D_h
            + sum_(d | gcd(h,g)) D_d C_(h/d,g/d).       (B4)

Proof. Decompose the child into the SCC containing its root and the
ordered words before and after that component. Their combined path count
d multiplies both the source and child counts. A singleton root component
gives the first term. A nontrivial root component gives precisely one
summand in the second term. This is a bijection on rooted isomorphism
classes, so it retains multiplicities as well as existence.

Every C_(h,g) is a finite polynomial. An elementary bound suffices:
a strong tournament of order m has a Hamiltonian path starting at every
specified vertex, hence H>=m. To see the first assertion, start a path at
that vertex. Strongness supplies an edge from its current vertex set to
an outside vertex. Insert that outside vertex immediately after the last
path vertex that beats it. Its successor, if any, is beaten by the inserted
vertex. The first vertex stays fixed, and iteration spans the tournament.
Thus degree C_(h,g)<=g-1. Busch's stronger bound greatly reduces this cutoff.

For a fully explicit rational description, let P_s(x) count nontrivial
strong tournaments with H=s. Let W_(d,L)(x) count ordered words of L such
components whose H values multiply to d:

    W_(1,0)=1;
    W_(d,L)=sum_(s|d, s>=3) P_s W_(d/s,L-1).

All other length-zero words vanish. Since each nontrivial component has
at least three paths, L<=floor(log_3 d). Splitting such a word between
its left and right external parts gives L+1 cuts, hence

    D_d(x)=sum_L (L+1) W_(d,L)(x)/(1-x)^(L+2).        (B5)

The denominator counts the independent singleton runs; at the cut there
are two distinct runs. Equations B4-B5 are a complete, effective family
classification for every fixed pair (h,g). All atomic alphabets are finite.
They do not assert an efficient uniform way to enumerate those alphabets
as h and g grow without bound.

Every nonzero B_(h,g) is eventually polynomial in source order, of degree
at most floor(log_3 gcd(h,g))+1. More precisely its degree is one plus the
largest realizable nontrivial-component word length in any contributing
exterior d (including d=h in the identity term). Positivity prevents
leading-pole cancellation. The bound is sharp for h=g=3^r: the exterior
word of r directed triangles gives degree r+1.

## 5. Coprime counts, ternary targets, and source information

For h!=g and gcd(h,g)=1, equations B4-B5 reduce to

    B_(h,g)(x)=C_(h,g)(x)/(1-x)^2.                     (B6)

Each primitive pointed strong child therefore generates exactly n-m+1
births at source order n>=m, where m is its source order. These correspond
to all distributions of n-m neutral singleton vertices before and after
the root component. Distinct distributions remain distinct rooted types.

The inherited J1/P1 criterion says any such birth must absorb every nontrivial source SCC:
the untouched exterior H divides both h and g, so it is one; H=1
means transitive, whose SCCs are singletons. This is especially useful for
g=3h+2 and g=3h+4 with odd h, because both gcds are one.

Thus all-order existence of a one-vertex ternary step is exactly
C_(h,3h+2)!=0 or C_(h,3h+4)!=0, respectively. This is a finite test for
each h, not a proof that it always succeeds. Nor do two nonzero kernels
select the same actual source. Multiplying scalar transition tables loses
the middle tournament and does not compose actual births.

The inherited [different-source induction](../decomposition16_odd_spectrum_20261110/README.md)
requires both steps for every realized odd h>=95, plus its finite head and
join operation. Nothing here proves that universal antecedent. The
[source-preserving frontier theorem](../chanpak_five_core_bridges_20261110/algebra/LOW_RATIO_FRONTIER.md)
also remains distinct: it classifies multi-birth targets containing a fixed
strong source below a ratio threshold. Our formulas count one-birth types
and retain all sources; neither theorem substitutes for the other.

## 6. Complete coefficients below 125 and explicit infinite examples

**CITED cutoff.** Busch, Corollary 1, gives f(10)=125 and an increasing
minimum thereafter. Therefore every strong tournament of H<125 has order
at most nine. The [primary paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf)
was checked; the repository route is
[THM-4744, interval-chains-attain-busch](../../../01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md).

**FINITE-EXACT, complete at every order in this H band.** The two counting
directions agree on 687 ordinary strong classes, 5,086 pointed strong
classes, and 872 nonzero C kernels. B4-B5 then give 946 nonzero B families
for odd 1<=h<=g<125. Their eventual degree histogram is:

| degree | 1 | 2 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|---:|
| families | 744 | 161 | 34 | 6 | 1 |

Pairs absent from the complete table have no one-vertex birth at any order.
This is a classification of a bounded path-count band across unbounded
graph sizes, not an order-nine-only statement.

In particular:

    C_(33,101)=10x^6+2x^7+4x^8;
    C_(33,103)= 6x^6+2x^7;
    C_(11,35) = 2x^6;
    C_(13,41) = 2x^6.

Hence B_(33,101)[n]=16n-90 for n>=8 and
B_(33,103)[n]=8n-42 for n>=7. At n=100 these are 1510 and 758
ordinary pointed births. They are not counts of unpointed target classes
or literal blue/black lines. All coefficients at smaller n are supplied by
B6, with negative-degree binomial terms interpreted as zero.

Examples of genuine all-order one-birth exclusions in the same certified
band are 3->11, 15->47, 17->55, 35->107, and 35->109.
These are exclusions of a construction step; the target values themselves
are realizable. They do not refute the inherited threshold-95 induction.

The table stores positive terms [shift,pole,mass], meaning
mass*x^shift/(1-x)^pole. At any n the contribution is
mass*binom(n-shift+pole-1,pole-1), or zero if n<shift.
All 5,086 primitive births retain actual source rows and the selected
attachment mask in primitive_births.tsv. The family query returns a literal
primitive witness plus its neutral padding for coprime positive queries.

## 7. Complete common-source complexes in the same target band

The inherited P1 theorem also lets us retain the coordinate which a scalar
kernel loses. For h>1, remove only the exterior singleton components of
each primitive source, keeping its middle core K, its actual cut, and its
left/right padding requirements (a,b). Do not remove internal singletons.
For each K and target g, Pareto-minimal pairs describe exactly the upper
set of padding (A,B) where that target is available. Canonicalization of K
retains an explicit permutation and transports the cut.

For a set of coprime targets, take coordinatewise maxima of one requirement
pair for each target and Pareto-minimize again. The set has a common source
iff at least one K survives. Thus the maximal target sets for each K give
the facets of the complete common-source compatibility complex in this
band. This is an application of P1, with new complete coefficient and
witness tables, not a replacement of its source condition by a scalar.

If a surviving K has m vertices and frontier pairs
(a_1,b_1),...,(a_t,b_t) with a increasing and b decreasing, its ordinary
source-class generating function is

    x^m [sum_i x^(a_i+b_i)
         -sum_(i<t) x^(a_(i+1)+b_i)]/(1-x)^2.         (B7)

Proof. Split the upper set into disjoint strips a_i<=A<a_(i+1),
B>=b_i, and the final infinite strip. Sum the two geometric series.
For h>1, stripping exterior singleton SCCs is unique, so different K
or different padding pairs do not identify source classes. Each nonempty
core family has eventual slope one.

The exact atlas contains 200 middle-core types across 36 nontrivial
source counts, for all coprime targets below 125. It stores facets,
incompatible target pairs, minimal padding, normalized cuts, and actual
core rows. The simultaneous ternary-source series are complete through h=39.
For h=33 and targets 101,103 the result is

    common_source_33(x)=(2x^7+x^8)/(1-x)^2.

There are exactly three middle-core types: one of order six with minimal
padding (0,1) or (1,0), and two of order seven requiring respectively
(1,0) and (0,1). Consequently there are no common order-six sources,
two of order seven, and 3n-19 of every order n>=8. In particular there
are eleven at order ten. This counts actual sources supporting BOTH
alternative additions, unlike the separate birth counts 16n-90 and 8n-42.

For h=11 and targets 35,37 the corresponding series is
(2x^6-x^7)/(1-x)^2: exactly n-4 common sources for n>=6.
Thus seven occurs at source order eleven in this particular sequence.
These numerical roles follow from the displayed counting maps and do
not identify seven with a possible full Hamiltonian-path count.

All normalized minimal witnesses and one common-source construction per
ternary core are checked literally. A separate comparison with the complete
source census verifies every eligible source's full target-availability
set, not only its pairwise margins. The atlas's h>1 restriction is explicit;
for h=1 all sources are transitive and any finitely many available targets
can share a sufficiently long transitive source.
