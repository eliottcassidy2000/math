# Odd-zeta moduli spaces, Petersen factors, and the K5 intersection lattice

**2026-10-19. Status: elementary PROVED extensions, with FINITE-EXACT checks.**
No new irrationality theorem, polytopality theorem, or priority claim is made.
The facts about stable-curve strata and their Chow rings are CITED; the explicit
maps, counting formulas, tensor calculations, and failure tests below are proved
here from those inputs. Independent review is still required before promotion to
repository canon.

**Inheritance.** This extends
[`odd_zeta_parallels_petersen_20261001.md`](odd_zeta_parallels_petersen_20261001.md),
§3.4. The recent graph input is
[THM-4701, Kuratowski graphs and finite unit-circle graphs](../../01-canon/theorems/THM-4701-kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs.md).
Its exclusion of K5 and K3,3 is a **subgraph** exclusion, not a minor exclusion.
The earlier note's corrected distinction between dinner arrangements and convergent
configurations is retained. The concept board is: stable-tree gluing; graph
products/minors; divisor intersection lattices; cellular-form valuations; period
products versus primitive odd-zeta targets.

## 1. Main conclusions

1. A stable-tree stratum in M̄0,N with r factors M̄0,5 must satisfy N ≥ 3r+2.
   The maximum r at N=8,10,12 is therefore **2,2,3**, not 2,3,4.
2. Pure Petersen-product strata exist at all three sizes:
   M̄0,5² in M̄0,8 and M̄0,10, and M̄0,5³ in M̄0,12.
   Their natural corner-adjacency graphs are P□P, P□P, and P□P□P.
3. There are three different graphs/complexes: Cartesian corner adjacency,
   strong-product corner intersection, and the join boundary dual complex.
   Conflating them loses dimensions and higher intersections.
4. The K5 connection can retain **the entire integral divisor intersection
   lattice**, rather than only graph nonplanarity. If B is the unsigned
   vertex-edge incidence matrix of K5 and A is Petersen adjacency, then

       Q = A - I10 = J10 - BᵀB = Bᵀ((1/4)J5 - I5)B.

5. Tensor products of this intersection lattice have rank 5^r and signature
   ((5^r+(-3)^r)/2, (5^r-(-3)^r)/2). They do not supply an irrationality criterion.
6. For Brown's actual odd-zeta configurations, the displayed Petersen-product
   strata lie on cell boundaries where the integrands have **zero logarithmic
   residues**. One cannot obtain the original odd-zeta integral by simply taking
   residues onto those products.

## 2. How many Petersen factors fit in a stable tree?

Let T be a stable tree with N labelled leaves and k internal vertices. Write d_v
for the number of incident flags at v, including leaves, so d_v≥3. Its closed
stratum is

    X_T ≅ ∏_v M̄0,d_v,
    codim(X_T in M̄0,N) = k-1,
    dim X_T = Σ_v(d_v-3).

These are the standard gluing descriptions of genus-zero boundary strata. Summing
valences and using the k-1 internal edges gives the exact identity

    N = 2 + Σ_v(d_v-2).                                    (1)

If r vertices have valence five, each contributes three to the right side and
every other vertex contributes at least one. Hence

    N ≥ 3r+2; equality holds exactly when every internal vertex has valence five.

For any r≥1 the bound is attained: take a path of r five-valent vertices and add
leaves to fill their remaining flags. Any surplus N-(3r+2) can be supplied by
inserting trivalent vertices into edges and putting one extra marked leaf on each.
Thus floor((N-2)/3) is the **sharp maximum number of five-valent factors in the
canonical stable-tree factorization**. This is not an exclusion of every possible
abstract graph embedding or non-boundary subvariety.

Explicit labelled examples:

- **N=8:** two components with four marks each, joined at a node. One divisor,
  X_T=M̄0,5², dimension four, codimension one.
- **N=10:** a four-component chain of valences (5,3,3,5), with marking blocks
  {1,2,3,4}, {5}, {6}, {7,8,9,10}. Then X_T=M̄0,5², dimension four,
  codimension three. A shallower chain (5,4,5), with block sizes (4,2,4), gives
  M̄0,5²×M̄0,4 in codimension two.
- **N=12:** a central unmarked trivalent component with three five-valent tails,
  bearing {1,2,3,4}, {5,6,7,8}, {9,10,11,12}. Then X_T=M̄0,5³,
  dimension six, codimension three.

The dimensions N-3 are 5,7,9, the natural weight bounds for the odd-zeta spaces
under discussion. The tree bound concerns boundary factorization, not arithmetic
independence at those weights.

### Exact census by valence profile

Suppose m_d vertices have valence d, k=Σm_d, and N=2+Σ(d-2)m_d. The number of
stable trees with N individually labelled leaves and this internal valence
multiset is

    (N+k-2)! / ∏_d [m_d! ((d-1)!)^m_d].                    (2)

Proof: temporarily label all internal vertices. The degree-sequence form of the
Prüfer bijection counts (N+k-2)!/∏_v(d_v-1)! trees. Divide by the permutations of
internal vertices with equal valences. This action is free: an automorphism of a
stable tree fixing every labelled leaf fixes every vertex.

| N | Internal valences | Product | Codimension | Strata |
|---|---|---|---:|---:|
| 8 | 5,5 | M̄0,5² | 1 | 35 |
| 10 | 5,5,4 | M̄0,5²×M̄0,4 | 2 | 5,775 |
| 10 | 5,5,3,3 | M̄0,5² | 3 | 51,975 |
| 12 | 5,5,5,3 | M̄0,5³ | 3 | 525,525 |

The central-trivalent stars at N=12 form only 12!/(4!³3!)=5,775 of the last
525,525 strata. Counting those stars as all Petersen-cube strata would be wrong.
The script checks (2) independently by enumerating internal tree skeletons and
allocating labelled leaves.

## 3. Three products that must not be interchanged

Let S=M̄0,5 and label its ten boundary curves L_a by a∈binom([5],2). Distinct
L_a,L_b intersect in one point iff a∩b=∅; there are no triple intersections.
Their graph is P=KG(5,2), the Petersen graph.

In X=S^r, define a **corner cycle**

    C_a = L_a1 × ... × L_ar,       a=(a1,...,ar).

It has dimension r, half the dimension of X. For distinct a,b:

- C_a∩C_b is empty if some coordinate has unequal, nonadjacent labels;
- otherwise dim(C_a∩C_b)=r-h(a,b), where h counts unequal coordinates.

Consequently:

- adjacency by intersection dimension r-1 gives **P□r**;
- adjacency by any nonempty intersection gives **P⊠r**;
- the full divisor boundary of X has dual complex **P * ... * P**, a join.

The last complex has 10r vertices, not 10^r. Its faces choose, in each factor,
either no curve, one curve, or an intersecting pair. Its face polynomial is
(1+10t+15t²)^r. Since P is a connected graph with first Betti number six,
the join has reduced integral homology Z^(6^r) in degree 2r-1 and zero in other
degrees. Indeed it is homotopy equivalent to a wedge of 6^r spheres of that degree.

| r | Corner cycles | Cartesian edges | All intersection edges | Dual-complex vertices | Top homology rank |
|---:|---:|---:|---:|---:|---:|
| 1 | 10 | 15 | 15 | 10 | 6 |
| 2 | 100 | 300 | 750 | 20 | 36 |
| 3 | 1,000 | 4,500 | 31,500 | 30 | 216 |

For r=3 the 4,500 Cartesian adjacencies are intersections along **surfaces**, not
curves. There are also 13,500 pairs meeting along curves and 13,500 at points.
The boundary dual complex here is the boundary of the **product stratum** X; it
is not the entire dual complex of M̄0,N containing X.

## 4. A faithful K5–Petersen bridge: the intersection lattice

Petersen vertices are also the edges of K5. Let B be the 5×10 matrix with column
b_ij=e_i+e_j. For distinct edges of K5, the off-diagonal entries of BᵀB are one
when the edges meet and zero when they are disjoint; its diagonal entries are two.
Therefore

    A = J10 + I10 - BᵀB.                                 (3)

The boundary curve self-intersections are -1. This follows directly from a Keel
relation: L12+L34=L13+L24; intersection with L12 gives L12²+1=0.
Thus the actual divisor intersection Gram matrix is

    Q=A-I10=J10-BᵀB.

Because every column of B sums to two, this also gives

    Q=BᵀHB,     H=(1/4)J5-I5.                            (4)

B has rank five; H has signature (1,4). The familiar Petersen eigenvalues
3,1,-2, with multiplicities 1,5,4, become the intersection eigenvalues
2,0,-3. Hence ker Q=ker B, of rank five.

### Integral, not merely rational

The image of B is

    Λ={u∈Z5 : Σu_i is even},
    h(u,v)=(Σu_i)(Σv_i)/4 - Σu_iv_i.                    (5)

Every column has even sum. The five columns indexed by
12,13,23,14,15 have determinant -2, so they form a basis of Λ. In particular h
is integral, odd, unimodular, and of signature (1,4).

The Keel relations are precisely balanced four-edge relations in K5, e.g.

    [12]+[34]-[13]-[24]=0.

They generate ker B over Z. One elementary elimination is

    [24]=[14]+[23]-[13],   [34]=[14]+[23]-[12],
    [25]=[15]+[23]-[13],   [35]=[15]+[23]-[12],
    [45]=[14]+[15]+[23]-[13]-[12].

All five eliminated edges occur with coefficient one. Thus no unverified
saturation assumption is needed. Keel's divisor presentation now identifies
Pic(S), with its intersection form, with (Λ,h).

This is a connection with a precise contract:

- source: the ten labelled boundary divisor classes on S;
- target: even-sum integral vertex vectors of K5;
- map: [Lij]↦e_i+e_j;
- preserved: every integral linear divisor relation and every intersection number;
- lost by the bare graph: the self-intersection diagonal and which cycle/chain
  is used for an actual period integral;
- restoring data: Q and, for periods, polar/cell divisors, forms and integration
  chains.

### Five forgetful fibrations and the E4=A4 root system

The K5 incidence map is geometric. Let f_i:S→M̄0,4 forget marking i and let F_i
be its fibre class. For a partition {a,b}∪{c,d}=[5]\{i}, the corresponding
reducible fibre is L_ab+L_cd. The three pairings of the four labels give the same
class by Keel's relations. Under (5),

    F_i ↦ 1-e_i,       L_ab · F_i = 1 exactly when i∈{a,b}.

Thus B records the intersections of the ten boundary curves with the five actual
forgetful fibrations. Petersen's fifteen edges split into five sets of three,
colored by the omitted label; each set consists of the three reducible fibres
of one f_i. This retains algebraic information that a K5 minor alone does not.

The anticanonical class κ=-K maps to 2·1. To check this without choosing a
blowup model, adjunction gives (-K)·L_ab=1 for every boundary (-1)-curve, while
h(2·1,e_a+e_b)=1; the ten curves generate Pic(S) and its pairing is nondegenerate.
Consequently

    κ²=5,   h(u,κ)=(Σu_i)/2,
    κ⊥={u∈Z5:Σu_i=0}, with pairing -u·v.

The twenty roots are e_i-e_j, i≠j, and simple roots e_i-e_(i+1) have the negative
A4 Cartan matrix. This identifies the degree-five del Pezzo root system **E4=A4**
directly in the Petersen/K5 lattice. It does not identify this E4 notation with
an Eisenstein series. The degree-six/degree-four neighboring cases belong to
the separate E-system atlas.

### Tensor consequences

The r-dimensional corner cycles in S^r have Gram matrix Q^⊗r, by the product
intersection formula, including all repeated-coordinate self-intersections.
Their class span is Pic(S)^⊗r. Thus:

    rank = 5^r,   formal-coordinate nullity = 10^r-5^r,
    signature = ((5^r+(-3)^r)/2, (5^r-(-3)^r)/2).

At r=2 these are rank25, nullity75, signature(17,8); at r=3 they are rank125,
nullity875, signature(49,76). These signatures refer to the induced nondegenerate
form on the corner-class span. There are additional middle-dimensional classes:
the full middle Chow rank of S^r is the coefficient of t^r in (1+5t+t²)^r,
namely 27 at r=2 and 155 at r=3. The extra M̄0,4 factor in the N=10 codimension-two
stratum is not silently included in this middle-pairing statement.

## 5. K5 and K3,3 as explicit minors, not subgraphs

Use standard Petersen vertices a_i,b_i, i mod5, with edges a_i a_(i+1),
b_i b_(i+2), and a_i b_i.

- Contract each spoke a_i b_i. The outer pentagon and inner pentagram supply
  complementary edge sets on five resulting vertices, giving exactly K5.
- Delete a0 and suppress the degree-two vertices a1,a4,b0. The resulting K3,3
  has parts {a2,b3,b4} and {a3,b1,b2}. Explicit connected branch sets are

      {a2,a1}, {b3,b0}, {b4,a4}; {a3}, {b1}, {b2}.

P contains neither K5 nor K3,3 as a subgraph. The first is impossible by degree;
the second would use all three incident edges at each of six vertices and hence
be a proper connected component of the connected Petersen graph.

Cartesian products preserve these minor witnesses: take products of connected
branch sets. Thus P□r has both K5□r and K3,3□r as minors. The script explicitly
checks the K5-product quotient through r=3. These minors are combinatorial
contractions, not asserted morphisms of algebraic surfaces and not period maps.
They are compatible with THM-4701's subgraph exclusions in the residue graphs.

## 6. The missing sidecar: cellular valuations and a zero-residue test

Let δ,δ′ be cyclic orders and let r_δ(S) be the number of cyclic runs of a stable
subset S. Its complement has the same number of runs. Brown's valuation formulas
(3.9) and (3.12) become

    ord_DS f_(δ/δ′) = r_δ′(S)-r_δ(S),
    ord_DS ω_δ′ = r_δ′(S)-2,
    ord_DS(f_(δ/δ′)^k ω_δ′)
       = k(r_δ′(S)-r_δ(S))+r_δ′(S)-2.                   (6)

For verification, the number of same-side cyclic edges is N-2r_δ(S). Substituting
half that number into Brown's I_D gives (6) immediately. This supplies an exact
integer decoration on every boundary vertex; it is stronger than its unlabeled
intersection graph.

Take δ=(1,...,N) and the actual odd-zeta orders from Brown's Definition7.1:

    N=8:  (8,2,7,3,6,4,1,5),
    N=10: (10,2,9,3,8,4,7,5,1,6),
    N=12: (12,2,11,3,10,4,9,5,8,6,1,7).

Use the explicit strata of §2. Their defining splits are:

- N=8: {1,2,3,4};
- N=10: {1,2,3,4}, {1,...,5}, {1,...,6};
- N=12: {1,2,3,4}, {5,6,7,8}, {9,10,11,12}.

Their respective orders are 2k+1; (3k+2,3k+2,3k+2); and
(3k+2,k,3k+2). All are δ-cell facets, and all these orders are nonnegative for k≥0.
Every ordinary logarithmic residue along these nodes is consequently zero.
This directly refutes the proposed shortcut “the Petersen product boundary
factors the odd-zeta cellular integral by taking its boundary residue.” One may
study different polar strata, regularized degenerations, or different forms,
but those require additional arguments.

The script checks all stable splits: 119 at N8, 501 at N10, and 2,035 at N12;
it verifies convergence on every cell facet and prints the exact affine orders
on the displayed strata. This is a finite check of these specified orders, not
an arithmetic approximation estimate.

## 7. Product periods are a different operation

Nodal gluing of r five-pointed components needs 3r+2 marks. By contrast, the
product of r independent two-dimensional periods can be described by forgetting
marks with three common anchors: it uses 2r+3 marks. The preimage of the product
of real cells is generally a union of shuffle cells. Integration over that union
gives the product period. A single cellular-integral identity needs Brown's
additional multiplication hypotheses.

For example, products ζ(2)^r occur in this way, and ζ(3)ζ(2)^r has natural
ambient size N=2r+6. Thus the same ambient sizes N8,N10,N12 can support the
product periods ζ(3)ζ(2), ζ(3)ζ(2)², ζ(3)ζ(2)³. These are not identities with
ζ(5),ζ(7),ζ(9). Graph products, nodal gluing, and period multiplication have
different maps and different marking budgets.

The compact intersection lattices above are algebraic, pure-Tate data. Odd-zeta
periods require open/relative data, including extension and comparison
information. An irrationality argument still needs a sequence of nonzero linear
forms with controlled rational coefficients, a legitimate denominator-clearing
bound, and analytic decay that beats that bound. Neither a Kuratowski minor nor
the signature of Q^⊗r supplies those conditions. This note establishes no
irrationality result for any of ζ(5),ζ(7),ζ(9), individually or collectively.

## 8. Reproduction and the next valid computation

The companion [`odd_zeta_boundary_graphs_20261019.py`](../../04-computation/experiments/odd_zeta_boundary_graphs_20261019.py) uses Python3 and SymPy;
its matching `.out` ends with `ALL CHECKS PASSED`. It tests both positive
constructions and hostile distinctions: three graph products; integral relation
saturation; extra middle classes; total-profile versus star-only counts;
subgraphs versus minors; and zero residues on the selected strata.

Run from the repository root:

    python3 04-computation/experiments/odd_zeta_boundary_graphs_20261019.py

A defensible continuation is to attach the complete valuation vector (6), then
actual residue/restriction maps, to selected Brown configurations. Compute the
resulting relative-cohomology weight pieces before any claimed elimination of
unwanted zeta coefficients. Only after obtaining a useful linear-form family
should one compare decay and denominator growth. The current exact geometry
reduces bookkeeping and prevents false operations; it does not replace that
arithmetic step.

## Primary sources

- Sean Keel, *Intersection Theory of Moduli Space of Stable N-Pointed Curves of
  Genus Zero*, Transactions AMS330(1992),545–574: boundary product description,
  p548; Facts1–2, pp551–552; Chow-ring relations, pp550 and571ff.
  [Primary article PDF](https://www.math.colostate.edu/~renzo/teaching/Topics10/keel.pdf).
- Francis Brown, *Irrationality proofs for zeta values, moduli spaces and dinner
  parties*, [arXiv:1412.6508](https://arxiv.org/abs/1412.6508),
  [primary PDF](https://arxiv.org/pdf/1412.6508): §2.2 for strata; equations(3.9),
  (3.12) for valuations; §6 for period-product maps; Definition7.1 and
  Proposition7.2 for the odd configurations; §§8 and11 for the relative motive
  and the distinction between geometry and irrationality. No claim of current
  best irrationality bounds is imported from this historical paper.
