# Sharp order ten stability from five vertex incidence

Status: **PROVED**, with an exhaustive integer graph certificate and a
same-agent algebra audit. An independent researcher audit remains desirable.

## The theorem

Let A be any real weighted tournament of order ten: Aii=0, Aij+Aji=1,
0<=Aij<=1. Let H be its weighted directed Hamilton-path sum, u=det(I+A),
v=(-1)^10 det(A), Delta=H-u-|v|, q_ij=Aij Aji, and V=sum q_ij.
Write Cj for the weighted sum of directed j-cycles, counted once up to
cyclic rotation, tau=C3, and M_ab... for unordered vertex-disjoint cycle
collections. A two-cycle has weight q. Then

    Delta >= (39/64)V + (8849/18432)tau
             + C4/32 + C5 + C6/8 + C7 + C8/2 + C9.          (1)

The coefficient 39/64 is optimal. Equality in Delta>=39V/64 with V>0
holds exactly for the ordinal sum of five fair two-vertex tournaments,
up to relabeling. With V=0, equality holds exactly for transitive tournaments.

Together with the inherited proofs at n=3..9 and n>=42, this leaves
**11..41 OPEN**. This is a weighted stability theorem, not a census of
order-ten tournament isomorphism classes.

Dependencies: the exact cycle expansions, triangle-variance lemma, and
five-vertex identity in
[SHARP_LOW_ORDER_STABILITY.md](../session138_stability_and_actual_contexts_20261102/stability/SHARP_LOW_ORDER_STABILITY.md);
cycle-complement grouping and the self-contained graph augmentation method in
[SHARP_NINE_GRAPH_CERTIFICATE.md](../session138_stability_and_actual_contexts_20261102/stability/SHARP_NINE_GRAPH_CERTIFICATE.md);
uniform margins u_r+v_r>=2^(1-r), u_r-v_r>=r*2^(1-r), also recorded in
[MULTISTAR_STABILITY.md](../session138_stability_and_actual_contexts_20261102/stability/MULTISTAR_STABILITY.md).
These mechanisms are inherited; the order-ten grouping and certificate are
the extension here.

## Exact signed ledger

Put G_s=H-u+s*v for s=+1,-1, and

    P=V-M22+M222-M2222.

For every cycle partition lambda into parts at least two, sum(lambda)<=10:
the coefficient in u is (-1)^(number of even parts); the coefficient in H
is 2^length(lambda) if all parts are odd, and zero otherwise. Empty
collections contribute one to both. The coefficient in v is zero unless
sum(lambda)=10, when it is (-1)^length(lambda).

These rules specify all 42 collections, including the empty one. They
avoid assuming a sign for det(A). The compiler in [audit.py](audit.py)
subtracts every group below from each exact ledger and checks that every
unused coefficient is nonnegative. The remaining terms are printed in
[audit.json](audit.json).

## Four cycles with six vertex complements

Sum over a selected directed four-cycle C, and let R be its complementary
six-set. The complete terms involving a four-cycle are

    E4,- = sum_C w(C)[u_R+v_R+C4_R/2-M24_R],
    E4,+ = sum_C w(C)[u_R-v_R+C4_R/2].                     (2)

Explicitly,

    E4,- = C4-M24+M34-M44+M45-2M46+M224-M234
           +2M244-2M2224+2M334,
    E4,+ = C4-M24+M34-M44+M45+M224-M234.

The factors matter: selecting a four-cycle in M44 or M244 counts twice.
For R of order six, M24_R<=C4_R/4: every selected four-cycle leaves one
possible two-cycle of weight at most 1/4. Thus E4,->=C4/32 and
E4,+>=3C4/16. The common bound is C4/32.

## Six and eight cycles

After removing (2), the six-cycle terms are

    E6,- = C6-M26+M36+2M226,
    E6,+ = C6-M26+M36.

On the complementary four-set R their brackets are respectively
1-V_R+tau_R+2M22_R and 1-V_R+tau_R. Since

    u_R=1-V_R+tau_R-C4_R+M22_R,
    u_R>=5/16,  M22_R<=3/16,

we get E6,->=5C6/16 and E6,+>=C6/8. The u_R bound is the average of
the two uniform four-vertex margins. The common eight-cycle group is
C8-2M28>=C8/2; inserting the negative term also into G_+ only weakens it.

## The five vertex sum of squares pays the five cycle coupling

For a five-vertex remainder R the inherited exact identity is

    M22_R = V_R/8 + M23_R - S_R/8,
    S_R = sum_e q_e (Kij+Kjk+Kki)^2 >=0,                (3)

where {i,j,k}=R minus the endpoints of e and K=A-A^T. Also
M23_R<=tau_R/4, because each triangle leaves a single two-cycle.
Both signed ledgers contain

    E5=C5+M25+3M35-M225
      =sum_C5 w(C)[1+V_R+3tau_R-M22_R]
      >=C5+(7/8)M25+(11/4)M35 >= C5.                 (4)

This is where five-vertex stability contributes an exact joint incidence
identity, stronger than simply bounding M22_R by its box maximum 15/16.
No term used in (4) was spent in (2) or the six-cycle groups.

## A positive triangle reserve

The common triangle group can be taken as

    E3=tau+3M33-M223-2M2233.

For each pair of triangles, the complementary four-set has M22<=3/16,
so M2233<=3M33/16. Consequently

    E3 >= tau+(21/8)M33-M223
       =sum_C3 w(C)[1+(21/16)tau_R-M22_R],             (5)

where R now has seven vertices. On such an R, let W count pairs of
incident distinct edges with product q_e q_f. The inherited local lemma
gives tau_R>=W_R, while Cauchy-Schwarz on vertex variance degrees gives

    W_R >= (2/7)V_R^2-V_R/4,
    M22_R <= V_R^2/2-W_R.

Therefore

    M22_R-(21/16)tau_R
      <= -(9/56)V_R^2+(37/64)V_R
      <= 9583/18432.

The last inequality is an exact completed square. Equation (5) gives
E3>=(8849/18432)tau. The unused M233, M2223 and 7M333 terms are positive;
their presence and all full-cover signs are checked by the ledger compiler.

Combining the groups, with the unused positive C7 and C9, yields

    Delta >= P+M23+(8849/18432)tau
             +C4/32+C5+C6/8+C7+C8/2+C9.              (6)

## The integer graph certificate

Let P21(q) count a two-edge wedge and a vertex-disjoint edge. The triangle
lemma on the complementary triples gives M23>=P21. Define

    F(q)=P+P21-(39/64)V.

Every monomial uses distinct edges, so F is separately affine on
0<=q_e<=1/4. At a box vertex corresponding to an undirected graph G,
write m,m2,m3,m4 for its edge and matching counts, and p21 for a wedge
and disjoint edge. Then

    256F = 25m-16m2+4m3-m4+4p21.                    (7)

[graph10.cpp](graph10.cpp) verifies (7) on every graph of order ten up
to relabeling, using no external graph census:

1. Explicitly partition all 2^21 labeled seven-vertex graphs into 1,044
   permutation orbits, checking coverage of every one of the 2,097,152 graphs.
2. Attach all eighth-vertex neighbourhoods and canonicalize to 12,346 graphs.
3. Attach all ninth-vertex neighbourhoods and canonicalize to 274,668 graphs.
4. Attach all 512 tenth-vertex neighbourhoods to each, checking 140,630,016
   augmentations. This is exhaustive coverage, not an assertion of that many
   nonisomorphic ten-vertex graphs.

Canonicalization first refines canonical vertex colours, then enumerates
EVERY permutation within the remaining equal-colour classes. Refinement
alone is not used as an isomorphism test. Induction by deleting the last
vertex proves coverage at each augmentation stage. The minimum of (7) is
zero; its only zero types are the empty graph and a perfect matching.
See [graph10.json](graph10.json).

To make this feasible, the update on a new neighbourhood S is quadratic.
For a nine-vertex base G with degrees d_i, adjacency indicators Aij and
wedge count W, let m_j(G-i) be matching counts after deleting i. Put

    Qij = m-d_i-d_j+Aij,
    ell_i = W(G-i)+d_i*(m-d_i+1)-sum_(j adjacent i) d_j,
    Li = 25-16m_1(G-i)+4m_2(G-i)-m_3(G-i)+4ell_i.

Then the integer expression in (7) on G+S is

    Fint(G+S)=Fint(G)+sum_(i in S) Li+4 sum_({i,j} in S) Qij.  (8)

Here Qij is the nonnegative number of edges outside {i,j}. The linear
part of p21 records the new vertex as an endpoint of the disjoint edge
or as a leaf of the wedge; its quadratic part records the new vertex as
the wedge centre. Matching updates use
m_j(G+S)=m_j(G)+sum_(i in S)m_(j-1)(G-i).
The program compares (8) with independent direct subset matching counts
and complementary-degree wedge counts on 6,200 augmentations.

The compiler has an all-size version. For arbitrary coefficients a_j,b,c
and f(G)=sum_(j>=1) a_j m_j(G)+b*p21(G)+c*W(G), its new-neighbourhood
polynomial has linear coefficient

    sum_j a_j m_(j-1)(G-i)+b*ell_i+c*d_i

and quadratic coefficient b*Qij+c. No higher degree occurs. Given the
deletion matching responses, all 2^n extensions are evaluated in
O(n*2^n) arithmetic operations. This is a certificate-polynomial compiler,
not a canonical augmentation algorithm for tournament isomorphism classes.

Separate affinity extends (7) to the entire real q-box. Substituting it
into (6) proves (1).

## Equality and sharpness

Equality in the unreserved bound forces tau=0, hence W=0. All positive-q
edges are therefore a matching. On r such edges the endpoint values of
256F for r=0..5 are 0,25,34,31,19,0. The multilinear interpolation of F
is an expectation over independently retained matching edges. If V>0,
zero expectation requires exactly five edges, each retained with
probability 4q_e=1. Thus all five pairs are fair.

With no cyclic triple, these pairs form transitive ordinal blocks. Indeed
an outside vertex must see both ends of a fair pair in the same direction,
or one positive-weight cyclic triangle occurs; the block quotient has no
cyclic triple. Conversely five fair ordinal pairs have H=1,
u=(3/4)^5, |v|=(1/4)^5, V=5/4, attaining 39V/64.

## The next order is a different graph problem

At eleven, the analogous bare matching polynomial includes +M22222.
Using the proposed sharp coefficient 781/1280, the endpoint expression
for P+P21-kappa V is

    5120F11=499m-320m2+80m3-20m4+5m5+80p21.

It is -56 on P3 disjoint union 4K2, so the order-ten graph reduction
does not extend unchanged. Adding alpha W requires alpha>=7/40 on this
witness. This is an obstruction to that intermediate certificate, not
to the order-eleven weighted stability conjecture. Its triangle reserve
could still pay the extra wedge term; the complete order-eleven signed
ledger and a universal graph certificate remain OPEN.
