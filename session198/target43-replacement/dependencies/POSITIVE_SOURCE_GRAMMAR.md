# A finite positive source grammar and a physical halving descent

Session 198. **PROVED:** every nontrivial start-one source has a three-letter
positive arithmetic presentation. The complete C2 and C3 bilinear forms become
strictly positive matrices in these coordinates. Removing at most two source
letters halves a witnessed C2 value; removing at most three letters halves a
witnessed C3 value, apart from finitely many short-source boundary cases.

**OPEN:** these are descents of existing representations. They do not prove
that every prescribed odd integer has a representation. In particular, no
universal odd-coverage conclusion is made.

## 1. Sources and the missing distinction

We inherit the complete physical source grammar and two-source reductions from
`dependencies/BILINEAR_C2.md` and `dependencies/C3_TWO_SOURCE.md`.
A source is a nonempty tournament V with a root c having exactly one Hamilton
path starting at c. Its response is

    (D,E,U,d) = (H(V), e_c(V), R_c(V), H(V-c)).

Here R counts ordered two-path covers whose second path is nonempty and starts
at c; the first path may be empty. H(empty)=1. A source is nontrivial if it has
at least two vertices. For a nontrivial source the standard injections give

    E <= d,   D-E >= d,   U-D >= d-E.                 (1)

For completeness, deletion of c injects endpoint-c paths into paths of V-c.
Insert c before the first vertex it beats in any path of V-c to prove the second
inequality. Such a vertex exists because c starts a spanning path and V is
nontrivial. Cutting each Hamilton path before c gives D distinct reserved
covers. Of the d covers whose reserved path is singleton c, precisely E are in
that image. Thus there are at least d-E additional covers, proving the third.

Set

    y=(a,f,e,p)=(D-E-d, U-D+E-d, E, d-E).            (2)

All four coordinates are nonnegative. The inverse is

    D=a+2e+p, E=e, U=a+f+2e+2p, d=e+p.              (3)

The change is unimodular: (2)-(3) are inverse integer linear maps. The
nontrivial two-vertex source is t=(1,0,2,1), with y(t)=(0,0,0,1).
The singleton has y=(-1,0,1,0), so it is deliberately excluded from the
nonnegative cone. Nonnegative coordinates do not certify physical membership.

## 2. Three one-vertex physical operations

From a nontrivial source (V,c), use these actual operations:

* O: append a universal sink and retain c.
* Z: append a vertex beating only c, and retain c.
* J: append the same vertex beating only c, and mark that newborn as root.

All preserve the condition of having one root-start Hamilton path. For Z, a
path starting at the old root cannot use the newborn's sole outgoing edge,
so the newborn is last. For J, the path must begin with newborn,c and then
follow the old unique root-start Hamilton path. The ordinary path-cover split
gives the response maps

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).                  (4)

For Z, covers with the newborn on the unreserved path or the reserved path
give two copies of the old U covers, with exactly the d singleton-root covers
overcounted; this gives 2U-d. The other identities are the same one-seam split
used in the inherited source theorem. They are also checked on every labelled
rooted source through order five by direct subset enumeration.

In coordinates (2), (4) becomes the nonnegative linear grammar

    O(a,f,e,p)=(a+e, a+2f+e+2p, 0, e+p),
    Z(a,f,e,p)=(2a+f+e+p, f, e+p, 0),
    J(a,f,e,p)=(f+p, 0, a+e+p, e).                 (5)

**Theorem 1.** The orbit of (0,0,0,1) under (5) is exactly the set of
nontrivial physical source responses, with an actual source recipe for every
orbit point.

Soundness follows from the constructions. For completeness, the inherited
grammar starts at the singleton and uses O and indexed ordinary gadgets M_L.
An M_L gadget consists of a new root and a transitive tail. Construct its root
by J, then append the tail in its transitive order. Each tail vertex either
loses to the current root (use O) or beats that root and no other earlier
vertex (use Z), retaining the current root throughout. These are exactly the
gadget's prescribed arcs. Thus every M_L is a finite J/O/Z word. The first
vertex extension of the singleton is t, whether it is an O or a J. Every
nontrivial source therefore comes from t. This is a graph construction, not an
assumption that the integer cone is inverse-closed.

In particular every physical point lies on at least one of e=0, p=0, f=0,
according to its final operation. This is necessary, not sufficient.

## 3. An exact finite inverse for the positive source grammar

For a candidate y=(a,f,e,p), reject negative coordinates, d=e+p<1, U<2,
or even D. These are necessary physical conditions (Hamilton-path counts in
tournaments are odd, also preserved directly by the source grammar).
Accept t separately.
The following lists every possible arithmetic predecessor; retain only
nonnegative candidates with d at least one and U at least two, and recursively require physical
membership.

* O is possible only if e=0 and f-a is even. For integers q with
  0<=q<=min(a,p), its predecessors are

      (a-q, (f-a-2p)/2+q, q, p-q).

* Z is possible only if p=0 and a-f-e is nonnegative and even. For
  0<=q<=e, its predecessors are

      ((a-f-e)/2, f, q, e-q).

* J is possible only if f=0. For 0<=q<=min(a,e-p), its predecessors are

      (e-p-q, a-q, p, q).

These formulas simply solve (5). In particular a predecessor's fourth
response is not invented or silently forgotten. Multiple candidates are
necessary because a one-vertex move can discard a coordinate.

Every accepted nontrivial predecessor has strictly smaller U. In the old
coordinates the increments are U, U-d and D for O,Z,J, respectively, and
U>=2d>=2. Hence this recursion terminates and, by Theorem 1, is a sound and
complete physical source membership test with witness reconstruction.

There is a logarithmic depth bound. Put u=(1,1,2,2), so U=u y. For any two
successive letters X then Y, direct multiplication gives u Y X >= 2u
entrywise. The nine difference rows uYX-2u, in order OO,OZ,OJ,ZO,ZZ,ZJ,JO,JZ,JJ,
are

    (2,2,4,4), (2,2,3,3), (1,0,2,1),
    (2,2,2,2), (2,2,1,1), (2,1,2,2),
    (2,0,4,2), (1,0,2,1), (2,0,3,2).

Thus a word of length k has U>=2^(1+floor(k/2)). This bounds the depth, not
the branching count or bit complexity. No polynomial-bit-time claim follows.

The hostile formal tuple (D,E,U,d)=(1,0,3,1) has y=(0,1,0,1).
It satisfies (1), but O fails its parity test and neither Z nor J is possible.
The inverse rejects it. More generally D=1 forces a transitive source, whose
root is first and U=2^(|V|-1); hence arbitrary formal U at D=1 is invalid.

## 4. Positive C2 forms and sharp one-step growth

Let y,z be two actual nontrivial source vectors. The complete C2 forms become

    N=y^T P_N z,  M=y^T P_M z,

    P_N=[[6,4,10,9],[4,3,7,6],[10,7,15,14],[9,6,14,13]],
    P_M=[[11,8,20,18],[8,6,14,13],[20,14,34,31],[18,13,31,29]]. (6)

Both are symmetric and every entry is positive. These are exact substitutions
of (3) into the inherited bilinear identities.

**Theorem 2.** For F=N or M and X in {O,Z,J},

    F(Xy,z) >= (4/3) F(y,z),
    F(y,Xz) >= (4/3) F(y,z).                        (7)

Proof: 3 X^T P_F-4 P_F is entrywise nonnegative for each of the six cases.
This is an exact 4-by-4 integer calculation; the verifier prints all six
matrices. Multiply by nonnegative y,z. Symmetry gives the right-source case.

The constant 4/3 is sharp even on physical sources. Let u,v be powers of two
at least two. The transitive sources t_u=(1,0,u,1), t_v=(1,0,v,1) are actual.
Their exact values are

    N(t_u,t_v)=3uv+1,
    N(Jt_u,t_v)=4uv+u+3v-1,
    M(t_u,t_v)=6uv+u+v+1,
    M(Jt_u,t_v)=8uv+2u+6v+1.

Let u=v tend through powers of two. Both ratios tend to 4/3. No relaxed cone
example is being substituted for a physical sharpness witness.

There is a stronger two-step conclusion: for every two-letter word W,
W^T P_N >= (15/7)P_N and W^T P_M >= (7/3)P_M. Their minimum entry ratios,
over the nine words, are respectively 15/7 and 7/3, both at W=JJ. In
particular each two-letter extension more than doubles either C2 value.

## 5. The same coordinates give positive C3 forms

For two nontrivial sources the eight complete C3 forms have the following
matrices, with the family names of the inherited C3 theorem:

    A=[[14,10,21,20],[10,7,16,15],[21,16,33,31],[20,15,31,29]],
    B=[[14,10,26,23],[11,8,20,18],[20,15,35,32],[19,14,34,31]],
    C=[[23,18,44,40],[18,14,34,31],[44,34,80,73],[40,31,73,67]],
    T=B^T,
    J0=[[12,12,24,24],[12,12,21,21],[15,12,30,27],[15,12,30,27]],
    J1=[[12,12,21,21],[12,12,18,18],[18,15,32,30],[18,15,31,29]],
    J2=J1^T, J3=J0^T.                              (8)

Again every entry is positive. Thus the earlier zero entries in the Type-II
positive cone were caused by a weaker choice of coordinates, not by a zero
physical leading direction in this cone.

For every three-letter word W and every matrix P in (8),

    W^T P >= 2P,   P W >= 2P                         (9)

entrywise. This is a fixed finite integer certificate: form the 27 products
of the three displayed 4-by-4 generator matrices and check both inequalities.
The minimum entry ratios over all words and both sides are

    A:27/8, B:61/20, C:113/31, T:61/20,
    J0:2, J1:9/4, J2:9/4, J3:2.

The code independently reconstructs (8) from the original response formulas
before checking all 432 matrix inequalities. No sampled source bank is a
premise for (9).

The C3 theorem permits singleton source boundaries. They do not belong to
the cone (2). Direct evaluation gives the two positive linear readouts

    C(o,w)=21a+16f+36e+33p,
    T(o,w)=12a+9f+15e+15p.

For Type II the two allowed singleton boundaries are scalar converses of
J0(v,o)=3(D-E)+9U-6d, which is also the second displayed readout.
The same 27 products each multiply both rows entrywise by at least two.
Their exact minimum ratios are recorded by the verifier. Thus they have the
same three-letter halving descent, without using negative singleton
coordinates in (9).

## 6. What the descent proves, and the remaining coverage obligation

Given an actual C2 representation, remove a final source letter. By (7) its
value drops by a factor at least 4/3 while remaining a physical C2 value in
the same N or M family. Continue to N(t,t)=13 or M(t,t)=29. The regular-five
value 15 is a separate seed, exactly as in the complete C2 theorem.

Removing two final letters from a source, whenever available, drops a C2
value by more than half. For C3, remove three final letters from either source
and use (9); singleton boundaries use the explicit linear readout. If no
source has three letters, each is one of the finitely many source words of
length at most two. These are the finite terminal cases. The representation
descent therefore has logarithmic depth and never leaves actual source
membership.

The direction matters. This says

    existing representation -> smaller existing representation.

Universal coverage still requires the reverse existential statement:
for every remaining odd target n, at least one of the physical inverse
branches of the positive bilinear system exists. A finite bank covering a
large interval, a nonnegative integer solution of (6) or (8), and the
existence of the terminating source inverse do not establish this statement.

For example, the rejected formal tuple (1,0,3,1) paired with t gives the
formal N value 19. The exact physical C2 inverse rejects 19. Thus relaxing
physical membership can create false scalar coverage even at a small target.

## 7. A short-certificate reformulation of the unresolved range problem

Let b=(0,0,0,1)^T and write G_w for the product of source matrices associated
with a word w. The exact nonexceptional C2 scalar language is

    { b^T G_u^T P G_v b : u,v in {O,Z,J}*, P=P_N or P_M }. (10)

For C3 use the eight matrices (8), together with the two positive singleton
readouts. Thus both complete spectra are ranges of explicit nonnegative
four-state matrix products, with one designated middle matrix. The two words
are freely chosen, but every word is physical. No membership test of a
relaxed numerical cone occurs in (10).

As an elementary complexity corollary, binary membership in either complete
class has certificates of length O(log n) bits, verifiable in polynomial
time in log n. For C2, iterating (7) in the two source words gives

    n >= 13 (4/3)^(|u|+|v|)

in N, and the analogous bound with 29 in M. For C3, remove triples from the
two source words using (9). Positivity of the terminal integer yields

    floor(|u|/3)+floor(|v|/3) <= floor(log_2 n),
    |u|+|v| <= 3 floor(log_2 n)+4.                  (11)

The singleton boundaries have the same argument with just one word. A
certificate specifies the family and its words, or the C2 exceptional seed
15. Replay (5), then evaluate the appropriate fixed matrix or row. Positive
readout coefficients bound every source coordinate by n at the end, and
U strictly increases along each word, so all intermediate integers have
O(log n) bits. A verifier may reject once U>n. A fixed number of additions
per letter and the final fixed-size bilinear evaluation therefore suffice.
This establishes NP membership for these two restricted scalar classes;
it says nothing about NP membership of the unrestricted tournament spectrum.

The remaining uniform-coverage question is exactly a range-universality
problem for (10) and its C3 analogues. A source-word shortening certificate
answers membership for a supplied positive instance, not whether every
sufficiently large odd integer is a positive instance.

## 8. Reproduction and scope

Run `python3 verify.py` and `python3 -O verify.py` in this packet. The verifier
uses explicit exception gates, not removable assertions. It compares the
new three-letter grammar with the independent indexed grammar, checks the
arithmetic inverse and hostile cases, tests literal one-vertex graph moves,
reconstructs every bilinear matrix from the original formulas, and verifies
the exact matrix growth certificates and physical sharpness identities.
The frozen dependencies are named and hashed in the receipt.

This packet proves a complete positive source language and physical descent
of represented values. C1/C2/C3 universal odd coverage remains open here.
