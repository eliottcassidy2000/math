# A bounded central word and an explicit arithmetic coverage lemma

Status: **PROVED, independently audited in the stated scope** for the physical compiler,
finite inverse, normalization, and conditional implication below. The uniform
arithmetic occupancy assertion in Section 7 is **OPEN**, not a proved theorem.
Session 201. No source-pair existence is inferred from a height bound or from
surjectivity of a signed scalar.

## 1. Inheritance and the actual remaining question

The reviewed physical adjoint gives a complete inverse for a target already
known to be a C1 path count. The sharper adjacent-cut theorem in
`dependencies/PHYSICAL_EXISTENCE.md` improves that inverse to source height
floor(sqrt(3(n-1))). It does not show that an inverse pair exists. Likewise,
the signed-positive-section result identifies an actual physical slice on
which a signed value becomes a positive count; it does not show that every
signed fiber meets that slice. Append-only repair of an arbitrary signed
witness is now known to fail, even at represented targets 3 and 33.

The older `dependencies/TWO_INDEX_40.md` already isolates an OPEN two-index
arithmetic lemma using forty C2/C3 graft templates. The present route is
different: both free parameters are explicit unit-source indices, there is
only one Hamilton cycle, and a bounded finite menu supplies the entire
middle of the source word. Its arithmetic lemma is a stronger sufficient
condition than unrestricted C1 existence, not a reformulation proved
equivalent to universal tournament coverage.

Live concepts: physical adjunction; explicit binary unit sources; bounded
central words; arithmetic progression inversion; signed-fiber occupancy.
The hostile boundary is target 17383, which cannot use a central word of
length at most seven even though it has a length-eight construction.

## 2. The physical unit source and chronological conventions

Use column profiles s=(D,E,U,d)^t and chronological words: the first letter
acts first. The seed is t=(1,0,2,1). The reviewed operations are

    O(s)=(D,0,2U,d),
    Z(s)=(D-E+U,d,2U-d,d),
    J(s)=(D-E+U,D-E,D+U,D).

Each adds one vertex to an actual source. O and Z retain the old root; J
retains the new root. To construct its adjacency matrix, begin with the
two-vertex source 0->1, root 0. A new O vertex loses to every old vertex.
A new Z or J vertex beats the old root and loses to every other old vertex.
For a chronological word w=c1...ck, put M_w=M_ck...M_c1.

Define g(1)=1 and

    g(2h)=g(h)                 (h>=1),
    g(2h-1)=g(h)+h-(h mod 2)  (h>=2).

For u>=2 set tau(u)=(g(u),u mod 2,u,1)^t. Write u-1 in binary, remove its
initial 1, and replace each remaining 0 by Z and each 1 by O. Call the
result B(u). Induction using the displayed operations gives

    M_B(u)t=tau(u),
    O tau(u)=tau(2u),       Z tau(u)=tau(2u-1).

Thus tau(u) is an actual nontrivial source, not a formally chosen profile.
Its source order is 1+ceil(log2 u), and |B(u)|=ceil(log2 u)-1. No singleton
parameter u=1 is admitted in this family.

Let

    Q = [0  0  0  1]
        [0  1  0 -1]
        [0  0  1 -1]
        [1 -1 -1  1],          K(v,w)=v^t Q w.

The dagger reverses a word and exchanges O and Z, fixing J. The exact
reviewed identity is K(M_w v,x)=K(v,M_wdagger x). Also
K(t,s)=F(s)=D-E+U, and appending J to any nontrivial source word constructs
a strong tournament with exactly F(s) Hamilton paths and one Hamilton cycle.

## 3. A finite arithmetic family of actual C1 graphs

For m>=0 let W_m consist of the empty word; J if m>=1; and all words of
length 2 through m that begin and end in J. Define

    Psi_W(u,v)=tau(u)^t Q M_W tau(v),  u,v>=2, W in W_m.      (1)

There are 1095 words in W_8. Expansion of their integer matrices Q M_W
leaves 592 distinct matrices. Keeping the first occurrence by increasing
length preserves the least central length for each such matrix.

Every value in (1) is the path count of the actual C1 graph whose complete
chronological construction word is

    B(v) W dagger(B(u)) J.                                  (2)

Indeed adjunction turns (1) into F(M_dagger(B(u)) M_W M_B(v)t).
Its exact order is

    1+ceil(log2 u)+ceil(log2 v)+|W|.                         (3)

There are no independent parity or g coordinates: both are determined by
their indices. The matrix list is generated from three displayed linear
operations, so (1) is a finite explicit list of bivariate integer formulas,
not an existential quantifier over arbitrary physical source profiles.

Conversely, any C1 source word has an initial O/Z block, a final O/Z block,
and either an empty middle or a middle that begins and ends in J. Absorbing
the outer blocks into the two unit indices gives (1), but with no a priori
bound on |W|. This converse does NOT prove that m=8 suffices.

## 4. Complete finite inverse, with all rounding explicit

Use the physical nonnegative coordinates

    y=(a,f,e,p)=(D-E-d,U-D+E-d,E,d-E),
    s=P y,     U=u0^t y,     d=e+p>=1,

where

    P=[1 0 2 1]            u0=(1,1,2,2)^t.
      [0 0 1 0]
      [1 1 2 2]
      [0 0 1 1]

For each word W put A_W=Q M_W, C_W=P^t A_W P, and

    b_W=min_(i,j) 4 C_W[i,j]/(u0[i] u0[j]).                 (4)

This is an integer because each denominator divides four. The coefficient
comparison gives the exact rigorous bound

    Psi_W(u,v) >= (b_W/4) u v.                              (5)

In particular b_W>=2. One all-word justification is that
P^t Q P >= (1/2)u0 u0^t entrywise, every operation preserves the positive
cone, and every operation increases or preserves the source cover. The
finite checker also directly verifies (4) for each matrix in W_8.

The family W_m is closed under dagger, and A_W^t=A_Wdagger. Therefore the
whole family is symmetric in its two indices; restricting to u<=v loses
no value. For a target n, a complete inverse domain is

    2<=u<=isqrt(floor(4n/b_W)),
    u<=v<=floor(4n/(b_W u)).                                (6)

These bounds are necessary for a representation of n and do not promise
one. They prove that a failed exhaustive inverse is an all-height
nonmembership certificate for this finite-word family.

There is a second exact inverse, independent of enumeration of every v.
For fixed W,u write

    tau(u)^t A_W=(alpha,gamma,beta,delta),
    n=alpha g(v)+beta v+gamma(v mod 2)+delta.                (7)

For all W in W_8 and every u>=2, alpha,beta>0. The finite matrix certificate
is that columns 0 and 2 of P^t A_W are entrywise nonnegative and strictly
positive in rows e,p. Since e+p=1 for a unit source, this proves the claim
for all indices, not just the tested indices.

Writing v=2t+q for q=0,1, any solution of (7) satisfies

    2 beta t = n-delta-(beta+gamma)q  (mod alpha).            (8)

Divide by gcd(2 beta,alpha), provided it divides the right side, invert the
remaining coefficient modulo alpha/gcd, and test the resulting progression
only inside (6). The modulus-one case is all integers t. Testing (7) after
the congruence filter is essential: the congruence alone is not sufficient.

## 5. Boundary absorption and the 235-word normal menu

The global physical word equality JZ=OJ is exact. For each of the finitely
many words in W_8, form its finite equivalence component under this local
rewrite in either direction; length is preserved, so this process terminates.
If the component contains a word with an O/Z boundary letter, remove all
initial and final O/Z blocks. The resulting middle has strictly smaller
length. Both removed blocks are absorbed into actual unit indices using
O tau(u)=tau(2u), Z tau(u)=tau(2u-1), and adjunction on the left.

This proves containment of that component's value set in the value set of
a shorter central word. Since every shorter central word was already in
W_8, discarding such a component does not alter the union of values. If no
boundary can be exposed, keep one representative. The explicit finite
rewrite certificates leave 235 retained words and 376 discarded components.
Their retained counts by lengths 0,...,8 are

    1,1,1,1,3,8,21,55,144.

This is a union-of-value-sets equality, not a graph-isomorphism claim and
not an assertion that every local rewrite equivalence is known. Nor is it
an all-length confluence theorem. In particular one must not impose u<=v
on the reduced list merely because the original full list was transpose
closed. The author inventory uses the full 592-matrix transpose-closed
list, so it has no such representative-selection defect.

## 6. Graph size and the necessary length-eight hostile control

For lengths 0,...,8 the exact minimum b_W values in (4) are

    2,4,4,8,12,20,36,60,100.

Combining (3), (5), and
ceil(log2 u)+ceil(log2 v)<=floor(log2(uv))+2 shows that every represented
target in W_8 has an actual C1 graph of order at most floor(log2 n)+7.
This is an upper bound for this compiler, not a lower bound on optimal
tournament order. For example, at length eight, uv<=n/25 gives
floor(log2(uv))<=floor(log2 n)-4 and then (3) gives the stated bound.
The other eight lengths follow from the displayed exact minima.

The complete arithmetic inverses establish:

    1877=Psi_JJJJJOJJ(3,5),
    17383=Psi_JJJJJOJJ(3,55).

Neither target has any representation with |W|<=7. At 17383 a second
length-eight matrix representative JJOJJJJJ also works with indices 3,55.
Literal graph counting verifies respectively 14 and 17 vertices and one
Hamilton cycle. Thus the assertion "every odd n>2647 uses W_7" is false.
This is a restriction of this compiler, not an obstruction to other
positive source constructions.

## 7. The precise missing arithmetic lemma and its implication

**OPEN bounded-central-word occupancy lemma.** For every odd integer n>2647
there exist integers u,v>=2 and a word W in W_8 such that

    (g(u),u mod 2,u,1) Q M_W (g(v),v mod 2,v,1)^t=n.         (9)

One may equivalently use the 235-word normal menu, allowing both index
orders. No source orbit, signed fiber, or arbitrary admissible profile is
hidden in (9). Its complete target-dependent search domain is (6).

If (9) holds, every odd n>2647 is an actual C1 path count, with the explicit
graph-size bound in Section 6. Exhaustive arithmetic enumeration below
that threshold, together with the 27 literal repair graphs in
`dependencies/REPAIR_GRAPHS.json`, supplies every remaining positive odd
target except 7 and 21. The target 1 is supplied by a transitive tournament.
This conditional conclusion is positive existence. The general exclusion
of 7 and 21 is a separate inherited theorem, not proved by this packet.

The 29 targets absent from the arithmetic family through 5001 are

    7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,
    299,421,449,503,533,569,583,637,821,1067,1339,1861,2647.

The author checker also compares against the complete physical source
inventory U<=5000; F>=U+1 makes this complete at these target values.
Thus these fixed small gaps are not artifacts of the central-length cap.
The 27 repair graphs for the other gaps have two or three Hamilton cycles.

The larger C++ inventory, the Python progression inverse, and literal
Hamilton path/cycle counts are recorded separately. Finite survival of
(9), however far checked, is not a proof. In particular this packet does
not assert that the smaller index may be bounded by a fixed constant,
that coverage obeys a simple binary induction, or that a signed witness
can be appended into the family. The all-target occupancy step remains
the mathematical problem to solve.
