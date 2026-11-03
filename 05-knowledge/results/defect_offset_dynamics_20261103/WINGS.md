# Two integer wings, a four-parameter class family, and a ramified clock

**PROVED**, with exact symbolic run coefficients and independent literal counts.

## W1. Actual source and composition law

The six-vertex tournament Q has bit rows

    (2,4,9,19,39,15),

where bit j in row i means i->j. It is strong, prime under homogeneous
substitution, rigid, and is taken to its converse by i->5-i. These finite
properties are checked exhaustively in `wing.py`.

Replace vertices (0,1,2,3,4,5) by transitive blocks of sizes

    (m,l,1,1,k,n),   m,l,k,n>=1.                           (W1)

The four variable vertices induce the transitive order 4->5->0->1, so R1
applies. Define the actual left-wing response

    x=2^m(2^l-1)-2^(l-1)+2,  b=2^(l-1)-1,

and the right response (y,c) by replacing (m,l) with (n,k). Then

    H=2(xy+bc)-1.                                         (W2)

Proof certificate: the complete run polynomial in `wing.json` has sixteen
positive-index coefficients. Substituting r!S(t,r), cancelling, and grouping
the four dyadic variables gives exactly W2. `wing.py` performs the symbolic
identity; sixteen independent subset-DP checks cover all block sizes 1 or2.
The general proof is the compression/lift bijection in RUNS.md, not those
sixteen samples.

This defines two attachments by their joint response to the fixed two-vertex
middle, and keeps the dot-product coupling. It does not assert that every
positive integer pair is physically realizable. The precise guards are

    b+1 is a power of two;
    (x+b-1)/(2b+1) is an integral power of two at least2.   (W3)

These conditions are necessary and sufficient for a wing. The recovered
exponents recover the actual transitive blocks and their incidences.

The growth laws on a wing are

    m+=1: (x,b) -> (2x+b-1,b),
    l+=1: (x,b) -> (2x+(x+b-1)/(2b+1)-2,2b+1).            (W4)

The integer division in W4 is licensed by W3. It cannot be performed on an
arbitrary residue pair after forgetting its realization witness.

## W2. The two- and three-parameter source formulas

Setting l=k=1 gives

    H(T_mn)=2(2^m+1)(2^n+1)-1.                            (W5)

This is the incoming Session164 complementary-offset family in canonical
block coordinates; its formula and protected-port proof are credited in
[INCOMING.md](INCOMING.md). The present enlargement is W1-W4 and W6.

Thus (m,n)=(1,1),(3,1),(3,3) give 17,53,161. Increasing a single block
implements 3H+2 exactly when its size changes 1->3: solving
2^r+1=3(2^s+1), r>s>=1, gives s=1,r=3 by the 2-adic valuation.
No +4 move can stay anywhere in W5, since all its values are1 modulo4.

Setting l=1 instead gives, with A=2^m, B=2^n, C=2^k,

    H+1=(A+1)[(2B-1)C-2B+4].                             (W6)

The operation m:1->3 therefore supplies 3H+2 on an entire two-parameter
bank, for every n,k>=1. In particular n=1 gives the uniform defect-four
return 9*2^k-1 -> 27*2^k-1 found independently in RUNS.md. The integer m=3
output is retained as a consumed response coordinate; it is not reset to1.

## W3. Why independent affine wing updates cannot supply the other offset

The actual wing pairs affinely span R^2; for example (3,0),(5,0),(6,1) do.
Suppose fixed affine maps v->Av+a and w->Bw+b were to implement 3H+4 for
every pair of wings, at the response level. W2 would require

    (Av+a) dot (Bw+b) = 3 v dot w +1.

Affine spanning, independently in each variable, forces A^T B=3I,
A^T b=0, B^T a=0, and a dot b=1. The first makes A,B invertible; the next
two give a=b=0, contradicting the last. Hence no such independent affine
updates exist, even before imposing W3. The +2 lane instead requires just
the factor3 in the bilinear form. Actual realizability still has to be checked.

This locates the missing resource: a +4 construction must couple the two
wings, change the middle relation, add an observed constant coordinate, or
use a more general nonlinear operation. A formal extra coordinate is not a
free physical path. The literal complement escapes in ESCAPES.md change the
old relation and therefore lie outside the obstruction's hypotheses.

## W4. A closed modulo-seven clock has nine states per wing, not seven

For m,l modulo3, put M=2^m, L=2^l in F7. There are nine possible pairs.
The current residue response is

    b=L/2-1,  x=M(L-1)-L/2+2.                            (W7)

If L!=1, x,b recover M,L. If L=1, all three M values collapse to the same
point (x,b)=(5,3). One legal l-growth sends x to M+1 and separates all three.
Thus the seven current wing observations require exactly nine states under
the two actual independent growth operations. Nine suffice by retaining
(m mod3,l mod3); fewer cannot pass the current and next-l observations.

For two wings, there are49 current joint observations and81 phase states.
Even observing only H modulo7 after arbitrary independent block growths
separates all81: `wing.py` lists all81 future-observation signatures over the
81 growth residues and finds them distinct. Conversely those81 phases compute
every such response, so this is minimal for this exact operation language.
The equality to81 in another tournament automorphism problem is not a map.

The exact mechanism is a disappearing division denominator, not random
information loss. At L=1 the denominator 2b+1 vanishes modulo7; the integer
source still contains M, and the next growth makes it observable again.
This is the same type of sidecar requirement as the earlier transport work,
now in a fully realized attachment family.

## W5. Golden/F49 transport stays separate from the wing quotient

In W6 the fixed vertex1 remains a named processing root. Let
D=[(2B-1)C-2B+4]/2. The guarded polynomials give

    E_1(1)=2D-1,  S_1(2)=H,
    v=(H,X,Y)=(2(A+1)D-1, AD, (A+2)D-1).

In the inherited chart F49=F7[theta]/(theta^2-4theta-2),
alpha=a+b_theta*theta and g=H-Y, this is

    (a,b_theta,g)=(5AD,(A+2)D-1,AD) modulo7.              (W8)

The symbol b_theta is not the wing coordinate b. All five literal Q4
presentations are retained, using their inherited distinct integer matrices.
The normal-basis Frobenius is still
(a,b_theta,g)->(a+4b_theta,-b_theta,g); it is not automatically a realizable
block growth. The full relation, ordered ports and exponents are retained
when a phase leaves this source chart. The 81-state claim in W4 is for the
specified Hamilton-response clock, not for every joint port tensor or Q4
observer. F49 is not arithmetic modulo49.

## W6. Exact class and path-marked enumeration inside this family

The prime quotient Q recovers its six substitution blocks. Its automorphism
group is trivial, and each transitive block is rigid. Consequently the four
positive sizes give distinct ordinary classes, all rigid. Their series is

    A(z)=z^6/(1-z)^4.

For clarity, any proper homogeneous set meeting two different blocks would
project to a nonsingleton homogeneous set of Q. Primality forces that
projection to be all Q; but an omitted vertex in any block then sees both
an in-neighbour block and an out-neighbour block, contradicting homogeneity.
Thus the full blocks are exactly the maximal proper homogeneous sets and
are intrinsically recovered, including singleton blocks.

Converse sends (m,l,k,n) to (n,k,l,m), so self-converse classes have series
S(z)=z^6/(1-z^2)^2 and converse-merged classes have series (A+S)/2.
Root-marked classes have series z A'(z).

Because automorphisms are trivial, each actual Hamilton path gives a different
path-marked class. Summing W2 gives

    P(z)=z^6(17-76z+148z^2-144z^3+64z^4)
         /((1-z)^4(1-2z)^4).                             (W9)

Its first coefficients, at orders6 onward, are 17,128,630,2540,9067,29756,91752.
For example set B=z/(1-z), A=2z/(1-2z). The single-wing sums are
X=A^2-3AB/2+2B^2 and Y=AB/2-B^2; then P=z^2[2(X^2+Y^2)-B^4].
Both poles in W9 have order4, giving minimal eventual scalar recurrence
characteristic (t-1)^4(t-2)^4, of order8 over Q.

These are exact subclass contributions and marked-presentation recurrences.
They do not extend the known values of A000568 or classify the full metagraph.
Actual complements can leave this family, as the literal certificates do.
