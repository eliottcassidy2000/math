# Tiling wedges, modular centers, and paired arithmetic families

Date: 2026-10-03 (America/Denver). Owner proposal, with two owner-confirmed
corrections: the first triple is (1,1,3), and q returns (O,E).

**Status:** PROVED for the elementary identities proved below; FINITE-EXACT
for the enumerations in the companion output. The interpretation of p's
two numbers as family indices is an explicit assumption, not a confirmed
definition. No new theorem ID or prize-problem consequence is claimed.

## Inheritance and concept board

- Closest proved mechanism: the fixed-path tiling convention in
  `01-canon/definitions.md`, the presentation/class distinction in
  [THM-1430](../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md),
  and the loop/path discrete derivative in
  [THM-4524](../../01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md).
- Hostile: several fixed-path tilings can represent one isomorphism class;
  in this note the complete order-five census gives 64 presentations and 12
  ordinary classes. Also, a boundary cell has fewer than eight symmetry copies.
- Corrected near miss: the tiling/self-converse count correction in
  `01-canon/MISTAKES.md` (search "SC(n)=A000568(n-1)") warns against reading a
  coordinate count as an isomorphism-class count. We do not use that identity.
- Relevant sidecar: retain coordinate addresses, reflection signs, chosen
  path, and loop bits. Bare modular residues do not determine tournament bits.
- Cards used: "Search the statement before the method" and "Correct the
  object before sharpening the technique" in `00-navigation/META-PATTERNS.md`.
  No new meta-pattern promotion is warranted by this single thread.

| Lane/object | Operation and invariant | Loss/boundary and decisive test |
|---|---|---|
| Anchor: six free tiles inside fifteen positions | Add path and loop positions | Presentation versus class; enumerate all 64 tilings |
| Niche: modular multiplication wedge | Signed reflections and transpose | Keep sign and zero axes; reconstruct all cells |
| Odd centers / paired atoms | Invert 4 modulo an odd modulus | Specify the nonzero-index table; test both mod-4 families |
| Wildcard: p | Addition and odd multiplication | Output typing; discriminate and invert the quadratic |
| q | Affine lattice embedding | Multiple-of-four and strict inequality constraints; invert |

The new common coordinates connect the first two lanes by an explicit
position bijection, the third by parity, and p/q by an exact mixed-family
coordinate change. None supplies a binary comparison of modular residues
that is invariant under tournament isomorphism.

## 1. The fifteen-cell wedge and its nine-cell L

For modulus m >= 2 set h=floor(m/2). Define

    D_h = {(a,b): 1 <= a <= b <= h},  M_m(a,b)=ab mod m.

Every nonzero residue r can be written r=epsilon*a (mod m), where
a=min(r,m-r) and epsilon in {+1,-1}; choose +1 when r=m/2.
For two residues, fold both coordinates and sort them. Then

    M_m(r,s) = epsilon_r epsilon_s M_m(min(a,b),max(a,b)) mod m.

The zero row and column are known separately. This proves that the wedge
reconstructs the entire multiplication table. The actions are the two
coordinate sign changes and transposition; a one-coordinate reflection
negates the entry, while a double reflection and transposition preserve it.

For m=10 the wedge is:

    a\b  1  2  3  4  5
     1   1  2  3  4  5
     2      4  6  8  0
     3         9  2  5
     4            6  0
     5               5

The arm b=5, traversed from a=5 to a=1, followed by a=1 from b=4 to b=1,
is exactly 5,0,5,0,5,4,3,2,1. Removing these positions leaves
2<=a<=b<=4, six entries arranged in lengths 3,2,1.

For m=2n, the same boundary is always elementary: M(1,b)=b, while
M(a,n)=n for odd a and 0 for even a. It contains 2n-1 positions;
the interior contains (n-1)(n-2)/2 positions. The boundary's starting value
depends on the parity of n.

"One eighth" means a symmetry fundamental region, with stabilizers on its
boundary. For m=10, the 15 nonzero-coordinate orbits have sizes
6*8+8*4+1=81; the 19 zero-axis cells complete the 100-cell table.
For odd m=2h+1 the corresponding count is C(h,2)*8+h*4=(m-1)^2.
For even m=2h it is C(h-1,2)*8+2(h-1)*4+1=(m-1)^2.

## 2. Exact coordinate match to loops, path arcs, and free arcs

Fix vertices 1,...,n and path n -> n-1 -> ... -> 1. For 1<=v<=u<=n,
the following is a bijection from unordered pairs with loops to D_n:

    (u,v) -> (v,n)       if u=v           [loop]
             (1,v)       if u=v+1         [path arc]
             (v+1,u-1)   if u>=v+2        [free arc].

These images are respectively b=n, a=1 with b<n, and the interior
2<=a<=b<=n-1. The three regions are disjoint and exhaust D_n.
For n=5 this realizes 15 = 5 loops + 4 path arcs + 6 free arcs.

The map preserves position and the three roles; it has not been shown to
transport graph incidence or an arithmetic invariant of isomorphism classes.
Modular entries are labels, while tournament choices are binary orientations.
The two free arcs (4,1) and (5,3) both get residue 6 modulo 10, so residue
alone loses even the pair address.

THM-4524 supplies a separate genuine graph relation: loop bits l_v switch
arc uv by l_u XOR l_v, and the path orientations are the adjacent discrete
derivative of l. Global complement of all loop bits gives the same result.
At n=5, six free bits and five loop bits give 2048 selfies mapping exactly
two-to-one onto 1024 labeled tournaments. The four path entries are fixed in
the tiling gauge, not four further independent bits. If one chooses to read
the alternating 5/0 arm as alternating 1/0 loop bits, that switching reverses
every path arc; this is an additional binary interpretation, not automatic.

## 3. The two odd families are the two explicit inverses of 4

Let F_n=(n,ceil(n/2),2n+1), n>=1. Write

    B_N = F_(2N-1) = (2N-1,N,4N-1),
    A_N = F_(2N)   = (2N,N,4N+1).

Thus the list alternates B,A,B,A,..., and the owner's atom is precisely

    N <-> ((2N-1,2N),(4N-1,4N+1)),

with the first pair giving x-coordinates and the second pair z-coordinates
of the two triples sharing y=N. Ordered pairs retain the intended matching.

Use the nonzero-index multiplication table, with rows/columns 1,...,m-1.
For odd m=2n+1 its central rows and columns are n,n+1. Put
d=n^2 mod m and c=n(n+1) mod m, using least nonnegative residues. Then

    center = [[d,c],[c,d]],   4d=1 (mod m),   c=-d (mod m).

Proof: n+1=-n modulo m and 2n=-1 modulo m. Since 4 is invertible and d is
nonzero, d+c=m. Resolving the inverse of 4 gives, for every N>=1:

    m=4N-1 (B):  [[N,3N-1],[3N-1,N]],
    m=4N+1 (A):  [[3N+1,N],[N,3N+1]].

These are all odd moduli >=3, without a primality assumption. In particular
mod 11 gives [[3,8],[8,3]], and mod 9 gives [[7,2],[2,7]]. The smallest
cases m=3 and m=5 satisfy the same formulas.

For even m=2k, the nonzero table has one central entry k^2 mod 2k:
it is 0 when k is even and k when k is odd. For k>=2, the like-side diagonal neighbors
(k-1)^2 and (k+1)^2 are center+1 modulo 2k; the crossed neighbors
(k-1)(k+1) are center-1. At m=10 these are center 5 and neighbors 6 and 4.

## 4. What p preserves, under the family-index interpretation

As written, p returns two numbers, not two triples. Assume its intended
lift is p(F_a,F_b)=(F_s,F_t), with

    P(a,b)=(s,t)=(a+b,(2a+1)b),   a,b>=1.

If instead s,t are meant as coordinates of one triple, a coordinate choice
and a different closure test are needed. The following scalar identities
hold regardless of that semantic choice.

P is injective on positive ordered integer pairs. Indeed b solves

    2b^2-(2s+1)b+t=0,
    D=(2s+1)^2-8t=(2a-2b+1)^2.

The two potential roots sum to (2s+1)/2, a half-integer, so at most one is
an integer. Exactly one is b for data in the image; then a=s-b. Explicitly,
choose the sign for which b=(2s+1 +/- sqrt(D))/4 is integral. The complete
image test is: s>=2, D an odd nonnegative square, and that integer root
satisfies 1<=b<=s-1.

An equivalent triangular-number identity, T_k=k(k+1)/2 for integer k, is

    t = T_s - T_(a-b).

At fixed s>=2 the maximum is T_s, achieved uniquely at
a=floor(s/2), b=ceil(s/2). The inverse needs both outputs, not family letters
alone. Under the index lift, parity yields:

| Input families | Output families |
|---|---|
| A,A | A,A |
| A,B | B,B |
| B,A | B,A |
| B,B | A,B |

There is also an actual multiplication law on family indices:

    a star b = a+b+2ab = a+t,
    z(F_(a star b)) = z(F_a) z(F_b).

It is commutative and associative, because z identifies it with ordinary
multiplication of odd integers >=3. An identity would require adjoining
index 0, corresponding to odd integer 1; that is outside the original list.

## 5. q is an injective affine chart with a restricted image

For positive integers 1<=I<J the confirmed definition is

    q(I,J)=(O,E)=(2(I+J)-1,4J).

Its inverse is J=E/4 and I=(2O+2-E)/4. Consequently its image is exactly

    E divisible by 4, E>=8, O odd, E/2+1 <= O <= E-3.

The inequalities are precisely I>=1 and I<J. For each fixed J there are
J-1 outputs. Merely requiring E even does not suffice: (5,6) is impossible.
The excluded endpoint (7,8) requires I=J=2. Allowing I=0 or I=J changes the
image and must be stated separately.

There is a precise link to the atoms: for B_I=F_(2I-1) and A_J=F_(2J),

    O=x(B_I)+x(A_J),  E=2x(A_J).

Thus q records the sum and the doubled second coordinate of this mixed pair.
One recovers x(A_J)=E/2 and x(B_I)=O-E/2, and hence

    P(2I-1,2J) = (O, E(2O-E+1)/2).

This is an invertible coordinate change between q and P on the stated
mixed-family domain I<J. For example q(1,2)=(5,8) corresponds to
P(1,4)=(5,12). q is therefore related to p by an explicit map, rather than
only sharing the alternating family pattern. Under the index lift the
outputs of q themselves select a B member and an A member.

## Reproduction and stopping boundary

Run `python3 04-computation/experiments/tiling_modular_atoms_20261003.py`.
Output: [tiling_modular_atoms_20261003.out](tiling_modular_atoms_20261003.out).
All arithmetic is exact; no external tables, random sampling, or inherited
filters are used. The script reconstructs 348550 cells for moduli 2..101,
checks both central families for N=1..500, 10000 ordered p inputs and 4950
q inputs, the position bijection through n=30, and independently canonicalizes
all 1024 order-five tournaments against all 64 fixed-path presentations.
It separately checks the 2048-to-1024 switching map. Positive controls are
the owner's mod-9/10/11 formulas; hostile controls retain boundary orbit
sizes, q's excluded endpoints, and the presentation/class distinction.

The all-parameter claims rely on the proofs above, not extrapolation from
the computation. The independent algebra and wedge reviews agreed with the
proofs. The open semantic point is the intended output typing of p. A further
tournament/arithmetic theorem would need a specified binary observable and
a proof of what survives a change of Hamiltonian path or vertex labeling.
