# Exact source-cover minima and the capacity of a quadratic source norm

Session 199. **PROVED:** a sharp all-order minimum for the reserved-cover
count of start-one sources; an exact quadratic-norm presentation of diagonal
C2 counts; and density-zero support for finitely many constructions using two
fixed extensions of one common source. **OPEN:** independent-source C2/C3
coverage and the full odd spectrum. This is not a new forbidden target claim.

## 1. Recovery, definitions, and the proposed bridge

The unrelated source was the represented-norm theme of THM-3333,
`01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md`,
and THM-3334,
`01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md`.
Their surrounding derivations distinguish an integer norm from its represented
Gaussian factors and distinguish a norm-preserving mutation from a physical
operation. The norm-65 collision is the canonical hostile. The idoneal discussion
in `05-knowledge/results/digit_idoneal_difference_20261018b.md`, Section 5,
further distinguishes a complete genus quotient from an unproved lift; no
idoneal completeness or class-field theorem is used here. The Stern--Brocot
source THM-4059 was read with its actual inverse/word-reversal sidecars rather
than treated as a source reachability theorem.

The target is the complete physical source grammar from Session 198
`bilinear-descent/PROOF.md`, frozen as `dependencies/POSITIVE_SOURCE_GRAMMAR.md`.
A source is a nonempty tournament with a distinguished vertex having exactly
one Hamilton path starting there. Write

    (D,E,U,d)=(H(V), e_c(V), R_c(V), H(V-c)),

where R counts ordered spanning two-path covers with the reserved nonempty
path starting at c and the other path allowed empty. H(empty)=1. For every
nontrivial source the positive coordinates and their inverse are

    v=(a,f,e,p)=(D-E-d, U-D+E-d, E, d-E),
    (D,E,U,d)=(a+2e+p,e,a+f+2e+2p,e+p).

All coordinates of v are nonnegative; e+p>=1. The seed for order two is
b=(0,0,0,1), and every actual source of order m>=2 has an O/Z/J word of
exactly m-2 letters. Every such word is physical. The matrices are

    O=[[1,0,1,0],[1,2,1,2],[0,0,0,0],[0,0,1,1]],
    Z=[[2,1,1,1],[0,1,0,0],[0,0,1,1],[0,0,0,0]],
    J=[[0,1,0,1],[0,0,0,0],[1,0,1,1],[0,0,1,0]].       (1)

O adds a sink; Z adds a vertex beating only the retained root; J does the
same and marks the newborn. Column products act rightmost first. Put
u=(1,1,2,2), so U=u v.

The concept board is **represented norm / physical chamber / Pareto rows /
source-class entropy / independent source choices**. The map to a quadratic
norm preserves the diagonal C2 scalar, but not arbitrary norm-preserving
moves. The complete physical source word remains the membership sidecar.
The cheapest decisive test is the order-two source under the explicit
norm involution in Section 4.

## 2. An exact finite recurrence proves all-order cover minima

For a finite set of nonnegative row vectors, let Min retain exactly its
componentwise-minimal distinct rows. Let R_k be the minimal rows among
u G_1 ... G_k for all choices G_i in {O,Z,J}. This is a row-product ordering;
its reversed word is the chronological physical source word. All length-k
words occur, so no source is lost by this convention.

Because every matrix in (1) is nonnegative, a dominated row remains dominated
after multiplication by any one of them. Consequently

    R_0={u},
    R_(k+1)=Min { r G : r in R_k, G in {O,Z,J} }.       (2)

This is an exact recurrence, not a relaxation. Every retained row is an
actual word product, and every discarded row dominates some retained row.
The operation on the right of (2) commutes with multiplication of all rows
by the same positive scalar.

Here is the complete initial calculation; entries within braces are rows.

    R_1={(2,2,3,3),(2,1,4,3)};
    R_2={(3,2,6,5),(4,4,5,5)};
    R_3={(5,4,10,9),(6,5,9,9),(6,3,11,9)};
    R_4={(9,6,18,15),(10,9,15,15),(10,5,19,15)};
    R_5={(15,10,30,25),(18,15,27,27),
         (18,9,33,27),(20,19,25,25)};
    R_6={(25,20,50,45),(27,18,54,45),
         (30,25,45,45),(30,15,55,45)};
    R_7={(45,30,90,75),(50,45,75,75),(50,25,95,75)}
       =5 R_4.                                          (3)

All rows in (3), including every removal by dominance, follow by multiplying
the at most four preceding rows by the three displayed matrices. Thus (2)
and the exact equality R_7=5R_4 imply by induction

    R_(k+3)=5 R_k for every k>=4.                         (4)

This is not extrapolation from a stabilized finite census: equality of the
full recurrence states, plus scale-equivariance of the transition, proves
all future steps.

**Theorem 1 (sharp cover minimum).** For every j>=0 the minimum of U over
all rooted start-one sources of the indicated order is

    order 3j+3: 3*5^j,
    order 3j+4: 5*5^j,
    order 3j+5: 9*5^j.                                  (5)

At order one the singleton has U=1; at order two the minimum is U=2.

Proof. At order m>=2, minimize u G_1...G_(m-2) b. This is exactly the
minimum fourth coordinate among R_(m-2). Equations (3)-(4), and R_0, give
(5) including its first two cycles of boundary values. Every retained row
is an actual word, so these are minima over physical sources, not merely
lower bounds on a numerical cone.

There are also simple uniform witnesses. Starting at the edge, the words
J, JZ, JZZ give balanced responses (h,e,h,e) with (h,e)=(3,1),(5,1),(9,1).
For any balanced source the chronological three-letter word JZZ has response

    (h,e,h,e) -> (5h,h,5h,h).                            (6)

This follows by three substitutions in the physical transfer laws. Repeating
(6) on the three seeds realizes every equality in (5), with exactly three
new vertices per repetition. QED.

In particular every source of order m>=2 satisfies

    U^3 >= 5^(m-1).                                    (7)

For the three residues, the necessary scalar comparisons are 27>=25,
125=125, and 729>=625. Order two gives 8>=5. This sharp source proof does
not cite or derive the unrestricted strong-tournament minimum theorem.
It independently produces the same growth scale within this source class.
It does not classify all equality sources.

## 3. Source-counting capacity, not a physical cone theorem

The independently proved rooted class count, frozen as
`dependencies/FIBONACCI_CLASS.md`, together with the unique weak singleton-tail
decomposition gives the generating function for all rooted sources

    W(z)=z(1-2z)/(1-3z+z^2).

Thus there is one class at order one, and at order m>=2 there are exactly
F_(2m-3) rooted source graph isomorphism classes. With

    lambda=(3+sqrt(5))/2,
    beta=3 log(lambda)/log(5),
    rho=beta/2=0.8969811534981604... ,                    (8)

the number of source graph classes having U<=B is O(B^beta).

Indeed (7) bounds their order by 1+3 log(B)/log(5), and summing the geometric
upper bound F_(2m-3)<=lambda^m over those orders proves the assertion.
Different graph classes may have identical response vectors; that only makes
the number of numerical profiles smaller. No source is obtained by assuming
a vector with this bound is reachable. All asymptotic constants are absolute.

The strict inequality rho<1 is exact: lambda^3=9+4sqrt(5)<25. This matters
in the application; a merely polynomial state bound would not suffice.

## 4. The C2 diagonal is a norm, but its arithmetic involution is not physical

For a nontrivial actual source define

    x=D-E+U,       y=(U+E)/2,       z=U-d.                (9)

The inherited parity U=E modulo 2 makes y integral. The inequalities
D-E>=d>=E and U>=2d show x odd, y>=1 and z>=1. Moreover

    x-2y=D-2E>=1,
    4y-x=U+3E-D>=0.

The latter is even minus odd and hence is at least one. Therefore every
physical point belongs to the strict chamber

    2y < x < 4y.                                        (10)

For two sources v,w, direct expansion of the complete C2 N form yields

    N(v,w)=x_v x_w+2(y_v z_w+z_v y_w).

In particular

    N(v,v)=x^2+4yz,                                    (11)

so every diagonal N value is 1 modulo 4. These are actual C2 tournaments
by the inherited graft theorem. The coordinates are not an arbitrary
parametrization of the integer triples in (10). In fact their relation to
the previously injective physical signature (A,B,C) is

    A=x+2y-2z, B=2y-z, C=z.

Thus (9) itself is injective on physical response profiles by
`dependencies/SOURCE_SIGNATURE.md`, but deciding which triples are physical
still requires that theorem's source-membership sidecar.

Consider the standard positive-triple norm move on the chamber x>2y:

    (x,y,z) -> (x-2y, x-y+z, y).                        (12)

Direct expansion preserves x^2+4yz. Its inverse on the image chamber
x'<y'-z' is (x',y',z')->(x'+2z', z', y'-x'-z'). These are the two outer
branches of the usual Zagier involution; no sum-of-two-squares theorem is
needed here. Every positive input in x>2y has a positive output, but

    x'-2y'=-x-2z<0.                                   (13)

Thus (12) sends **every** physical source triple outside the necessary
physical chamber. The smallest example is the transitive edge:

    (D,E,U,d)=(1,0,2,1), (x,y,z)=(3,1,1), N=13;
    (3,1,1) -> (1,3,1).

The norm 13 survives, actual source membership does not. The same loss
occurs before questions about genus labels, prime factors, or idoneal
completeness could help. Other target-dependent returns or other norm moves
are not ruled out by this one-step obstruction.

## 5. A density-zero obstruction for common-source constructions

Let F be any of the two complete C2 forms N,M or the eight complete C3
forms A,B,C,T,J0,J1,J2,J3. In this section **both arguments are nontrivial
physical sources**. All these forms have the physical product bound

    F(v,w)>=3 U_v U_w.                                (14)

For C2 this is weaker than the inherited sharp product bounds. For C3 it
follows from the positive matrices in the source-grammar dependency; the
verifier also checks P_F>=3 u^T u entrywise for all ten matrices.
Every O/Z/J step increases U, so U_(Xv)>=U_v for every word X.

**Theorem 2 (finite common-source capacity).** Fix any finite collection of
triples (F,X,Y), where F is one of those ten forms and X,Y are fixed finite
source words, possibly empty. The number of distinct counts at most N in

    { F(Xv,Yv): v any nontrivial physical source,
                   (F,X,Y) in the fixed collection }

is O(N^rho), with rho from (8).

Proof. By (14), any contributing common source has U_v<=sqrt(N/3).
Section 3 bounds the number of such source graph classes by O(N^(beta/2)).
Each class gives at most one scalar per fixed triple. Taking a finite union
proves the result. QED.

This includes every diagonal form F(v,v) and every finite bank of fixed
physical perturbations of a diagonal. It excludes the allowed singleton
boundary families of C3, whose one-source readout is linear rather than
quadratic. It also excludes two independently chosen general sources.

**Growing-word version.** If X,Y may vary with the target but each has at
most K(N) letters, the number of choices of their words is O(3^(2K(N))).
The same proof gives

    O( N^rho * 3^(2K(N)) ).                            (15)

In particular K(N)=o(log N) still gives N^(rho+o(1))=o(N). Any scheme of
this exact common-prefix type with positive-density scalar coverage must
therefore pay logarithmic independent suffix complexity. More precisely,
a uniform K(N)<=delta log N still forces density zero whenever
rho+2 delta log(3)<1. This is a statement about the declared presentations,
not a lower bound on unrestricted tournament complexity.

### Arbitrary products do not repair the fixed-bank deficit

Let A be the distinct values from a fixed collection in Theorem 2, together
with finitely many additional integers >=2. For each s>rho, the estimate
#(A intersect [1,N])=O(N^rho) gives sum_(a in A) a^(-s)<infinity by dyadic
summation. Hence the product

    product_(a in A) (1-a^(-s))^(-1)

converges and bounds the Dirichlet sum over the multiplicative closure of A.
Consequently that closure has O_epsilon(N^(rho+epsilon)) values at most N
for every epsilon>0. Taking 0<epsilon<1-rho proves density zero even after
arbitrary ordered tournament joins and singleton padding. Joins multiply
Hamilton-path counts; they need not remain strong or have two/three cycles.
This is the same Dirichlet mechanism recovered from Session 197
`dependencies/ENTROPY_BUDGET.md` (recovered from
`session197-bridges/PROOF.md`), applied to a different physical class.

## 6. Scope, validation, and remaining positive route

The sharp minimum is an all-order theorem from the exact recurrence state
cycle (3)-(4), not a fitted periodic minimum. The entropy bound counts graph
classes and therefore remains valid despite nonunique O/Z/J histories. The
norm identity uses actual source responses, and its failed involution is an
explicit operation-domain failure rather than a putative missing Hamilton
count.

Run `python3 verify.py` and `python3 -O verify.py`. The standalone checker
uses only the standard library. It compares the Pareto recursion with complete
word-row enumeration through depth ten; checks the exact recurrence-state
cycle and long iterations; checks all ten product-bound matrices; verifies
the algebraic norm identity, the outer norm move and physical chamber on a
complete U-bounded source bank; and independently counts literal source
path/cover inventories and C2 diagonal graph counts for small witnesses.
Normal and optimized JSON receipts must agree. Every test gate raises an
explicit exception. The proof's all-order claims do not use census saturation.

The useful positive lesson is precise: an exact all-odd construction through
the bilinear forms cannot be compressed into a fixed or sublogarithmically
thickened common-source norm family, even with arbitrary products afterwards
in the fixed-bank case. Independent source choices, or some genuinely different
physical scalar operation, remain necessary for this proposed route. This
neither disproves nor proves that the full independent-source C1/C2/C3 union
covers every odd target other than 7 and 21.
