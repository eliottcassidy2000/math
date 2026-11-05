# Exact defect and exceptional physical fibers for singleton reassociation

Session 198. Status: **PROVED** algebraic identity and physical fiber theorem;
**FINITE-EXACT** minimal hostile control. This diagnoses a proposed scalar
reassociation, not all target-dependent source-pair constructions.

Let s=(D,E,U,d) be an actual nontrivial start-one source and let J=M_1 be its
physical singleton birth, with its newborn root retained:

    J(s)=(D-E+U, D-E, D+U, D).

The reviewed symmetric C2 forms N and M accept two nontrivial sources. It is
tempting to combine scalar symmetry with source-word decomposition and transfer
a singleton birth from the right source to the left without changing the
output. That implication is false. The exact lost scalar is a rank-two wedge
for N and a rank-four alternating form for M.

## 1. Exact commutation defect

Define two positive source observables

    p(s)=D+d,        rho(s)=D-E+U-d.

For any four-vectors v,w, not merely physical sources, one has

    N(v,Jw)-N(Jv,w)=rho(v)p(w)-p(v)rho(w).          (1)

For the second C2 form,

    M(v,Jw)-M(Jv,w)=rho(v)p(w)-p(v)rho(w)
                         +2(U_v D_w-D_v U_w).    (2)

These are exact integer identities. For physical nontrivial sources p>0, so
(1) vanishes if and only if rho(v)/p(v)=rho(w)/p(w).

To verify the algebra without any graph inference, write B_N,B_M for the
symmetric matrices of N,M in coordinates (D,E,U,d), and write A for J:

    B_N=[[1,-1,1,0],[-1,1,0,-1],[1,0,3,-1],[0,-1,-1,0]],
    B_M=[[1,0,2,0],[0,0,-1,-1],[2,-1,6,-1],[0,-1,-1,0]],
    A=[[1,-1,1,0],[1,-1,0,0],[1,0,1,0],[1,0,0,0]].

Then B_N A-A^t B_N is

    [[0,1,-1,2],[-1,0,0,-1],[1,0,0,1],[-2,1,-1,0]],

which is rho^t p-p^t rho and has rank two. The M defect differs only in its
(D,U) entry -3 instead of -1, and its (U,D) entry 3 instead of 1. Its determinant
is four, so it has rank four. Matrix multiplication proves (1)-(2) at every
source size. No abstract cone point is used to claim physical realizability.

## 2. The complete physical zero-defect fiber at the smallest source

Let t=(1,0,2,1), the transitive two-vertex source. Since p(t)=rho(t)=2,

    N(t,Jw)=N(Jt,w)  iff  U_w-E_w=2d_w.            (3)

**Theorem.** Among actual nontrivial sources, the right side of (3) holds
precisely for t and for

    w=(3d,d,3d,d)=M_2(z),                         (4)

where z is an actual balanced source, meaning H(z)=R_root(z)=d. The singleton
z is allowed and yields the triangle. The assertion does not claim that every
odd integer d has a balanced source.

Proof. Every weak source has E=0 and is a strong start-one core or singleton
followed by a nonempty transitive tail. The universal physical inequality
U>=D+d-E and (3) give D<=d. Since d<=D, this forces D=d. A nontrivial strong core
has E_core>=1 and D-E_core>=d, which contradict D=d. Thus the core is a singleton,
D=d=1, and (3) gives U=2. This is exactly t.

For a strong source use the reviewed exact last-gadget alternatives. In a
nonsingleton-last presentation, d=E and D>=3E; in a singleton-last presentation,
d=U-D+E and D<=3E. The triangle may be placed in the first case.

In the nonsingleton case, (3) gives U=3E. Since U>=D>=3E, all are equal. In the
singleton case, (3) and d=U-D+E give U=2D-3E. The inequalities U>=D and D<=3E
again force D=U=3E and d=E. Thus every strong member has the form (4).

For D>3, the reviewed last-step theorem at D=3E,U=D forces the last ordinary
gadget to be M_2. Its predecessor has H=E=d and cover response
(U-E)/2=d, hence is balanced. The triangle has the singleton predecessor under
M_2. Conversely, applying M_2 to any actual balanced source
(d,E_old,d,d_old) gives exactly (3d,d,3d,d), regardless of its older endpoint or
deletion data. These are actual strong sources and satisfy (3). QED.

The zero-defect fiber is unbounded but highly constrained. For instance,
repeated M_2 moves supply d=3^j. That family does not justify commuting J across
an arbitrary source pair or arbitrary word.

## 3. Minimal physical failure

Take actual transitive sources

    v=t=(1,0,2,1),       w=(1,0,4,1).

They have respectively two and three vertices, rooted at their first vertices.
Then Jv=(3,1,3,1) is the triangle and Jw=(5,1,5,1). The physical C2 N grafts give

    N(v,Jw)=41,          N(Jv,w)=45.               (5)

Both graphs have seven vertices and exactly two Hamilton cycles. Formula (1)
gives their difference -4: the source ratios are 1 and 2, not equal.

For M the corresponding eight-vertex physical grafts give

    M(v,Jw)=85,          M(Jv,w)=93.               (6)

Here the difference in (2) is -8.

These are minimal by total order of the two original nontrivial sources: order
four forces both to be the unique rooted two-vertex source and gives zero
defect. At total order five the only rooted source responses at order three
are the transitive (1,0,4,1) and the triangle (3,1,3,1). The first supplies (5),
while the triangle lies on the N zero-defect fiber. The verifier independently
enumerates all rooted labelled tournaments through order three to check this
finite statement, then constructs and counts the four actual graft graphs.

## 4. Scope and continuing obligation

This rules out the specific scalar-preserving move
N(v,Jw) -> N(Jv,w) outside its exact ratio fiber, and likewise diagnoses M.
Symmetry N(v,w)=N(w,v) does not imply gadget commutation. The invariant required
by that proposal is rho/p for N; M needs the extra alternating (D,U) response.

It does not rule out replacing a failed transfer by another physical source,
a different gate, a longer word, or a target-dependent construction. In
particular it neither proves nor refutes the conjectural equality between the
full C2 scalar image and the generalized one-source families, and does not
settle coverage of the remaining odd integers.

Run `python3 verify.py` and `python3 -O verify.py`. The verifier has no removable
assertions and uses literal adjacency constructors plus independent subset path
and cycle counts for the hostile graphs. It also checks the algebra and physical
zero-defect fiber throughout its explicitly bounded forward source inventory.
The proof, not that finite test, establishes the all-order identities and fiber.
