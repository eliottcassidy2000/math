# Exactly which rooted contexts transport the neutral repair

**PROVED from an explicit path-cover bijection and the inherited weighted-cut
identity; FINITE-EXACT independent controls.** The result concerns one-sided
uniform cross edges. Arbitrary cross matrices and more general repairs are
outside its scope.

## 1. The response of a one-sided rooted amalgamation

Let A and B share only their named root r. Keep their internal arcs and
orient every edge from X=A-r to Y=B-r. Write

    z_A=(h_A,s_A,q_A,F_A),       z_B=(h_B,s_B,q_B,F_B),
    e_A=# Hamilton paths of A ending at r,
    G_A=# ordered two-path covers whose first path starts at r,
    P_X=P_(A-r)(2).

Assume |A|≥2. The exact rooted profile of their amalgamation T is

    H(T)=(G_A-q_A)(F_B-q_B)+(h_A-e_A)q_B
                              +q_A(h_B-s_B)+e_A s_B,
    s_T=s_A q_B,
    q_T=q_A q_B,
    F_T=P_X(F_B-q_B)+F_A q_B.                       (1)

**Proof.** A Hamilton path can cross from Y back to X only through r.
Partition by whether it has a Y vertex before r and an X vertex after r.
If both occur, its form is X1,Y1,r,X2,Y2, with Y1 and X2 nonempty;
restriction gives exactly(G_A-q_A)(F_B-q_B) choices. If only X2 occurs,
the A path does not end at r and all of Y follows, giving(h_A-e_A)q_B.
If only Y1 occurs, all of X precedes a B path not starting at r, giving
q_A(h_B-s_B). If neither occurs, concatenate an A path ending at r and
a B path starting there, giving e_A s_B. Empty X1,Y2 are included.

For a path starting at r, all of X must precede Y, proving the s formula;
deleting r gives X→Y, proving q. For F, if its first bin contains a Y
vertex, the two X parts form any ordered two-bin cover and the B part
has a vertex before r, giving P_X(F_B-q_B). Otherwise all Y lies in the
second bin, giving F_Aq_B. Restriction reverses every construction. QED.

These formulas retain boundary responses rather than only H(A),H(B).
They are the precise gluing analogue of the older HN boundary-count method.
The extra e_A,G_A,P_X are needed on the context side.

## 2. A sharp classification for a pure deletion repair

Suppose B and B' have the same h_B,s_B,F_B, but q_B'=q_B+w with w≠0.
Let T,T' be their one-sided amalgamations with the same A. Put

    delta_out(A,r)=G_A-h_A-q_A+e_A.

The inherited [weighted-cut identity](../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md#3-the-defect-detects-indegree-and-the-legal-inverse),
applied to the converse of A, gives

    delta_out≥d_out(r)(d_out(r)-1),
    delta_out=0 iff d_out(r)≤1.                     (2)

Each pair weight is strictly positive; for degree0 or1 the sum is empty.
Subtracting (1) gives the exact transported difference

    Delta(H,s,q,F)
       =w(-delta_out, s_A, q_A, F_A-P_X).            (3)

**Theorem.** Both the current count H and the first root-child target
H+F-s remain unchanged **if and only if r is a sink of A**.

Proof: preserving H forces delta_out=0, hence degree0 or1. At degree0,
s_A=0 and F_A=P_X, so(3) preserves both counts. This is precisely the
ordinal-prefix context X→B used in REPAIR23.

At degree1 let u be the sole out-neighbour and U=A-r. Then

    s_A=s_U(u),       F_A=P_U(2)-F_U(u).

The first identity deletes the initial r. For the second, append r to
the first bin of an arbitrary ordered cover of U, including an empty
first bin; the forbidden choices are exactly those whose first bin ends
at u. Thus the first-birth target changes by

    -w(F_U(u)+s_U(u)),                              (4)

which is nonzero because F_U(u)≥H(U)>0. Degree≥2 already changed H.
This proves the iff. The singleton A boundary simply gives T=B and also
has r as a sink, so the theorem extends to that boundary separately. QED.

## 3. A minimal context hostile

Take A to be the directed edge r→x and insert the23-core pair from
[REPAIR23](REPAIR23.md). Their amalgamations have profiles

    (81,9,9,111) and(81,11,11,109).

The current count81 survives, but the next targets are183 and179.
This is the ordinary R visibility of the inherited pure-q repair in its
smallest context. The theorem extends it to **every** one-sided context:
no clever choice of a nontrivial feedback block in A can make this pure-q
change neutral for both readings. Selective cross edges or another repair
direction are necessary to escape this obstruction.

## 4. Relative integer geometry and the physical lifting obligation

For fixed odd current/next targets h,M, the formal profile fibre has

    F-s=M-h,       Delta=(0,v,w,v).

Actual deletion counts are odd, so physical differences have w even.
For chain-entry coordinates c=h-s+4q-F and d=F-q, the fibre maps to

    Delta(c,d)=(4w-2v,v-w).

Writing w=2t, the integer matrix has determinant-4. Its image is exactly
the affine congruence

    c+2d≡M+2 (mod4).                               (5)

Indeed q=(c+2d-M)/2 must be odd, and
s=(c+4d+2h-3M)/2 then recovers an integer profile. This is an exact
lattice calculation; it imposes no positivity, graph realization or
degree constraints. The one-arc H9 repair realizes(v,w)=(1,2), and the
23-core repair realizes(0,2). Their numerical directions generate the
whole physical parity lattice, but they live at different actual bases.
Their integer span does **not** let them be composed at an arbitrary target.

This is the useful distinction from the planar-Jacobian response-fibre
theorem: there, a proved shear action fills the whole nonempty fibre under
its hypotheses. Here(1)–(4) prove an obstruction to transporting even one
generator through a specified class of contexts. Integral compatibility
is still only a necessary condition for a physical lift, as in the older
Hodge and cycle-displacement discussions.

Reproduce: `python3 -B 05-knowledge/results/surplus_selector_20261113/context_transport.py`.
