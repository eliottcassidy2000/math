# Self-converse C2 tournaments: a diagonal source classification

Session 197. **PROVED from the named reviewed full-C2 normal form and source enumeration.** This is a separate corollary, not a change to those frozen proofs. It classifies self-converse graphs at every order and gives their exact scalar image and a complete bounded inverse. It does not assert cofinite scalar coverage.

## 1. Named input and physical conventions

The reviewed `BILINEAR_C2.md` proves that every strong tournament with exactly two directed Hamiltonian cycles, except the regular five-tournament R5, has a unique presentation by a type N or M and an **ordered pair of rooted graph isomorphism classes** (V,c),(W,a). Each source is nontrivial and has exactly one Hamiltonian path starting at its root. Uniqueness is at graph-class level, not at the level of response quadruples.

For a source write

    v=(D,E,U,d)=(H(V), e_c(V), R_c(V), H(V-c)).

Here R counts ordered spanning two-path covers (P,Q) with Q nonempty starting at c and P possibly empty. Let S* be the set of these actual nontrivial source quadruples. The source oracle and its physical recipes are proved in the named normal form.

The graph constructions, recalled to specify every arc, are as follows. Keep the original V vertices first. In type N, add one vertex b which beats only c among V. In type M, first add such a b, then add u which beats only b among V+b. Finally add the converse of W, with marked vertex r corresponding to a. All arcs from the left graph to this right graph point rightward except:

- type N: r beats b and c;
- type M: r beats u and c.

These are the original physical N/M constructions. Both sources have order at least two. The exceptional R5 has rows [6,12,24,17,3], H=15, and C=2.

## 2. Converse is source exchange

**Lemma 1.** For either type F in {N,M},

    converse(F(V,W)) is isomorphic to F(W,V).

**Proof.** Identify each vertex of the old V with its corresponding vertex in the new right-hand copy converse(V), and each old right-hand vertex from converse(W) with its corresponding vertex in the new left-hand W. For N, map the central vertex b to the new central vertex. For M, exchange the two central positions: old b maps to the new u, and old u maps to the new b.

The maps reverse every within-source edge because the opposite occurrence is a converse copy. All cross-source edges are uniform except r->c, which maps to the same exceptional root-to-root orientation after reversal. In N, b->c is paired with r->b. Every other incidence at b is uniform and reverses as required. In M, b->c is paired with r->u, V->u is paired with b->converse(W), and the edge u->b maps to itself with its two endpoints exchanged. All remaining central incidences are uniform. This checks every pair of vertices.

For exact array indexing, put n=|V|, m=|W| and k=1 or2. The old V index i maps to m+k+i; an old right index n+k+j maps to j. In N, n maps to m. In M, n maps to m+1 and n+1 maps to m. QED.

**Theorem 2 (graph classification).** A nonregular strong C2 tournament F(V,W) is self-converse if and only if (V,c) and (W,a) are rooted-isomorphic.

**Proof.** Equal rooted source classes suffice by Lemma1. Conversely, a self-converse F(V,W) is isomorphic to F(W,V). Uniqueness of the ordered graph presentation in the reviewed normal form forces the first rooted class to equal the second. No claim that numerical response equality implies rooted isomorphism is used. QED.

For an actual equal labelled source on both sides, Lemma1 gives an explicit anti-involution. It exchanges the two source copies. Type N fixes precisely its central vertex; type M exchanges its two central vertices and has no fixed vertex. Every nonregular C2 graph is rigid by the reviewed complete classification. Hence this anti-isomorphism is unique whenever it exists: two anti-isomorphisms would compose to an automorphism. R5 separately has five anti-isomorphisms.

## 3. Exact all-order enumeration

The reviewed rooted-source enumeration, including the unique weak singleton tail, gives

    W(z)=z(1-2z)/(1-3z+z²).

Remove the singleton source. Theorem2 permits exactly one source graph class, used on both sides. Type N adds one central vertex and type M adds two. Therefore the complete unlabelled self-converse strong-C2 generating function is

    SC2(z)=z^5+(z+z²)(W(z²)-z²)
          =z^5+z^5(1+z)(1-z²)/(1-3z²+z^4).       (1)

This is an ordinary class-counting substitution: it doubles the size of one chosen rooted source, rather than taking an unordered pair of arbitrary sources. The ordered-pair uniqueness already removes every possible coincidence between the two types and excludes R5.

Let F0=0,F1=1 be the Fibonacci numbers. For every k>=2, the number of source classes of order k is F_(2k-3). Consequently the numbers of nonregular self-converse C2 classes at orders 2k+1 and 2k+2 are both F_(2k-3). Add the R5 class at order five. The total counts at orders five through ten are

    2,1,2,2,5,5.

Since every nonregular class is rigid, its labelled contribution is n! per class. R5 contributes24 labelled graphs instead of120. These statements concern graphs up to isomorphism, not the number of distinct scalar counts.

## 4. Exact diagonal scalar image and inverse

The symmetric bilinear forms in the reviewed theorem become

    N(v,v)=(D-E)^2+2DU+3U²-2d(E+U),
    M(v,v)=D²+4DU+6U²-2EU-2d(E+U).                 (2)

**Theorem 3.** The Hamiltonian-path spectrum of self-converse strong C2 tournaments is exactly

    {15} union {N(v,v), M(v,v): v in S*}.            (3)

Every value in the right side has a physical self-converse graph: choose any one actual source realizing v and use that same source twice. The converse implication follows from Theorem2 and the established scalar formulas. Different graph classes can share a quadruple, but this does not harm either direction of the existential scalar statement.

The positive source bounds in the normal form specialize to

    N(v,v)>=3U²+1,
    M(v,v)>=6U²+2U+1.                              (4)

For an exact positive target H, handle15 and enumerate actual source states with U<=floor(sqrt((H-1)/3)); evaluate the two diagonal forms and retain equality. The source cap oracle is complete, terminates, and records a physical recipe. Thus this is a complete numerical-target inverse for (3), not a finite-order census and not a promise that every target in a residue class occurs. The sharper M bound may prune further. The crude inherited O(B^4) source-state bound already gives finite polynomial work in numeric H; it is not polynomial time in log H.

## 5. A residue obstruction at every order

**Lemma 4.** Every actual source quadruple satisfies D odd and U congruent to E modulo two.

**Proof.** Use the complete physical source oracle in the normal form. Its singleton seed (1,1,1,1) has these properties. A sink step sends (D,E,U,d) to (D,0,2U,d). For an ordinary gadget indexed by L, the constants g(L) and Q(L) are odd; the recurrence preserves this, since each increment in g at an odd step is even, and each added power in Q is even. When L>=2, its B is even, so the new E and U are both odd and the new D is odd because U-E was even. For L=1, the new E=D-E and U=D+U have the same parity, and the new D=D-E+U is odd. These cases prove the induction. QED.

**Corollary 5.** Every nonregular self-converse strong C2 tournament has

    H congruent to1 modulo4.                       (5)

Indeed E+U is even. The first formula in(2) is congruent to (D-E)^2+2DU+3U² modulo4. If E,U are even, this is1 because the terms 2D(U-E) vanish modulo4; if they are odd, (D-E)^2 is0 and the remaining two terms are2+3=1. In the second formula, D² is1, 4DU vanishes, and 6U²-2EU=2U(3U-E) is divisible by4 because U-E is even. This proves(5). R5, with H15, is the single exceptional self-converse C2 graph outside this residue class.

Thus a self-converse C2 construction cannot cover targets3 modulo4 other than15. This is a structural obstruction for that class, not an obstruction for unrestricted tournaments or for the full C2 class.

## 6. Equal responses do not make a particular graph self-converse

Two actual rooted sources of order seven have rows

    V=[124,65,122,114,98,2,32], root5,
    W=[124,113,122,66,104,8,32], root5.

Both have quadruple (27,9,27,9), but their unique Hamilton paths from the marked root give different canonical rooted adjacency arrays. They are not rooted-isomorphic. The physical graphs N(V,W) and N(V,V) both have H3321 and C2. The former is not self-converse; the latter is. Literal cycle-based anti-isomorphism checks in the verifier confirm both statements without using the source factorization. This prevents replacing Theorem2's rooted graph condition by equality of four numbers.

## 7. Verification and dependencies

The verifier constructs the rooted source classes from their canonical initial-SCC and singleton-tail decomposition, checks start counts and unique-root-path canonical labels, and builds every ordered N/M presentation through graph order ten. It checks the literal converse permutations, exact H/C counts, source equality versus all anti-isomorphisms found by Hamilton-cycle matching, the generating-function coefficients and the residue claim. A separate order-fifteen response-collision hostile tests the graph/scalar distinction in Section6. Bounded diagonal support is included only as a finite control of the exact inverse.

Run `python3 verify.py` and `python3 -O verify.py` in a temporary copy; the verifier writes VALIDATION.json. Every check raises an explicit exception rather than relying on Python assertions. The all-order conclusions follow from the proofs and named classifications, not from stabilization of the finite graph list.

Frozen named dependencies:

- BILINEAR_C2.md, SHA256 `5eab9deb2871819d8a1b5684c4987f2f4da807adc55063a30675d531c5ad3813`, published at commit f01df1766885c807aa4acae17685c2f4777f3ead under session197/bilinear-c2-inverse/PROOF.md.
- C2_COMPLETE.md, SHA256 `4cc0bf11efc435d7a2746daf6d8f1a58c49c94e57e608911b24464facdd1057d`, for the complete class and rigidity.
- FIBONACCI_CLASS.md, SHA256 `d95a133e360b96589009ebd8a1c8449525c4a2847a88b33a4b3f2b839739e63e`, for rooted strong-source enumeration; the weak-tail factor is explicitly proved in the bilinear dependency.

The converse-swap observation was independently suggested and checked by the low-cycle worker; its separate Hamilton-cycle isomorphism census through order ten agrees with(1). This packet supplies the explicit all-order maps, diagonal classification, exact inverse and residue consequence, with no literature-priority claim.
