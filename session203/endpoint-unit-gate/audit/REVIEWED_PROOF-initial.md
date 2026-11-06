# A mod-four obstruction to a fixed rooted unit-correction gate

**PROOF CANDIDATE, pending independent audit.** The theorem below is all-order in its stated physical context class. It proves no missing unrestricted Hamilton-path count and no eventual spectrum threshold.

## 1. Physical contexts and the target

A rooted-module context consists of a fixed nonempty tournament K, a distinguished vertex q, and a set M contained in K minus q. On input a rooted tournament (U,v), retain K and U minus v, identify v with q, keep every input edge inside U, and orient every x in U minus v toward precisely the core vertices in M. Every core vertex outside M and q beats every such x. Write T_K,M(U,v) for this actual tournament. No core order bound is imposed.

Two contexts may use wholly different cores, orders, edges and masks. Thus old core edits, vertex deletions and replacements are permitted between the parent context P and the child context C. The restriction is that each context copies the input root once and every other core vertex is uniform toward U minus v. The cores and masks are fixed across the input family. One may also complement every copied root incidence in either context; the argument is unchanged.

**Theorem.** Fix any two such contexts and any even integer k. There is no constant c such that

    H(C(U,v)) = k H(P(U,v)) + c

for every strong rooted tournament U having start_v(U)=1. More strongly, four explicit five-vertex inputs already refute the identity, regardless of the two core sizes. There is no uniform bound on their output orders.

In particular a fixed rooted-module pair cannot produce a unit correction H(child)=2H(parent)+1 even after restricting the source to an actual, physically realized unit-start class. The conclusion does not assume that a unit endpoint is freely attainable for a prescribed total count.

The same conclusion holds if c is replaced by any affine function of the two free input bits used below, since its mixed difference vanishes.

## 2. Exact twin curvature

Let W be a tournament, let x,y be distinct vertices with identical incidences to every vertex of W minus {x,y}, and adjoin an apex q. Hold its other incidences fixed; put q->x iff b_x=1 and q->y iff b_y=1. Let F(b_x,b_y) be the resulting Hamilton-path count. Contract the pair {x,y} to a single vertex with their common outside incidences, obtaining the tournament Q.

**Lemma.** For either orientation of x-y,

    F(1,1)-F(1,0)-F(0,1)+F(0,0) = -2 H(Q).        (1)

**Proof.** A Hamilton word whose apex is an endpoint depends on at most one of b_x,b_y. The same is true of an internal apex unless its two neighbors are exactly x,y. Let E_xy count words of W with x immediately followed by y and every step except that exposed step directed correctly. Inserting q there contributes

    E_xy (1-b_x)b_y.

The opposite order contributes E_yx(1-b_y)b_x. Hence the mixed coefficient is -(E_xy+E_yx).

Because x,y are twins outside their pair, contracting the exposed block is a bijection from either class of words to Hamilton words of Q. The internal orientation is irrelevant because that one step is exposed. Thus E_xy=E_yx=H(Q), proving (1).

This is the twin equality case of the inherited exposed-path capacity formula in `session203-complement-gadget/PROOF.md`, not a new general cut identity. Complementing both bits leaves the mixed difference unchanged. Q is an actual tournament, not a formal contraction with inconsistent incidences.

## 3. Four strong unit-start sources force the obstruction

Let U_ab have vertices v,x_1,x_2,x_3,x_4. The four x vertices form the transitive chain in that order. The root incidences are

    v->x_1,    v->x_2 iff a=1,    v->x_3 iff b=1,    x_4->v.

These are four explicit physical sources, with (a,b) in {0,1}^2. Every one is strong: v,x_1,x_2,x_3,x_4,v is a directed spanning cycle. Every one has exactly one Hamilton path starting at v, namely v,x_1,x_2,x_3,x_4, because deleting v leaves a transitive tournament with only that Hamilton order. Their ordinary H values, in order 00,01,10,11, are 9,9,11,9. In particular the two freely varied bits do not manufacture four arbitrary source profiles.

In either output context, remove q. The two adjacent chain vertices x_2,x_3 then have identical incidences to every remaining vertex: x_1 beats both, both beat x_4, and every fixed-core vertex treats all input nonroot vertices uniformly. They are therefore a literal twin pair. Let Q_P,Q_C be the corresponding contractions after deleting q. Lemma (1) gives

    Delta H(P(U_ab)) = -2H(Q_P),
    Delta H(C(U_ab)) = -2H(Q_C).

Consequently

    Delta[H(C(U_ab))-kH(P(U_ab))]
      = -2[H(Q_C)-kH(Q_P)] = 2 mod 4              (2)

when k is even. The final equality uses only the inherited parity theorem that every tournament has an odd number of Hamilton paths: THM-001, `01-canon/theorems/THM-001-redei.md`. A constant or affine correction has mixed difference zero, contradicting (2). This proves the theorem.

The two quotient tournaments may have different orders and different counts. Their oddness, not any equality or structural identification, is what prevents cancellation.

## 4. A positive parity-compatible endpoint identity

The obstruction has a useful sharp contrast. For an arbitrary rooted nonempty U, write h=H(U), s=start_v(U), and let L count ordered spanning two-path covers (P,Q) with P nonempty and ending at v, while Q may be empty.

Construct P(U) by adding a vertex a with v->a and a->(U minus v). Construct C(U) by adding two vertices a,b with

    a->v->b->a,
    a,b -> (U minus v).

All edges of U remain unchanged. Then

    H(P(U)) = h+L-s,
    H(C(U)) = h+2L,
    H(C(U))-2H(P(U)) = 2s-h.                    (3)

These are elementary specializations of the inherited ordered-root-amalgamation formula in `session197-odd-cycle-path-polynomials/ROOTED_ORDERED_AMALGAMATION.md`; no novelty claim is made for that calculus or these immediate specializations. For a direct interpretation, amalgamate on the left either the rooted edge v->a or the rooted triangle a->v->b->a, with every new nonroot vertex beating U minus v. The edge has rooted data (h,e,s,R,d)=(1,0,1,2,1), and the triangle has (3,1,1,3,1). Substitution gives (3), including singleton U.

For actual unit-start inputs the correction is 2-H(U), not 1. This is an endpoint-sensitive physical identity, but the source Hamilton count remains coupled. For instance a rooted directed triangle gives parent H=5 and child H=9, with correction -1. The four U_ab above give a nonzero mixed difference for 2-H(U_ab), exactly as the theorem requires.

The larger family-C rooted lift discovered concurrently has the same universal correction 2s-h. Its independent exact tensor check is recorded in `KNOWN_C_TENSOR.json` and `FAMILY_C_AUDIT.md`. It is a connection to the existing local repair, not a unit extractor.

## 5. Scope and stopping boundary

This theorem forbids a fixed one-copied-root module interface for an even-multiplier recurrence with a constant correction. It allows arbitrarily large, unrelated parent and child cores and does not assume induced preservation of one core inside the other.

It does not exclude a target-dependent menu, two or more sensitive copied vertices, incidence rules depending on the positions of input vertices, changing the input's internal arcs, nonlinear encodings, or recurrences restricted to a smaller source family avoiding this four-input configuration. None of these exits is proved sufficient. In particular failure of a constructor does not establish nonexistence of any scalar target in the unrestricted tournament spectrum.

The inheritance board is: physical unit-start source; exposed-path capacity; literal twin contraction; parity modulo four; and correlated root-amalgamation responses. The closest proved mechanism is strict pair curvature. The corrected near miss is replacing a physical endpoint coordinate by a freely selectable integer. The underused sidecar is the parity of the actual contracted tournament.

## 6. Reproduction

    python3 verify.py > NORMAL.json
    python3 -O verify.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

The standalone verifier checks the twin identity by literal Hamilton-path DP for every labelled core through order4 and every uniform module mask, using all four strong unit-start inputs and both copied-root polarities. It also verifies the small positive constructors on every rooted labelled tournament through order5, independently counting L by subset path covers. Optional `tensor_search.cpp` and `check_known.py` are exploratory/algebraic companion checks, not dependencies of the all-order proof.
