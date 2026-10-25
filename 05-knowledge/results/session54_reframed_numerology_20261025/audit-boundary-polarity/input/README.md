# The native Paley-7 boundary is multiplicative F4 projectivization

Status: PROVED explicit commutative square, FINITE-EXACT verification of every blocker pair and symmetry. This supplement refines the frozen63 bridge without changing that packet. It uses its independently verified63 blocker masks and reconstructs every map.

## Adjusted field coding

Keep the Paley tournament on F7, G={x->ax+b:a in {1,2,4}}, and converse iota:(u->v)->(-v->-u). Identify F64× with exponents e modulo63 using a primitive element g. Change only the original source/sink phase choice from1 to10:

- Star at c: e=9c.
- Common-out center c: e=9c+21.
- Common-in center c: e=9c+42.
- Source pair ordered u->v, v-u=2^k: e=9u+10*16^k.
- Sink pair ordered u->v, v-u=2^k: e=9v-10*16^k.

Call this bijection Psi. The translation, multiplier, and converse actions remain e->e+9, e->16e, and e->-e. The same proof as the original phase choice applies: residues modulo9 separate the five families, translations cover the seven lifts, and the source/sink formulas exchange under converse. The checker verifies all63 values, every affine map, and converse. The old baseline1 code is not fiber-compatible; an explicit witness is recorded rather than treating the two codes as identical.

## Native boundary and the commuting square

Let partial(B) be the XOR of the two endpoint indicators over all deleted arcs. It is a nonempty even subset of the seven vertices. Define a map rho from Z/21 to its21 boundary images by, for c in F7,

    rho(9c)   = V\{c},
    rho(9c+7) = {c} union N^-(c),
    rho(9c-7) = {c} union N^+(c).

These are disjoint seven-element families and exhaust the21 residues, since the three displayed forms occupy the three residues modulo3. Then

    partial(B) = rho(e(B) mod21).

Consequently partial(B)=partial(B') iff e(B)=e(B') modulo21. This is checked for all3969 ordered pairs, not merely by counting fibers.

For a direct proof, the three centered seven-families all have boundary V\{c}; their exponents differ by21. A source pair (u,v) with d=v-u has boundary the closed in-neighborhood at c=u+5d. Check the base pair(0,1), whose boundary is {1,3,4,5}, and transport by G. Similarly a sink pair has boundary the closed out-neighborhood at c=u+3d. For d=2^k,

    10*16^k = 3*2^k+7 modulo21,

which gives exactly the two rho formulas. Thus the square commutes by formulas, not only by finite evidence.

Since F4×=<g^21>, equality modulo21 means multiplication by an element of F4×. Therefore the native boundary is literally the quotient

    F64× -> F64×/F4×

after the stated adjusted coding and rho. The21 boundary images can consequently be identified with the21 points of PG(2,4), because F64 is three-dimensional over F4.

## What this does not identify

A fiber aF4× has three elements, but it represents ONE projective point, not a projective line. An actual PG(2,4) line has five points. No native blocker-incidence relation is asserted to be its projective line relation. Native arc XOR and Hamming distances remain different from the transported field operations.

Converse descends to the point permutation [x]->[x^-1], but it is NOT a collineation of this transported PG(2,4). The field checker constructs all21 projective lines from two-dimensional F4 subspaces. None of their inverse images is a line. One explicit witness in exponent-mod21 coordinates is

    line:          {0,1,6,8,18},
    inverse image: {0,3,13,15,20}, not a line.

This is an independently tested limitation, not an assumed incidence preservation.

## Inversion's exceptional derivative quartet

Extend inversion to F64 by f(0)=0. For each nonzero a, the equation

    f(x+a)+f(x)=b

has four solutions at b=1/a, two solutions for30 other b, and no solutions for33 b. The four exceptional solutions are exactly aF4. Hence its three nonzero solutions are exactly the multiplicative fiber recovered above.

Proof: x=0,a contribute precisely when b=1/a. For x outside those two values, put y=x/a. The equation becomes y^2+y=1/(ab). A quadratic y^2+y=c has two roots exactly when the absolute trace of c is zero. Among63 nonzero c,31 have trace zero. At ab=1, trace(1)=0 because the extension degree is six; the two roots are the nontrivial elements of F4, and together with0,1 they form F4. The other30 trace-zero cases each give two roots; all remaining b give none. The code checks all63 nonzero a and all64 inputs, tallying every one of the64 possible b via exact polynomial arithmetic, with modulus X^6+X+1 and primitive element2.

This derivative statement explains exactly the same scalar fiber. It supplies no cryptographic conclusion or additional graph-incidence claim.

## Evidence and reproduction

verify_base.py and blockers.json are copied, pinned dependencies from the previously audited63 packet; the full census is not silently recomputed by this supplement. verify.py validates the new bijection, all affine/converse equations, the quotient decoder, all-pairs fiber equivalence, and a baseline1 hostile. algebra.py independently constructs F64/F4, all projective lines, the inverse-line hostile, and the full derivative census.

Run replay.py --self-test from any directory. It checks the exact manifest file set and digests, normal/-O byte-identical results, and file-add/content-mutation hostiles. This packet does not change the original63 or P11 packets.

## Typed incidence repair

TWO-SHEET.md proves and checks a separate positive repair: two typed copies of the21 boundary labels form PG(2,4) point/line incidence, with converse as trace polarity rather than a point collineation. Its lift through the actual three-element boundary fibers has126 vertices and945 edges. This relation is explicitly new, and an equal-native-intersection hostile prevents conflating it with native arc incidence.
