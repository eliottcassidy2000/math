# A proposed precise lambda-cubed route: Frobenius over F5

**Status:** elementary identities PROVED from cited standard Frobenius facts;
the finite model, extension-count, magic-arrangement and quotient statements
are FINITE-EXACT. This is a **proposed precise interpretation** of the unresolved
lambda-cube cue, not confirmation of the user's intended referent, not a new
theory of elliptic curves, and not a solution of the magic-square-of-squares
problem. No global originality claim is made.

## 1. Inheritance and the precise new connection

At repository revision `74e81ed37c9175d94e8fad57fc5809cf78cee572`, Session143's
`magic_cycle/FULL_CUBIC_CLASSIFICATION.md` establishes the following classification:
for all real A,B,C,D, the indexed values of

    A t^3 + B t^2 + C t + D,    t = -4,-3,...,4,

are magic-arrangeable if and only if B=0 and either A=0 or C=-15A. Repetitions
are included. The nonaffine family is therefore exactly A(t^3-15t)+D.
This is the closest proved mechanism; its existing source snapshot is retained.

The same packet's `lambda_octonions/PROOF.md` leaves the intended lambda-cube
referent UNRESOLVED. Its signed octonion cube is not reinterpreted here.
Its `magic_cycle/sources/S57_synthesis.md`, section 6, already records a
**different classical elliptic bridge**: rational square-magic reduces to a
guarded arithmetic progression of eligible doubled x-coordinates on one
congruent-number curve. We do not rebrand that construction, nor replace its
simultaneous rational-square conditions by finite-field point counts.

There is also an older exact lexical hit in
`04-computation/experiments/sixth_clock_branches_20261004.py`, function
`sextic()` (lines 200-205 at the pinned revision). Its check named "lambda
cube" concerns alpha=1+3*phi in Q(phi), phi^2=phi+1, with alpha^2=5*phi^4 and
norm(alpha)=-5. That negative-norm construction is not identified with the
positive determinant-5 Frobenius pair below. The source is retained for
provenance, not rerun as part of this bounded finite-field task. Neither hit
confirms what the user intended by the phrase.

This bounded follow-up connects the exceptional polynomial t^3-15t to the
third power of Frobenius. The source is a nonsingular short Weierstrass model
over F5; the map is

    E -> n=#E(F5) -> t=6-n -> Tr(Frob_5^3)=t^3-15t
      -> #E(F125)=126-(t^3-15t).

All n in {2,...,10} occur. On these nine scalar states the final map is
bijective, and its values are magic-arrangeable. It forgets the equation,
isomorphism class, group structure and individual points. A chosen model is a
sidecar, not a canonical section. This transfers a **numerical magic** predicate,
not squarehood, octonion multiplication, or a curve-level cycle action.

## 2. Why the cube has coefficient -15

Let pi be the q-power Frobenius of an elliptic curve over Fq. The standard
facts used here are: pi has degree q; its dual satisfies pi*hat(pi)=[q];
pi+hat(pi)=[t]; and #E(Fq)=q+1-t. One way to get the last equality is that
1-pi is separable, its kernel is E(Fq), and

    [deg(1-pi)]=(1-pi)(1-hat(pi))=[1+q-t].

Multiplication by pi in pi+hat(pi)=[t] now gives

    pi^2 - [t] pi + [q] = 0.

Equivalently its characteristic polynomial on the two-dimensional rational
ell-adic Tate module, ell distinct from the characteristic, is X^2-tX+q.
The standard trace/determinant and extension-field statements are given in
[Snowden, Lecture 3, Point counting](https://public.websites.umich.edu/asnowden/teaching/2013/679/L03.html);
the separable-kernel argument is also in
[Sutherland, Lecture 7, Theorem 7.3](https://math.mit.edu/classes/18.783/2025/LectureNotes7.pdf).
These geometric inputs are cited, not independently formalized here.

Let lambda,mu be the two characteristic roots. Then lambda+mu=t and
lambda*mu=q. Consequently

    lambda^3+mu^3 = (lambda+mu)^3-3 lambda mu (lambda+mu)
                 = t^3-3qt.

More generally s0=2, s1=t, and sr=t*s(r-1)-q*s(r-2). Applying the point-count
formula to pi^3 gives #E(F(q^3))=q^3+1-s3. At q=5 this is exactly

    #E(F125) = 126-t^3+15t.

The checker verifies the characteristic identity and cubic trace symbolically
for the companion matrix [[0,-5],[1,t]] over Z[t], not merely at nine samples.

## 3. Complete finite census and a magic square

There are 25 equations y^2=x^3+ax+b with a,b in F5. Exactly five fail
4a^3+27b^2 != 0: (0,0),(2,2),(2,3),(3,1),(3,4). Thus 20 are nonsingular.
The complete direct count, including the unique point at infinity, is:

| a | b and #E(F5), for nonsingular models |
|---|---|
| 0 | 1:6, 2:6, 3:6, 4:6 |
| 1 | 0:4, 1:9, 2:4, 3:4, 4:9 |
| 2 | 0:2, 1:7, 4:7 |
| 3 | 0:10, 2:5, 3:5 |
| 4 | 0:8, 1:8, 2:3, 3:3, 4:8 |

The checker saves every F5 affine point. It independently counts all 20
curves over F125=F5[z]/(z^3+z+1), both by a square-value histogram and by all
125^2 coordinate pairs. The modulus has values 1,3,1,1,4 on F5, hence no root
and, being cubic, is irreducible. Field inverses, 125th powers and Euler's
criterion are checked as arithmetic controls. Every direct extension count
agrees with the Frobenius formula.

| t | -4 | -3 | -2 | -1 | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|---|---|---|---|
| #E(F125) | 130 | 108 | 104 | 112 | 126 | 140 | 148 | 144 | 122 |

These nine distinct integers have the arrangement

    108 148 122
    140 126 112
    130 104 144

All eight lines sum to 378. The trace indices in the corresponding cells are
[-3,2,4; 1,0,-1; -4,-2,3]. There are exactly eight arrangements, the dihedral
variants: the checker uses the exhaustive zero-centered form

    u     -u-v    v
    v-u    0     u-v
    -v     u+v   -u.

Only 144 among these entries is an integer square. In particular this is
**not a magic square of squares**.

For the user's cycle, the numerical map sends

    (4,8,3,7,2,6,10,5,9,4)
    -> (148,104,144,112,122,126,130,140,108,148).

This is a relabeling of an abstract nine-cycle, not a proved elliptic-curve
operation realizing its successor. Complement n->12-n corresponds to t->-t
and sends the extension count N to 252-N.

**Quadratic hostile control.** Direct counting in F25=F5[w]/(w^2-2) verifies
#E(F25)=36-t^2 for every model. Its nine indexed values have sum 264, so the
required magic center would be 264/9=88/3, not an integer. This multiset cannot
be arranged as an integer magic square. Extensions beyond degrees 1,2,3 are
outside this follow-up.

## 4. The important Hecke boundary

The third symmetric-power trace is a different expression:

    Tr(Sym^3(pi)) = lambda^3+lambda^2 mu+lambda mu^2+mu^3
                 = (t^3-3qt)+qt = t^3-2qt.

At q=5 it is t^3-10t, not t^3-15t. Likewise, for a normalized weight-two
eigenform with trivial character, at a prime p not dividing its level, the
local Euler factor gives a_(p^r)=a_p*a_(p^(r-1))-p*a_(p^(r-2)), with a_1=1.
Thus a_(p^3)=a_p^3-2p*a_p; the literal scalar cube a_p^3 is a third quantity.
See [Sutherland, Lecture 24, sections 24.4-24.5, especially Theorem 24.27](https://math.mit.edu/classes/18.783/2025/LectureNotes24.pdf).
The power-trace recurrence instead starts with s0=2: equal recurrence
coefficients do not imply equal sequences.

The checker verifies both symbolic polynomials, and an independent exhaustive
arrangement test finds no magic square for either {t^3-10t} or {t^3} on
t=-4,...,4. The inherited all-real cubic classification explains this
failure structurally. No weight-two level-13 space, particular modular form,
or new Hecke correspondence is assumed to exist.

**Formal-operator corollary (why 5 is forced).** Fix a real q and nine real
or complex 2-by-2 matrices M_t with trace t=-4,...,4 and common determinant q.
Such a family always exists formally, for example M_t=[[0,-q],[1,t]]. The
characteristic identity gives Tr(M_t^3)=t^3-3qt and
Tr(Sym^3(M_t))=t^3-2qt. Applying the inherited all-real cubic classification,
with leading coefficient 1, shows that the first indexed multiset is
magic-arrangeable **if and only if q=5**, whereas the second is
magic-arrangeable **if and only if q=15/2**. In particular no integer common
determinant gives the latter property on these nine trace states. This is a
formal matrix-family corollary, not a claim that all nine traces occur for
elliptic curves over arbitrary finite fields. The F5 census is the separate
realization step that makes the q=5 case arithmetic.

## 5. Matching small counts does not identify tournaments and curves

In characteristic distinct from 2 and 3, an isomorphism of short Weierstrass
models fixing infinity is a scaling: comparing the y and x^2 coefficients
in a general Weierstrass change of variables forces all translations and
shears to vanish. Hence the equivalence is

    (a,b) -> (u^4 a,u^6 b),   u in F5*.

Here u^4=1 and u^6 is either 1 or -1. The 20 models yield 12 classes: two
classes with a=0, four with b=0, and six with ab nonzero. Their trace fibers
have sizes 2 at t=-2,0,2 and 1 at each of the other six traces.

The independently regenerated census of all 1024 labeled five-vertex
tournaments instead gives 12 isomorphism classes and 9 adjacency spectra,
with fiber sizes 3,2,1,1,1,1,1,1,1. Thus **no bijection of the 12-object sets
can identify equality of curve traces with equality of tournament spectra**:
the equivalence relations have different block sizes. This excludes a
quotient-preserving identification, not arbitrary set bijections or every
possible relation between curves and tournaments.

The nonsquare twist d=2 is (a,b)->(4a,3b). It negates t by the quadratic
character point-count formula and acts freely on the 12 classes. If a!=0,
it changes a to -a, which scaling cannot undo; if a=0, the two b classes are
swapped. It therefore has six orbits. Tournament converse instead fixes
eight classes and has ten orbits. No bijection can intertwine these two
involutions. The retained theorem is the scalar Frobenius-cube magic map,
not a curve/tournament or octonion equivalence.

## Reproduction and stopping boundary

From this directory:

```sh
python3 check_bridge.py > /tmp/frobenius.normal.json
python3 -O check_bridge.py > /tmp/frobenius.optimized.json
cmp /tmp/frobenius.normal.json result.json
cmp /tmp/frobenius.optimized.json result.optimized.json
python3 test_controls.py
python3 verify_manifest.py
```

The finite checker uses only exact standard-library arithmetic and explicit
exception gates. The source pins record bounded provenance inspection rather
than claiming a repository-wide novelty search. The intended user referent,
squarehood transport and curve-level cycle operation remain unresolved. No
higher-extension search or broader conjectural bridge is part of this packet.
