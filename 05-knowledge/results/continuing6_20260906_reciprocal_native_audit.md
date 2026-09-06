# Independent audit of the native reciprocal packet theorem

**Verdict: PASS for the analytic iff, exact counting, finite controls and
three conditional-on-primality infinite constructions. No repair requested.**
The infinitude of prime terms remains CITED Dirichlet. The result concerns
only the exchange of outputs0 and1 in the p-point completed reciprocal
graph. It does not claim a successful infinite family, a classification of
arbitrary output swaps, or a 2p-point construction.

This audit read the complete producer report and source, but imports neither
the producer nor a repository implementation. The independent exact engine
uses the Python standard library only.

## 1. Primitive native packets: both directions and multiplicity

The repaired graph is invariant under transposition. Three untouched points
cannot be collinear because a line meets the modular conic xy=1 in at most
two points. A triple containing both moved points P=(0,1),Q=(1,0) would
require an untouched point with x+y=1, which is impossible in the standard
integer box. It follows that every triple has exactly one moved point,
and transposition pairs the P and Q families disjointly.

Every line from P through an untouched point has a unique positive primitive
slope a/b, with equation b(y-1)=ax. Its two untouched points are therefore
(bk,ak+1),(bl,al+1), where positive integers k<l are unique. Modulo p their
distinct y coordinates solve b y^2-b y-a=0 and sum to1. The only possible
standard sum is p+1, so a(k+l)=p-1. The bounds used in this step are native
integer bounds, not a choice of modular representatives.

Writing h=k+l, the reciprocal identity gives bkl=h+cp. The estimate
bkl-h> -p excludes negative c. If c=0 then b=1/k+1/l is not a positive
integer, both when k=1 and when k>=2. Thus c>=1. Finally bl<=p-1 gives
h+cp<=k(p-1), so c<k. Hence k>=2, l>=3, h>=5. This proves all of the
producer's necessary packet conditions, including the short-cofactor bound.

Conversely, those conditions place both coordinates in the standard box.
For example al+1=p-ak<=p-2 and bl<=p-1. The equation

    bk(ak+1)=-abkl=-ah=1 modulo p

holds because ah=p-1 and bkl=h+cp; the same argument applies at l. Both
points are therefore untouched reciprocal points, and their integer line
through P is the specified one. Primitive slope and k<l reconstruct the
packet uniquely. The conic permits at most two untouched intersections on
that line, so one packet corresponds to one triple, not a larger collection
of choices from a long line.

This proves the bijection and X=2N in the claimed sense. The bound
b<2a follows from bh=bk+bl<2(p-1)=2ah. The p61 packet attains h=5 and
therefore verifies sharpness of that particular necessary bound.

## 2. Quadratic splitting is insufficient without standard lifts

In the allowed slope range a,b are nonzero modulo p. The roots of
ab k^2+bk-1 have sum h modulo p and are nonzero. Their discriminant is
Delta=b(b+4a). A nonzero square gives two distinct roots; zero discriminant
must not be counted as a pair.

Set H=min(h-1,floor((p-1)/b)). If both standard roots lie in [1,H], their
positive sum is less than2h, so the sum congruent to h modulo p is exactly
h. They supply the native coordinates and all packet conditions above.
Conversely every native packet puts both standard roots in that interval.

The one-root formulation also checks out. If a root k lies in [h-H,H],
then h-k lies in the same positive interval, below p, and is the other
standard root. This reasoning still requires nonzero discriminant. It is
not valid for an arbitrary modular representative, or after multiplying a
root by b without reducing and rechecking its coordinates.

At p13,a2,b1 the roots7,12 of the split quadratic both miss H5. Their
associated modular collinearity has actual determinant65, while the
repaired graph is successful. This is a decisive native-lift hostile to
a character-only criterion. At p5 the character-five polynomial has a
single repeated root, correctly giving no triple.

## 3. Character-five obstruction and the three progressions

For p!=5 with 5 a square, the two roots of x^2+x-1 are distinct and avoid
0,1,p-1. Their standard sum is p-1, and the reciprocal points (x,x+1)
lie on the actual integer line y=x+1 through P. Thus the character condition
is sufficient for failure. The converse is refuted at p37 by the stated
primitive packet, and the complete smaller-prime bank establishes the
bounded minimality claim.

An independent algebraic route confirms the special reciprocity statement
used for the progressions. For an odd p!=5, let zeta be a primitive fifth
root of unity over the algebraic closure of F_p. Then t=zeta+zeta^-1 solves
t^2+t-1=0. Frobenius fixes t exactly when zeta^p is zeta or zeta^-1:
these are the two roots of X^2-tX+1. Therefore 5 is a square in F_p exactly
when p is1 or4 modulo5. All three listed progressions are instead2 or3
modulo5, so their prime terms have character -1. No empirical character
pattern is being extrapolated.

The template proof does not need primitive a,b for existence. Its affine
coordinates satisfy exact degree-two reciprocal-product and line identities
for every integer parameter r>=0, with coefficientwise native bounds. I
independently reconstructed all coordinates and quotients. The two product
quotient pairs for the three progressions are respectively

    (11+96r, 17+150r),
    (3+4r, 18+25r),
    (8+9r, 22+25r).

The exact identities were verified at three parameter values after proving
the degree bound two. This is polynomial identity verification, not sampling
of the infinite prime family. Linear coordinate bounds have nonnegative
constant and slope coefficients and therefore hold for all r>=0.

All three residue/modulus pairs are coprime. I read the opening theorem in
the primary English translation of Dirichlet's1837 paper directly; it
supplies infinitely many prime terms for each of these progressions.
[Dirichlet's original arithmetic-progression theorem, translated by Ralf
Stephan](https://arxiv.org/pdf/0808.1408v2). This is the sole external
infinitude input. It does not replace the explicit integer-coordinate proof.

Primitive normalization is needed for counting. If the displayed a,b have
gcd d, divide both by d; h,k,l,c then all multiply by d while the geometric
points remain unchanged. At p113 this produces exactly (4,3,28,8,20,4).
The same graph has four triples, so the two triples forced by a template
must not be advertised as its total defect count.

## 4. Independent exact evidence

The independent source covers the exact declared17 primes: all odd primes
through43, then61,97,113,197. It enumerates all **1715293** unordered triples
of their literal repaired integer graphs. From each triple through P it
extracts the primitive slope, spacings and carry, recovering every saved
packet without using the producer's anchored-line grouping.

A separate complete scan covers all **303** eligible primitive slopes. It
finds modular roots by directly evaluating the quadratic at every residue,
without the producer's discriminant square-root formula or h,k,c packet
loops. Every bounded-root decision, the equivalent one-root interval,
every primitive packet and the doubled native triple counts agree.
There are17 accepted packets across these controls.

The referee checks the entire successful-control list, the first negative-
character failure, the h5 sharpness example, the p13 lift hostile, the p5
ramification case, the p113 primitive-gauge change and additional triples,
and all coefficients or degree-bounded evaluations of the three templates.
No random or floating computations occur.

## 5. Frozen reproduction and pins

    python continuing6_20260906_reciprocal_native_audit.py
    python -O continuing6_20260906_reciprocal_native_audit.py

The optional `--certificate PATH` overrides the input. Defaults support
adjacent outside artifacts and a filed 04-computation source with results
in 05-knowledge/results. Always-active `require` gates are used throughout.
Stdout is explicitly LF; normal and optimized stdout were captured as raw
bytes, are identical and contain no carriage returns.

**998 exact gates PASS in both modes.**

Frozen producer pins checked:

- Source `d0981913961389520200a8215b1c2f8611ea7090bd379ba2024ceb9a6b345c11`.
- Certificate `1da5e84836de4aee72e0636e792a43715c9cb85eac9e621203e215ac7e228380`.

Independent referee pins:

- Source `d1ab92b90becba019037f16634be910dbbaf1b6a47b67e7a9d821093f9b941de`.
- Normal/optimized output `e24e96ef9757fba9687acf5e4aec16f478f2eefaa49bdba4c65f77fb75711456`.

No producer, repository file, maintained navigation or Git state was changed
by this audit. The stated analytic and finite scopes are accepted without
an outstanding repair.
