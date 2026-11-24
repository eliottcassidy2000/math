# Determinant-seven Farey transfer for the H=23 two-slot kernel

The H=23 two-module calculation gives the ordered nonempty cover
kernel K=[[4,11],[9,23]]. Fix the second marked module B (vertex 4 in the stated H=23 labeling)
with h Hamilton paths and c ordered nonempty two-path covers. The two
coefficients supplied to the first marked module A are M(h,c)^t,
where M=[[11,4],[23,9]] (take A's two-cover coordinate first). In particular det M=7. The following is
a determinant-seven extension of the rational-address perspective in
THM-4057 and THM-4059, not a strengthening of their separate depth
packet character estimates.

## Exact local arithmetic and the positive cone

The corresponding increasing projective map is
F(t)=(11t+4)/(23t+9). It maps [0,∞] continuously and bijectively onto
[4/9,11/23], and is bijective on the rational points of these intervals.
Indeed the inverse on integer column vectors is projectively given by
J=adj M=[[9,-4],[-23,11]], whose two coordinates are nonnegative exactly
on the image cone.

If v=(x,y) is primitive, the content g(v) of Mv divides 7 since
J M v=7v. It is exactly 7 if x+y≡0 mod 7 and 1 otherwise.
The kernel line has slope −1, while the image lattice consists exactly
of the pairs (X,Y) with Y≡4X mod 7 (as follows from J or index 7).
Thus for a primitive output vector w in the image cone, the unique
primitive antecedent is Jw/gcd(Jw); its unnormalized image is w on the
image line mod 7, and 7w otherwise. This distinction cannot be removed
by left or right unimodular integer changes of basis.

For any primitive Farey neighbors v,v' (absolute determinant 1), their
primitive output distance in the Farey determinant is
7/(g(v)g(v')). The two input vectors cannot both lie on the kernel line
mod 7. Consequently the output remains a Farey edge precisely when one
endpoint lies on that line; otherwise the transformed determinant is 7.
The analogous inverse rule exchanges the kernel and image lines. No
global neighbor-preserving SL₂ realization of this projective map
exists: two invertible real matrices of the same projective action are
scalar multiples, and a rational scalar cannot change determinant 7 to 1.

## Positive-address transducer

Use L=[[1,1],[0,1]] and R=[[1,0],[1,1]], with addresses evaluated as
products in the written order acting on (1,1). They enumerate all
positive primitive columns by the subtractive Euclidean algorithm.
For a nonnegative determinant-seven matrix, repeatedly subtract the
second row from the first if the first dominates componentwise
(emitting L), or the first from the second under reverse dominance
(emitting R). Each subtraction strictly reduces the sum of entries,
and determinant rules out ambiguity. The final reduced matrix is of
the form

    A=[[c+r,b],[c,b+s]],  b,c>=0, r,s>=1,
    c s + b r + r s = 7.

Indeed the rows are incomparable, and the opposite diagonal ordering
would give a negative determinant. Thus this is a finite list of
nonnegative integer matrices.

Start by reducing M, which emits RRLLLL and leaves ((7,0),(1,1)).
For each input letter, multiply the retained matrix by that letter
on the right, reduce as above and append the emitted address. At the
end, apply the residual matrix to (1,1), divide its content, take its
Euclidean address, and append that address. Every step is an exact
matrix identity, proving that the final address is the primitive
projective image. The companion exact enumeration runs over all
b,c in 0..7 and r,s in 1..7 in the displayed equation: all 25 obtained
states are reachable, and each step emits at most seven letters,
with equality at some state. In particular output address length is
at most seven times input length plus an absolute constant. This is
a positive integral word transfer, not merely a real change of chart.

## Why the native state does not fill the real interval

For a tournament A of order n≥1, each Hamilton path gives n−1 distinct
ordered nonempty two-path covers by cutting at one edge. The map from
(path,cut) is injective, because the ordered cover recovers the path
by concatenation and the cut by its first path length. Hence
c/H≥n−1. The normalized physical reading
ρ=(11H+4c)/(23H+9c) therefore satisfies

    0 < ρ−4/9 = 7/[9(23+9c/H)] ≤ 7/[9(9n+14)].

Every closed subinterval away from 4/9 contains only finitely many
native ratios: the inequality bounds n and there are finitely many
tournaments of bounded order. The endpoint 4/9 is not attained, but
is approached by transitive modules of increasing order (their H=1
and c=2^n−2). Accordingly the closure of the physically attainable
ratios is exactly their discrete set together with the unique
accumulation point 4/9. The projective real interval is full only
after replacing actual module states by arbitrary nonnegative input
rays. In particular continued-fraction density in the completed
chart is not an existence theorem for tournament modules.

## Two physical probes recover the pair, including exact digit loss

Use the probing box A first as TT₁ and then as TT₂, whose (h,c)
coordinates are (1,0) and (1,2). If their height readings are x,y,
then Y=x and X=(y−x)/2 in M(h,c)=(X,Y), so

    h=(9y−17x)/14,   c=(45x−23y)/14.

Thus this pair of physical probes is complete: equality on them is
equivalent to equality of (h,c), and so to equality in every other
probing box. At any odd modulus the two probes detect precisely the
congruences detectable by arbitrary probing boxes. At a power 7^j,
the ambiguity in the formal underlying pair is exactly

    Δh≡Δc≡0 mod 7^(j−1),   Δh+Δc≡0 mod 7^j.

Indeed congruence of the two readings means MΔ=0 modulo 7^j.
Multiplying by the adjugate forces the first two divisibilities;
after dividing them out, the kernel line yields the last condition.
The converse is immediate. This is a cyclic order-seven ambiguity
of the formal residue states at every depth, not a claim that every
representative of each pair residue is achieved by a native module.

Replay: run farey-seven-check.py with Python 3.

## A sharp failure of native iteration

The real lift \(M/\sqrt7\) belongs to \(\mathrm{SL}_2(\mathbb R)\). It is hyperbolic (trace \(20/\sqrt7\)), with positive fixed slope \((1+\sqrt{93})/23\) and hyperbolic translation length \(2\,\operatorname{arcosh}(10/\sqrt7)\). These are statements about the real completed chart. There is an even stronger boundary than a mere lack of density: applying \(F\) to any *physical output slope* never gives another physical output slope. In fact
\[
 F([4/9,11/23])=[80/173,213/460]\subset(23/50,19/41).
\]
Modules of order at least four give slopes at most \(23/50\) by the cut injection. The singleton gives \(11/23\); the two-vertex module and three-cycle give \(19/41\); and the transitive triple gives \(35/77<23/50\). These exhaust smaller orders. Thus there is no second native step on the output state set, even though signed real iteration and abstract Farey transfer are well-defined.
