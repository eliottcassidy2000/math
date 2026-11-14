# The smallest independently movable three-point window has h=9

Here a full-profile window means several physical marked punctures
with exactly the same (h,f,P,Q,S,T), so arbitrary independent wings
add the same amount and only the ordinary augmented count K changes.
The [known smallest two-point switch](../chanpak_f1_physical_carry_20261114/MINIMAL_SWITCH.md)
is at h=5. It was not known from that result whether arbitrarily
long transitive padding around a count-five core could also give a triple.

**Theorem.** At h=5, the diameter of K is at most two even at
fixed *unsigned* coordinates (h,f,P,Q), without fixing S,T.
Thus no independent three-point window exists at h=5 at any
order. The minimum h for a full-profile three-point window is
exactly **9**, by the [literal h=9 triple](MOVABLE_TRIPLES.md).
The proof uses a finite carry automaton and an exact eight-coordinate
linear certificate, not a cutoff on the length of either wing.

## Strip the common parts before comparing any two views

The old puncture of count five has one nonsingleton strong
component, the unique count-five four-vertex core (8,1,3,6).
For completeness odd-cycle packing makes five equal to the empty
packing plus two weights of two. There can be only two directed
triangles and they cannot be disjoint. Every vertex in a strong
component lies on a triangle: if a vertex has both outneighbors
and inneighbors but lies on no triangle, every inneighbor
beats every outneighbor, preventing strongness. Thus the
nontrivial component has at most five vertices, and a strong
five-vertex tournament has at least three triangles. For the
last assertion its outdegrees lie in {1,2,3}, sum to ten, and
have equally many 1s and 3s, so
10−Σ binom(outdegree,2)≥3. The factorization of the ordinary
count across strong components excludes additional nontrivial
components at h=5. All other components are singletons.

Thus f=32·2^D, where D is the total number of singleton
vertices. Compare two punctures with equal P,Q at this f.
Suppose the first has its four-core to the left of where the
second puts its core. The shared initial singleton word is
fixed by Q modulo 2^a, where a is its length: the contribution
is 5 times its reverse-complement address, and 5 is a unit.
The shared terminal word is determined by P modulo 2^b in the
same way. Remove these identical two outer words. Equal P,Q for
the full graph now mean equal P,Q for the middle (the ordinal
transfer multiplies their differences by nonzero factors).
If the positions coincide, the sixteen intrinsic (P,Q) are
distinct, so the same mask was used and K agrees.

If the positions differ by d>0, the two remaining views are
C→Z and X→C', where C,C' are the four-core with different
allowable apex masks and X,Z are length-d transitive words.
Here X,Z denote the ordinary addresses, not reverse-complements.
Let p,q and p',q' be the intrinsic (P,Q) of C and C'. With
L=2^d and reversal declared at width d, the equalities are

    32 X − 5 Z = p L−p',
    32 rev(Z) − 5 rev(X) = (27−q')L+q−27.

These are exact, so no residue matching is substituted for
the actual interface.

## Carry automaton and certified exact count

Read digits x_i,z_i of X,Z from least to most significant.
Use a state (c,e) starting at (p',27−q') and ending at (p,27−q).
A digit pair acts by

    c'=(c+32 x_i−5 z_i)/2,     e'=2e−32 z_i+5 x_i.

The first carry must divide integrally. Both carries are in
[−5,32]: the ordinary division process preserves this interval,
and the e-process is the reverse of the ordinary division
process for the reversed equality. The exact endpoint
conditions give the two displayed integer equations, so this
finite graph accepts precisely the genuine shifts.

To check K on those shifts, track first and last digit bits in
the discrete state, and the following vector with linear
integer updates: the constant 1, current power L, both
addresses X,Z, both reversed addresses, and both ordinary
counts g(X),g(Z). When prepending bit x to a read word, its
update for g is

    g(new) = g(old)+x(L−rev(old)−first(old)).

Set first(empty)=1 for this update. The formula follows from
the weighted inversion count; the other coordinates update
by ordinary binary appending and reversal.

For every productive mask pair and augmented discrete state,
retain an exact rational basis of the span of the reachable
vectors there, adding a vector only if it is independent.
Propagate retained reachable vectors by each legal transition.
There are at most eight independent vectors in each state.
Linearity proves that the final span includes **every** reachable
vector, by induction on path length. This finite proof therefore
certifies all word lengths, even on productive cycles.

At an accepting state the signed count change before restoring
the outer wings is the linear expression

    K'−K + 5(g(X)−g(Z)) + X q' − p (L−1−rev(Z))
         + p(X) T' − S(1−first(Z)).

Restoring common outer wings adds (T_middle'−T_middle)
times the left address parity and (S_middle'−S_middle)
times the right address parity. These are binary choices
(even including an empty wing, with parity zero):
T_middle'=5(1−first(X)), T_middle=T;
S_middle'=S', S_middle=5 p(Z).
The script checks the exact output for all four exterior choices.
In each accepting discrete state the functional is a fixed
constant times the 1-coordinate on its entire span, and each
constant is −2, 0 or 2. This certificate has eight positive-width
mask pairs and 108 productive boundary states. The older h=5
unit switch attains two, so the upper bound is sharp.

There is no h<5 switch even unsigned by the
[triangle support classification](../chanpak_f1_physical_carry_20261114/SUPPORT_DECODER.md);
h=7 is excluded by the scalar 7 theorem. This proves the
claimed minimum at nine. It does not constrain windows that
change their exposed responses, such as the five-point
tilted families which balance a linear combination instead.

Run **node replay_h5_shifts.js** for the exact span certificate.
No floating-point rank computation or finite length bound is used.
