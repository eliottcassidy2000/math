# Unit signed shield in a two-apex reverse bridge

The word "factor" here has two precise meanings. The two-apex deletion
is a transitive tournament (Hamiltonian count 1). Additionally, a
primitive signed linear form on the right module's two counts is exactly
1 at every height. We do not assert that the entire transfer operator
acts as the identity or that a one-apex factor-1 interface can hide a
gap at fixed slope.

## An exact physical repair of +2 at arbitrary heights

For a binary word w, put its letters on a forward transitive spine,
and orient apex→letter i iff w[i]=1. Write h(w) for the resulting
Hamiltonian count, A(w) for the word read normally in binary, A(ε)=0,
and p for its last bit, initially 0. The one-apex recurrence is

    1: (h,A,p) ← (h,2A+1,1)
    0: (h,A,p) ← (h+A+p,2A,0),   initial (1,0,0).

The strong words start with 1 and end with 0. Let rc reverse and
complement a word. Reflection of the spine after reversing all arcs
gives h(rc(w))=h(w).

Take the reverse bridge from the
[previous two-apex construction](../../chanpak_target_ribbon_audit_20261116/TWO_APEX.md):
join strong one-apex modules U(u), U(v) by directing all cross arcs from
U to V except that the mutual apex arc is oriented y→x, from the apex of
V to the apex of U. Let P=A(u)+1 and Q=A(rc(v))+1. Then

    H = h(u) h(v) + P Q − 1.                             (1)

For completeness, x-only and y-only odd cycles stay inside their
modules, contributing (including their disjoint pairs and the empty
packing) h(u)h(v) in the odd-cycle packing formula. A crossing cycle
uses y→x and is otherwise a rooted path starting at x, a forward cut
arc, and a rooted path ending at y. Their unsigned counts are P,Q;
each signed even-minus-odd count is −1, as its strong word ends in 0
(in the right module, pass to the converse). There are therefore
(PQ−1)/2 odd crossing cycles, of weight 2. These contain both apices and
cannot be in a packing with any other odd cycle. The modules are
individually strong, and the reverse bridge makes the whole tournament
strong. Deleting both apices leaves a transitive spine, while two
disjoint triangles certify that its feedback vertex number is exactly 2.

Now use the same order-ten spine for either left word

    u0 = 1101011000         (h0,P0) = (769, 857)
    u1 = 1111111110         (h1,P1) = (513,1023)

so (Δh,ΔP)=2(−128,83). These differ by four apex–spine arcs at
zero-based locations 2,4,7,8; every other arc stays fixed.
Let t>=1, D=2^t. Choose v=rc(w) with w from either of the following
families:

| w | h(w) | A(w)+1 |
|---|---:|---:|
| 11 0^t 1011010 | 249 D+59 | 384 D+91 |
| 101 0^t 1011010 | 415 D+59 | 640 D+91 |

Both satisfy the exact unit equation

    −128 h(v) + 83 Q = −128·59 + 83·91 = 1.

Consequently the whole physical edit changes H by exactly 2 at
every height in both families. Baselines are respectively

    H0 = 520569 D + 123357;
    H0 = 867615 D + 123357.

The first pairs are (1164495,1164497) and (1858587,1858589).
This is a different seam from the variable binary reader: the
increased path count nets a loss of 256 in h against a gain of 166
in the cut endpoint count. Four left-apex arc reversals are optimal
for this particular profile change: 166 is not a sum of at most
three distinct signed powers of two. A largest exponent at least 10
would give magnitude at least 256, and the remaining finite signed
sums can be checked directly.

## Inverse: exactly two strong right rays for this shield

The unit equation forces A(w)≡90 (mod 128); hence a word of
length at least seven ends in y=1011010. A shorter strong word
cannot encode 90. Write w=xy. Let C be the binary sum giving the
reverse address of the zeros in y, i.e. ∑ 2^i over zeros at positions
i from the left, starting at zero. Here C=82, and h(y)=61.
The recurrence gives

    A(xy) = 128 A(x)+90,
    h(xy) = h(x)+82A(x)+60,

because the first bit of y is 1. Substituting in the unit
equation gives exactly A(x)+1−h(x)=2; thus x is nonempty and
ends in 0 (parity), and starts in 1 (strongness of w).

Set δ=A+1−h. A 1 increases δ by A+1, and a 0 decreases it by p.
After an initial run of r ones and its first zero, δ=2^r−2.
Additional zeros leave δ unchanged. Any later run of s≥1 ones,
completed by a zero, increases δ by

    (2^s−1)(A+1)−1 ≥ A,

where A is the address just before that run. Therefore δ=2 either
comes from the initial run 11 alone, or starts at δ=0 with 10
and then takes precisely one more 1 immediately followed by a
zero. No further 1-run fits. These are exactly x=11 0^t and
x=101 0^t, t>=1. This also derives the tables and rules out hidden
dense strong guards within the very same primitive projector.

## The general dyadic primitive projector and its limitation

Take k≥1, an odd p≥1, and ask for strong right words w satisfying

    −2^k h(w)+p(A(w)+1)=1.                              (2)

Up to the finite words of length less than k, put Q0=2^k.
The congruence forces a unique suffix y by its numerical address
a = p^(-1)−1 mod Q0. Its last bit must be 0, as w is strong.
Put r=rev_k(a), C=Q0−1−r, and M=(p(a+1)−1)/Q0. If x is the
nonempty prefix, a first zero in y carries precisely the last bit
p_x of x into the h recurrence. Thus equation (2) is exactly

    h(x)+p_x [y begins in 0] = (p−C) A(x)+1+M−h(y).      (3)

Only the choice p−C=1 can have infinitely many strong guards:
* For p−C≥2, use h(x)≤A(x)+1 (follows from δ≥0 for all
  leading-one prefixes), which bounds A(x). For p−C≤0, note$
  that h(y)−(1+M)>−1 when p=C, hence h(y)≥1+M there. In fact
  h(y)=1+a−J(y) by the recurrence and, when p=C, 
  h(y)−M = E(y)+(a+r+2)/Q0 > 0, where
  E(y)=ar/Q0−J(y)≥0 is the [signed reversal defect](../../harmonic_reader_residual_20261116/RESIDUAL.md).
  This strict lower bound persists for p<C as M decreases.
  Thus 1+M−h(y)≤0 and (3) has no prefix solution.
* For p−C=1, p odd forces y to start in 1. The only allowed
  prefix defect is δ=h(y)−M=E(y)+(r+1)/Q0, a positive integer.
  Its parity must be even for a nonempty strong prefix ending in
  0. In fact evenness of δ forces the last bit of x to be 0,
  since δ≡A(x) mod 2. For any fixed δ, trimming all but the
  first trailing zero from x leaves finitely many shapes: an
  initial 1-run is bounded by 2^r≤δ+2; any subsequent 1-run
  starts at A≤δ and ends with its first zero at A'≤4δ+2 by
  the increment formula above. There are also finitely many
  prefixes ending in 1 for fixed δ: the last incomplete
  1-run increases δ by (2^s−1)(A+1), so if such a run follows
  a zero then both A and s are bounded, and initial runs
  alone have δ=2^s−1.

So the only infinite locus for each fixed (k,p) is a finite
union of trailing-zero rays before the known suffix (possibly
empty). Equivalently, an infinite projector must come from a
suffix obeying both p=Q0−r and (a+1)r≡−1 mod Q0, with the
required even prefix defect accessible. This classification
does not bound the union when the projector itself varies
with the target. It shows that one fixed signed shield cannot
by itself give universal absolute placement.

## The one-apex rigid factor versus the new seam

For a single distinguished root z, if g=H(T−z)=1 then T−z
must be transitive (the odd-cycle packing identity detects
any directed triangle). Its weighted rooted paths count
R+S=1+A, since every complementary induced subtournament is
transitive. Thus Z=R+S−g=A fixes the entire orientation word
at a fixed order. Leading zero padding adds only dominating
sources and does not alter h,R,S. Therefore the usual one-root
requirement of identical (g,Z) but distinct constant I cannot
yield a factor-1 +2 pair. The reverse bridge avoids this
rigidity by retaining two profile coordinates and matching
a primitive signed combination rather than freezing the
ordinary root mask.
