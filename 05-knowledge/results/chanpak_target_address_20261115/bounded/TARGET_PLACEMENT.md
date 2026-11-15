# Target-first placement to O(N^(3/8)) with one apex

**Proved.** Every sufficiently large integer target N lies within
O(N^(3/8)) of the exact Hamilton count of an explicitly constructible
strong one-apex tournament on O(log N) vertices. The construction
uses integer square roots and elementary arithmetic; it does not
search through tournaments or enumerate host words. It is a placement
result, not an exact selector.

## 1. The inherited bounded-error product

Given a transitive spine with bits w=b0…b_(m−1), add an apex beating
the 1-sites. If the first bit is 1 and the last 0, the tournament is
strong: the apex reaches each spine vertex through the first one,
and every spine vertex returns via the last one. Let A be the binary
integer of w and B that of its reverse complement, in exactly m bits.

The [one-apex spectral law](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md)
gives H=AB/2^m+O(m). Explicitly, forward 1-before-0 pairs give the
inversion weights in AB/2^m; an adjacent 10 contributes one extra
unit and backward contributions sum to less than m/2. Including the
constant 1 gives absolute error at most m/2+1.

The new step here is to place this product **from a prescribed
target**. Two ordinary independent dyadic prefixes leave uncontrolled
reverse-complement tails of order sqrt(N); holding **both** leading halves on a common arithmetic grid keeps
the main piece of those tails **exactly fixed**. Treat the fixed
piece as part of the product baseline and then balance around a
near-square. Placing both addresses on an a priori grid instead
would create a larger, unnecessary midpoint error.

## 2. Pick the scale and align the tails

Choose m so that 2/5 < N/2^m ≤ 4/5 (for example the least integer m
with 5N≤4·2^m). Set m=2r+e with e in {0,1}, L=2^r, U=N/2^e. Then

~~~text
sqrt(2/5)L < sqrt(U) ≤ sqrt(4/5)L.
~~~

Set g=2^floor(r/4), so g is comparable to L^(1/4). Let a0 be the
least integer no smaller than sqrt(U). Require both leading
addresses a,b≡a0 mod g. Write p(a),p(b) for reverse complements
in exactly r bits. Their leading log_2(g) bits are identical,
fixed by a0 mod g. If v is that common leading integer, then

~~~text
p(a)/L = v/g + O(1/g),    p(b)/L = v/g + O(1/g).
~~~

There is no search for the residue. It is specified directly by
the target's rounded square root. In particular either tail can
be large; only its uncertainty is small. Put λ=v/(2^e g).
If e=0, take w to be a's bits followed by p(b)'s bits. If e=1 insert
a zero in the middle. Thus

~~~text
A = 2^e L a + p(b),
B = 2^e L b + e L + p(a).
~~~

On a,b≤L with the aligned tails, the exact product satisfies

~~~text
AB/2^m = 2^e (a+λ)(b+λ+e/2) + O(L/g).
~~~

This follows by factoring 2^e L from the exact A and B.
Both truncated tail fractions and the middle bit in B are
kept in the main product; replacing them by unstructured errors
would lose the improvement.

## 3. Choose a and b by descent from the target

Put μ=a0+λ+e/4. Since a0 was already the ceiling of sqrt(U),
we have 0≤μ−sqrt(U)<2, and hence

~~~text
0 ≤ μ²−U = O(L).
~~~

Choose a nonnegative integer t nearest to
(sqrt(μ²−U)+e/4)/g (rounding down is sufficient). Use

~~~text
a=a0+tg, b=a0−tg,
(a+λ)(b+λ+e/2)=μ²−(tg−e/4)².
~~~

Since tg=O(sqrt(L)+g)=o(L), both a,b stay strictly within
[L/2,L) for large N, using the fixed margins from the scale choice.
Their common residue retains the same λ. Rounding the square
root incurs an error of size at most g in that coordinate, so

~~~text
|(a+λ)(b+λ+e/2)−U| = O(g sqrt(L)+g²).
~~~

Thus the literal strong word constructed above obeys

~~~text
|H−N| = O(m + L/g + g sqrt(L)+g²)
      = O(L^(3/4))
      = O(N^(3/8)).
~~~

All decisions can use integer arithmetic (placement.py). For
example express μ=M/D with D=4·2^e g, then take
t=floor((isqrt(M²−D²U)+eD/4)/(gD)). The integer (D²U) is
built exactly even when e=1. A check on H uses the weighted
inversion recurrence; writing the full adjacency matrix is
only O((log N)^2) work. The initial a0 rounding and the
second square root are the only target-scale decisions.

## What this does and does not close

This is a direct baseline placement improvement: no probabilistic
equidistribution, assumed graph profiles, or matching of unrelated
residues. It leaves an absolute repair interval of order N^(3/8),
not a bounded interval, and changes the actual interface as a and
b move. A separate repair must work on those same manufactured
states; a generic unit move after an unrelated prefix is not enough.
The bounded-W theorem in this packet also explains why a fixed
finite source library with a tiny marked-closer correction cannot
on its own supply such a repair.
