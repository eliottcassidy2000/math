# Reflect the source, then solve the exact target line

Let h_A=H(A), f_A=F₂(A), where F₂ counts ordered two-slot path
covers (empty slots allowed). For nonempty tournaments A,B make the
literal cyclic join A→B→c→A with one new vertex c. This is strong:
inter-block arcs alone give a return through all three parts.

## Two physical readers, one signed reflection

The ordinary Hamiltonian count N of the cyclic join is exactly

    N = (f_A−h_A)(f_B−h_B) + 2h_A h_B.
    3N − 2 f_A f_B = (f_A−3h_A)(f_B−3h_B).

For the first equation use the [odd-cycle packing expansion](https://arxiv.org/src/2307.05569v1/anc/redeiberge-long.pdf#page=14).
Packings not using c contribute h_A h_B. Any other packing has
exactly one cycle through c. Its other vertices are a nonempty path
in A followed by a nonempty path in B. The cycle is odd iff those
two subsets have the same parity. The factor 2 per cycle cancels
the half in the parity projector. The unsigned choices are
(f_A−h_A)(f_B−h_B); the signed choices are (−h_A)(−h_B).
The latter identity follows by coloring each odd cycle or isolated
vertex into either slot: the signed weights cancel on any nonempty
tournament, and removing the empty subset subtracts h_A or h_B.

For any positive integer q the same proof gives, writing F=F_q,

    F_q(join) = q [F(A)F(B)
          + (F_{q+1}(A)−F(A))(F_{q+1}(B)−F(B))
          + (F_{q−1}(A)−F(A))(F_{q−1}(B)−F(B))].

There are q choices of slot for the marked cycle. A selected path
and a complementary q-slot cover give F_{q+1}; signing the selected
slot gives F_{q−1} since its minus sign cancels one of the q other
choices per odd packing atom. Removing the empty selected path
subtracts F_q. This formula includes q=1 with F₀=0 on nonempty sets.

In particular a transitive reader B of order k≥1 gives

    N = (2^k−1)f_A +(3−2^k)h_A
      = 2h_A+(2^k−1)(f_A−h_A).

For k=1 and k=2 this gives f_A+h_A and 3f_A−h_A,
a literal reflection about 2f_A with the source unchanged.
There is no subtracting one tournament's total from another as
an operation; both outputs are separate strong tournaments.

Example: a four-vertex count-five source, followed by six
transitive singleton vertices inside A, has (h_A,f_A)=(5,2048).
The two readers give exactly the prime **2053** and **6139=7·877**
with prime cofactor. The [replay](replay_absolute_reflection.js)
counts both full tournaments. The following actual cores permit
more changing-source choices at the same dyadic baseline.

| h | f | Outgoing adjacency masks |
|---:|---:|---|
| 5 | 32 | (12,1,2,6) |
| 19 | 256 | (44,49,2,38,45,4) |
| 217 | 4096 | (236,113,18,70,73,220,4,94) |
| 223 | 4096 | (236,241,82,6,201,220,8,76) |
| 599 | 8192 | (44,241,50,38,73,208,141,29) |

Both h and f multiply across ordinal sums. For example for any
L≥5a+8b+12c+12d+13e, combine the listed cores with multiplicities
a,b,c,d,e and use singleton padding to reach f_A=2^L.
Then both 2^L+h_A and 3·2^L−h_A are attained exactly, with
h_A=5^a19^b217^c223^d599^e. The upper and lower signs can thus
respond to the same physical core change in opposite directions.
This is an infinite explicit palette, not a complete fiber.

## Exact inverse when both color counts are dyadic

For prescribed odd N>1, set L=ceil(log₂ N). If f_A and f_B
are powers of two, their product **must** be 2^L.
Indeed f≥2h for every nonempty block (the two empty-slot choices).
With α=h_A/f_A and β=h_B/f_B, the relative count is
1−α−β+3αβ. On 0<α,β≤1/2 this is strictly between 1/2 and 1
(check the corners, noting α,β>0).

Hence the concrete inverse is: factor 3N−2^(L+1), try both signs
of its divisors, and try all splits s+t=L with s,t≥1. For each
factor pair u,w solve

    h_A=(2^s−u)/3,    h_B=(2^t−w)/3.

Reject unless the counts are positive integers and belong to the
actual fibers (h_A,2^s), (h_B,2^t). Divisors include all possibilities
because the two factors are odd and nonzero. This is a necessary
and sufficient absolute test with an exact construction upon
acceptance: f≥2^|A|, so both fiber checks are finite. Under a
transitive reader of order k, instead try the unique
h_A=(2^L−N−2^(L−k))/(2^k−3) for each k=1,…,L−1.
For general non-dyadic inputs one can vary both f as well:
f_A,f_B≤2(N−2) follows from f−h≥f/2, f−h≥1 and h≥1.
A prime target forces gcd(h,f)=1 for each input, since any shared
divisor divides N but is at most h≤N/3.

A further exact pruning rule uses the same cut inequality proved
below: f≥(n+1)h. The factor f−3h is negative only for a singleton,
where it is −1. The factor +1 occurs only for the two-vertex
transitive block: for n≥3, f−3h≥h, and equality at +1 would
force h=1 and f=4, contrary to f≥2^n. Therefore a positive
prime 3N−2^(L+1) forces one input to be the two-vertex
transitive reader even in the arbitrary two-input dyadic class.
A negative transformed target forces a singleton reader.

## A sharp upper window, even with all dyadic sources and readers

At fixed s+t=L, relax the possible h to 1≤h_A≤2^(s−1) and
1≤h_B≤2^(t−1). Since N is bilinear, its maximum is at a corner.
The corner h_A=h_B=1 gives 2^L−2^s−2^t+3. The other three
are at most 3·2^(L−2), because one ratio is 1/2 and the other
is at most 1/2. The first expression is maximized by balancing
s,t. For L≥6 this balanced value exceeds 3·2^(L−2)
(the two powers sum to at most a quarter of 2^L for L≥6,
with equality possible at L=6); for L=3,4,5 it is the
largest odd value not exceeding that bound. At L=2 it is 3.
Both h=1 corners are transitive, so the exact maximum is

    2^L − (2^floor(L/2) + 2^ceil(L/2) − 3).

This leaves a window of order sqrt(N) before every dyadic upper
boundary. Enlarging the dyadic source bank alone cannot close it.
The next positive target is to move the product f_A f_B itself
while supplying both actual reflected factors; the factorization
unambiguously tests success, but does not guarantee a choice for
every admissible odd N. Chan–Pak is not proved here.

## Zero density even with every dyadic fiber available

If T has n vertices, cutting each full path at any of the n+1
positions injects (full path, cut) into the ordered two-slot
covers. In the other direction each ordered cover injects into
a permutation plus a cut (concatenate the two paths). Thus

    (n+1) H(T) ≤ F₂(T) ≤ (n+1)n!.

Consequently for F₂=2^s, the ratio H/F₂ is O(log(s+2)/(s+1)).
Indeed the factorial bound forces n+1 to be at least a constant
times (s+1)/log(s+2). This is an absolute constraint, not a
statement that the ratio is typically small.

At height 2^(L−1)<N<2^L=M, put J=floor((log₂ L)/2).
If both exponents s,t exceed J, the two ratios are
O(log(J+2)/(J+1)), and therefore N lies within
O(M log(J+2)/(J+1)) of M. If one exponent s is at most J,
the other is at least L−J. For each such small-side ratio α,
the exact formula N/M=1−α+(3α−1)β puts N within
O(M log L/L) of M(1−α). Across the small exponents there
are at most O(2^J)=O(sqrt L) possible ratios, even allowing
every integral h. Thus the second case covers at most
O(M log L/sqrt L) possible integers. The total per octave is
o(M); summing the octaves gives o(X) outputs through X.
So no bank or refinement restricted to *dyadic* color counts
can cover the scalar spectrum by this cyclic join, even when
both of its input tournaments vary over all orders.

## Move the product to reach targets outside the dyadic window

For two inputs of order at least two, both f's are divisible by
four: in the [odd-cycle expansion](https://arxiv.org/src/2307.05569v1/anc/redeiberge-long.pdf#page=14)
at two colors, each term has at least two factors of two
(either two free vertices or the color and weight of an odd cycle).
Both shifted factors u=f−3h are positive. Hence a nondyadic
factor-first decoder can try multiples M of 16 with N<M<3N/2,
factor 3N−2M=uw, and choose the two physical sources with
those shifted factors and product f_A f_B exactly M. The
one-singleton branch is the separate test f+h=N.

The shifted fibers are finite independently of M. For n≥3,
the cut bound gives u≥(n−2)f/(n+1)≥2^(n−2), so an input
for u has order at most 2+floor(log₂ u). For instance +1
belongs only to TT₂; +3 only to the directed triangle;
+5 only to TT₃. The next small candidates +7,+9,+11 do not
occur: order at least five gives u≥16, and orders at most
four have u from −1,1,3,5,13,15,17. This replaces arbitrary
source choice at a specified factor by a finite literal search.

Two exact window repairs are certified in the replay:

* The prime **8147** exceeds the entire dyadic-class cutoff
  8003 in its octave. Use (h_A,f_A)=(5,32) and (h_B,f_B)=(29,320)
  with B outgoing masks (12,17,2,38,13,23). Then M=10240,
  3N−2M=3961=17·233, and both factors are physical.
* The target **3983=7·569** exceeds its cutoff 3971. Use TT₂
  and (h_B,f_B)=(109,1364) with outgoing masks
  (124,9,42,80,6,90,22). Here M=5456 and the physical factors
  are 1 and 1037.

These are prescribed absolute outputs beyond the dyadic
restriction, not residue hits. An all-target proof still
needs to ensure one successful M and factor pair, or another
family when the cyclic join fails.

The singleton-cycle readout was derived in parallel with the
bilateral reader work. The two-input exact divisor test and
sharp upper window follow by allowing both color fibers to move.
Run **node replay_absolute_reflection.js** for literal color counts,
the core palette and the prime/7p reflection example.

## A second physical route: change the response without changing H

The [complementary-reader construction](COMPLEMENTARY_READER.md)
solves a different exact equality. For an explicit unbounded
family it produces two strong tournaments of different orders
with exactly the same H, one with R₊ asymptotic to H and the
other with R₊ asymptotic to 2H. At H=831 the responses are 881
and 1567. This is a concrete way to retune the reader without
an ordinary-count error, not a realization of H=881 by that move.
The construction comes with an exact inverse for its own family
and a separate literal replay.
