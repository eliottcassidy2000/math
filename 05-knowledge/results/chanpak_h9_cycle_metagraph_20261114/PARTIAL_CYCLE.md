# Partial cycle-marked tournament classes: rotation and converse Burnside

Let (n=k+m), (kge 3). We count pairs ((T,C)) up to
isomorphism, where (C) is a distinguished directed simple
(k)-cycle, unrooted. There is no induced restriction apart
from the compulsory cyclic edges; chords remain free. We also
count pairs fixed by taking converse and reversing the cycle.
This extends the full-cycle case in [PROOF](MIXED_MOMENTS.md) and solves
the nonrigid cycle-scaffold variant named by the
[marked-path packet](../marked_path_tournament_series_20261021/session23-symmetry/README.md).
Converse-merged classes then number ((A+B)/2).

For an integer partition (lambda) put
(z_lambda=prod_ell ell^{m_ell}m_ell!) with (m_ell) its
multiplicities. For odd partitions define the familiar tournament
exponent for pair orbits
$$
 E(lambda)=sum_i(lambda_i-1)/2+
               sum_{i<j}gcd(lambda_i,lambda_j).
$$
The undecorated order-(m) tournament enumeration uses
(sum_{lambdadash m, mathrm{odd}}2^{E(lambda)}/z_lambda)
(the empty sum convention is one at (m=0), from the empty partition).

## Ordinary cycles

The exact cycle-marked count is
$$
 A(n,k)=rac1k
 sum_{substack{dmid k\d	ext{ odd}}}
 arphi(d),2^{k(k-3)/(2d)}
 sum_{substack{lambdadash m\lambda	ext{ odd}}}
   rac{2^{E(lambda)+(k/d)sum_jgcd(d,lambda_j)}}{z_lambda}.
$$
Indeed, fix the cycle on the first (k) vertices. Its orientation
forces the group (C_k	imes S_m) when relabelling. A cycle
rotation of order (d) permutes the scaffold in (k/d) cycles
of size (d); even (d) or any even outside permutation cycle
is impossible in a tournament automorphism, by exchanging a pair.
There are (k(k-3)/(2d)) free internal chord orbits, and the
edges in the unmarked set have (E(lambda)) free orbits. Between
vertex cycles of sizes (d) and (lambda_j) there are exactly
(gcd(d,lambda_j)) orbits of pairs. All of those are free.
Independence of orbit choices and counting each permutation type
by (m!/z_lambda) gives the result. The (d=1) inside term is
exactly the marked-path packet's coloured Burnside channel;
the (d>1) terms retain what path rigidity suppresses.

## Self-converse cycles

Write (f=1) if (k) is odd and (f=0) otherwise, and
(u=(k-f)/2). Let
$$
 I(k)=
 egin{cases}
 (k-3)(k+1)/4&k	ext{ odd},\
 k^2/4-k/2-1&k	ext{ even}.
 end{cases}
$$
Let (Lambda) contain the partitions of (m) whose parts are
all (2mod 4), except for at most (1-f) singleton parts.
For these write
$$
 E^-(lambda)=sum_{lambda_i>1}lambda_i/2
                   +sum_{i<j}gcd(lambda_i,lambda_j).
$$
Then the number of cycle-marked classes fixed by converse is
$$
 B(n,k)=
  2^{f-1},2^{I(k)}
  sum_{lambdainLambda}
    rac{2^{E^-(lambda)+f,|lambda|
          +usum_jgcd(2,lambda_j)}}{z_lambda},
$$
where (|lambda|) in this last exponent counts parts, not their
sum. In particular there are no such classes when (k) is odd
and (m) is odd.

For proof, use the twisted coset version of Burnside with an
orientation reversal and a reflection of the scaffold. Each
allowed scaffold reflection fixes (f) vertices and exchanges
the other (u) pairs; all (k) reflections work internally for
odd (k), but only the (k/2) edge-axis reflections for even (k).
Their internal free exponent is (I(k)), as proved in the
full-cycle packet. For the outside anti-permutation, an odd
cycle of size greater than one has a pair orbit with odd length
and unswapped final endpoints, impossible under edge negation.
For an even cycle of length (ell), the antipodal pairs have
orbit length (ell/2) and swapped final endpoints, so
(ell/2) must be odd. Two fixed vertices anywhere also obstruct
the sign. These checks are sufficient: all surviving non-antipodal
internal orbits and all cross-cycle orbits have even length
and no endpoint swap, except for the single allowable fixed
vertex which contributes no pair with itself. Thus each even
outside cycle supplies (ell/2) free internal edge orbits,
and exterior-cross and scaffold-cross free orbit counts are
exactly the gcd sums displayed. Twisted Burnside divides by
(k m!), giving the prefactor (2^{f-1}).

The first-moment count here is of *cycle-marked isomorphism
classes*, not simply the total raw cycle counts summed over
unmarked classes. An unmarked tournament contributes one such
class per orbit of its directed cycle set under its automorphisms.
The full-cycle rotor formula specializes at (m=0).
