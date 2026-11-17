# Prime rotation, orbit phases, and what the diagonal forgets

This is a different use of the cyclic group ring from the Hadamard
orthogonality equation. It does not impose autocorrelation flatness.

Let a tournament have a fixed automorphism g of odd prime order
ell. Partition its vertices into a fixed set F and m moved
ell-orbits. Label each moved orbit by F_ell so that g adds 1;
for distinct orbit indices i,j, write D_ij for the possible
increment b−a of an arc (i,a)→(j,b). Then
D_ji is the complement of −D_ij. Write d_ij=|D_ij| mod ell.

Reduce G mod ell. Every nonfixed odd-cycle pack cancels in an
ell-orbit. A fixed pack splits into a pack on F and atoms on
moved orbits: a fixed cycle cannot mix a fixed vertex with a
moved one, since a nontrivial rotation of a directed cycle
fixes no vertex. A nonfixed cycle mixing them would intersect
a translate at that fixed vertex, so cannot belong to a fixed
pack either.

On the moved orbits there are three atom types:

* An unused orbit gives y^ell.
* An internal invariant cycle of one orbit is a constant
  step in one of (ell−1)/2 allowed directions. It gives
  x. Together with the unused choice this weighs
  U=y^ell−x/2 modulo ell.
* An odd simple orbit cycle of length r≥3 chooses one
  increment on each edge. If their total is zero, the lift
  consists of ell disjoint odd cycles and has weight x^ell;
  if nonzero, it is one invariant cycle and has weight x.
  These are all possibilities: a g-invariant cycle is a
  rotation with distinct orbit positions in its quotient;
  a free orbit of disjoint cycles can use each moved orbit
  at most once per cycle, otherwise its translates overlap.
  Since ell is odd the projected length must itself be odd.

Fix one orientation C of a simple unoriented odd cycle.
The degree of closure is the constant coefficient of
the product of the corresponding indicator words in
F_ell[C_ell]. Denote the all-ones group word by J.
Here J²=0 and fJ=f(1)J. The reverse indicators are J−f*;
expanding the reverse product and adding (the cycle length
is odd) cancels its zero-J term against the original when
taking constant coefficients. Exactly one-J terms survive.
The total zero-phase count for both orientations is therefore

    e_(r−1)(d_e : e in C) = sum over e of product over f≠e d_f,

while the total count of all increment choices is zero mod ell,
since reversing changes every d_e to −d_e. Thus the combined
weight of this unoriented cycle is

    (x^ell−x) e_(r−1)(d_e : e in C).

Sum over disjoint sets of unoriented odd cycles on the moved
orbit indices, with a factor U for each unused index. Multiply
by G of the fixed set. This gives the exact polynomial G mod ell.
It shows why winding matters for the integer read but can be
invisible mod ell. Along x=2q, each U becomes q^ell−q and each
cycle factor contains the same divisor. If any vertex moves,
P as a polynomial mod ell is divisible by q^ell−q (at least
its square for a positive even number of moved orbits, by
parity of the projected cycle lengths). At one orbit this
reduces to G≡y^ell−x/2 and P≡q^ell−q. It agrees with the
free action on ordered Hamiltonian paths at q=1.

## Why exact autocorrelation alone is not the whole interface

Take three orbits of size three with internal +1 cycles,
and for each i<j allow an arc i→j precisely at the increment
delta_ij. Set delta_01=delta_12=0 and delta_02=z, where z is
0,1,2. Every individual interface has identical periodic
autocorrelation in all three cases (each is a translate
of a singleton word, and the reverse interface is its
translated reverse complement). Scores and simultaneous
cyclic symmetry also agree. But the H values are
2697,2625,2721 by literal packing count. The loop phase
delta_01+delta_12−delta_02 cannot be removed by relabelling
the origins of the three orbits. This is a block-cyclic
example, not a circulant Hadamard counterexample.

A separate elementary carry obstruction shows why one
should not copy ordinal product constructions wholesale.
If a sign row is the lexicographic product c_(im+j)=a_i b_j
with both block lengths n,m>1, then its shift um+v gives
C_a(u) times the noncarrying partial correlation of b and
C_a(u+1) times the carrying part. Flatness at v=0 would force
C_a(u)=0 for all nonzero shifts. At u=0,v=m−1 the first
part is n b_0 b_(m−1), nonzero, and the carrying part is now
zero. No such nontrivial product is a circulant Hadamard word.
This is a restricted no-go only, not a proof of the full
order classification.
