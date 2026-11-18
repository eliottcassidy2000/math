# A fixed-sum carry makes distinct tape products overlap

**Status:** all-depth algebraic carry, infinite physical overlapping
families for every depth m≥65, and a restricted exact target inverse.
There is no complete Chan–Pak selector here.

## What to retain from the summand and Zeckendorf work

The repaired [summand graph](../arithmetic_braids_20260917_summand.md)
uses pairs of *distinct* positive parents. As n arrives, it has
floor((n−1)/2) parent pairs: zero for 1,2, one for 3,4, and multiple
from 5 onward. These are three birth types, not three modulo classes.
The missing ascending arcs are the doubling diagonal n→2n;
restoring equal parents changes the graph. In the present construction
the fixed-sum fibre has two distinct ordered-pair realizations, and
the transfer is an edge within that fibre. Forgetting the two
parents and keeping the sum would erase exactly the change we need.

The [three-colour Zeckendorf audit](../modular_fourier_zeckendorf_routes_20261003.md)
distinguishes the actual marked extra-unit representations from
other colourings. Their marker is separate from the integer; the
three nonzero Fibonacci charges also have a neutral fourth state.
In particular, fixed-sum arithmetic uses an integer carry that
cannot be inferred from a colour alone. This motivates retaining
the whole arithmetic state here: **balanced sum, changed product,
binary carry debit and repair residue**. The step below uses binary
contact words, not Fibonacci expansion of a tape; the transfer
from Zeckendorf is the exact carry plus its marked realization.

## A product change, not an address swap

Use the [folded protected core](FIGURATE_OVERLAP.md) with
G=1105, g=17²G^(2m), B=448²(224·976)^(2m), unsigned root
reads U=Ubar, and signed reads
(V,Vbar)=17G^(2m)·(−17,−15). Put
```
P=2·17²G^(2m−1),  M=P/2,  g=G M.
```
An outer tape is a real transitive-deletion arm; its address
determines its contact bits. With two tape addresses x,y the exact
center, before a repair, is
```
h + Bxy + U(x+y) + g[C(x)+C(y)] − Vbar (x mod2) − V (y mod2).
```
The checksum satisfies C(2a+1)=C(a) and
C(2a)=C(a)+a+(a mod2), with C(0)=0. The retained repair choices
add P times 17j+130k for 0≤j,k≤2m, including
[2193,294m−2193] for m≥65.

Compare the two physical tape pairs
```
E = (4a,4b),       O = (4a+1,4b−1),
b=a+L=2t+1,       t=2z+eta, eta in {0,1}.
```
They have the same sum, but O's product exceeds E's by 4L−1.
Both have the same core and the same available repair palette.
The four checksum evaluations give
```
C(4a+1)+C(4b−1)−C(4a)−C(4b)
    = C(b−1)−C(b)−2a−3b−1
    = −18z+2L−6−8eta.
```
For the second equality use C(b−1)−C(b)=t+(t mod2).
The two low bits of each tape change from 00 to odd, so the
signed boundary term adds −V−Vbar=1040P. The difference
between the actual centers is therefore
```
B(4L−1) + g(−18z+2L−6−8eta) + 1040P.
```
Merely balancing the tape sum does not balance the product.
The checksum is the corrective term and must not be dropped.

## Solve the exact cancellation at every depth

First set 4L≡1 mod M. Let L0 be the positive inverse of 4 mod M,
put k0=(4L0−1)/M, and write L=L0+Ms. Since B is even, the
resulting normalized center difference (in units P) is the integer
```
c_eta + (2B+g)s − 9G z,
c_eta=(B/2)k0 + G(L0−3−4eta) + 1040.
```
The first congruence selects the product gap, and the second
selects the absolute debit. Because B and G are coprime,
d=gcd(2B+g,9G) is a divisor of 9. The actual residue schedule
is particularly small. Modulo 9,
B=4·7^m and g=4^m, so d=9 when 3 divides m and d=3 otherwise.
Also M=4g mod9. If m≡0 mod3, then k0≡L0+2 mod9 and
c_0≡6 mod9; in the other two cases k0≡L0−1 mod3 and
c_0≡0 mod3. Switching eta to 1 subtracts 4G≡1 mod9.
Taking the centered residue modulo d gives:

| Depth mod 3 | eta=0: center shift / P | eta=1: center shift / P |
|---:|---:|---:|
| 0 | −3 | −4 |
| 1 | 0 | −1 |
| 2 | 0 | −1 |

For each entry delta, solve the linear congruence
(2B+g)s≡delta−c_eta mod 9G by division by d. Then z is
the exact quotient (c_eta+(2B+g)s−delta)/(9G).
Advancing s along its solution progression increases z
linearly. The slope of b−2L per unit s is
[4(2B+g)/(9G)]−2M=(8B−14g)/(9G)>0.
Thus discard a finite prefix to ensure b≥2L, so
a=b−L≥L≥1. Both tape pairs are then positive, ordered and
balanced (larger at most twice smaller). Products increase
without bound along the remaining progression. This proves
infinitely many genuine distinct-product matches at each m.
For eta=0 and 3∤m the two centers even coincide exactly.

For eta=1 the shifts are nonzero. At m≥65 the full certified
progressions in the two representations overlap, with at most
four normalized steps lost at the edge. Relative to the lower
of the two centers, the rising lattice has three regions:
lower-only, common, upper-only. At depth 65 the common
part consists of 14,724 steps inclusive, and the extreme
regions each have one step. As with the summand graph,
the useful object is the representation fibre rather than a
single label assigned to the visible count. No missing target
outside this progression is declared absent from the spectrum.

## Target inverse and the two different continuation states

At fixed m and eta=1, both a,b are affine functions of the
nonnegative index n along the solution progression. Let
E(n)=(x,y) and use the lower envelope
Q(n)=h+Bxy+U(x+y). Both x and y grow by at least four per
step. The checksum bound 0≤C(i)≤i gives an upper error
g(x+y)+294mP, and 294mP<B/2 by the old folded-core
bound. The other center differs by at most 4P and lies above
Q(n), since on E the checksum alone contributes at least
3g(a+b)≥6g. The next quadratic envelope starts higher:
Q(n+1)−Q(n)≥4B(x+y). Thus the two-address windows at different
indices are disjoint in absolute height. Given N, find the last
Q(n)≤N by exponential and binary search and try *both* exact
debits. A successful return constructs an actual tournament;
a negative answer rules out only these progressions.

A shared H here has **two distinct marked continuations**,
unlike the earlier fixed-address cap reflection. The tape
widths do not change: x is a positive multiple of four and
x+1 has the same width; y=4b with odd b≥3 is never a power
of two, so y−1 also keeps its width. For right width w the
weighted start read of the glued graph is
```
U_final = 2^w (Bx+U)+g U_right.
```
The product-changing step increases x by one, whereas
|Delta U_right|≤2^w−1. Therefore the two weighted starts
differ by more than 2^w(B−g)>0. The six-read interfaces
agree *within* each repair window, not between the two
tape pairs. Appending a sink keeps the common H and sets
the ordinary root end count to zero; a subsequent one-back-arc
singleton then makes the new counts differ by twice the
original difference in weighted root starts. Thus an overlap
can switch which continuation is available at one fixed
absolute count. This is a concrete use of the Zeckendorf
lesson of keeping the marked representation, not just H.

**Scope and next attack.** This construction resolves the
different-products question constructively; the solution
progression still leaves gaps between its quadratic envelopes.
It does not prove that these gaps cannot be filled by other
tape pairs. A next proof test is to find independent carries
that reach different solution progressions while allowing the
same target to change the continuation state. Merely adding
more outputs on this one repair lattice cannot settle absolute
coverage.

Replay: `python3 distinct_products.py`. It checks an independent
small literal graph calculation, the three depth classes,
inverse checks at indices through 10^30, and compiles both
tournaments at one shared depth-65 target of 2101 digits.
Counts of the large compiled graphs follow from the all-order
transfer; they are not asserted to have been counted by subset
enumeration.
