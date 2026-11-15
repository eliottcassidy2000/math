# An exact observable quotient, and two terminal probes that undo it

The [previous physical window](../chanpak_collatz_norm_clock_20261115/README.md)
has counts
```
H_j(x,y)=432xy+129x+258y+15[x odd]+30[y odd]+45W+135+2j
W=p(x)+p(y).
```
The construction realizes three literal counts at increments +2. The
subsequent results are about what a reader can detect, not about
producing every absolute count.

## 1. A complementary word has a hidden reflection

Use $A=x+1$, $B=y$, $r=[y\text{ odd}]$, and $K=W-B+r$, and keep complementary input
parity. Set $P=AB$, $D=A-B$, and $T=A+B$. Appending the bit $1-t$ to x and
$t$ to y gives the exact transition
```
A' = 2A − t, B' = 2B + t, r' = t,
P' = 4P + 2tD − t², D' = 2D − 2t,
K' = K + t(D − 2r), T' = 2T.
```
For example, when t=1 the checksum increments by x+[x odd]; subtract
the new B and add the new parity to obtain the displayed K law.
The observed count is
```
H_j = 432P + 129D + 45K + 21 − 30r + 2j.
```
Consequently no complementary continuation, of any length, distinguishes
two states of the same parity with the same (P,D,K). In particular the
involution $J(A,B,K,r)=(-B,-A,K,r)$ commutes with every complementary
letter and leaves every ordinary reading unchanged. Its lost coordinate
is the sign of T, where $T^2=D^2+4P$. This is a symmetry of the formal
transition system and its finite reductions, not permission to use
negative addresses in a tournament.

Over the integers, from an odd-x, even-y state, the three readings at
suffixes empty, 1, 11 (x bits) recover all of (P,D,K): the three
coefficient rows, after dividing columns by 432,129,45, are
(1,1,1),(4,2,1),(16,4,1), with determinant −6. The positive real branch
of T then determines A,B uniquely. Over a field with prime characteristic
not 2,3,5,43, the same observations recover precisely (P,D,K); any two
distinct indistinguishable states on this parity sheet are related by J.
On composite residue rings these reflection choices can vary by CRT
component; do not assert only two possibilities there.

At prime 43 the coefficient 129 vanishes, but complementary observations
still recover (P,D,K). Read before/after one x-bit 1 to recover P and K
(432,45 and 4−1 are all invertible). Repeat the two readings after
appending x-bits 01, with y-bits 10. The new product is
$P'=16P+8D-4$, which recovers D without the coefficient 129.
Thus for every prime other than 2,3,5 this is exactly the observable
quotient, not merely a sufficient list of states.

## 2. Two noncomplementary terminal probes remove the reflection

Start again with r=0. Probe U appends 1 to both addresses, so the new
addresses are (2A−1,2B+1), with unchanged W=K+B. Its exact readout is
```
U_j = 1728P + 1122A − 303B + 45K − 123 + 2j.
```
Probe V follows U by appending x-bit 1, y-bit 0. The addresses are now
(4A−1,4B+2) and the checksum is K+3B+2. Its readout is
```
V_j = 6912P + 3972A − 561B + 45K − 237 + 2j.
```
Therefore $U_j(s)-U_j(Js)=819T$ and $V_j(s)-V_j(Js)=3411T$.
Both probes are physical appends on positive states; both still admit
the full center bank at increments +2. Since gcd(819,3411)=9, for every
prime except 2,3,5, ordinary complementary probes plus these terminal
probes recover the entire prefix state (A,B,K) modulo that prime.
Indeed one terminal gives the missing sign on each field with T nonzero;
T=0 is already fixed by J. Equivalently, once P,D,K are known, the two
terminal equations have T coefficients 819/2 and 3411/2. This reasoning
also recovers the state over any odd modulus coprime to 15 by CRT and
the linear equations, rather than by assuming a composite ring has only
two square roots. For prime powers of 43 the direct two-reading method
still works because it uses units independent of 129.

There is a meaningful arithmetic exception: the single both-one probe
still does not see this reflection modulo 7 or 13, even though it breaks
it over the integers and dyadic rings; the second probe removes that
particular blindness. The dyadic behavior is subtler than an
odd-field inversion because division by 2 is not available. This note
does not assert full dyadic observability. The earlier mod-4 escape
and physical modular clock remain valid and independent.

## 3. When an analogy transfers

The complementary system, with its physical suffix composition and
its ordinary reader, really does factor through (P,D,K). This is why
the pattern persists under every repetition, not merely for a list of
observed words. Changing a reader changes the quotient, and the
breakout coefficients specify exactly at which primes the analogy
still works. Compare the [older marked Bott–Collatz note](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_bott_marked_periodicity_20261007d.md):
repetition of the whole marked representation is much stronger than
repetition of its spectrum or its abstract similarity class.
Cayley–Dickson gives the complementary warning: the quadratic length
persists after doubling even when the multiplication identity that
made it a norm stops transferring. Here the test is concrete: specify
what is preserved by every admissible composition and change the reader
if it does not distinguish the state needed by the next move.

This is a probe theorem, not an exact-target placement theorem. It
identifies what continuation data the previous repair family really
uses and two ways to recover data previously invisible to its reader.

## 4. An exact periodic-word clock, and its price

A full repeated motif of h complementary bits on the y side, ending
in 0 and starting on r=0, has the form

    D' = L D − 2v,       K' = K + ρD + c,       L=2^h,

where v is the motif as a binary integer, ρ its bit reversal, and c
depends only on the motif. This follows by summing the elementary
transitions; the coefficient of the incoming D at letter i (numbered
from zero) is t_i·2^i. Put D*=2v/(L−1) and
Z=K−ρ(D−D*)/(L−1). These changes have odd denominators and are legal
at every dyadic precision. One period acts exactly as

    D−D* ↦ L(D−D*),     Z ↦ Z + (ρD*+c).

Thus successive readings of a periodic word need not return to the
same state. The transverse coordinate contracts in the 2-adic
metric, while the retained coordinate translates. This is a
literal theorem about the tournament clock, not an application
of Bott's period-eight theorem.

For the y-word 10, L=4, v=2, ρ=1, c=0, so the drift is 4/3.
After a fixed transient depending on k, varying the repeat count n
modulo 2^(k−2) sweeps a full checksum class modulo 2^k: in
K_n=K+4n/3+(D−4/3)(4^n−1)/3 the changing geometric term vanishes at
that precision once 2n≥k. Pad all choices to the same length with
y-zeros, which leave K unchanged; at least k padding zeros also reset
the endpoint y modulo 2^k. Complementary padding on x preserves the
actual address sum across the choices. This shows a second mechanism
for the old clock: iterable drift, rather than independent digits.
Its length is exponential in k, so it does not beat the quadratic
guard. The period of the local readout and the cost of placing it
are two different quantities.
