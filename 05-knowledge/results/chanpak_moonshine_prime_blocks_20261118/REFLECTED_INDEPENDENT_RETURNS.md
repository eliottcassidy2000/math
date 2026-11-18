# Independent reflected continuations with a repeatable exact return

**Status:** literal small witnesses, all-order paired return,
and a terminating inverse for the restricted count language.
This is not universal Chan–Pak coverage.

## Reflection preserves the debit, not the address

In [Riemann's functional equation](https://dlmf.nist.gov/25.4),
the completed function pairs s and 1−s; the symmetry alone
does not determine the location of its zeros. Here the precise
analogue is not an analytic assertion: an involution preserves
one integer part of the tournament transfer while changing
another. This lets us give two independent physical words the
same absolute output. It avoids the obstruction to identical
suffixes in [the ordered phase carry](ORDERED_PHASE_CARRY.md).

Recall the [retained-root alphabet](ANS_INVERSE.md). At profile
(H,e,t,d), a binary word of width k read low bit first has
numeric address W and ends at
```
H' = H + t W − d K(W) − e b_0,
t' = 2^k t − d rev_k(W),
e' = d b_(k−1),   d'=d,
```
where K is the sum of 2^(j−i−1)+[j=i+1] over pairs of ones at
i<j, and rev_k reverses the k bits. Reversal preserves the
distances between pairs of ones, so K(W)=K(rev_k(W)).
Leading or trailing zeros contribute no terms to K.

Suppose two literal, rooted tournaments have profiles
(H,e,a u,d) and (H,e,b u,d), with a,b distinct positive,
coprime odd integers. Find an odd k-bit W, let R=rev_k(W),
and require
```
a W=b R;       W=b v, R=a v.
```
The other endpoint bit is also one. Then append W to the
first tournament and R to the second, never identifying
their roots or requiring them to have the same arcs. Both
absolute counts agree because the leading terms, checksum
debits and endpoint debits agree. Moreover their next
weighted starts are respectively a(2^k u−dv) and
b(2^k u−dv). **The ratio a:b survives the return.**
The two ordinary end counts are now d in both tournaments.

More freedom is available: precede both reflected words
by j common zero letters, with j independently chosen at
each return. The appended addresses are 2^j W and 2^j R;
their reversals at width k+j are still R and W. If C=aW and
K=K(W), the all-order reduced return is
```
H' = H + C 2^j u − d K − e[j=0],
u' = 2^(k+j) u − dv,
e'=d.
```
No factorization theorem, root approximation or modulo-only
heuristic enters this equality.

## A finite test for the reflection itself

This is the trimmed reverse-multiple problem, not true for
every ratio. At bit position i the equation has coefficient
a b_i−b b_(k−1−i). The running ordinary binary carry c_i
satisfies
```
c_0=c_k=0,     c_(i+1)=(c_i+a b_i−b b_(k−1−i))/2.
```
Division must be integral at every step; forward carries are
in the fixed interval [−b,a] by induction. Read from both
ends toward the middle. A pair of carries (u0,v0) and an
outside bit pair (x,y) advance to
```
((u0+a x−b y)/2, 2v0−a y+b x).
```
The even middle closes when the two carries coincide. The
odd middle tests one common bit. For odd a,b the exterior
bits must agree by parity; if both are zero they can be
stripped without changing the ratio. Hence it is enough to
start with (1,1) and exhaust the finite carry automaton.
This decides every allowed width for the trimmed problem.
For example 7:19 has no such word; the ratio 11:13 takes
only four bits. This finite exact obstruction does not
rule out more general non-reflected pairs of continuations.

## Two physical returns

- The [1063 pair](ORDERED_PHASE_CARRY.md) has e=d=1 and
  a:b=1847:1967. The automaton finds
  W=51214312118786901 and R=48089900601626541
  (56 bits). Their shared debit K is 34179481178302758.
  The literal graphs after the first return have order 68
  and the common count **94558655002221104451**, exact by
  the transfer. Their start ratio is unchanged, so arbitrary
  further returns are allowed. The two order-68 arc tables
  are built and checked; their counts are not claimed to
  have been recomputed by full subset enumeration.
- A substantially shorter physical pair begins on two
  order-seven rooted tournaments with (H,e,d)=(29,0,3).
  Their outgoing row masks and zero-based roots are
  (0,41,67,5,15,29,59), root 2, and
  (0,65,83,7,11,31,57), root 6.
  Their weighted starts are 66 and 78, ratio 11:13.
  Use W=13 (1101 in ordinary notation), R=11 (1011), and
  K=8. The first return yields **863**, a prime, in two
  distinct marked states (863,3,1023,3) and (863,3,1209,3)
  at order 11. All four profiles here are confirmed by
  direct subset counting, not inferred only from the two
  formulas. The two words in physical order are
  (1,0,1,1) and (1,1,0,1).

## A Collatz-style inverse that really descends

From the 863 pair take u=93, d=e=3 and k=4. Every next u
is odd, and the exact predecessor is
```
s = v2(u'+3),   j=s−4,   u=(u'+3)/2^s.
```
The old H is recovered by subtracting its prescribed return
increment. Require s≥4, both states positive, and eventual
entry into the exact seed (H,u)=(863,93). Each legal inverse
above the seed strictly decreases u. Thus there is no room
for a cycle or divergence *on this restricted state language*.
This is why the Collatz analogy is useful: valuation chooses
an actual inverse branch. It is an odd-part map of u+3 on
certified interfaces, not a claim about 3u+1 or about arbitrary
input states being physically available.

There is an exact scalar test as well. Let T be the total
number of appended letters after the 863 seed, so T is a
sum of parts k+j≥4. If T_i are the successive partial totals,
then u_i/2^T_i=93−3 sum_(l≤i) 2^(−T_l), between 92.8 and 93.
The contribution of the last return is thus at least
143·92.8·2^(T−4)−27. The sum of all positive contributions
is at most 143·93·2^(T−4)·16/15, since successive partial
totals differ by at least four. In particular, with
A_T=143·93·2^(T−4), every nonempty output satisfies
```
863 + 0.99 A_T <= H <= 863 + 1.07 A_T.
```
These absolute bands for different T are disjoint.
For a prescribed count choose its only possible T (or reject)
and test all compositions of T into parts at least four
against the exact return, not against a rounded formula.
A positive hit constructs both tournaments. A negative hit
excludes only this language.

The number of compositions of T with parts at least four
satisfies f(T)=f(T−1)+f(T−4) for T≥4, with f(0)=1 and
f(1)=f(2)=f(3)=0: split by whether the first part is four
or can be decreased by one. This is the slope-three
recurrence from the earlier
[figurate stagger atlas](../../../01-canon/theorems/THM-1360-figurate-stagger-atlas-pisot-ladder-proth-mersenne-relational.md),
now counting actual decoder programs rather than a loose
similarity in growth. It bounds the number of distinct
absolute counts up to X by O(X^(log_2 lambda)), with
lambda>1 solving lambda^4=lambda^3+1. The exponent is
strictly less than one. Indefinite return does not imply
a universal selector.

**Next question:** what physical operation changes the ratio
or initial height while entering another accepted carry
automaton, and can these changing entries close the gaps in
the scalar bands? A successful bridge needs a physical entry
and an exact height, not a functional equation by itself.

Replay: `python3 reflected_return.py`. The outside-to-inside
carry test is in `reflection_search.py`; the replay exercises
both physical constructors, the marked inverse, and the exact
restricted scalar decoder.
