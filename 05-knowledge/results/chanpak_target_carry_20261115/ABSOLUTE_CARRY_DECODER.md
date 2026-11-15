# The target chooses the high sum: no all-ones rounding

**All-order theorem with exact finite certificates.** For each
M in {2,3,4,5,6,8}, every sufficiently large odd integer N has two
literal strong tournament counts bracketing it, both congruent to
N modulo 2M, at absolute distance O_M(N^(3/8)). This improves the
previous O_M(N^(2/5)) congruent exponent **for these moduli**.
No exact universal placement or extension to all moduli at exponent
3/8 is asserted. Modulus 1 already follows from the unconstrained
two-tape result.

## 1. Exact finite compiler for arbitrary prescribed high sum

Use p from the [paired-tape transfer](../chanpak_absolute_steering_20261115/TWO_TAPE_STEERING.md)
and write c(x)=p(x)/2. Its appending recurrences are
```
c(2x+1)=c(x),       c(2x)=c(x)+(x+(x mod 2))/2.
```
Consequently, for a fixed word b of length L appended to u, the
increment c(2^L u+b)−c(u) modulo M depends only on u modulo 2M:
each step's input and parity determine its drift, and appending
preserves equality modulo 2M.

Fix the exact desired sum K=U+V. Put D=2^L, s=K mod D, and try
one aligned block U=Dj+b with 0<=b<D. If b<=s, use low V digit s−b
and high prefix k=floor(K/D)−j. If b>s, use low digit s+D−b and
high prefix k=floor(K/D)−j−1. Keep the larger of the two segments
(ties choose the first). Thus the **target** chooses the branch
using its actual s; neither K nor s is moved to an all-ones word.
When the two prefixes are nonnegative the entire selected segment
gives physical high addresses.

Here is the exact finite certificate. For the indicated L, for
every s, both high prefixes modulo 2M, and each delta in {0,1},
the selected segment attains every pair of residues
```
  U mod oddpart(M),  [c(U)+c(V)+delta(U mod 2)] mod M.
```
The checked rows are
```
M       2    3    4     5     6    8
L       5    8    6    10     9    8
```
The [certificate](replay_target_carry.py) computes c by the displayed
recurrence on 0<=x<2MD. It checks all s and all (2M)^2 representatives
with a bitset of missing pairs, for both deltas. Replacing a representative
with any actual prefix in its class adds a constant to the second output
and leaves the first unchanged, by the drift observation. Thus the finite
certificate proves the displayed formula at **every** prefix height, not
merely up to a numerical cutoff. In the actual decoder, only the required
syndrome is tested: take the least b that meets it, without generating a
menu of desired output residues.

## 2. Absolute placement with that certificate

Use one fixed cyclic carrier and the physical formula
`F(x,y)=B+Qxy+Tx+Sy+45[p(x)+p(y)]` from the earlier note, with
(B,Q,T,S)=(2855,5760,4034,3914). Choose an integer a within O_M(1)
above sqrt((N−B)/Q), initially floor(sqrt((N−B)/Q))+2. Take q=2^r
of order a^(1/4); at large N it is divisible by 2^(v2(M)+1).
Write a=qu+v, 0<=v<q, and K=2u. If necessary increase a by O_M(1)
to achieve
```
  B+2Sq u+(S+T)v ≡ N mod d,   d=gcd(45,M).
```
Indeed (S+T) is a unit mod d, and as v ranges across the q-cell
each residue occurs in d steps; crossing a cell boundary adds
at most another d steps. This only moves the center by O_M(1).

Use x=qU+v and y=qV+v with U+V=K. In the arbitrary-K compiler
specify U=0 mod m=oddpart(M). Therefore the smooth terms of F are
fixed mod 2M (q handles the dyadic component). Since K is even,
the binary split gives
```
 p(x)+p(y)=2p(v)+(q−1−rev_r(v))K
           +2[c(U)+c(V)+delta(U mod 2)],
 delta=1−bit_(r−1)(v).
```
The first two summands are fixed and even. The rightmost term
times 45 supplies an arbitrary multiple of 90 mod 2M; the displayed
condition modulo d and oddness of the counts give exactly the
necessary solvability condition modulo 2d. This is why the finite
certificate checks *both* delta channels; dropping the last adjacent
digit would give a wrong congruence.

Now program one actual U,V pair per aligned block, from the first
block wholly on or above U=u through the penultimate block before
U=K. At the first block, x,y=a+O_M(q), the product loss off the
midpoint is O_M(q²), and the initial extra margin above the prescribed
target is of order sqrt N. All the other scalar terms in F are
nonnegative, so this count is above N. At the final block, y=O_M(q),
so that count is below N. Bisect on the **actual programmed counts**
to adjacent blocks on opposite sides; monotonicity of the full
counts is not assumed.

Write x=a+t and y=a−t. The smooth part is a constant minus Qt²
plus (T−S)t; its peak differs from N by O_M(a). The unfrozen
surcharge is between 0 and 45(K+2)=O(a/q). The upper witness
therefore has |t|=O_M(sqrt a). The adjoining block changes t by at
most 2qD=O_M(q), and the two actual outputs differ by at most
```
 O_M(q sqrt a+q²+a/q)=O_M(N^(3/8)).
```
This closes the previous alignment loss without making a numerical
choice of checksum independent of the word. It uses O(log N)
block probes, each testing at most D fixed-size options after
the two high-prefix charges have been evaluated.

## 3. Exact successor graph and a test for further moduli

There is a finite-state test for whether the same larger-segment
certificate ever holds at a given M. Read the bits of s from the
most significant end. Keep the residues of the two prefix addresses
mod 2M, charge mod M, and the carry c_(i+1) into the next higher
column. For an input column s_i, enumerate binary digit pairs and
c_i with x_i+y_i=s_i+2c_(i+1)−c_i. Append the chosen digits, advance
the two exact modular drifts, and at the last column require c_0=0.
The largest-segment choice (leading bit and the tie word 011...1)
is a regular condition; use it to choose c_L. For each pair of
starting residues, each required output pair and each delta, the
successful input words form a regular language. Intersect these
finitely many languages. It is decidable whether this intersection
contains **all** words of some positive length: determinize, and
iterate the set of reachable states under either next symbol.
The sets eventually repeat; a length works iff the reachable set
at that length contains only accepting states. This is a decision
procedure for this certificate at a given M, not a bound on the
required length, nor a proof that each M passes.

For a particular target, finite suffix reachability marks the
viable states; choose the lexicographically first viable bit pair
at each step. On viable nonterminal states this has precisely one
outgoing edge. It preserves the exact sum and the required residue,
and the number of remaining columns proves termination. The
functional-graph comparison with Collatz is literal but limited:
a deterministic outgoing edge does not supply viable states for a
prescribed final **integer**, and the Collatz forward graph has no
remaining-column rank. The final placement error still need not
vanish.

## 4. Direct Rule 30-style scalar control at the first output bit

The [weighted cut law](CARRY_CUT_LAW.md) shows what the first
certificate row is doing. Take L=6 and keep the same larger
segment chosen by the target. For each of the 64 values s, there
is a b and a free position j in {1,2,3}, counted from the right,
where the number of other free positions at distance 1 or 2 is
odd. The free positions are the ones of b XOR (c−b), with c=s or
s+64 according to the branch. This is the 64-case test in
[replay_first_bit.py](replay_first_bit.py). No high-prefix
residues need be enumerated: all neighbors of j lie inside the
six columns, and the bit at the physical bottom is unchanged.

Swap the owner of this one free digit. The exact sum and the
entire carry word are unchanged. The cut law changes the
checksum by an odd integer, so the target's residue forces one
of the two actions. The possible adjacent-digit correction
delta(U mod 2) is also unchanged. This is stronger than finding
both parities somewhere in the block: both choices are in the
same carry cell and differ at one specified site. The next
obligation is to share guards for successive sites without
forcing large amounts of inaccessible alignment.

This is inspired by [Rule 30's sideways reconstruction](https://writings.stephenwolfram.com/2019/10/announcing-the-rule-30-prizes/),
which solves a single binary input from its output and the
other two inputs. The identities are not the same automaton.
The tournament checksum has larger, precision-dependent read
radius, consistent with the repo's
[Rule 30 strip and lost-current work](../nextforest_20260926_automata.md).
