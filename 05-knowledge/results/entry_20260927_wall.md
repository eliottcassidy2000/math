# A genuine carry wall, its exact division clock, and why it is not yet a temporal wall

**PROVED elementary identities / FINITE-EXACT controls / OPEN temporal
amortization, 2026-09-27. No novelty claim.** This note transfers an actual
absorbing mechanism from the Gilbreath discussion, with the spatial and
temporal clocks kept separate. It does not prove Collatz or Gilbreath.

## 1. Inheritance and the live board

The closest proved mechanism is [THM-4511, the size-four Gilbreath
wall](../../01-canon/theorems/THM-4511-gilbreath-size-four-wall-theorem.md).
For a row beginning with1 and then values in{0,2,4}, an initial2 before
the first4 protects the leading1 forever, including all re-emissions.
Each subsequent column has a{0,4} phase followed by a{0,2} phase. This
uses the closed alphabet, directional dependence, and the absence of a
later return from2 to4. Allowing6 breaches that argument:

    [1,2,6] -> [1,4] -> [3].                       (W0)

The relevant Collatz mechanism is the exact binary carry scan in
[the excursion forest](forest_20260926_excursions.md), supplemented by
[the macro-carry synchronization theorem](nextforest_20260926_automata.md),
section4. The canonical temporal hostile is the growing family containing27
from [the swap-lift note](reset_20260926_swaplift.md), section6. The corrected
near miss is the earlier claim that additive difference dynamics forbid a
Collatz edge encoding: [the signed-difference audit](forest_20260926_gilbreath.md)
shows that a trace can be encoded, but its computation can be hidden in an
infinite seed. Preserving the finite arithmetic source is indispensable.

The least-used sidecar here is the *image of the entire carry bank*, rather
than the one known physical carry. It measures sensitivity to a boundary
condition, not the size of the actual orbit. The live board has five items:

| Object | Predicate | Missing coordinate or cheapest hostile |
|---|---|---|
| Gilbreath column | no4 after entering its{0,2} phase |6 defeats alphabet closure |
| Single carry scan | boundary dependence has disappeared | Does the next row preserve this? |
| First repeated source-bit pair | exact position of the first output1 | retain the unread integer tail |
| Ordinary finite binary source | canonical digits end | the next row can have a longer word |
| Recursive descent certificates | selected iterate is smaller | inherited27 family forbids a fixed horizon |

## 2. A true irreversible phase in one carry scan

For a source bit b in{0,1}, tripling with carry uses

    f_b(c)=floor((3b+c)/2), c in{0,1,2},
    emitted bit=(3b+c) mod2.                        (W1)

Thus f_0 maps(0,1,2) to(0,0,1), and f_1 maps it to(1,2,2).
Let S initially be{0,1,2} and update S to f_b(S) while scanning from low
bits to high bits. The cardinality of S cannot increase under any input.

**Single-scan wall theorem.** A nonempty binary word synchronizes all
three incoming carries if and only if it contains00 or11. Before that
first repeated adjacent pair, S has exactly two elements; after it, S
has exactly one element permanently, under every continuation.

Proof. After the first bit, S is{0,1} when that bit is0 and{1,2} when
it is1. An alternating next bit interchanges those two sets. A repeated0
sends{0,1} to{0}; a repeated1 sends{1,2} to{2}. A function maps a singleton
to a singleton. This proves both directions and the permanence. Hence the
only nonsynchronizing length-L words, for L>=1, are the two alternating
words; the reset count is2^L-2, agreeing with the m=3 specialization of
the inherited macro-carry census. The empty word has rank3.

What disappears is dependence of the *future carry and future emitted
bits* on the entering carry. Already emitted low bits can differ and are
not erased from the arithmetic output. For the actual integer there is
already one specified carry c_0=1; the three-state bank is a controlled
comparison of three boundary conditions.

This is an exact finite absorbing phase, as requested. Its clock is the
digit position within a single multiplication, not the Collatz iteration.

## 3. The first wall is precisely the Collatz division exponent

Let n be positive odd, let b_i be its binary digits padded by zeros beyond
the canonical most significant bit, and put

    s=min{i>=1: b_i=b_(i-1)},
    a=v2(3n+1),     U(n)=(3n+1)/2^a.

**Division-clock theorem.** The wall is always finite and s=a. If a is
even it is a00 wall with outgoing carry0; if a is odd it is an11 wall
with outgoing carry2. In the physical scan c_0=1, all emissions before
position a are0 and the emission at position a is1.

Proof. Since n is odd, b_0=1. Until the first repeated pair the input
is1,0,1,0,... . Reading1 with carry1 emits0 and leaves carry2;
reading0 with carry2 emits0 and leaves carry1. If the first repeat is
at even position a, its bit is0 with entering carry1, which emits1 and
leaves0. At odd position a its bit is1 with entering carry2, which emits1
and leaves2. Finite zero padding ensures some repeat occurs. The first1
in the output3n+1 is therefore at exactly a=s.

This identifies the new spatial wall with an existing exact arithmetic
budget. It does not supply an additional independent contraction.

There is a useful complete coding of the remaining state. Define

    epsilon_a=1 if a is even, and5 if a is odd,
    r_a=(epsilon_a*2^a-1)/3,
    h=floor(n/2^(a+1)).

Then, for every positive odd n,

    n=r_a+2^(a+1)h,       U(n)=6h+epsilon_a.         (W2)

Conversely every a>=1 and integer h>=0 in(W2) give a positive odd source
whose first wall is exactly a. Indeed0<r_a<2^(a+1), and
3n+1=2^a(epsilon_a+6h) with an odd parenthesis. This proves the converse
as well as the forward formula. The parity of a specifies exactly which
nonzero odd residue modulo6 the successor occupies; the nonnegative
integer h is unrestricted. Thus the unread tail is, up to a fixed affine
map, the *entire next odd iterate*. Losing it loses the target arithmetic.

For example the a=1 cylinder is n=3+4h and U(n)=5+6h. The a=2 cylinder
is n=1+8h and U(n)=1+6h. The initial synchronizing event is identical
within either cylinder, while the successor values are unbounded.

**Finite-support endpoint.** The canonical binary word of a positive odd
n is strictly alternating if and only if

    n=(4^k-1)/3 for some k>=1,                      (W3)

and these are precisely the sources with U(n)=1. The canonical word is
101...1 of odd length2k-1, so its first repeated pair appears only in the
zero padding at position2k. This also follows from(W2): h=0 and even a
give the direct basin, and an odd wall cannot give successor1.

Consequently every ordinary odd input either contains an internal carry
wall or already lies in the direct basin of1. This is a complete local
dichotomy, not a termination proof: having an internal wall merely invokes
the unrestricted tail transformation(W2).

## 4. Why the phase does not persist under temporal evolution

The smallest odd source with an internal wall is3, whose canonical word11
is synchronizing. Its next odd iterate is5, whose canonical word101 is
nonsynchronizing. The actual first-wall positions, after zero padding,
are1 and4. Thus “an internal wall has appeared” can become false on the
very next physical Collatz row. This is a minimal positive-integer
counterexample to importing Gilbreath's temporal phase order verbatim.

Even a large number of spatial walls gives no bounded first-descent
horizon. The inherited family

    n_k=4*8^k-5, k>=1,
    27,251,2043,16379,...                           (W4)

has canonical length3k+2, first wall at position1, and exactly3k-1
equal adjacent pairs within its canonical digits. Its low bits are
1,1,0,1,1,...,1; only the two pairs adjacent to its single0 are unequal.
Nevertheless the swap-lift theorem proves

    U^(2j)(n_k)=4*9^j*8^(k-j)-5, 0<=j<=k,
    U^(2j+1)(n_k)=6*9^j*8^(k-j)-7, 0<=j<k.          (W5)

Every one of the first2k iterates is strictly above n_k. Therefore the
first-descent time, allowing infinity, is greater than2k despite the
fixed earliest wall and the unbounded number of other walls. This is a
source-varying family, not an infinite growth claim for one integer.

There is still a sharp positive local conclusion: for n>1,

    U(n)<n if and only if a>=2.                     (W6)

For a=1, U(n)=(3n+1)/2>n; for a>=2, U(n)<=(3n+1)/4<n.
Thus the wall position exactly decides *one-step* descent. The issue is
compensating earlier growth on an actual orbit, not whether a wall is
present or whether a single descent has occurred.

The inherited macro-carry theorem gives another sharp scale distinction.
After a parity prefix with j odd steps, the unread-tail multiplication
has3^j possible carries and no synchronizing word shorter than
ceil(j log2(3)). A two-bit wall for one multiplication does not remain
a two-bit wall for the macrostep. This statement is about the complete
residual machine, not a complexity lower bound for one fixed integer.

## 5. What transfers, and what does not

| Source -> target | Map and preserved predicate | Lost data / required sidecar / hostile |
|---|---|---|
| Boundary uncertainty -> carry uncertainty | S -> f_b(S); singleton images are absorbing under every continuation | already emitted low bits; retain their ordered word |
| Actual odd input -> first carry wall | scan from the low bit with c_0=1; first repeat equals v2(3n+1) | unread h; retain(W2) exactly |
| Canonical finite support -> padded scan | append zeros; every scan eventually has a wall | the next row's end marker;3->5 reverses the internal-wall predicate |
| Many internal walls -> proposed descent clock | no valid monotone map is established |(W4)--(W5) refute a uniform horizon even with first wall fixed |
| Gilbreath phase -> proposed Collatz temporal phase | both have an absorbing comparison state within their proved operation | there is no map commuting the two time evolutions;3->5 is the cheapest test |

The connection to [HYP-9162, the signed Hadamard tournament
tower](../hypotheses/HYP-9162-sierpinski-tournament-tower-frobenius-21.md)
remains a separate construction. Its skew-Hadamard doubling is proved;
the claimed automorphism equality at all levels remains OPEN. The
Mersenne sizes and zero-triangle counts in a single-seed XOR sea do not
identify a pairwise arithmetic observable on carry states. The carry
maps here are many-to-one transformations, with literal ties, so replacing
them by a tournament would discard the very coalescence that proves(W1).
No tournament is introduced, and no count of XOR zero regions is used as
an orbit distribution law. The shared useful operation is loss of boundary
dependence, not a coincidence of graph sizes.

The strongest next question is consequently narrow: can the exact
replacement h ->6h+epsilon_a be coupled to a source-preserving certificate
whose successful phase survives subsequent replacements? A wall flag
alone cannot do this. The existing swap-lift bank gives selected descent
certificates with exact tail bounds; the remaining entry problem needs
either a persistent arithmetic invariant or an amortization of regenerated
tail complexity. Equation(W2) makes the coordinate that must be paid for
explicit, while preserving the fixed ordinary source and its finite input.

## 6. Reproduction and scope

Run from the repository root:

    python -X utf8 -B 04-computation/experiments/entry_20260927_wall.py
    python -X utf8 -B -O 04-computation/experiments/entry_20260927_wall.py

The [script](../../04-computation/experiments/entry_20260927_wall.py) and
[saved output](entry_20260927_wall.out) check all131070 binary words of
lengths1..16; every131072 positive odd source below2^18; every29523
Gilbreath tail over{0,2,4} of lengths1..9, including every visible column's
phase order; both minimal hostile examples;128 members of(W4); and16384
independent(a,h) cylinder constructions. All checks use explicit failures
that remain active under Python optimization. Finite controls support
the all-height algebraic proofs; no tested cap is taken as an orbit bound.
