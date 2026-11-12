# A complete marked-spine modular grammar and an infinite prime-birth inventory

**Status: PROVED elementary; small literal controls FINITE-EXACT.**
This is a complete grammar for a specified infinite class, not all tournaments.
The point is to display exactly what “knowing a family intimately at every
size” can mean, and where the next family parameter must be admitted.

## 0. Inheritance and information ledger

The binary presentation and Hamilton-path transfer for a tournament with a
marked vertex whose deletion is transitive are inherited from
[`feedback_word_hamilton_spectrum_20261020/proofs.md`](../../feedback_word_hamilton_spectrum_20261020/proofs.md).
That source already proves the global omission of the attainable count13
from this construction class. We do not call the binary presentation new.

The prime-extension partition is inherited from
[THM-4169, prime-parent augmentation](../../../../01-canon/theorems/THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer.md).
The usual W critical family is classical, and the repository already uses
its criticality. Our elementary proof below needs no external classification
of all critical tournaments. The added package is the **complete module
interval dictionary, induced-prime inventory, and recursive birth reading**.

The closest hostile is the prime-four-vertex gap: there is no prime tournament
on four vertices, although prime tournaments recur at every size≥5. The
corrected near miss is treating an infinite family of primes as a finite
list of substitution atoms. The least-used sidecar here is the multiplicity
of each induced prime core, not merely its presence or Hamilton-path count.

`Source → target`: a pointed tournament with transitive deletion → its binary
word → module intervals and induced-prime inventory. The word loses no pointed
isomorphism data. The inventory loses run positions and lengths, generally
loses the tournament, and is not asserted to determine future arbitrary cuts.
Forgetting the marked vertex introduces further isomorphism identifications.

## 1. Exact objects and complete module dictionary

Let the spine be `v_0→v_1→...→v_(m−1)` with every arc directed from lower to
higher index, and let `x→v_i` iff the letter `b_i=1`. Write `T(b)`.
The empty word is allowed and gives the singleton apex.

The mapping from binary words of length m to **apex-marked** isomorphism classes
is a bijection. Any isomorphism preserving x must preserve the unique linear
order on its transitive deletion. In particular there are exactly `2^m`
such marked classes. This is not an unrooted tournament-class census.

A module is a vertex set seen uniformly by every vertex outside it. The
following list gives **all proper nontrivial modules** of T(b):

1. A spine interval `[a,b]` of length≥2 on which all letters agree.
2. `{x}∪[a,b]`, with a nonempty spine interval, where every omitted prefix
   letter is0 and every omitted suffix letter is1; exclude the whole graph.

**Proof.** A spine-only module cannot omit a point between its extreme indices:
that point beats one part and loses to the other. It is therefore an interval.
The only remaining outside vertex not automatically uniform is x, giving the
constant-letter condition. For a module containing x and a nonempty spine
set, the same between-point argument forces an interval. An omitted prefix
vertex beats all included spine vertices, so it must beat x, meaning bit0.
An omitted suffix vertex loses to all included spine vertices, so it must lose
to x, meaning bit1. Conversely these conditions make every outside vertex
uniform. The singleton{x} is trivial and the empty spine interval contributes
nothing else. QED.

This is already a complete all-size decision rule for every module of every
object in this class, with no exceptional sizes needing discovery.

### Count without listing modules

For m≥1 let `r_1,...,r_s` be the constant-run lengths, p the initial zero-run
length (zero if the first bit is1), and q the terminal one-run length. Then

`number of proper nontrivial modules = Σ_i binom(r_i,2)
                                      +(p+1)(q+1)−1−[p+q=m]`.       (1)

The first term counts constant spine intervals. For the second choose how
many leading zeroes and trailing ones to omit. Remove the whole-graph choice;
when the word is `0^p1^q`, remove also the unique choice omitting all spine
vertices. For m=0 the module count is0 separately.

## 2. Canonical pointed modular normal form

Strip the initial zeroes and final ones. If nothing remains, T(b) is transitive.
Otherwise the reduced word begins1 and ends0. It has an even number2r of
nonempty constant runs, alternately1 and0. Replace each run by a single
vertex. The quotient, including the apex, is

`W_(2r+1) = T((10)^r)`.

Thus the full marked object is an ordinal prefix, followed by W_(2r+1) with
its spine vertices replaced by transitive chains of the run lengths and its
apex left singleton, followed by an ordinal suffix. This representation is
unique with the marking retained. It constructs the object, rather than
merely counts it.

The complete module dictionary implies that for n≥3 the object is prime
(no proper nontrivial modules) iff its word is `(10)^r`, r≥1. In particular,
W_(2r+1) is prime. Every one-apex prime has this form. This conclusion is
about the one-apex class, not all prime tournaments.

Appending a bit has a fully specified structural action: extend the final
run if it agrees, or append a run if it differs, then apply the same prefix/
suffix rule. The run lengths are the exact sidecar needed to reconstruct the
new tournament. The prime quotient grows when a new completed10 alternation
is exposed; a longer homogeneous run enlarges a module rather than a prime
quotient. The marked extension tree itself is exactly the binary word tree.
More precisely r equals the number of adjacent `10` transitions in the whole
word. Appending0 after1 increments r by one; every other append leaves r
unchanged. This is an exact prime-arity birth clock. It predicts the largest
prime quotient, but not its multiplicities; those require the next section.

## 3. All induced prime types and their multiplicities

For j≥1 let `p_j(b)` be the number of vertex subsets inducing a prime
tournament of order2j+1. Every induced prime of order≥3 contains x, since
spine-only subtournaments are transitive. By §2 its selected spine subsequence
must be exactly `(10)^j`. Consequently

`p_j(b) = number of subsequences (10)^j in b`,                       (2)

and there are **no other induced prime types of order≥3**.
These are vertex-subset multiplicities, not isomorphism-class counts: all
subsets counted in p_j have the same prime type W_(2j+1).

Set `O(t)=1+Σ_(j≥1) p_j t^j`. Keep also `E(t)`, whose t^j coefficient counts
subsequences `(10)^j1`. Initially `(O,E)=(1,0)`. Appending one letter gives

`1: (O,E) ↦ (O,E+O)`,
`0: (O,E) ↦ (O+tE,E)`.                                           (3)

Each update simply either skips the new vertex or uses it as the last
required letter. This two-polynomial transfer is exact at every size.
For an eventually periodic word, its matrix powers give rational generating
functions in the number of repeated periods. No bounded-word experiment is
needed for this conclusion.

If the reduced word has2r runs, its largest induced prime has order2r+1;
every W_3,W_5,...,W_(2r+1) occurs and none larger occurs. Selecting one vertex
from each run constructs the largest one. An alternating subsequence cannot
use more than one vertex from a constant run. A truncated choice constructs
every smaller member. Thus the first-order structural question has an exact
answer for the entire class.

### The critical tower's closed inventory

Write `P_k(t)=O_(10)^k(t)`. Then

`P_k(t)=Σ_(j=0)^k binom(k+j,2j)t^j`,                               (4)
`P_0=1, P_1=1+t`,
`P_k=(2+t)P_(k−1)−P_(k−2)`,                                      (5)
`Σ_(k≥0) P_k(t)z^k = (1−z)/(1−(2+t)z+z²)`.                       (6)

**Counting proof of (4).** Write selected one-indices as2a_i and zero-indices
as2b_i+1. The conditions are `0≤a_1≤b_1<a_2≤b_2<...<a_j≤b_j≤k−1`.
The2j numbers
`a_i+i−1, b_i+i` (in their displayed alternating order) are strictly
increasing in `{0,...,k+j−1}`. Subtracting the shifts gives the inverse.
Hence there are binom(k+j,2j) choices. Equations(5) and(6) follow equivalently
from the period matrix `[[1+t,t],[1,1]]`, determinant1 and trace2+t.

At t=1, P_k(1) is the odd-index Fibonacci number F_(2k+1), with F_1=F_2=1.
The total number of prime induced subsets of order≥3 is P_k(1)−1.
This is a literal combinatorial bridge to a classical recurrence, not a
numerical resemblance. In particular W_(2j+1) first appears in this tower
at size2j+1, occurs binom(k+j,2j) times at size2k+1, and gives exactly2k−1
prime induced cards of size2k−1 in W_(2k+1), for k≥2.

For example, formula(4), not an induced-subgraph enumeration, gives:

| Family order | Induced W3 subsets | Induced W5 subsets | All prime subsets of order≥3 |
|---:|---:|---:|---:|
|51|325|17,550|20,365,011,073|
|101|1,275|270,725|573,147,844,013,817,084,100|
|201|5,050|4,249,575|453,973,694,165,307,953,197,296,969,697,410,619,233,825|

These are exact inventories of **one family member** at each displayed order,
not A000568 or a census of all tournaments at that order.

The new prime **arity** appears at every odd order, but the family rule does
not become more mysterious: the same tiny transfer predicts every future
inventory. This answers “new pattern or old family at a new scale?” differently
at the quotient-arity and parametric-family levels.

## 4. Rigidity, even prime births, and the bounded-seed limit

For k≥2, W_(2k+1) has trivial automorphism group. Spine vertices2r and2r+1
both have score2k−1−2r, while the apex has score k. If k is even, the apex
is the unique even-score vertex, each remaining score class is an ordered
pair, and every vertex is fixed. If k is odd≥3, only the score-k class is
larger: it contains x,k−1,k. The top-score pair0,1 is still fixed pointwise,
and0 loses to x but beats the other two members of that triple. This fixes x,
and then every ordered score pair. W_3 instead has three automorphisms.

Prime tournaments exist at every order n≥5, not only odd orders. For odd
n use W_n. For even n=2k+2, k≥2, start with W_(2k+1) and add y beating exactly
x and v_0. By THM-4169 a prime parent has only uniform and clone attachments
as nonprime extensions. This is nonuniform. It cannot clone any spine vertex
v_i with i>0 because y beats v_0 and v_i loses to v_0. It cannot clone v_0
because y beats x and v_0 loses to x. It cannot clone x because x beats v_2
but y loses to v_2. Hence the extension is prime.

There is no prime tournament of order4: a strong order4 tournament contains
a triangle, and every triangle attachment is nonprime by the q=3 case of
THM-4169; a nonstrong tournament of order≥3 is nonprime. Thus the prime-arity
birth set is `{3}∪{n:n≥5}` (orders1,2 are conventionally trivial primes and
are excluded from this statement).

**Bounded-seed obstruction.** A substitution grammar whose only prime
quotients have order≤B cannot generate W_n for odd n>B, or the above even
prime examples for n>B. A nontrivial substituted block is a module of its
output; a prime output therefore cannot be a genuinely enlarged substitution
of a smaller quotient. This remains true when ordinal sums and arbitrarily
long transitive inflation blocks are allowed. It does **not** exclude a
finite *parametric* grammar: the one formula W_(2k+1) already includes infinitely
many prime arities. Nor does it rule out arbitrary-cut vertex extension,
which is a complete rule but retains exponentially many choices.

## 5. Chan–Pak boundary and reproducibility

A complete structural grammar can have an incomplete numerical image. The
inherited one-apex all-order exclusion H≠13 already establishes that boundary
for this class, even though its modules and prime profiles are now explicit.
No claim that all tournament structures, or all allowed path counts, are
covered is made here.

`verify.py` has no nonstandard dependencies and imports no repository
mathematics. It checks the module dictionary and module-number formula on
all1,023 binary words through length9 against literal subset modules. It
counts all induced prime subsets independently on all255 words through
length7, verifies the binomial/transfer recurrence through k=100, and checks
rigidity and the even prime extension literally through k=7. The normal and
`python3 -O` outputs match exactly. Checks use explicit exceptions, not
assertions removed by optimization. The computations test the implementation;
the all-order statements are established above.
