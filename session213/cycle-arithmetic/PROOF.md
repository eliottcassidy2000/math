# What the universal cycle compiler does to total path counts

**Status: PROVED elementary compiler identities and scoped obstructions, with
FINITE-EXACT implementation checks. This does not prove the Chan–Pak
Hamiltonian-path spectrum conjecture. No literature-priority claim.**

## 1. Inheritance and the question

Let `H(T)` count directed Hamiltonian vertex words and let `HC(T)` count
directed Hamiltonian cycles modulo rotation. The question is whether the
Session212 realization of every positive cycle count `c` can construct a
tournament with exactly `2c+1` paths. All odd values at least 23 would follow
from such a construction for every `c>=11`. The results below show exactly
why two direct approaches do not supply that construction.

The source packet is `session212/cycle-bridge/PROOF.md` and
`session212/cycle-bridge/ENDPOINT_BRIDGE.md`. Its endpoint core `T(h,m)` has
`s_1(T)=h`, and adding a sole-exit vertex gives `HC=h`. Both identities are
inherited, not new here. The one-apex path formula and its binary recurrence
are also inherited: THM-4876, `deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak`,
statement 6; THM-4875, `old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21`,
part D; and THM-4882, `signed-one-apex-responses-slope-and-parity-classification-session198-clock-merged-positive-count-barrier`.

The closure inequality `H>=n HC` is inherited from THM-4729,
`closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement`,
L2, and THM-4526,
`shaved-tournaments-path-plus-chord-and-the-largest-unavoidable-graphs`.
The no-`+2` insertion theorem is inherited from THM-4876, statement 2.
The source manifest pins exact local source bytes and repository paths.

The contribution is a direct application of those mechanisms to the
Session212 binary compiler: an explicit correction tournament, its exact
zero-correction classification, and the resulting operation boundaries.

Anchor: unrestricted path-spectrum construction. Niche: the residual
one-apex correction retained by an endpoint compiler. Wildcard: which
operations must change old arcs rather than merely add vertices.
Concept board: cycle closure / marked endpoints / binary code /
one-apex correction / induced monotonicity / arithmetic recursion.
The hostile example is `h=3`: a valid three-cycle compiler does not yield
seven paths. The least-used sidecar is the whole residual bitword, rather
than its marked endpoint count alone.

## 2. A self-contained one-apex calculation

For a word `b=(b_0,...,b_(m-1))`, let `A(b)` consist of a transitive chain
`y_0,...,y_(m-1)` and an apex `z`, with `y_i -> z` exactly when `b_i=1`.
The empty word gives a single vertex. Write

    q = sum b_i 2^i,
    N(b) = sum_{i<j, b_i=0, b_j=1} w(j-i),
    w(1)=1,       w(d)=2^(d-2) for d>=2.

**Lemma 1 (inherited formula, direct proof).** `H(A(b))=1+2N(b)`.

**Proof.** Delete the apex from a Hamiltonian path. The vertices before
and after it form increasing chain subsequences. If the two subsequences
do not cross in chain order, the path uses a prefix followed by a suffix.
There are

    [b_0=0] + [b_(m-1)=1] + #{i: b_i=1,b_(i+1)=0}
      = 1 + #{i: b_i=0,b_(i+1)=1}

such paths, for `m>=1`; the identity is the telescoping transition count.
If the last vertex before the apex is `y_j` and the first after it is
`y_i`, with `i<j`, we need `b_i=0,b_j=1`. Each intermediate chain vertex
chooses a side, yielding `2^(j-i-1)` paths. Summing gives twice `w(j-i)`:
for adjacent pairs the extra one comes from the prefix/suffix count, and
for longer pairs it comes directly from that power of two. The empty
word has one path and no pair. QED.

Appending zero bits at the high-index end changes neither `q` nor `N`.
Thus define `D(q)=H(A(b))` for any binary representation in this bit order.
In particular `D(0)=1`.

## 3. The exact correction in Session212's endpoint core

Take `m>=0`, `1<=h<=2^m`, and write `h-1=sum b_i 2^i`.
The Session212 core has a chain `x_0,...,x_(m+1)` and apex `z`, with

    z -> x_0,      x_1 -> z,
    x_(i+2) -> z iff b_i=1.

Its marked root is `x_1`. This is exactly the one-apex word `(0,1,b)`.

**Theorem 2 (compiler correction).** For every allowed `(h,m)`,

    s_(x_1)(T(h,m)) = h,
    H(T(h,m)) = 2h + D(h-1) = 2h+1+2N(b).                 (1)

**Proof.** The first equality is Session212's endpoint bijection: before
`z` a path starting at `x_1` either stops immediately (one choice), or
ends at a bit vertex `x_(i+2)` with `b_i=1` (exactly `2^i` choices).
After `z` the complementary increasing chain starts at `x_0`.

For the total count apply Lemma 1 to `(0,1,b)`. The pair `(x_0,x_1)`
contributes one to `N`. For each `b_i=1`, the pair from `x_0` to
`x_(i+2)` contributes `w(i+2)=2^i`. The vertex `x_1`, whose bit is one,
creates no further pair. Pairs between bit vertices contribute `N(b)`.
Consequently `N(0,1,b)=1+(h-1)+N(b)=h+N(b)`, proving (1). QED.

**Corollary 3 (exact equality locus).** The following are equivalent:

- `H(T(h,m))=2h+1`;
- `D(h-1)=1`;
- the bitword has no zero before a later one;
- `h-1=2^k-1` for some `0<=k<=m`;
- `h` is a power of two.

All summands of `N` are positive, so the equivalence follows at once.
The word must be a run of ones followed by a run of zero padding. The
case `k=0` is the all-zero word, including the empty word at `m=0`.
Thus the exact target values in this particular `2h+1` attempt are
`2^(k+1)+1`. Its successful `h` values have density zero. This says
nothing about the density of the family's full path-count image, or
about other realizers of `2h+1`.

For example `h=3` gives `D(2)=3` and `H=9`, not 7; `h=10` gives
`D(9)=5` and `H=25`, not 21. Even beyond those globally forbidden
values, `h=11` gives `D(10)=9` and `H=31`, not 23. Other tournaments
can and do realize 23; this obstruction is compiler-specific.

## 4. What happens to the binary growth law

**Lemma 4 (inherited recurrence, applied here).**

    D(2q+1)=D(q),
    D(2q)=D(q)+q+(q mod 2).                                (2)

Prepending a low-index one creates no new zero-before-one pair.
Prepending a zero creates
`sum_j b_j w(j+1)=(q+(q mod 2))/2` such weighted pairs. Lemma 1 gives (2).

In the `tau=1` cycle compiler the two cycle updates are `h -> 2h`
and `h -> 2h-1`. Let `K=H(T(h,m))` and `D=D(h-1)`. Their path updates are

    h' = 2h:     D'=D,                         K'=2K-D,
    h' = 2h-1:   D'=D+h-1+((h-1) mod 2),
                 K'=2K-D+h-3+((h-1) mod 2).                (3)

These are a joint-state recursion, not a scalar path-doubling rule.
For power-of-two `h`, the first branch has `D=1` and gives `K'=2K-1`,
which simply remains on the same sparse exact-equality ray.
Equation (3) does not rule out another family with useful doubling laws.

## 5. Induced-extension obstructions

**Lemma 5 (inherited closure bound).** Every tournament on `n>=3`
vertices satisfies `H(G)>=n HC(G)`.

Cut each Hamiltonian cycle at any of its `n` edges. The resulting paths
are distinct: closing a path recovers its unique cycle and cut edge.
No converse assertion that every path closes is used.

Hamiltonian-path counts are monotone under induced containment. Indeed,
given an old path and a new vertex, insert it before the first vertex if
it beats that vertex, after the last vertex if that vertex beats it, or
at a transition from an old vertex beating it to one it beats. At least
one position exists. Choosing one position deterministically makes an
injection, because deletion recovers the old path. Iterate for more new
vertices.

**Corollary 6.** If a construction retains an input cycle tournament `G`
as an induced subtournament and `HC(G)=c`, its output has at least `nc`
paths. Therefore it cannot output `2c+1` when `nc>2c+1`. In particular,
no additive-vertex gadget retaining a Session212 comb of order at least
five and positive cycle count can do this.

Deleting the comb's sole-exit vertex leaves `T(h,m)` and escapes that
specific closure bound. But any subsequent induced extension of that
core still has at least `2h+D(h-1)` paths. It therefore cannot achieve
`2h+1` for any non-power-of-two `h`. A successful repair must change
more of the old object: deleting further vertices, reorienting old arcs,
or constructing a different tournament. No claim that any such repair
is impossible is made.

**Corollary 7 (multi-vertex `+2` obstruction).** If `X` is an induced
subtournament of `Y` and `H(X)>=5`, then `H(Y)-H(X)` is not two.

This is an elementary extension of the inherited one-vertex theorem in
THM-4876(2): a one-vertex extension with increment two has old count
one or three. Add the new vertices one at a time. Every increment is
nonnegative by induced monotonicity and even by Rédei parity. If their
sum were two, one step would have increment two, while its old count
would be at least five. Contradiction. The one-vertex structural theorem
is a dependency, not newly reproved or inferred from finite tests here.

## 6. Verification and remaining problem

Run `python3 verify.py` and `python3 -O verify.py`. The verifier uses
explicit exceptions, a separately coded subset dynamic program, and
complete permutation checks on the stated small universe. It checks
the one-apex formula, endpoint correction, equality locus, both binary
updates, and actual sole-exit cycle counts. It also checks closure and
monotonicity on every labeled tournament of orders at most five. These
finite controls verify the implementation, not the all-order proofs.

The universal Hamiltonian-cycle spectrum has not closed the unrestricted
Hamiltonian-path spectrum. This packet neither proves nor refutes that
every odd integer at least 23 is attained. Its output is an explicit
correction coordinate and precise failure boundaries for the most direct
cycle-to-path arithmetic proposals.
