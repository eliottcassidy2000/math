# Paid adaptive controllers with lossless recursive receipts

2026-10-04. **PROVED:** every nonempty legal episode of the controller below
strictly pays its immutable source, including arbitrary pattern changes and
early stops; no infinite legal episode exists. **PROVED:** an exact affine
interface test and a recursively composable certificate grammar; the incoming
four-slot family funds one credit with sharp contraction 15/16.
**FINITE-EXACT:** the stated experiments. **OPEN:** guard coverage for arbitrary
positive integers and a supplied root proof for arbitrary exit states.

Artifacts: [program](../../04-computation/experiments/adaptive_credit_potential_20261004.py)
and [saved output](adaptive_credit_potential_20261004.out).

## 1. Inheritance and the actual payment target

The closest proved mechanism is the
[155 mod 2048 common-future dependency](collatz_checkpoint_reroute_20261004.md),
with its [exact repeat fuel](collatz_recursive_dependency_kernel_20261004.md).
The [paid portrait controllers, sections 6–7](collatz_paid_portrait_controllers_20261004.md)
already pay across a specified anchor switch while preserving the original
source, and show that a refill cofactor can produce arbitrarily much new
binary precision. That refill observation alone does not pay a size debt.
The [root rank](collatz_branch_toll_rank_20261004.md) is inherited.

Canonical hostile: an actual expanding pattern after a shrinking dependency
can overspend the saved decrease. The corrected near miss is resetting the
comparison source at the pattern switch. The useful sidecar is an explicit
credit account tied to one immutable source. The live board is **native
guard / immutable source / ordered carry / credit / recursive interface /
terminal root obligation**.

The session uses the existing META-PATTERNS card **Audit and close sections
under their next native operation**: arithmetic carry and geometric port
frames pass composition tests; slope-only, current-only, and uncompensated
chart quotients fail explicit hostiles.

Anchor: prove payment despite adaptive switching. Niche: a lossless ordered
ternary-tree receipt. Wildcard: the same tree has an exact geometric carrier
in [adaptive triangle refinement](small_port_refinement_20261004.md).

Write U(x)=(3x+1)/2^(v_2(3x+1)) on positive odd integers. The two operations are

    H(x)=(729x+669)/1024,   native guard x=155 mod 2048,
    G(x)=(9x+5)/8,         native guard x=11 mod 16.

G is the actual valuation word (1,2). H is a common-future dependency:
with v=(1,2,1,1,1,2), one has F_v(x)=4H(x)+1 and
U(F_v(x))=U(H(x)). Thus each operation preserves the existence of a route
to 1, but H is not itself a forward Collatz edge or word. This type
distinction is retained by the compiler.

Every legal H has x>=155 and 111<=H(x)<x. Every legal G strictly increases
x. Both affine maps have positive coefficients.

## 2. The credit rule and an exact potential

An episode starts with original source N, current value x=N, and credit k=0.
An H operation adds two credits. A G operation spends one credit, and is
allowed only when k>=1. Both operations also require their native guards.
The controller may choose adaptively between available operations and may
stop at any time. Neither source N nor the accounting convention changes at
a pattern switch.
Payment is asserted at completed controller-action boundaries. Actual Collatz
paths used inside an action may rise above the original source.

Define the positive rational potential

    E(x,k)=(x+5)(9/8)^k.

**Neutral debit.** Since G(x)+5=(9/8)(x+5), a G operation preserves E exactly.
Its actual numerical growth is paid by the consumed credit.

**Strict credit minting.** An H operation has

    E(H(x),k+2)/E(x,k)
      =(81/64)[729/1024+(67/32)/(x+5)]
      <= 2349/2560 =: kappa < 1.

The ratio decreases with x; the least legal input is 155, and equality is
attained there because H(155)=111. This is a uniform, sharp bound for the
individual H transition with the stated credit award.

**Adaptive payment theorem.** Let a nonempty finite H/G word be executed
from a positive odd N, with every native guard satisfied and no negative
credit prefix. If it has h letters H and g letters G, then

    k=2h-g >= 0,
    116 <= E(x,k) <= (N+5) kappa^h < N+5,
    x<N.

Proof: the first letter must be H. Every subsequent state is at least 111:
H outputs satisfy this bound and G grows its input. The transition identities
give the upper bound by induction. Nonnegative credit gives x+5<=E(x,k),
proving the strict endpoint inequality for every nonempty prefix. The lower
bound is x+5>=116. No choice rule, final valuation, eventual root certificate,
or bound on the number of switches is assumed.

This also pays the inherited rank. For positive odd z>1 write
z-1=2^K t with t odd, and set R(z)=(3^K t^2,K), lexicographically;
R(1)=(0,0). Here N=155 mod 2048 has K=1. For every positive odd x<N,

    3^K(x) t(x)^2 <= 3(x-1)^2/4 < 3(N-1)^2/4.

Hence R(x)<R(N), independently of the valuation at the exit.

**Finite episode theorem.** Every legal episode contains only finitely many
operations. After h H operations,

    116 * 2560^h <= (N+5) * 2349^h,       g<=2h.

The first inequality bounds h by an explicitly computable integer, and the
second bounds the remaining operations. The program computes that integer
with exact multiplication, without floating logarithms. An infinite sequence
would require infinitely many H operations, contradicting the lower bound
on E.

Termination here means that this guarded episode stops. It can stop at an
unrooted proof obligation because no next operation is admitted. A new
available pattern or a root certificate still has to be supplied.

## 3. Exhaustion, nested reuse, and the sharp third-credit hostile

The theorem includes all interleavings, not only H^m followed by G^(2m).
The [grouped-switch package](adaptive_switch_payment_20261004.md) supplies
explicit all-height source families where the H repeat budget and the G
repeat budget are both exhausted. It also constructs completed root families
at every depth. Those are additional terminal witnesses, not a premise of
the payment theorem.

At zero credits, the controller can change from G to another legal H without
resetting N. It then earns two more credits under the same potential. A stop
with unused credits is also already paid. A refill of binary precision never
creates spendable credits by itself.

Three credits per H would fail for this size-payment account. Direct
composition gives

    G^3(H(x))=(531441x+1598741)/524288 > x

for every positive x. The word HGGG is genuinely legal on the source cylinder
110747 mod 1048576, so this is an arithmetic hostile, not just a formal map.
Its least instance is 110747 -> 112261. That particular endpoint still lowers
the inherited graph rank: size growth alone does not prove rank failure.
The next source in the same cylinder gives a genuine rank hostile:

    1159323 -> 1175143,
    R(source)=(1008020624763,1),
    R(endpoint)=(1035719040123,1).

Both source and endpoint are 3 mod 4. Thus two is the greatest universal
integer credit award per H for this G-denominated size account; some
three-block routes can still pay a different rank.

**Recursive replacement principle.** A separately proved, legal subroutine
that sends a positive odd x to y<x may be inserted as a zero-credit action:
E(y,k)<E(x,k). It can itself be represented by a previous paid receipt. This
preserves payment against the same original N. If such subroutines are atomic
and terminating, the extended controller still cannot execute infinitely:
E>=6 bounds the number of H operations, the credit count bounds G operations,
and between those finitely many operations only strictly decreasing positive
integer values remain. Rank-only decreases that enlarge x do not receive
this interface automatically.

## 4. An exact test for importing another pattern

Let a candidate pattern have legal affine action F(x)=a*x+b with a>0,
positive output on its domain, a justified source lower bound x>=M, and a
proposed integer credit change d. Keep the same potential. Its ratio is

    E(F(x),k+d)/E(x,k)
      =(9/8)^d [a+(b+5-5a)/(x+5)].

This is monotone or constant on x>=M. Its exact supremum is therefore

    rho=(9/8)^d * max{a, (F(M)+5)/(M+5)}.

The inequality rho<=1 is necessary and sufficient for nonincrease on the
whole real half-line. It is a sufficient test on any arithmetic subdomain
there; when that subdomain contains M and is unbounded, the same supremum
makes it necessary as well. A finite guard or a nonattained loose lower bound
does not inherit that converse. Actual path/dependency legality and
nonnegative resulting credits are separate premises.

Thus a controller can import another pattern after its current budget is
exhausted by proving this interface inequality. It retains the old source
and accumulated budget. A strict inequality at some earlier transition
continues to imply strict payment after neutral interfaces. Payment and
termination remain different obligations: an arbitrary library of neutral
zero-cost loops would not inherit the finite episode theorem.

For a mid-episode entry, the starting register must already satisfy
E(current,k)<=N+5, with either strict entry or a later strict transition for
strict payment. Resetting N to a larger current value proves a different
statement. The [lossless carrier package](lossless_payment_carriers_20261004.md)
retains N explicitly and audits this boundary.

### 4.1. Concurrent four-slot patterns can themselves fund new credits

The incoming [four-slot payment theorem, section 6](collatz_four_slot_compression_20261004.md)
proves strict descent for every actual positive odd word starting with 1
and permuting the multiset {1,1,2,a}, a>=3. Retaining its ordered carry gives
a stronger interface here: **each such pattern can earn one credit**.

Indeed its affine map is (81x+B)/Q, where Q=2^(a+4) and
B<=45+14*2^a. Its native source is 3 mod 4 and cannot be 3, whose first
four valuations are (1,4,2,2). Hence x>=7. The half-line bound gives

    E(F(x),k+1)/E(x,k)
      <= (9/8) max{81/Q, 47/96+51/Q}
       = (9/8)(47/96+51/Q)
      <= 1023/1024 < 1.

The equality between the two expressions follows from Q>=128. Retaining
the native source sharpens this to the exact uniform factor **15/16**.
At source 7 the only actual four-slot word is (1,1,2,3), ending at 5;
the ratio is (9/8)(10/12)=15/16. Every other source is at least 11, where
the same carry bound gives

    E(F(x),k+1)/E(x,k)
      <= (9/8)[47/128+117/2^(a+5)]
      <= 1899/2048 < 15/16.

The slope limit is smaller than the displayed endpoint bound, as before.
This proves sharpness over every a, ordering, and legal source. One is also
the largest uniform integer award that keeps E nonincreasing at every P step:
the native word (1,1,2,3) sends 7 to 5, and awarding two credits would give
E(5,2)/E(7,0)=135/128>1.
This is a failed potential interface, not a size-payment counterexample:
5 remains below 7, and the example alone supplies no legal overspending continuation.

**Extended adaptive theorem.** Allow H to earn two credits, any such four-slot
pattern P to earn one, and G to spend one. Keep every native guard and the
immutable source. If h and p count the two funding types, then every nonempty
legal episode obeys

    6 <= E(x,k) <= (N+5)(15/16)^(h+p) < N+5,
    g <= 2h+p.

Here the H factor 2349/2560 is also below 15/16. Thus the controller pays
size and the inherited original-source rank, and cannot continue
infinitely. The first funding action guarantees N=3 mod 4. First-hit root
certificates stop at 1; a formal four-slot word with a remaining valuation 2
is truncated there rather than padding the root. Including arbitrary
terminating zero-credit strict descents after this funded entry preserves
finiteness as in section 3. This bounds executed legal actions, not the runtime
of an unspecified search that may keep looking for an unavailable pattern.

Two concrete imports are (1,2,1,3) on 187 mod 256, with map
(81x+85)/128, and (1,1,2,3) on 7 mod 256, with map (81x+73)/128.
At zero credit their exact interface bounds are 31/48 and 5/6 respectively.
For the H/P/G library alone, its entry domain is precisely the union of the H guard and the
four-slot bank's guards: the first action must fund credit, and each such
funding action works alone. This remains a restricted domain.

The balanced receipt grammar now permits binary as well as ternary nodes:

    S = empty | H S G S G S | P S G S.

P retains its actual valuation word and source guard as a label. Its one
credit is the net pending-slot change at a binary node. The H/G affine
decoder below is proved only for that two-letter alphabet; enlarging the
library does not automatically retain injectivity. Preserve the typed
expression tree until a decoder for the enlarged library is established.
Geometric binary bisections can represent these nodes, but arbitrary local
bisections need boundary matching or neighbor refinement before the resulting
partition is claimed to be a conforming triangulation. A mixed tree with h
ternary and p binary internal nodes has 2h+p+1 leaves; each leaf's normalized
area is 3^(-h_path)2^(-p_path). Keep this typed tree if a geometric conformity
repair introduces additional cells: those cells are not new arithmetic steps.
The triangle package gives an exact collision of two distinct mixed chart
histories of length five. Thus even a full final matrix with marked ports
does not automatically recover the mixed history; the retained tree is
essential unless an additional decoding theorem is supplied.

### 4.2. A source-sensitive controller can earn a second credit

The native guard supports a stronger optional award: give P one credit at
x=7 and two credits at every other legal source. Sources 11 and 15 have
actual prefixes (1,2,3,4) and (1,1,1,5), so neither is in this bank.
Consequently every legal source above 7 is at least 19. At 19 the unique
P word is (1,3,1,2), ending at 13, and its two-credit ratio is 243/256.
For all remaining sources x>=23, the carry bound gives

    E(F(x),k+2)/E(x,k)
      <= (81/64)[47/224+477/(7*2^(a+4))]
      <= 7695/8192 < 243/256.

Together with the one-credit source 7, this proves the sharp uniform
contraction 243/256 for the adaptive award. H contracts more strongly.
The extended payment and finite-episode theorem therefore holds with
factor 243/256 and g<=2(h+p). This permits more growth actions while
keeping exactly the same arithmetic entry guards.

The one-credit grammar in section 4.1 remains a useful uniform subclass.
For source-sensitive receipts, each P node must retain its entry condition
and its actual award: a fixed syntactic letter count no longer certifies
the entire budget. The state-dependent credit check is part of the receipt.

## 5. Exact source cylinders and actual proof receipts

Compose the triples (P,Q,B) representing (P*x+B)/Q in chronological order.
For H and G, P and B are odd and Q is a power of 2. For every nonempty word,
the exact source cylinder is

    n=(Q-B)P^(-1) mod 2Q.

Necessity is final oddness. For sufficiency, reduce Pn+B=Q mod 2Q at the
first dyadic cut. The remaining odd multiplier forces the first division to
be integral. If a nonempty tail remains, its odd carry and even denominator
force that first endpoint to be odd; then recurse. The last endpoint is odd
by the original congruence. Hence every primitive native guard holds.
Positive sources and positive coefficients give positive intermediate values.
Every word has infinitely many positive sources on its cylinder, including
words that fail the separate credit condition.

The full affine carrier also retains the order: a
[3-adic carry decoder](lossless_payment_carriers_20261004.md) reconstructs the
H/G word. This is a new application of the earlier ordered-carry mechanism;
H/G are not silently reclassified as ordinary Collatz valuation letters.

To produce an actual common-future receipt, let y be the controller exit,
a=v_2(3y+1), and start with the one-step word (a) from y to U(y). Traverse the
controller word backward. A G prepends (1,2). An H replaces a current suffix
(b,tail) by (v,b+2,tail), with v as in section 1. The resulting source word
and the terminal word (a) have exactly the same endpoint U(y). Its lengths are

    source odd depth = 6h+2g+1,
    source halving cost = 10h+3g+a.

No prefix pads a root: all controller states exceed 1, G's intermediate state
grows, and every nonempty prefix of v grows its guarded source. If the common
future equals 1, it is the first hit. Otherwise the receipt is a paid
dependency on y and still needs a root proof for y.

## 6. Lossless ternary receipts and their geometric carrier

The credit increment is +2 for H and -1 for G. A balanced word has no negative
prefix and total zero. Such words have the unique recursive grammar

    S = empty | H S G S G S.

This is exactly an ordered full ternary tree. In preorder an internal node
creates three pending child slots, replacing one slot and giving net +2;
a leaf consumes one. The balanced word omits the final leaf, so credit equals
pending slots minus one. This identifies the combinatorial roles without
calling a G operation a terminal root proof.

With m internal nodes the word has m H and 2m G, and the number of trees is

    binomial(3m,m)/(2m+1).

For completeness, append the omitted leaf. Among all words with m increments
+2 and 2m+1 increments -1, exactly one of each 3m+1 cyclic shifts has
nonnegative proper partial sums and final sum -1. The usual last-minimum
argument proves that uniqueness; total sum -1 prevents a smaller period.
Dividing binomial(3m+1,m) by 3m+1 gives the displayed count.

Every such word is legal on its own nonempty source cylinder. Its composite
has P=3^(10m), Q=2^(16m), and modulus 2^(16m+1), regardless of tree shape.
Its carry records that shape. The cylinders of different trees need not be
disjoint; this is a receipt language, not a partition of the integers.

A prefix, including an unbalanced funded one, has a small budget summary
(delta,minimum), with delta=2h-g and minimum the least prefix balance including
zero. Concatenation obeys

    (d_u,m_u) * (d_v,m_v)
      = (d_u+d_v, min(m_u,d_u+m_v)).

Starting with credit k is valid precisely when k+minimum>=0. This summary
retains that predicate, not the whole word: HGG and GGH both have total zero
but different minima. The full affine carrier or the ordered expression
retains the history. A computed summary is a cache, not a self-authenticating
certificate supplied without its provenance.

The [triangle package](small_port_refinement_20261004.md) maps the same ordered
ternary tree to centroid subdivision. Exact chart matrices and parent/port
labels recover every node and child address. This gives a compositional
geometric store for the receipt. Native arithmetic guards, credit annotations,
and terminal root proofs remain labels; a triangulation alone supplies none
of them. Six barycentric flag charts and seven first-refinement vertices are
related geometric objects with different retained information.

The small sizes have the following exact roles:

| Size | Specified structure and retained information |
|---|---|
| 3 | Homogeneous monitor (current, original, 1); separately, the three child slots forced by the two-credit grammar |
| 5 | The affine-closed observer (x^2,N^2,x,N,1), minimal when x,N,1 and x^2-N^2 are all required |
| 6 | The full quadratic observer adds xN; separately, six ordered triangle flags retain relative port order |
| 7 | The three vertex, three edge, and one face barycenters; their ranks prevent extra bare-graph symmetries |

The triangle's seven coarse support/equality strata are another quotient,
not those seven vertices. All seven first coexist on an integer triangle
grid at degree 12. The degree-six grid lacks the two-equal-large interior
type, a useful hostile to inferring all allowed types from one small sample.
The six flag cells refine the three centroid-fan cells by a genuine two-way
edge bisection. These are explicit maps with declared information losses.

## 7. Reproduction and the remaining coverage obligation

Run:

    python -X utf8 -B 04-computation/experiments/adaptive_credit_potential_20261004.py
    python -X utf8 -B -O 04-computation/experiments/adaptive_credit_potential_20261004.py

The script uses exact integers/rationals and explicit exceptions. Its universe
contains all 510 nonempty H/G words of length at most 8 at source lifts 0,1,7:
1,530 independent native-guard/cylinder checks. The 645 funded instances each
check every prefix's potential, size/rank payment, finite bound, affine
endpoint, interface bound, and compiled actual common-future paths. It also
checks all 345 ternary trees through five internal nodes, 225 summary joins,
24 compressed repetitions, and eight malformed/type/guard/budget hostiles.
Three-credit examples separately test size failure and actual rank failure.
The imported four-slot bank is independently checked for a=3,...,14,31,64,
all six orderings, and lifts 0,1,7: 252 native source cases, with first-hit
truncation and a hostile to awarding two credits throughout the bank. The
uniform theorem over every a is the explicit inequality in section 4.1.
Mixed H/G/(1,2,1,3)/(1,1,2,3) programs through length four and three deeper
switching examples give 741 exact source replays at lifts 0,1,7, including
its original-source payment and exact native guards.
The same universe is checked again with the source-sensitive credit rule;
the saved output records the enlarged admitted count. Both sharp witnesses
7->5 and 19->13 are replayed as actual first-hit-safe words.

The switch theorem proves a uniform adaptive controller on its declared
domain. The union of sources admitting a nonempty funded program is exactly
155 mod 2048: every such program starts with H, and H alone works throughout
that cylinder. Section 4.1 adds the independently guarded four-slot library.
Neither statement proves that
every possible exit has an available next paid controller. For a library of
strictly paid episodes, successive original-source ranks cannot decrease
forever; proving that their guards cover every nonroot source is the missing
global step. Lossless storage preserves that obligation rather than erasing it.
