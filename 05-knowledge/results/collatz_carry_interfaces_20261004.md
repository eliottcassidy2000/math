# Carry interfaces connect paid controllers and Fourier phases

2026-10-04 local / 2026-10-05 UTC. **PROVED:** lossless translation coding
for the five-operation library below; exact composition of its native
arithmetic progressions and their finite admissibility automaton;
the resulting golden-ratio sublanguage and reciprocal trace law;
the source-address Fourier identity and its uniform
exponentially small error; the binary carry-square law and its state-count
boundary; marked-centroid coding on the two pure refinement alphabets.
**FINITE-EXACT:** the independent controls in the accompanying program.
**OPEN:** H1, universal paid coverage, and universal Collatz completion.
No novelty or priority claim is made for the elementary coding mechanisms.

Artifacts: [program](../../04-computation/experiments/collatz_carry_interfaces_20261004.py),
[exact output](collatz_carry_interfaces_20261004.out), and
[structured witnesses](collatz_carry_interfaces_20261004.json).

The reusable move is to identify a small compositional object, specify what
its quotient forgets, and recover only the coordinates needed by the next
operation. Here that move produces an exact bridge: the affine carry that
records a controller's order also determines its binary guard and its
Fourier phase. Keeping the native guard remains essential after cancellation.

## 1. Inheritance and incoming work

The closest mechanism is the ordered carry of
[four-slot compression](collatz_four_slot_compression_20261004.md).
The canonical hostile is equal clocks with unequal carries. The corrected
near miss is assuming that endpoint integrality reconstructs a guard after
algebraic cancellation. The least-used sidecar is a marked input/output
progression, including its actual period rather than just the reduced
denominator. Our board is **carry / native guard / marked seed / phase /
immutable source / termination**.

Anchor: lossless composition of paid Collatz controllers. Niche: the
source-address meaning of the concurrent Fourier calculation. Wildcard:
the same disjoint-image coding lemma in triangle refinement.

This session read the mathematical packages in the incoming batches,
including their corrections, rather than treating commit titles as evidence:

- Through `0392961fe6`: [adaptive credit](adaptive_credit_potential_20261004.md),
  [exhausted switches](adaptive_switch_payment_20261004.md),
  [lossless carriers](lossless_payment_carriers_20261004.md),
  [triangle ports](small_port_refinement_20261004.md), the Fourier update,
  and the coalescence audit. The carrier package supplied the H/G decoder;
  the triangle package supplied pure-alphabet injectivity and a mixed collision.
  The Fourier update supplied the exact cross-level root law and calibrated
  controls. Its retracted root-order rate law is not a dependency here.
- Through `37bf82ebfd`: [guard budgets](paid_guard_budget_20261005.md),
  [residue cover tree](paid_guard_cover_tree_20261005.md),
  [thresholded exits](paid_guard_27_boundary_20261005.md),
  [the 5/8/9 transfer](five_eight_nine_transfer_20261005.md), and the Fourier
  ensemble and Pascal-tower additions. The new L operation forced us to
  replace the initially sufficient denominator-only guard by a native
  progression. Its distinct ternary tag then extended the decoder.
- Through `9ba31b1b3b`: the exact full/odd-Haar moment proofs, integer-start
  ensembles, and q=5 heavy-tail caution in the Fourier package. These led to
  [the two-clock source-address theorem](source_address_two_clocks_20261004.md):
  a complete binary character partition and an exact adjacent-level energy.

The incoming coverage improvements belong to those packages; this note
does not add their densities again. We use the existing META-PATTERNS
cards on closing sections under native operations and retaining coordinates
at the first collision. The connection table in section 7 records the
actual transfer and its failure boundary.

## 2. A small disjoint-image coding lemma

Let injective maps f_i:X->X have pairwise disjoint images, and let a marked
seed o lie outside their union. Then every finite word has a distinct
value f_(i_k)...f_(i_1)(o). Its last letter is identified by the image
containing the value; applying that inverse removes the letter. The empty
word is recognized by o. Induction proves injectivity.

This lemma does not say that every element of X has a finite address.
Recognition needs either a decreasing cost or a supplied depth bound.
Disjointness must hold on the domain actually used, not merely on a sample.

For affine maps over the 3-adic integers, it suffices that each slope is
divisible by 9 and that the translations occupy distinct nonzero classes
modulo 9. Their images then lie in disjoint residue balls; o=0 works.
All arithmetic below can instead be read in exact rationals with denominator
a power of 2, where reduction modulo 9 is numerator times inverse denominator.

### Five operations whose translation records the whole word

Read words chronologically. A and B below are names of operations, not the
carry B_w used later.

| Letter | Affine action | Native source guard | Uniform credit | Translation mod 9 |
|---|---|---|---|---|
| H | (729x+669)/1024 | 155 mod 2048 | +2 | 3 |
| G | (9x+5)/8 | 11 mod 16 | -1 | 4 |
| A | (81x+85)/128 | 187 mod 256 | +1 | 2 |
| B | (81x+73)/128 | 7 mod 256 | +1 | 5 |
| L | (9x-3)/16 | 219 mod 256 | +4 | 6 |

G, A, B are the actual valuation words 12, 1213, 1123. H and L are the
proved common-future dependencies from the incoming notes. They are not
silently relabelled as forward Collatz words. The stronger source-sensitive
credit awards for A/B are optional and not used in this fixed grammar.

**Translation-code theorem.** The exact rational translation F_w(0) alone
determines every finite word in this five-letter alphabet, including words
whose native guard is empty. Consequently it determines the slope as well,
once decoded. This conclusion is specific to the declared alphabet.

Proof: all five slopes are in 9 Z_3, and the five translation residues in
the table are distinct and nonzero. Apply the lemma. Concretely, for the
last letter (p x+b)/q, recover the previous translation as

    t_old = (q*t-b)/p.

There is also finite recognition. In a generated nonempty word the
unreduced carry numerator is odd, including when it is negative, and the
denominator is the product of the letter denominators. Thus its reduced
dyadic denominator records the total cost. Strip the named denominator,
check divisibility of the remaining numerator by p, and continue. The cost
strictly decreases. Reject invalid residues, insufficient cost, nonintegral
carry predecessors, or a non-dyadic input. Re-encoding authenticates a
supplied full affine carrier; its slope is not accepted unchecked.

The representation is not a bounded-bit memory claim. Numerator and
denominator grow with the word; an explicit decoder takes time proportional
to its output length. Chosen sharing of repeated expressions is extra data.

Clocks alone fail especially sharply here:

    H, GA, GB, AG, BG all have slope 729/1024,
    but carries 669, 1085, 989, 1405, 1297 respectively.

H earns two credits, while each of the other four words has net credit zero.
The same multiplier therefore loses both order and the budget interpretation.

If K=(81x+53)/128 is added as a separate primitive, the alphabet is no longer
free: the incoming guarded identity is K=LG. Normalize K to LG, or retain an
expression modulo this proved relation. A decoder for one alphabet is not
automatically a decoder after adding a redundant generator.

## 3. Native progressions survive cancellation and compose exactly

Store an operation as a pair of aligned progressions

    [c,m ; d,l]:  c+m*t -> d+l*t,    t integer,

with 0<c<m, 0<d<l, and the operation's proof type. Positive sources then
correspond to t>=0. Its affine formula is d+(l/m)(x-c), but its native guard
is x=c mod m. The five rows above are

    H [155,2048 ; 111,1458],    G [11,16 ; 13,18],
    A [187,256 ; 119,162],      B [7,256 ; 5,162],
    L [219,256 ; 123,144].

The L formula has odd outputs already on 27 mod 32. That enlarged cylinder
is not its native dependency guard. The literal witness 27->15 fails the
required child prefix 12. The period 256 retains the three discarded bits.

For u=[c,m;d,l] and v=[e,n;f,k], put g=gcd(l,n). If g does not divide e-d,
the composite guard is empty. Otherwise choose the least solution

    t0 = ((e-d)/g)*(l/g)^(-1) mod(n/g),
    s0 = (d+l*t0-e)/n.

Then their chronological composite is exactly

    [c+m*t0, m*n/g ; f+k*s0, l*k/g].                 (1)

Proof: the interface requires d+l*t=e mod n. Its solutions are
t=t0+(n/g)j, and substitution gives (1). The displayed new representatives
remain in the interiors of their periods: 0<=s0<l/g. Associativity follows
from the same three compatible input/output equations, including empty
intersections. This is an arithmetic fibre product, not a slope product.

The old saturated carrier is the special case m=2Q, l=2P for
F=(Px+B_w)/Q. The new carrier permits a larger native period when cancellation
has reduced Q. It also stores the relative source/target placement needed
to test payment, since

    source - target = (c-d)+(m-l)t.                  (2)

Examples expose the semantic boundary. LG has exactly the K guard and map,
[219,256;139,162]. LB has empty guard: L outputs are 11 mod 16, whereas B
requires 7 mod 16. The formal balanced program LBGGGGG therefore has a
perfectly decodable translation and a nonnegative credit history, but no
native source. A tree's syntax is not its arithmetic realization.

On legal funded words, the incoming potential E(x,k)=(x+5)(9/8)^k applies.
G is neutral; H, A/B and L contract by at most 2349/2560, 15/16 and
6561/7168 respectively. Thus every funding step contracts by at most 15/16,
and every nonempty funded episode pays its immutable source. This is an
application of the incoming payment proof, not a new coverage theorem.
The decoded word permits regeneration of every native guard, budget prefix,
and typed local dependency. A terminal root proof remains a separate input.

With h ternary H nodes, p binary A/B nodes, and l five-child L nodes, the
formal balanced language has 2h+p+4l+1 leaves and

    2^p (3h+2p+5l)! / [h! p! l! (2h+p+4l+1)!]

words. The usual last-minimum cyclic-shift argument gives this count.
It counts receipt shapes; some have empty native guards, as above.

### The exact guard language recovers a reciprocal quadratic recursion

**Native-language theorem.** A finite word in H,G,A,B,L has a nonempty
native source guard if and only if it contains no adjacent pair LB.
The claim concerns the existence of sources, not legality at every source.

Proof: keep only the dyadic part of the target progression in (1).
If a legal prefix ends in L, its target period has valuation four and
its target is 11 modulo 16. Every other legal prefix, including the empty
prefix, has target period of valuation one and an odd target. All native
input periods have valuation at least four; the L input period has
valuation eight. The gcd test in (1) therefore gives these exact rules:

- From the ordinary state, every letter is possible. A following L leaves
  target valuation four and residue 11 modulo 16; every other letter
  leaves valuation one.
- From the L state, the next native input must be 11 modulo 16. H,G,A,L
  satisfy this and B does not. A legal following letter resets the state
  according to the same rule as above.

These statements follow by the formula
v2(l_new)=v2(l_old)+v2(l_letter)-min(v2(l_old),v2(m_letter)).
No higher bits obstruct existence because all input moduli are powers
of two. A nonempty progression supplies infinitely many positive sources.
Induction proves both directions.

There are two live states and one rejecting state. This deterministic
automaton is minimal: a B suffix distinguishes the two live states, and
the empty suffix distinguishes each from rejection. Its live counting
matrix, with rows and columns ordered ordinary,L, is

    M5 = [4 1]
         [3 1],       det M5=1, trace M5=5.           (2a)

Thus the number a_n of native-realizable words of length n satisfies
a_0=1, a_1=5, a_n=5a_(n-1)-a_(n-2):

    1, 5, 24, 115, 551, 2640, 12649, ... .

This makes the reciprocal quadratic structure native to the guard rules.
The eigenvalues are lambda=(5+sqrt(21))/2 and lambda^(-1), and

    trace(M5^(2n))=trace(M5^n)^2-2.                   (2b)

The correction has a boundary interpretation. At length two, linear
words exclude LB and number 24. Cyclic words also exclude BL, whose closing
boundary is LB, and number 23=5^2-2. Matrix trace counts cyclic words;
the sum from the initial ordinary state counts linear words. These
observables must not be interchanged.
Nor are cyclic label words arithmetic cycles. For example L repeated forever
passes this finite-state language, but L's fixed point is the rational
-3/7, not a positive integer. The native 2-adic guard contains that fixed
point; all finite L prefixes have positive realizations, without any one
positive integer realizing the infinite repetition.

There is an exact golden-ratio sublibrary. Restrict to G,B,L, retaining
the same native rules. The two-state matrix becomes

    M3 = [2 1] = [1 1]^2,
         [1 1]   [1 0]

with eigenvalues phi^2 and phi^(-2). Its native-word counts are
F_(2n+2)=1,3,8,21,55,..., whereas its two-letter cyclic count is 7.
Here the 7/8 difference and the golden ratio come from the actual forbidden
interface, not a matching decimal. These are counts of nonempty guards;
they are neither densities of their overlapping source cylinders nor
counts of grounded root certificates.
The incoming unit-clock transfer obstruction remains intact: M3 acts on
guard-admissibility states, while multiplication by 9/8 acts on arithmetic
values. No nonconstant quotient from the golden unit clock to that
contracting ternary action has been constructed.

More generally, weighting the five letters by h,g,a,b,l gives the matrix
[[h+g+a+b,l],[h+g+a,l]], with trace h+g+a+b+l and determinant b*l.
The trace-square defect is exactly twice the weight of the forbidden pair.
The determinant forgets which orientation is forbidden; the named
transition interface restores it. This is the same two-mixed-term
mechanism as the earlier reciprocal-product expansion, with a new
independently specified carrier.

The guard automaton can also be joined to the credit grammar. In its ordered
tree, LB occurs precisely when an L node's first child is B-rooted. Let T
count realizable balanced receipts, let z count executed letters, and let N
count those receipts whose root is not B. The nonnegative recursive system is

    T=1+z^3*T^3+2z^2*T^2+z^5*N*T^4,
    N=1+z^3*T^3+ z^2*T^2+z^5*N*T^4.

Equivalently T=1+z^3*T^3+2z^2*T^2+z^5*T^5-z^7*T^6. This counts exactly the
balanced programs with nonempty guards. It still does not test a particular
source, whose unbounded address is retained by the full progression.

## 4. The Fourier phase is the other end of the same progression

This section concerns ordinary valuation words, not H/L dependency letters.
Fix an odd multiplier q>=3 and a word w=(a1,...,am), all ai>=1. Put

    A=sum ai, Q=2^A, P=q^m,
    B_w=sum_(i=1..m) q^(m-i) 2^(a1+...+a_(i-1)).

Its affine action is (P x+B_w)/Q. Its exact odd source guard is

    c=(Q-B_w)*P^(-1) mod 2Q,    0<c<2Q,
    d=(P*c+B_w)/Q.                                   (3)

Final oddness reconstructs all preceding exact valuations. This follows
inductively by reducing at the first dyadic denominator and using the odd
carry of the remaining tail. The same proof holds on negative odd inputs.
Positive guarded sources stay positive, and negative guarded sources stay
negative because q>=3. Apply this at c and c-2Q to obtain

    0<d<2P,  d odd,  gcd(d,q)=1.

Hence the complete native progression is [c,2Q;d,2P], and

    B_w=Q*d-P*c,
    d/P = c/Q + B_w/(P*Q).                           (4)

Let e(z)=exp(2*pi*i*z). The Fourier summand for the q-adic Syracuse law is

    eta_w(u)=e(u*(B_w*Q^(-1) mod P)/P)=e(u*d/P)
            =e(u*c/Q)*e(u*B_w/(P*Q)).                 (5)

This is the promised bridge from the exact binary source address to the
terminal Fourier phase. At frequency u=1, or more generally gcd(u,q)=1,
the phase retains the canonical odd target d: among j and j+P exactly
one is odd. Nonunit frequencies can lose more. The phase does not retain c, Q, the word, or
the source's height parameter. For example 17 and 71 have the same length
and cost, and both target ports are 17, but their source ports are 483 and
469. Here 17 and 71 denote valuation words, not integer sources.

### An exponentially accurate replacement of the Fourier observable

The correction in (4) has the uniform bound

    0 < B_w/(P*Q)
      = sum_(i=1..m) q^(-i) 2^(-(ai+...+am))
      <= (2^(-m)-q^(-m))/(q-2).                      (6)

Every ai>=1, so replacing each suffix sum by its length proves the
inequality by a geometric series. Equality holds exactly when every ai=1.

Give each word weight 2^(-A); these weights sum to one over all words of
length m. Let mu_m be the law of B_w*Q^(-1) modulo q^m, and let nu_m be
the circle-valued law of c_w/Q_w modulo one with the same weights. Then

    |mu_hat_m(u)-nu_hat_m(u)|
      <= 2*pi*|u|*(2^(-m)-q^(-m))/(q-2).              (7)

Indeed |e(s)-e(t)|<=2*pi*|s-t|, so average (5)-(6). Absolute summability
justifies the infinite word sum. For any fixed frequency u and rate
rho>1/2, an O(rho^m) bound for one transform is equivalent to such a bound
for the other. At q=3,u=1 the desired rate lies between 1/2 and
log2(3)-1. This reformulates H1 as cancellation among the
normalized legal binary source addresses, with an explicitly smaller
O(2^(-m)) discrepancy. The cancellation bound itself is still OPEN.
These are prefix-word laws with the usual continuation at 1, not
distributions of first-hit certificates.

The source address is a classical inverse-parity object:

    c congruent to -sum_(i=1..m) q^(-i)2^(a1+...+a_(i-1)) (mod Q).

Equation (3) retains the additional endpoint-oddness bit modulo 2Q.

For q=3 this is the finite inverse-parity expansion in
[Bernstein and Lagarias, section 1, formula 1.6](https://websites.umich.edu/~lagarias/doc/bernstein.pdf).
Their 2-adic conjugacy supplies context, not an integer convergence theorem.
The explicit quantitative comparison (7) is proved above from our finite
carrier conventions; no external theorem is needed for that derivation.

The chronological direction also matters. Expanding the incoming Fourier
window recursion from level m downwards, a branch (a1,...,am), with prefix
sums sj, has rational phase coefficient

    sum_(j=1..m) q^(-(m-j+1))2^(-sj)
       = B_(reverse w)/(q^m 2^A).                    (8)

Reversal preserves the word weights, so it leaves the aggregate law
unchanged. Individual controller words must still retain their direction.

The subsequent [two-clock theorem](source_address_two_clocks_20261004.md)
uses this source coordinate to prove full-Haar orthogonality across odd-step
levels and an exact nearest-neighbor energy on odd seeds. It computes both
clock sums and isolates their unproved pointwise interchange at seed 1.

### Why a phase estimate does not certify every source

The exact word

    (4,1,1,1,1,2,2,1,2,1,1,2,1,1,1,2,3)

has P=129140163, Q=134217728, B_w=1106233681 and ports c=165,d=167.
Its least source grows, 165->167, although P<Q. Its next source pays:
268435621->258280493. Both have the same word, native residue guard,
and Fourier phase. The integer height in (2) decides the difference.
Both literal paths were independently checked without root padding.

Equation (7) transfers one averaged observable. It does not transport
universal source payment. A candidate pointwise proof must retain the
source column and show that its actual branch meets a paid wall or a
grounded terminal certificate. Density and phase cancellation cannot erase
the exceptional source obligation.

## 5. The root law has a finite carry tile and an unbounded tower

For an integer or 2-adic seed z, set

    R_(n,d)=z*q^(-n) mod 2^d,
    q*R_(n+1,d)=R_(n,d)+2^d*c_(n,d),  0<=c_(n,d)<q.

The next horizontal bit b_(n,d) satisfies
R_(n,d+1)=R_(n,d)+2^d*b_(n,d). Compatibility of the square gives

    c_(n,d)+q*b_(n+1,d)=b_(n,d)+2*c_(n,d+1).         (9)

The output bit and next carry are therefore uniquely determined by the
input bit and old carry:

    b'=(b-c) mod 2,    c'=(c+q*b'-b)/2.

This is an exact q-state binary transducer for division by q in Z_2,
starting with c=0 at depth zero. It retains the root branch omitted by a
bare equation e(theta')^q=e(theta). Its branch label is
c=-R_(n,d)*(2^d)^(-1) mod q.

At q=3 its six transitions (c,b)->(b',c') are

    (0,0)->(0,0), (0,1)->(1,1),
    (1,0)->(1,2), (1,1)->(0,0),
    (2,0)->(0,1), (2,1)->(1,2).

The concurrent Pascal tower expands q^m=(1+2)^m in phase coordinates.
Equation (9) is the same modular compatibility with its integer carry
retained. In both cases intermediate coordinates can be eliminated
exactly, but eliminating them does not make their information bounded.

**State boundary.** The composite division-by-q^m transducer has exactly
q^m distinct residual states in the deterministic, least-significant-bit
input/output model. At depth d its carry is

    C=-R*(2^d)^(-1) mod q^m.

For sufficiently large d, even restricting R to odd prefixes reaches every
one of the q^m values. From a residual carry C, an all-zero future input
has output -C/q^m in Z_2. Distinct carries give distinct output strings,
so no two states can be merged while preserving all continuations.
The q^m-state construction and this distinguishing test prove the claim.
Equivalently, one may retain the m local q-state tiles or an integer with
about m*log2(q) carry bits. A fixed finite alphabet is not fixed total memory.

This is not a lower bound on arbitrary encodings of a single fixed seed,
nor a no-go theorem for all symbolic algorithms. The row law (9) also holds
at z=0, where every phase is one and there is no cancellation. The actual
nonzero fixed seed and the weighted functional remain essential for H1.
Randomizing each level discards compatibility (9); the incoming hybrid
experiments test that loss, but their observed rates are not used as a proof.

## 6. The same coding lemma sharpens the geometric carrier

Let o=(1,1,1)/3 in the strict triangle interior. The incoming centroid-fan
charts A_c map that interior into the three disjoint regions with unique
minimum coordinate c. The six barycentric flag charts B_sigma map it into
the six disjoint strict coordinate-order chambers. The seed o lies in
neither union. Thus the lemma in section 2 gives:

**Marked-centroid theorem.** For either pure chart family, the single point
A_w o or B_w o determines every finite chart word w. Decode the first
minimum/order and invert its chart until the marked centroid is reached.
A matrix and its determinant are unnecessary for this restricted question.

This repairs a possible overreading of the incoming arbitrary-point
counterexample. A point with unspecified terminal local coordinates loses
history; the fixed marked seed imposes additional information. A terminal
port permutation still fixes o and is invisible, so retain it for local
coordinate composition. The chosen sharing of subexpressions is also absent.

Rational coordinates alone are not a termination certificate. The fan
inverse has the exact cycle

    (1,2,6)/9 -> (1,5,3)/9 -> (4,2,3)/9 -> (1,2,6)/9,

with branch labels 0,0,1. These are strictly interior triadic points, and
none is a finite marked-centroid address. A verifier uses a supplied depth
bound and rejects the cycle. Even a denominator that is a power of 3 does
not justify silently running the geometric decoder until it reaches o.

For the mixed fan/bisection alphabet the incoming matrix identity
A1 A0 A1 E0 E0=E0 E0 A0 A1 A1 remains a collision after evaluation at o;
the common point is (175,37,112)/324. Here the disjoint-image premise fails,
and the typed tree is still required. This gives both a successful transfer
of the coding lemma and a precise counterindication to extending it.

## 7. The resulting web and next obligations

| Source and target | Actual map | Preserved predicate | Loss and required coordinate |
|---|---|---|---|
| Divisor profile to controller | the prior P^2QR word profile, then ordered affine composition | operation multiplicities; with carry, the full action | chronology and native guards must be restored |
| Controller word to rational code | w -> F_w(0), using the five disjoint ternary tags | full word and therefore its typed local proofs | cost grows; redundant K must be normalized to LG |
| Native guards to a small automaton | target progression reduced to ordinary or last-L state | nonemptiness iff LB is absent; the G/B/L sublanguage has golden growth | actual source, period, height and terminal proof are lost |
| Native rules to a reusable boundary record | [c,m;d,l], composed by (1) | every legal source and endpoint, including empty domains | reducing the affine denominator alone loses guard bits |
| Collatz word to Fourier phase | [c,2Q;d,2P] -> e(u*d/P) | the exact summand; source-address replacement has bound (7) | source address, period and integer height are not recovered from the phase |
| Fourier tower to local tiles | the carry square (9) | cross-level branch compatibility | flattening m levels needs q^m distinguishable residual states |
| Pure refinement word to one point | w -> chart_w(o) | finite word with fixed interior seed | final port frame; arbitrary-point termination; mixed history |

A useful integer object is therefore a source N together with a typed
native progression, its recoverable operation code, the immutable payment
comparison, and either a smaller dependency or a supplied terminal proof.
Forgetting this structure returns N. The additional coordinates are
mathematical data with explicit decoders, not a promise of free information.

The most concrete new analytical target is (7): prove the required decay
directly for the source-address measure while retaining the compatibility
tiles. The most concrete controller target is to extend the native guard
cover, using the finite automaton to reject unrealizable switches and (1)
to recover their exact source cylinders before expensive searches.
Neither storing a finite receipt nor composing valid local receipts proves
that every integer has a completed receipt.

## 8. Reproduction and audit scope

Run from the repository root:

    python3 -B 04-computation/experiments/collatz_carry_interfaces_20261004.py
    python3 -O -B 04-computation/experiments/collatz_carry_interfaces_20261004.py

The program uses only exact integers and Fractions, with explicit exceptions
that remain active under optimization. Its saved output counts 1,936,562
checks. The principal universes are:

- All 5,461 H/G/A/B words through length six, two positive lifts and one
  negative control per lift; 10,921 compiled first-hit-safe common-future
  checks plus the identity at ROOT; 7,225 independent port-composition pairs.
- All 19,531 H/G/A/B/L words through length six. Translation decoding is
  injective on the entire formal universe; 3,546 native domains are empty.
  Nonempty domains are replayed at two positive lifts, with local dependency
  joins and every admitted credit prefix checked. All 29,791 triples of
  words of length at most two check native-port associativity, including
  empty domains. All 488,281 formal words through length eight independently
  check native-language and balanced-tree counts. The G/B/L sublanguage is
  enumerated through length eight; reciprocal trace identities are checked
  through n=20. The all-depth language proof is the two-state invariant.
- 13,888 ordinary-word cases: q in {3,5,7,9}, m=1,...,5, total valuation
  cost at most 14. Phases and correction bounds are exact rational
  identities. An independent probability recursion checks the complete
  finite alphabet a in {1,2,3}, through depth four for q=3,5,7.
- Carry tiles for five multipliers, seven signed/zero seeds, twenty levels
  and 32 bit depths; reachable and distinguishable composite states for
  q=3,5,7 and tower depths one through four.
- 3,280 fan centroid words through depth seven and 9,331 flag centroid
  words through depth five; the explicit triadic cycle and mixed collision.

Malformed codes, missing native bits, an empty balanced program, equal
phases with different source ports, and opposite height-payment outcomes
are hostile controls. The infinite conclusions follow from the proofs
above; these finite audits check their formulas and implementation.
