# Lossless source and payment carriers across controller switches

2026-10-04. **PROVED:** source-preserving affine composition; the precisely
scoped3/5/6 polynomial observer hierarchy; injectivity and exact source
guards for the H/G controller alphabet; and the prefix-credit composition
law. **FINITE-EXACT:** the controls declared below. This is a guarded
program compiler, not a proof that every source admits a paying program.

Artifacts: [program](../../04-computation/experiments/lossless_payment_carriers_20261004.py)
and [saved output](lossless_payment_carriers_20261004.out).

## 1. Inheritance and types

The closest mechanisms are the marked affine carry decoder in
[the four-channel kernel, sections5–6](reciprocal_four_channel_kernel_20261004.md),
the guarded dependency H in
[the recursive dependency kernel](collatz_recursive_dependency_kernel_20261004.md),
and the inherited valuation-dependent graph rank in
[the branch-toll rank note](collatz_branch_toll_rank_20261004.md).
The canonical hostile is equal slope with unequal carry/order. A corrected
near miss is resetting the source at a switch: local descent may leave
unpaid debt to the initial source. The least-used sidecar is that immutable
source coordinate. Our board is **source / affine carry / polynomial
observer / program order / exact guard / credit prefix**.

Use positive odd integers and

    H(x)=(729x+669)/1024,  guard x=155 mod2048;
    G(x)=(9x+5)/8,        guard x=11 mod16.

G is the actual valuation word(1,2). H is a common-future dependency,
not an actual forward Collatz word. Let V=(1,2,1,1,1,2). On H's guard,
F_V(x)=4H(x)+1. If a=v2(3F_V(x)+1), then a>=3 and the actual source
word V,(a) and child word(a−2) have the same endpoint. The script checks
these literal words, forbidding root padding. Composing H/G retains this
distinction; an H/G string alone is not an actual orbit word or a home
certificate. A separately supplied rooted terminal certificate can be
transported through the inherited common-future machinery.

## 2. The immutable source and the3/5/6 observer hierarchy

For F(x)=a x+b, keep state(x,n,1), where n is the original source. Its
matrix and an integer homogeneous version, for a=P/Q,b=B/Q,Q>0, are

    T3 = [a 0 b]              M3 = [P 0 B]
         [0 1 0]                   [0 Q 0]
         [0 0 1]                   [0 0 Q].

The row(1,−1,0) tests payment x<n. A homogeneous change of coordinates
must retain its transformed row and a positive affine-chart scale; the
guard is also transported. This is the interface to
[the triangle-chart carrier](small_port_refinement_20261004.md), not an
identification of three coordinates with three graph vertices.

If current x=(P n+B)/Q and the next program is(p x+c)/q, update

    (P,Q,B) -> (pP,qQ,pB+cQ).                              (1)

Consequently exact original-size payment is

    (P−Q)n+B<0.                                           (2)

Neither the current value alone nor a last-stage local decrease records
this comparison. For example the legal GGGH program gives

    212987 ->239611 ->269563 ->303259 ->215895.

The final H decreases its input, but215895>212987. Every stage guard and
common-future witness is checked independently in the program.

For quadratic observations, five coordinates

    (x²,n²,x,n,1)

are closed under all rational affine changes of x with n fixed. The
explicit matrix is

    [a² 0  2ab 0 b²]
    [0  1   0  0 0 ]
    [0  0   a  0 b ]
    [0  0   0  1 0 ]
    [0  0   0  0 1 ].

This five-dimensional polynomial space is minimal **if the required
starting observations are x,n,1 and x²−n², and arbitrary affine current
updates are allowed**. Comparing x²−n² with4x²−n² isolates x² and then
n². All five monomials are linearly independent. This is a linear
observable statement, not an information-theoretic lower bound on
arbitrary encodings. Without the original linear monitor, the orbit of
x²−n² spans only{x²,n²,x,1}, dimension4. Requiring closure alone does not
force5, and one fixed passive polynomial of n remains one coordinate.

The full degree-at-most-two space has six coordinates

    (x²,xn,n²,x,n,1)=Sym²(x,n,1).

The new mixed coordinate transforms as x'n=a xn+b n. Dropping it gives
the five-dimensional observer and commutes with every allowed update.
Keeping it permits arbitrary quadratic observables and the debt-square
(x−n)²=x²−2xn+n². These are the precise reasons for5 and6; the separate
seven-port geometry is a different object.

The inherited proper rank illustrates a useful five-coordinate observable
with an extra valuation sidecar. For odd n>1, write n−1=2^K t with t odd:

    E(n)=3^K t²=(3/4)^K(n−1)²,   rank=(E(n),K).

ROOT has rank(0,0). Given exact Kx,Kn, E(x)−E(n) is a linear functional
of the five coordinates, with coefficients

    (a,−b,−2a,2b,a−b), a=(3/4)^Kx, b=(3/4)^Kn.

The unbounded valuation has not disappeared. Numeric descent does not
universally lower this rank:15<17 but E(15)=147>E(17)=81. Conversely the
actual step171->257 lowers energy21675->6561. If the original n is3 mod4,
every smaller positive odd x has smaller energy, since E(x)<=3(x−1)²/4
and E(n)=3(n−1)²/4.

## 3. The full affine H/G map remembers its ordered program

Shift z=x+5. In this coordinate

    H:z->r z+d, r=729/1024, d=67/32;
    G:z->s z,   s=9/8.

For a chronological word with h copies of H and g copies of G, the
slope lambda=r^h s^g has

    u=v3(lambda)=6h+2g, v=−v2(lambda)=10h+3g,
    h=(2v−3u)/2,        g=5u−3v.                          (3)

Thus the slope already recovers the counts. Write its shifted translation
C=d T. Each H contributes one summand

    r^(number of later Hs) s^(number of later Gs)

to T. If the last H has k trailing Gs, its summand has3-valuation2k.
Every earlier H has valuation at least2k+6. The unique least term cannot
cancel, so v3(T)=2k. Strip those Gs, then the last H, by

    T <- (T/s^k−1)/r,

and decrement g by k and h by1. Repeat, ending with the remaining pure
G suffix of the reverse decoder. When h=0, T must be0. Exact re-encoding
checks all supplied coefficients. This proves injectivity of the full
affine representation of the free H/G word monoid, including the empty
word. In the original coordinate C=B/Q+5−5P/Q.

This decoder is particular to the marked H/G alphabet. It does not assert
injectivity of every finite affine generating set. The proof uses exact
3-adic separation of the carry contributions, not only a dimension count.
HG and GH have the same counts and slope but different carries, so counts
alone lose precisely information needed at a switch.

An independent unshifted decoder gives a shorter check. If the last
letter is G, its composite carry is9 B_old+5 Q_old, with3-valuation0.
If the last letter is H, it is729 B_old+669 Q_old, with3-valuation1.
Here Q_old is a power of2, and the empty prefix has B_old=0. Thus v3(B)
identifies the last letter directly. Divide P,Q by that letter's P,Q,
then recover B_old=(B−B_letter Q_old)/P_letter. Repeat to the identity.
This is independently implemented and checked against the shifted decoder.

## 4. A single exact congruence compiles all stage guards

For any nonempty H/G word let its unreduced integer carrier be(P,Q,B).
Here P is an odd power of3, Q is a power of2, and B is positive odd.
For a positive integer source n, the whole word is legal if and only if

    P n+B = Q mod2Q.                                     (4)

Equivalently its sources form the single odd cylinder

    n=(Q−B)P^(-1) mod2Q.

Each individual letter's guard is exactly the condition that its affine
endpoint is an odd integer; solving this gives155 mod2048 for H and11
mod16 for G. Necessity of(4) follows by composition. For sufficiency,
write the word as its first letter followed by a nonempty tail. Its
numerator has the form

    P_tail(P_first n+B_first)+Q_first B_tail.

Reducing(4) modulo Q_first forces the first endpoint to be an integer,
since P_tail is odd. Divide by Q_first; as Q_tail is even and B_tail is
odd, the first endpoint is odd. It therefore satisfies the first exact
guard. The divided congruence is the tail's odd-endpoint condition, so
induction finishes. The one-letter case starts the induction. Positive
slopes and carries preserve positivity; the exact first-letter cylinders
give the necessary minimum positive inputs. The empty word uses all odd
sources, represented by1 mod2.

Together, sections3 and4 say that the **full typed affine map plus exact
source** retains both program order and all H/G legality. A trace, slope,
modular observation, or unauthenticated summary does not. Decoding still
costs time proportional to the recovered word; this is an exact compact
representation, not a claim of constant-time expanded verification.

## 5. Composable prefix credit and payment scope

Give H increment+2 and G increment−1. A word has summary(m,Delta), where
Delta is its total credit change and m is the minimum prefix sum including
the empty prefix. For chronological concatenation uv,

    (m_u,Delta_u)*(m_v,Delta_v)
      =(min(m_u,Delta_u+m_v),Delta_u+Delta_v).             (5)

Splitting the prefixes into those in u and those extending into v proves
the law and associativity. Starting credit k is legal iff k+m>=0. This
summary can be cached on an authenticated composition tree without
expanding it again. It is not a self-authenticating certificate: the
leaves and their guards must be validated, and arbitrary claimed minima
cannot be trusted. HG and GH both end at credit1, but have minima0 and−1.
The full affine carry can recover that order; its slope cannot.

[The credit-potential theorem](adaptive_credit_potential_20261004.md)
proves original-size payment for every nonempty legal H/G word whose
prefix credit remains nonnegative, starting current=original and credit0.
Its potential is(x+5)(9/8)^k. The present program independently verifies
that payment over its stated finite universe, while proving the carrier
and summary interfaces for arbitrary words. An arbitrary checkpoint with
a reset source or unproved initial potential bound is outside that
theorem. Numeric payment, graph-rank payment, a finite guarded episode,
and a completed route to ROOT remain separate predicates.

## 6. Reproduction and audit scope

Run from the repository root:

    python -B 04-computation/experiments/lossless_payment_carriers_20261004.py
    python -B -O 04-computation/experiments/lossless_payment_carriers_20261004.py

The finite universe is all511 H/G strings of lengths0..8. Every full map
is decoded by two independent algorithms, and two positive members of every exact source cylinder are
independently replayed:1022 sources. Every cut checks affine and credit
composition; each of510 nonempty words also tests a single-bit guard
failure. All430 credit-valid nonempty source controls pay their original
source. There are3375 independently parenthesized summary triples using
all words of length<=3. For a,b,x,n from{−2,−1/2,0,1,3/2}, all three
observer dimensions give1875 exact state transports and1875 independent
matrix composition controls. Additional controls include the explicit
local/global hostile, rank boundaries, and eight invalid types/carriers.

These exact controls test the proved algebra and types. They do not show
that the guarded alphabet applies to every supplied source, pays every
unresolved state, or completes every Collatz orbit. No priority or external
literature claim is made.
