# Exact readiness at the three distinguished roots of the two-birth operation

**PROVED, independently audited by the root reviewer. FINITE-EXACT controls:**
all rooted labelled tournaments through order four, plus left/right padding
controls. Other old vertices are explicitly outside this classification.

## Operation and inheritance

For a nonempty tournament D with marked vertex u, write h=H(D),
q=H(D-u), and let e count Hamilton paths of D ending at u. Set H(empty)=1.
Let G=G(D,u) count ordered covers by a nonempty path starting at u and a
second path which may be empty, as in [the induced repair](../repairs/TERNARY.md).
Call the marked vertex ready when G-h=2.

Adjoin a sink a to D, then a vertex b which beats only u and loses to every
other old vertex and to a. Write D' for this tournament. The inherited
two-birth lemma gives H'=H(D')=h+2G. Thus a ready source gives 3h+4.
This note determines exactly what happens at u, a, and b; it does not
classify the remaining old vertices or forbid a later change of core.

## N1. Exact marked response

For every nonempty rooted D,

    G(D',a) = H' + h - e,
    G(D',b) = H',
    G(D',u) = 4G - q,
    H(D'-u) = q.

In particular the three readiness defects are respectively

    h-e,  0,  2G-h-q.

**Proof.** In D+b, Hamilton paths avoiding b->u are exactly the old paths
not ending at u, followed by b; there are h-e. Paths using b->u correspond
to a root-starting path and its complementary path, with b between them;
there are G, including the case of an empty complement. Thus
H(D+b)=h-e+G.

A path starting at a is either the singleton a, the path a,b, or a,b
followed by a nonempty old path starting at u. The weighted complements
in the three cases contribute h-e+G, h, and G. This proves the a formula.

A path starting at b is either the singleton b or b followed by an old
path starting at u; in the latter case it may additionally end at a.
The contributions are h, G, G, proving the b formula.

For a path starting at old u, distinguish which of a,b it contains.
Neither gives G; only a gives G; only b gives G-q because the singleton
old path u cannot be followed by b; both give G, since they must appear
at the end in the order a,b. Their sum is 4G-q. After deleting u, a,b
are two successive sinks, so the remaining count is q. QED.

## N2. Complete classification at the new vertex a

The vertex a is ready if and only if

    D = TT_s -> C3,  s>=0,  with u a vertex of the directed triangle.

In every such case h=G=3 and H'=9. In particular a can be ready only at
the resulting count nine, even when the old order is arbitrarily large.

**Proof.** Readiness means h-e=2 by N1. Hence u is not a sink. For every
Hamilton path of D-u, insert u immediately before its first vertex beaten
by u. All preceding vertices beat u, so this gives a Hamilton path of D
not ending at u. Deleting u recovers the old path, proving q<=h-e=2.
Tournament path counts are positive odd numbers, including q=1 for an
empty deletion, so q=1. A tournament with exactly one Hamilton path is
transitive; consequently e<=1. Since h is odd, h-e=2 forces h=3 and e=1.

The inherited cyclic-core bound K(3)=3 says that any tournament with h=3
is TT_s -> C3 -> TT_t. This follows directly from its SCC decomposition:
there are exactly three cyclic vertices, forming the unique nontrivial
component. Since D-u is transitive, u lies in that triangle. A nonempty
right padding would prevent a Hamilton path from ending at u; hence t=0.

Conversely, for TT_s -> C3 rooted in the triangle, exactly one Hamilton
path ends at u. A path starting at u stays in the triangle. The transitive
prefix is forced into its complementary path and does not change that
complement's count. There are exactly three root-starting paths in C3,
each with complementary count one. Thus h=G=3, e=1, and N1 gives readiness
at a and H'=9. QED.

The inherited facts about oddness, transitive H=1, and the cyclic-core
bound are used with their usual nonempty-tournament scope; the empty
deletion is handled by H(empty)=1. For the bound see
[Session 133 core reduction](../../session133_fractran_spectrum_refinement_20261101/core_reduction/THEORY.md)
and [the growth audit](../../strategy_fibres_absorption_20261113/GROWTH.md).

## N3. A ready source loses readiness at all three distinguished roots

If G=h+2, then none of u,a,b is ready in D'. At b the defect is zero.
At u it is h+4-q>=4: inserting u into an old Hamilton path gives the
standard injection q<=h. At a, N2 would force G=h=3, contradicting the
source preparation G-h=2.

This rules out automatic rearming by retaining the marked vertex or
choosing either newborn. Another old vertex, a count-preserving rewrite,
or a different source can still supply a return. The finite H79 two-arc
repair in [SOURCE_TRANSFER](../repairs/SOURCE_TRANSFER.md) is outside this
three-root obstruction and is not contradicted by it.

## Exact controls

`natural_roots_check.py` independently enumerates directed paths to compute
G, uses subset DP for complementary H counts, and checks whole Hamilton
counts/endpoints by literal permutations. Its complete small universe is
all 285 rooted labelled tournaments of orders 1 through 4, with no filter.
It checks all four N1 equalities, N2's iff, and N3 whenever its hypothesis
holds. Separate C3 controls use multiple left and right padding lengths,
including arbitrarily extendible positive left-padding and hostile
right-padding families. The finite checks corroborate the proof rather
than establishing its universal quantifiers.

Reproduce from the worktree root:

    python -B 05-knowledge/results/prime_terminal_20261113_1610/absorption/natural_roots_check.py

The output is retained in `natural_roots_receipt.json`.
