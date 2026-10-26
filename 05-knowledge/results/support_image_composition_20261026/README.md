# Recoverable frames, Hamilton-support images, and signed composition

**Status: audited lane proofs, FINITE-EXACT computations, and independent
integration review PASS.** [Integration review](audit_integration/ROOT_AUDIT.md). Research base `4950589d056a64f8807774cf7da893dcafcd63a0`.
Result labels below are local to their linked files. No priority claim is
made. The inherited reconstruction theorem remains a dependency, not a
new result of this packet.

The strongest closure is that **four markers recover every tournament
with two distinct ordered ports**, including nonstrong cores. The new
image criterion identifies exactly which oriented graphs can be Hamilton
supports. Signed marked deletion cards give an unrestricted augmentation
recurrence for A000568. Residue clocks supply a separate exact composition
interface, with a sharp distinction between finite lifting and one positive
integer realizing an infinite path.

## 1. From the singleton sandwich to arbitrary cores

Let `h=[C;s,t]`, `s≠t`, and let `*` be the rooted hinge product. The
inherited `e*h*e` decoder recovers strong cores. Now, with `A2=e*e`,

    E(h)=A2*h*A2=e*e*h*e*e

is injective on joint rooted isomorphism types for **every n≥2**.
Degree-one vertices recover the left markers; indegree-one vertices
recover the right markers. A sink entrance or source exit produces a
three-cycle ambiguity, whose rotations give the same rooted core type.
Precisely,

    Aut(E(h)) ≅ C3^epsilon × Aut(C;s,t),
    epsilon = [s is a sink] + [t is a source].

An exact quadratic image guard tests at most nine marker choices and
re-encodes the whole adjacency relation. Four added vertices are optimal
among **single fixed nonempty two-sided hinge contexts**; this is not a
lower bound for arbitrary encodings. Decode–hinge–encode transports the
associative product to unrooted encodings. It does not identify that product
with native Hamilton dynamics on the larger tournament.

See [sandwich S1–S4](sandwich/README.md) and its
[independent proof audit](audit_sandwich/REPORT.md). This closes only the
all-core `A2*h*A2` clause of [HYP-9416](../../hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md).

## 2. An exact image test, with its first false local candidate

Write U(D) for the arcs occurring in some Hamilton path of D. Let C(D)
be the partial reconstruction algorithm using ordered strongly connected
components and the inherited two-step rule. If C(D) is defined, let

    F_D(z) = sum over Hamilton paths P of C(D)
             of z^(number of P's arcs missing from D).

Then **D is the Hamilton support of a tournament if and only if**

    C(D) is defined, U(D)=D, and F_D(1)=F_D(0).

The polynomial has nonnegative coefficients: equality means completion
created no new Hamilton word. Splitting a word at its missing arcs gives
an exact expansion over ordered spanning D-path covers with their actual
endpoint connectors. This is an exponential exact test, not a new efficient
recognition bound.

The local decoder plus U(D)=D first fails at five vertices. In outgoing
bitmask notation,

    D=(10,20,16,6,1), C(D)=(14,20,16,22,1), F_D(z)=7+2z.

There are 480 such false positives among all 59,049 oriented partial
graphs of order five. Separately, the directed four-cycle gives a joint
connector obstruction `F=4+z_02*z_13`: each added connector alone is safe,
but the pair creates a path. That four-cycle does not pass the local
decoder; its role is to obstruct arbitrary independent connector tests.

The image is **not hereditary under induced deletion**: the Hamilton
support of a transitive triple is a directed path, whose endpoint-induced
pair is edgeless. The compatible deletion is decode–delete–resupport.
See [image I1–I3](image/README.md) and the
[independent audit](audit_image/AUDIT.md).

## 3. What this gives black-line connectivity

For two decoded tournaments on a common vertex set of size n≥3, they are complements
along a Hamilton path exactly when their common arcs form a spanning
directed path. The path is then unique. Its positional reversal tests
blue versus black. For ordinary isomorphism types, summing this indicator
over label transports and dividing by the target automorphism count gives
the exact black-path multiplicity; dividing again by the source
automorphism count makes the expression symmetric.

Two black paths in the same core with the same ordered ports have a common
strong singleton-sandwich hub after padding their complement states.
This produces a two-black-line walk between those padded source classes.
It joins this local operation to the inherited all-order black descent and
terminal-basin reduction, but **all-order black-core connectivity remains
OPEN**. Different port pairs and different terminal basins still need
realized connections. See [image I4–I5](image/README.md).

## 4. An unrestricted signed-card recurrence for A000568

“Even graph” means that every automorphism reverses an even number of
reference-oriented edges. It does not mean even vertex degrees. The
classical equality of unlabelled even-graph and tournament counts is due
to [Royle et al., Theorem 1.1](https://arxiv.org/pdf/2204.01947v2).

For every unlabelled graph G of order N−1 and every neighbourhood orbit S,
adjoin root v, obtaining X. Each marked isomorphism
`f:(G,S)→(X−w,N_X(w))` extends to an automorphism `f̂` taking v to w.
With `K=Aut(G)_S`, set

    Q(G,S) = (1/|K|) sum_w sum_f epsilon_X(f̂).
    A000568(N) = (1/N) sum_[G] sum_[S] Q(G,S).

Q is the root-orbit size if X is sign-even and zero otherwise. Thus the
recurrence uses order-N−1 graph states and marked-card transports, with
no order-N child canonicalization oracle. It still performs isomorphism
search and is **not a recurrence solely in earlier A000568 numbers**.
Odd parents are essential: deleting a leaf from the even graph P3 gives
the odd graph K2.

Root-moving signs multiply under composition. For a relatively-even root
in an odd graph, their cancellation lies in exactly one odd connected
component. This gives the further exact factorization `D(z)=O(z)E(z)`
for the rooted odd defect. See [even_graphs E1–E4](even_graphs/RESULT.md)
and its [independent audit](audit_even/PROOF_AUDIT.md). The inherited
Burnside formula and Hopf/Lie bijections remain prior work; this refines
their compositional interface.

## 5. Clocks, growing sieves, and the arithmetic realization boundary

For shortcut Collatz, use labelled residue arrows

    r --e--> s iff 2s = 3^e r + e (mod M), e∈{0,1}.

Every length-L path has exactly one source residue class modulo `M*2^L`.
Its positive realizations below a chosen height have an explicit floor
count. For M dividing N, reduction has unique **outgoing** path lifting
from every prescribed fine start exactly when N/M is odd. No two-sided
graph-covering assertion is made. Adding an odd sieve prime preserves this
property; adding temporal dyadic precision changes it.

The inherited growing sieve composes with these cylinders using its
prime-equality exceptions and a complete finite-height cutoff. But every
finite all-odd clock prefix has positive lifts while the full path forces
source −1. The strict twin wheels also retain −1; their different
Pythagorean map sends it to the degenerate signed triple (−1,0,1).
This is an exact common realization obstruction, not a Collatz–sieve
conjugacy.

The F8/F1/abc line retains the earlier boundaries: eight-point branch
parity is a Legendre determinant character over F7; F8 joint Frobenius
degree and abc radical/height are separate data. The incoming every-modulus
vector collision is reproved directly; the full transport matrix retains
the determinant at odd moduli M>1. No abc estimate or F1 equivalence follows.
A Turing-machine transport must also retain the fixed source, acceptance
guard and certified native time; a residue-clock cycle cannot replace a full recurring machine
configuration. See [arithmetic C1–C5](arithmetic/README.md).

## 6. The proposed unifying theory and its next obligations

[COMPOSITION.md](COMPOSITION.md) gives precise abstract criteria for
operation-relative recoverability. A quotient supports partial operations
exactly when legality and output observations are constant on its fibres.
An injective encoding with decoder and exact image guard transports the
already well-defined operation theory. Weighted composition sums over compatible joint
witnesses and their transports; associativity needs both parenthesizations
to describe the same admissible multi-input witnesses.

This explains the three retained coordinates: **joint incidence, transport,
and realization guards**. Nonnegative path defects detect new witnesses;
signed card sums intentionally cancel existing witnesses. Those two kinds
of zero must stay distinct. Finite path lifting additionally needs one
common source to support an infinite conclusion.

The next concrete targets are a smaller sufficient image/attachment state,
a tournament/even-graph correspondence compatible with attachments, and
realized black-line transports between every terminal basin. None is
silently inferred from equal counts or faithful encodings.

[BOARD.md](BOARD.md) records inheritance and the portfolio.
[AUDIT.md](AUDIT.md) records proof and finite-check scopes. Run
`python3 -I -B verify_packet.py --live` here for frozen sources, review
bindings, maintained edits and both normal/optimized receipts. Historical
replays should omit `--live` after maintained navigation has advanced.
