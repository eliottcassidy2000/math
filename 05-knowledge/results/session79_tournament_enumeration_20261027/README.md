# A000568 from re-rooting, prime composition, and observable rank

2026-10-27. **PROVED:** the exact operators and identities below, with the
classical decomposition dependency identified. **FINITE-EXACT:** the declared
independent checks. No claim of historical priority, polynomial complexity,
or a complete classified prime atlas.

## What this session adds

The repository already generates A000568 exactly, including an arithmetic
future-response implementation through order 200 and a Hamilton-relocation
component construction. This session does not present those as new discoveries.
It makes three different kinds of closure precise:

1. **Constructive closure.** A re-rooting operator identifies the child classes
   produced by all parent attachments. Its distinct rows generate the next
   representative bank, and its binary rank is A000568 at that order.
2. **Aggregate prime closure.** The entire prime-tournament cycle index can be
   recovered from the known all-tournament cycle index, without an unknown
   prime census as input. This counts and composes primes but does not list them.
3. **Observer-relative compression.** A terminal arithmetic interface has a
   strictly smaller linear representation than its individual-history quotient.
   Its exact rank also proves that this particular interface cannot have a
   uniformly bounded dimension.

These are different objects and different notions of generation. None is a
license to erase alignment, attachment bits, or future symmetry.

## 1. A concrete metagraph generator

See [the proof](ROOT_SWITCH.md) and [standalone checker](check_root_switch.py).
For each order N, first construct all pointed child states by attaching a vertex
to each order-(N−1) parent, keeping the attachment mask modulo the parent's
actual action. Define Q[i,j] to count the vertices at which rerooting state i
produces state j. Each lookup compares only deletion parents of order N−1,
with the attachment bits transported by the same relabeling.

The block belonging to one ordinary child tournament U is

    Q_U = 1 (m_1,...,m_r),

where the m_i are its vertex-orbit sizes. Consequently

    Q² = NQ,
    a_N = tr(Q)/N = rank_Q(Q) = rank_F2(Q).

Tournament automorphism groups have odd order, so every m_i is odd. This is
why binary rank works even at even N, when Q mod 2 is square-zero rather
than an idempotent. It is not the binary graph-Laplacian construction that
failed in earlier relocation examples.

Distinct Q rows have disjoint supports and correspond exactly to ordinary
child classes. Selecting one raw extension from each distinct row supplies
the next representative bank without child canonicalization. Thus the
construction is an all-order recursive generator, not just an identity that
assumes the child classes have already been found. Its storage and exact
parent comparisons are substantial; no speed advantage is established.

Over the rationals, Q/N is a stochastic projection: rerooting once already
gives the orbit-size-weighted distribution inside the underlying child class.
Uniformity is on actual vertices, not on pointed isomorphism types. This
explains the symmetry weights in the inherited attachment recurrence.

Independent checks cover every state and matrix entry through N=6, with an
additional complete labelled-universe oracle through five. The published
checker tests orders one through six. The all-order conclusion is proved.

## 2. Recover the aggregate prime inventory

See [the proof and exact bank](PRIME_CYCLE_INDEX.md),
[recovery](recover_prime_cycle_index.py), and
[small independent census](verify_prime_cycle_index.py).

Let B=Z_T be the classical all-tournament cycle index, F=B−1, and P=Z_prime
for modular-prime tournaments of order at least three. Canonical strong-module
decomposition and ordered strong components imply

    P ∘ F = 1 − 1/B − p_1,
    P = p_1/(1+p_1) − F^<−1>.

The inverse is **plethystic**, not ordinary univariate compositional inversion.
The explicit odd-permutation Burnside coefficients supply B; weighted-degree
induction then recovers every coefficient of P to any finite requested bound.

The exact bank through order 15 contains only 55 nonzero power-sum monomials.
It includes 796,876,429 prime tournament classes at order 11. This count is a
proved-formula output, not a new exhaustive order-11 classification. Independent
implementations agree coefficientwise through 15, and direct fixed-permutation
module tests agree through six.

There is also a small genuinely scalar boundary. Put f_n=2^((n−1)/2) for odd n.
If g_n counts prime tournaments fixed by one specified n-cycle, then g_1=0 and

    g_n = f_n − sum_(d|n, 1<d<n) g_d f_(n/d).

Thus g_9=12 and g_15=112, independently checked by direct circulant/module
enumeration. These are labelled fixed-point counts, not unlabelled circulant
classes. Every proper nontrivial cycle in a prime automorphism must share a
nontrivial gcd with another cycle; otherwise its vertices form a forbidden
proper module. This supplies an exact arithmetic sparsity test for the bank.

The recovery is not an independent numerical shortcut for A000568: B already
contains its Burnside information. Its gain is removing the unspecified prime
inventory from the aggregate composition theory. It does not recover individual
prime representatives, the coloured metagraph edges, or anti-automorphisms.

## 3. Smaller counting states, with a precise boundary

See [the rank theorem](TERMINAL_RANK.md) and
[standalone exact checker](check_terminal_rank.py).

At a fixed committed degree, suppose h distinct past cycle counts remain and
r further fixed vertices can be inserted. The linear response matrix has entries

    C[k,t] = 2^(k(k−1)/2) (2^t)^k / k!,   0 <= k <= r.

It is a row-scaled Vandermonde matrix. The exact minimal linear dimension for
all these residual outputs is min(h,r+1). Keeping moments can therefore be
strictly smaller than keeping every distinct future-word response.

For the degree-respecting terminal interface, the coordinate counts are
819 versus 722 at horizon 100, and 3302 versus 2872 at horizon 200. Their leading
terms are N²/12 and N²/14. This does not reduce the earlier peak of the full DP,
and a tested direct-tail variant did not provide a useful work reduction.
The theorem is a structural compression and a scoped unbounded-rank obstruction,
not an all-algorithm complexity lower bound or a practical speedup claim.

## Inheritance and the next precise question

The strongest inherited mechanisms remain [weighted-core transport](../weighted_core_transport_20261020/README.md),
[future-response arithmetic](../tournament_future_boundary_20261020/README.md),
and [packed integer enumeration](../tournament_enumeration_20261021a/README.md).
The [corrected modular composition work](../tournament_composition_classification_20261026/README.md)
keeps the order-12 cyclic-triple/automorphism counterexample in force: equal
unsigned cyclic triples do not preserve the symmetry needed for these counts.

The useful unifying principle is **retain exactly the responses required by
the next allowed operation, then quotient only what that operation cannot see**.
For re-rooting, the retained object is the aligned deletion-parent and mask.
For prime substitution it is the permutation-action cycle index.
For a fixed terminal arithmetic observer it can instead be a finite moment
vector. These quotients are not interchangeable.

The next open target is a smaller, computable structural representation of the
root-switch row supports that remains closed under every attachment, without
first materializing every pointed state. Current-parent symmetry alone is known
to fail that task because new symmetries can appear after attachment. No such
smaller unrestricted structural carrier is claimed here.

## Reproduce

From this directory, using Python's standard library:

    python3 check_root_switch.py
    python3 recover_prime_cycle_index.py 15
    python3 check_terminal_rank.py

Repeat with `python3 -O`; correctness checks remain active. The prime recovery
calls the companion independent census itself. No saved receipt, external
package, or pre-existing prime-class list is required.
