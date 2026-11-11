# Future response order and the arithmetic content it does not preserve

**PROVED in this note.** The exact q order compares future count sizes.
Its integral content is an additional resource for prescribed-target
coverage. The root resets in [PROOF.md](PROOF.md) remain valid; the two
notions of improvement must be distinguished when composing them.

## 1. The response content persists through the permitted grammar

For q=(H,U,U−d,H+U−E), set

    c(q)=gcd(H,U,U−d,H+U−E)=gcd(H,E,U,d).

Every common O/Z/J word acts by an integer matrix on q. A fixed right-root
amalgamation also acts by an integer linear map on (H,E,U,d): this follows
from the exact [amalgamation formulas](../../h13_paired_descent_20261111/modular_cycles/ADDITIVITY.md).
Consequently c(q) divides every coordinate of every subsequent state,
and in particular every subsequent H. Content can increase; it cannot be
removed by those operations alone.

For the small paired reset, the lower branch133 keeps content1, whereas
the upper branch changes

    (135,13,141,13), content1
      → (135,27,135,27), content27.

The latter is strictly smaller in the future q order, but every permitted
future count from that marked source is divisible by27. In particular a
3H+2 or3H+4 transition from that source is unavailable within this grammar.
The generic reroot formulas give its content explicitly:

    ε=0: gcd(H,d'),
    ε=1: gcd(H,d',3).

For the lower N(C) branch, H≡1 mod3, so ε=1 guarantees content1.
This is a positive primitive-response criterion, not just a bound on size.

## 2. There is an absolute q minimum at every attainable count

For every rooted tournament,

    U≥H,   U≥d,   E≤H.

The first inequality has a direct injection: split a Hamilton path just
before the root; its suffix is the distinguished root-start path and its
prefix is the complementary path. A Hamilton path already starting at the
root uses an empty complement. Concatenating complement and distinguished
path recovers the original word, so the map is injective. The root-alone
distinguished path contributes all d paths of the deletion tournament,
proving the second inequality. The third is immediate.

It follows that

    q≥(H,H,0,H).

This minimum is attained physically: append a global sink to any tournament
with count H and mark that sink. The resulting state is(H,H,H,H), so its
q vector is exactly(H,H,0,H). It preserves the count and has content H.

Thus a count-preserving decrease of q alone cannot be the missing
Chan–Pak induction. There is already a universal minimum, but it retains
only future counts divisible by H. The nontrivial value of the new
rebases is their realization inside specified small-core families, with
explicit root responses and sometimes primitive content, rather than an
unqualified claim that smaller q gives greater target coverage.

## 3. The next operation must carry both resources

A source-selection rule should retain the actual graph/root, its q vector,
and c(q), together with any promised core or port restrictions. A physical
ordinal factor can sometimes be removed to divide counts and responses;
formal division of an integer vector does not establish such a graph map.
The asymmetric453/399 pair restores common ports but has content3 on both
branches. Its exact factor structure must therefore accompany any attempt
to interpret that restoration as a reusable nonzero-digit transition.

This connects the3-versus1 pointer to an actual invariant: primitive
response content1 allows all residues in principle; content divisible by3
forces the entire permitted future into multiples of3. It is a necessary
condition, not a sufficiency theorem for reaching arbitrary targets.
