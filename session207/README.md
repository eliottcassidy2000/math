# Session 207: tournament counts, orders, and metagraph incidence

Let H(T) count directed Hamiltonian **vertex words**, S_n={H(T):|T|=n},
and m(h) be the least realizing order, or infinity. This is a navigation
summary of independently audited packets, not a replacement for their proofs.

## The answer, with its exact scope

**Every individual input is decidable by a terminating exact algorithm.**
For a fixed n, enumerate tournaments and compute their path counts to obtain
S_n. For a proposed h, the finite witness bound below decides attainability
and computes m(h), including a certificate of absence when the search fails.
Thus “no simple complete closed form” does not mean “no exact algorithm.”

**The uniform range conjecture remains open:** it is not proved that every
positive odd integer other than 7 and 21 occurs. Those two omissions, and
all nonpositive/even omissions, are proved. This session does not supply a
simple closed expression for unrestricted S_n or m(h) at every input.
The open-status and prior-source record is [here](../05-knowledge/results/rotation_owners_spectrum_20261020e/spectrum.md).

## 1. Exact forward and inverse answers

- [Forward spectrum and continuation state](size-spectrum/PROOF.md): an
  independent complete census through n=8, with all full/strong supports
  and 430 adjacency-code witnesses. At n=8 there are 320 counts, maximum
  661; the missing odds below that maximum are 7,21,611,615,617,619,623,
  625,635,647,655. This independently recovers known small-order data;
  the [inherited order-nine result](../01-canon/theorems/THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension.md)
  is not displaced or claimed new.
- [Inverse-size reduction](inverse-size-spectrum/README.md): every inverse
  set is exactly {n:n>=m(h)} or empty. Strong components multiply counts.
  If A_r is the strong count support, then S_n=union_r A_r S_(n-r), with
  S_0={1}. Minimum orders obey an exact divisor/min-plus recurrence using
  minimum strong-atom orders. Composite atoms must be retained: count 9
  first occurs on five vertices, earlier than two separate triangles.
- [Complete inverse prefix through 661](inverse-prefix661/PROOF.md): all
  329 allowed odd counts in [1,661] have certified exact minimum orders,
  recorded in [the full table](inverse-prefix661/minimum_orders.json).
  The nine previously missing allowed counts above have minimum order 9.
  Exhaustiveness through eight plus nine explicit witnesses proves this;
  it is not a new exhaustive order-nine census. Padding supplies every
  subsequent order.

## 2. A sharper all-input cutoff and a uniform restricted-family solution

[Base-eleven scalar compression](base11-compression/PROOF.md) proves that
any attained h>1 has a witness of order m satisfying

    3^5 * 11^m <= 11^3 * h^5,
    m(h) <= floor(3 + 5 log_11(h/3)).

This replaces components by different tournaments with the same scalar
count; it does not shrink the original SCC core unchanged. Together with
finite enumeration it gives a terminating exact membership/minimum-order
algorithm, still taking 2^O(L^2) time for L-bit h. The case h=1 is separate.
It also gives an exact forward stabilization window: if
11^3*h^5 < 3^5*11^(n+1), every globally attained h>1 already belongs to S_n.
This window does not assert that all its odd integers are attained.

[Moon interval chains](chain-partitions/PROOF.md) admit a **complete
all-order restricted-family solution**. A chain with block lengths w_i>=2
has order 1+sum w_i and count product(2^(w_i-1)+1). Using the known Moon
product law and cited Bang–Zsigmondy theorem, its only scalar relations are
permuting blocks and replacing a block of length 4 by two blocks of length 2.
Canonical partitions therefore exclude parts 1 and 4. With p(t) the ordinary
partition number and p(t)=0 for t<0, the exact number of distinct chain
counts at order n is

    p(n-1)-p(n-2)-p(n-5)+p(n-6).

After source padding, the exact constructed subfamily of S_n has
p(n-1)-p(n-5) distinct counts. Descending integer division recognizes this
infinite family and computes its unique unpadded chain order in polynomial
bit time, without factoring. That chain order need not equal unrestricted
m(h): count 15 has chain order 6 but unrestricted minimum 5.

## 3. The native metagraph needs more than count and size

Here the native move preserves a directed Hamiltonian path and reverses
**every other arc**. It is neither a single-arc flip nor a vertex-extension
edge. For n>=3 its labelled degree is H; ordinary class-line incidence uses
H/|Aut(T)|. Converse merging and paired-black conventions require the
additional corrections stated in their source packets.

[The complete H=3 incidence row](../05-knowledge/results/session207_degree_three_metagraph_20261106/PROOF.md)
classifies every source as Tr_a -> C3 -> Tr_b and computes its neighbour
count g(a,b) by tribonacci recurrences in both parameters. At order 5,
Tr_2 -> C3 has only count-9 neighbours, while Tr_1 -> C3 -> Tr_1 has only
count-13 neighbours. Both sources have H=3. Thus the H-partition is not
equitable, and the walk choosing a uniform Hamilton path is not strongly
lumpable under H alone. This is a fixed-order obstruction, not just a
comparison across sizes.

[Universal transport](../05-knowledge/results/session207_universal_metagraph_transfer_20261106/PROOF.md)
extends this to every fixed core T and chosen Hamilton path P. For
T_ab=Tr_a -> T -> Tr_b, the actual path-complement neighbour counts satisfy
tribonacci recurrences separately in a and b. Nine corner counts per path
orbit determine every padded count. After a,b>=3, the neighbour automorphism
group is the original neighbour's endpoint stabilizer; one additional group
order per orbit converts counts to ordinary fibre degrees. The resulting
finite signatures determine exact degree-labelled incidence measures throughout
the stabilized region a,b>=3; count rows are determined at every padding.
They do not determine the entire metagraph or identify distinct target classes.

In particular every attainable source count has actual neighbours of
unbounded degree along all sufficiently large consecutive sizes. This
strengthens the earlier [padding amplifier](../05-knowledge/results/session207_metagraph_fibre_amplification_20261106/RESULT.md).
No size-independent bound on neighbour degree can depend only on the
source count or its fixed ordinary fibre degree.

## 4. Infinite multiplicity and typical-degree families

[Mixed-prime multiplicities](mixed-prime-multiplicities/PROOF.md) determines
the first two coefficients of the eventual ordinary-class polynomial for
every fixed H=3^b5^c with b+c>=1, and its first-order deviation from prime-component
stacks. In the explicitly defined paired-black converse-merged metagraph,
a fraction 1-O(1/n) of that fixed-H stratum has degree 5^c. This graph's
conventions must not be substituted for the labelled graph above.

[Finite-defect universality](coefficient-defect/PROOF.md) goes further:
for fixed j and b+c>=max(1,j), the coefficient of n^(b+c-j), multiplied
by b!c!, is a rational
polynomial in b,c of total degree at most 2j. Only finite strong-atom
inventories of prime weight at most j+1 enter. These multiplicity results
do not assert a closed minimum-order formula for the same infinite family.

## What remains to do

The terminating algorithms, finite inverse table, exact restricted family,
and universal fixed-core transport are established at their stated scopes.
The global odd-except-7-and-21 existence conjecture, simpler unrestricted
minimum-order formulas, and all-core metagraph structure remain separate
obligations. All linked packets distinguish inherited ingredients, cited
theorems, finite-exact computations and their genuinely all-order deductions.
