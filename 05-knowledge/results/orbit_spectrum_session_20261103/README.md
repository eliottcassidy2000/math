# Exact tournament counts and constructive spectrum interfaces

**2026-11-03. PROVED scoped constructions and counting reductions /
FINITE-EXACT certificates / CITED inherited mechanisms;
the assertion that 7 and 21 are the only missing positive odd Hamiltonian
path counts remains OPEN.**

This session gives an actual infinite moving-port construction, a precise
arithmetic return guard, and improved ways to compute class-count data.
It also isolates an obstruction that the next induction must overcome.

## 1 An infinite blocked branch is repaired by moving an old mark

Start with the strong tournament T_m=x→y→TT_m→x, whose count is 2^m+1.
Keeping x,y as feedback vertices, **no suffix of any length** gives either
3H+2 or 3H+4 when m≥3. A positive clone-run cone proves this, with a small
nonrealizability argument closing its one numerical exception.

For m≥2, instead select the last old clone and x as the feedback pair.
Two specified births preserve every old arc and name and give B_m with

    H(B_m)=3*2^m+5=3H(T_m)+2.

The old physical vertex x has two-slot endpoint count 6*2^m−3. The inherited
two-birth reset therefore produces an actual balanced rooted profile
(b,t,t,b), with named root and child, where

    b=15*2^m−1,      t=3*2^m+5,
    gcd(b,t)=1 iff m is not 11 modulo 12.

Examples include primes 59 and 239, and 119=7*17, 959=7*137, 7679=7*1097.
When m=3j the quotient b/7 obeys q_(j+1)=8q_j+1, q_1=17, and the pair
is always primitive. The construction is valid for all m≥2, but primality
is claimed only for the verified examples. The m=1 input has its own
explicit plus-two base case. [Proof, named cuts and return data](coverage/NOTE.md).
The receipts retain full graph rows, physical old-name maps, all seventeen
two-apex continuation queries, and the literal balanced child graphs.

The numerical family 15*2^m−1 was already realized in
[prime-lane exchanges R2](../prime_lane_ports_five_groups_20261102/realization.md),
and the generic ready-root mechanism already supplied a plus-two count.
The added information here is the exact blocked-versus-repaired attachment
rule, physical mark transport, balanced child, and its arithmetic guard.

This is a reusable balanced interface, not a closed arbitrary-target
induction. The child is a different source family. Even after every
possible initial feedback re-marking, the plus-four branch fails on each
T_m with 3≤m≤10 under arbitrary suffix growth; these are complete
all-suffix results for eight specified sources, not an all-m theorem.

## 2 Compute exactly the requested part of A000568

[The certified-head calculator](counting/NOTE.md) applies an inherited
moved-support Burnside bound to a variable cutoff and checks the exact
distance to the next integer. It returns a rational certificate for the
binary length and requested leading bits, retaining only as many cycle
cores as that query requires.

For a(10000), it certifies **49,876,542 binary digits** and the first 256 bits
without constructing the complete integer. A complete
[a(256) checkpoint](counting/full256.json), with 30,957 bits, uses the
inherited exact engine. Its leading 256 bits and sixteen normalized low
bits are checked independently. The repository already had exact counts
through 300; this checkpoint is validation, not a new enumeration record.

The partial-query computation is not advertised as a faster algorithm
for the full integer. The measured small-order loss case is retained.
Incoming projective-response counting is credited to its older exact
quotient and integer gauges; its independent replay remains useful.

## 3 Exact class counting can retain the prime and seven-prime fibres

Let B_n(h) count fixed-Hamilton-spine presentations with H=h, and N_n(h)
ordinary isomorphism classes. The inherited prime asymmetry theorem gives

    N_n(p)=B_n(p)/p,                    p>3 prime.

For p≥11 prime,

    N_n(7p)=B_n(7p)/(7p)+(6/7)S_n(7p),

where S counts the symmetric classes. Summing N_n(h) over h is exactly
A000568(n); retaining these fibres preserves the target information lost
by the total class count.

The symmetric strong 7p fibres with 7p<19467 are now completely classified:

| H | Symmetric strong classes | Order |
|---:|---:|---:|
|4151|2|10|
|5789|2|10|
|14791|2|12|
|16079|4|12|
|17717|2|12|

There are no others in this range at any order. A positive marked-module
kernel and a complete quotient-order 7 basis prove the all-order cutoff.
For every nontrivial module M and strong rooted quotient Q of order≥7,

    H(Q[M])≥9(H(M)+F_M(2)),

where F_M(2) counts ordered two-slot path covers with empty slots allowed.
The quotient-order 7 bound is sharp. Neither this inventory nor the class
formula calculates all asymmetric fibres or proves target existence.
[Proof, ordinary padding series and exact symmetry correction](bridge/NOTE.md).

## 4 What the next induction must change

The seventeen-query two-apex observer is already complete for its forward
word operation. Seventeen **strong** sources span its full linear query
space. This proves that a fixed suffix giving aH+b on every strong source
must be only neutral sink padding. Adding more suffix length cannot repair
that universal template.

There is also a stronger source-specific boundary: every strong H15
tournament, after every two-apex marking with transitive deletion, misses
47 under all suffixes. Strong H15-neutral native resets followed by this
growth phase therefore cannot supply the needed branch. The exact missing
data are the operation and the ordered attachment cut, not just primitivity.

An interior insertion gives the positive control 112→10102, with H9→H31,
preserving both apices and every old arc. Its simple repeated-pattern
extension has the wrong exponential coefficient, so the isolated repair
is not promoted to a uniform law. [Strong-span proof and complete source menus](words/NOTE.md).

The remaining obligation is a target-dependent construction whose **actual
output source** returns to a family supporting the next required branch.
The present results use no exact-spectrum axiom. The known exclusions of 7
and 21 and the reduction to prime and7p realization remain inherited; proving
that every other positive odd value occurs remains OPEN.

See [recovery and concept board](RECOVERY.md), [audit and reproduction](AUDIT.md),
and [incoming connections](INCOMING.md). No literature-priority claim or
Lean formalization is made.
