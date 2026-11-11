# Complete reset dictionary for the neutral four-defect atlas

**PROVED all-k graph constructions and response laws; FINITE-EXACT complete
operator dictionary in the specified universe. Chan–Pak remains OPEN.**

The inherited [C4 atlas, A3–A4]
(../../ready_return_pentagon_20261103/ATLAS.md) classifies all neutral
four-defect roots by sixteen tails after a prefix11^k. We exhaust the
following operation, while preserving the prefix, root, both ports, and
all arcs between retained old vertices: delete ONE tail vertex, add two
global sinks in succession, and add one apex with a fixed cut that sees
the entire transitive prefix uniformly. We ask for the count3H+2 or3H+4
at every k>=1. Neutral new singleton padding may then be removed.

## 1. Five source families have exact uniform instructions

Write D=2^k, k>=1. Tail positions in the table are numbered from1.
The surviving cut descriptions use the retained tail and port labels;
a,b are the original ordered ports, B denotes the tail01 and A the
tail10. Every apex loses to every vertex not explicitly listed.

|Original tail|Original H|Delete|Add sinks|Apex beats|Final H|
|---|---:|---|---:|---|---:|
|00,01|9D−1|second tail vertex|2|b|27D−1=3H+2|
|00,10|11D+1|first tail vertex|1|a,b|33D+5=3H+2|
|01,01,10|34D+5|either01|2|retained B and first new sink|102D+17=3H+2|
|01,11,10|34D+5|middle11|2|B and first new sink|102D+17=3H+2|
|01,10,01|30D+5|first01|1|retained B,a,b|90D+19=3H+4|

All five final families are strong. The one-sink instructions arise
from exact hits in the declared two-sink universe where the second new
sink remains global after the apex is added. Removing this NEW singleton
preserves H and every retained old vertex and arc. No old-root or old-port
deletion is hidden in this simplification.

The last line is a new actual uniform+4 family. Its initial values are

    65→29→29→199,     125→55→55→379,
    245→107→107→739.

Its initial H=5(6D+1) is always composite, but the output90D+19 has no
fixed prime divisor, since gcd(90,19)=1. No claim of infinitely many
prime outputs follows. This branch supplies source-compatible+4 targets
outside the fixed-factor+2 construction; it does not cover every source
or every desired integer.

Every value90D+19 is3 modulo4 for k>=1, so this entire target family
is disjoint from the square-plus-one primes A002496.

The first output199 is excluded from BOTH strong1-/2-ready value images
by the complete ready inverse through425 in [READY_REBASE section3]
(READY_REBASE.md). Thus the new+4 operation is real, but returning its
first output to either standard ready interface at unchanged count is
impossible at every order. This is an explicit hostile boundary for
iterating this particular branch through the old ready recipes.

The first two new+2 traces begin17→5→5→5→53 and23→13→13→71.
The remaining two reproduce the earlier uniform141 reset. The first
deletion may leave a weak parent; all intermediate sink steps are
allowed weak tournaments, and the final sources are strong.

## 2. The five-core and defect transitions are explicit

In every target the original marked root still has precisely the
original in-ports a,b, and is neutral: s=q. The exact data are

|Target H|s=q|F−H=delta|Root-and-ports deleted residual|
|---:|---:|---:|---|
|27D−1|9|10|TT_(k+4)|
|33D+5|13|8|TT_(k+3)|
|102D+17|41|24|TT_k→S5a|
|90D+19|33|14|TT_(k+1)→C3|

Thus the operations physically move the neutral4 source into the
neutral8,10,14,24 strata. The +4 family enters the nontransitive residual
stratum with kernel C3; the original+2 family enters S5a. These are actual
residual graphs with retained ports, not an identification based only on
equal path counts. They connect directly to the current
[bounded-defect kernel classification](../defect/BOUNDED_TEMPLATES.md).

For the +4 target at k1 the rows are

    (222,204,80,116,98,199,128,28),

with old root5, ports3,4, and free prefix0. Deleting the root and ports
gives(30,28,8,16,4), visibly TT_2→C3. The other three distinct target
templates and all explicit old-label maps are retained in the receipt.

For each fixed target, remove the free prefix and old root to obtain W,
and let J be its two ports. The inherited all-k two-port identities are

    H=2^k E_J,W(2)+H(W)−C_J,W(2),
    s=q=H(W),       F−H=C_J,W(2).

Independent exact tables on W give(E,H(W),C) respectively
(27,9,10),(33,13,8),(102,41,24),(90,33,14). They prove every
displayed formula for all k. Strongness is certified by an explicit
spanning directed cycle in each k1 template; replacing its prefix
vertex by the ordered k-chain extends that cycle at every k.

## 3. Completeness of the specified operation dictionary

There are six two-letter tails and ten three-letter tails. At k1 the
declared operation examines exactly

    6·2·2^7 + 10·3·2^8 = 9216

literal final cuts, with no isomorphism filtering. Exactly77 cuts have
one of the two requested scalar target counts at k1. For EVERY such
cut the verifier computes the complete nonnegative module-run polynomial.
An all-k identity3(A·2^k+B)+c requires precisely

    A1=3(2A+B)+c,   A2=3A,   Aj=0 for j>=3.

Exactly six cut/deletion instructions pass: the six alternatives in the
five table rows (either01 gives two alternatives). Four have already
strong raw outputs; two have only the removable new global sink.
Every failed scalar candidate retains its complete coefficient list in
the receipt, recording the precise obstructing growth coordinate.

No other atlas tail has a uniform+2 or+4 reset of THIS form. In
particular the Fermat-compatible source16D+1 has none here. This does
not exclude more old-vertex deletions, different numbers or types of
births, nonuniform cuts along the prefix, parameter-dependent charts,
or other return programs. It classifies a bounded physical instruction
set, not all metagraph paths.

Reproduce with `python -B atlas_returns.py`. The verifier rebuilds every
source, all9216 cuts, all77 marked responses, every successful fixed
kernel table, all persistent arcs, and literal k1..6 controls.
