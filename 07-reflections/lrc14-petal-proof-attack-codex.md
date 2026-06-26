# LRC14 Petal Themes and Proof Attack

Codex session note, 2026-06-24.

This note extends the petal side of the LRC14 work and ties it back to the
marked Farey/tournament viewpoint.  The computation is in
[lrc14_petal_proof_attack_codex.py](/home/bigo/math/04-computation/lrc14_petal_proof_attack_codex.py:1),
with output at
[lrc14_petal_proof_attack_codex.out](/home/bigo/math/05-knowledge/results/lrc14_petal_proof_attack_codex.out:1).

The main point is not a new proof of LRC14.  It is a sharper proof target:

```text
Inside the minimal one-petal branch, mark-unital means AP or 12->24.
```

The exact terminal core in the `v<=300` single-replacement bank is:

```text
AP       -> M=1/14
12->24   -> M=1/14
10->20   -> M=2/27
8->16    -> M=2/23
```

So the local petal problem has become a clean marked-node problem.  The two
non-AP/GW terminal petals are not counterexamples; they are explicit migration
witnesses.

## Flower Denominators

The flower prompt gives a useful coordinate picture:

```text
1/pi = 0.318309886...
7/22 = 0.318181818...
```

After 22 petals at angle `1/pi` radians:

```text
22/pi = 7.002817496... radians.
```

Thus the 22-family is a radian-denominator effect.  It should not be confused
with full-circle closure, which uses:

```text
1/(2*pi^2).
```

The new petal ledger places that 22-family next to the LRC denominators:

```text
22  flower radian denominator
23  nonunit petal mark 2/23 from 8->16 and 9->18
27  pair-sum modulus 2*14-1 and unit petal mark 2/27
41  first nonplanar unit child 3/41 from 12->36
```

This is a coordinate system, not a proof.  The proof-bearing denominators are
still the marked Farey nodes and the divisor thresholds.

## Doubling Petals

For the minimal doubling petals `h -> 2h`, the exact marks are:

```text
7 -> 14   M=1/11   apex-blocked
8 -> 16   M=2/23   nonunit child
9 -> 18   M=2/23   nonunit child
10 -> 20  M=2/27   unit child
11 -> 22  M=2/25   nonunit child
12 -> 24  M=1/14   tight apex
13 -> 26  M=2/27   unit child
```

Only `12->24` is both marked at `1/14` and blocker-free in the
Goddyn-Wong window:

```text
[14-h, 27-2h].
```

For `h=12`, this window is `[2,3]`, and neither endpoint is coprime to 12.
Every other doubling petal has a coprime blocker.  This is now the cleanest
local mechanism for why the top petal is special.

The `n+2` recursion probe behaves differently.  The first valid same-parity
tails mostly fall into apex-blocked or loose-up rows.  It does not reproduce
the mark-preserving `12->24` event.  That suggests the `n*2` mode is the
actual petal mechanism in the local AP/GW branch, while `n+2` is better read as
a residue-lift stress test.

## Marked Farey Migration

The spectrum/Farey language keeps the failure modes unified:

```text
tight      -> marked node 1/14
loose-up   -> coarser divisor node, q(S)<14, or apex blocked by a multiple of 14
loose-down -> Farey child/descendant above 1/14
```

The named exemplars are:

```text
AP          -> 1/14
GW 12->24   -> 1/14
7->14       -> 1/11   apex blocked
12->26      -> 1/12   same-residue divisor liar
8->16       -> 2/23   nonunit child
10->20      -> 2/27   unit child
12->36      -> 3/41   unit child, first nonplanar child
```

The complete-bipartite packet readout is:

```text
p/q -> K_{p,q}.
```

Thus `2/27` is still bipartite-planar as `K_{2,27}`, while `3/41` is the first
unit child carrying a `K_{3,3}` obstruction through `K_{3,41}`.  This matches
the earlier Farey packet fact that `3/4 -> K_{3,4}` is the first reduced
ordinary Farey packet crossing the complete-bipartite nonplanarity threshold.

## Exact Atlas

The exact diagnostic atlas was kept intentionally rerunnable:

```text
unique AP single replacements, v<=80:
  rows = 872
  tight = 2
  below threshold = 0
```

The only tight rows in this exact atlas are AP and `12->24`.  The broad
`v<=300` bank is used for the cheap filter/terminal census rather than for an
all-row exact mark distribution.  That mirrors HYP-2920: exact marks matter
most at the terminal core, while the earlier filters are proof-obligation
gates.

The `v<=300` filter-strength tournament has no directed 3-cycles under the
pass-set inclusion orientation.  The terminal counts are:

```text
one_petal_or_ap       270
minimal_one_petal       4
top_petal_or_ap         2
```

The equal pass-set pair in this bank is:

```text
exact_odd_skeleton = even_dipole_or_zero.
```

That equality is local to this single-replacement bank, but it is a useful
warning: some filters that look conceptually distinct can collapse on a bounded
frontier.  A proof should keep the semantic labels, not only the row counts.

## Proof Attack

The current LRC14 route can be organized as follows.

First, use theorem-level divisor and unit gates:

```text
q=2..13 covered,
q=14 omitted,
unit residues mod 14 represented.
```

Failure here gives loose-up or apex-blocked migration.

Second, do not redo the AP-tail exact-measure work.  THM-542, THM-543, and
THM-544 already close the one-tail and two-tail below-second AP collar.  The
petal attack is about marked-node preservation, not safe-measure interval
arithmetic.

Third, prove the terminal petal lemma:

```text
If a primitive tight row enters the minimal one-petal branch, then the petal is
12->24, or the row is AP.
```

The exact escape witnesses to eliminate are already visible:

```text
8->16   -> M=2/23
10->20  -> M=2/27
```

Fourth, prove a global entrance-or-death statement:

```text
Every primitive tight survivor either enters the one-petal branch, or is killed
by a multi-dipole/far-tail theorem that preserves the marked Farey node.
```

This is where the multi-tail, summand/multiplicand, unital averaging, and
tournament metagraph work should connect.  The proof quotient must be unital in
the algebraic sense: it must preserve the marked identity node `1/14`.  Any
quotient that identifies AP/GW with `12->26`, `10->20`, or `12->36` without
remembering the mark is non-unital for LRC14.

