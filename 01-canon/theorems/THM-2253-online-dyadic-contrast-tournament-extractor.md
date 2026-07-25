---
id: THM-2253
title: "Online dyadic-contrast tournament extraction"
status: >
  PROVED + VERIFIED-EXACT. For any almost-surely nonconstant exchangeable
  source on a finite alphabet, orient the distinct alphabet pairs at each
  dyadic node by any tournament. Stop at the first aligned dyadic interval
  whose two halves are constant and distinct, and output their tournament
  orientation. Swapping the homogeneous halves is a mass-preserving,
  output-reversing involution of terminal cylinders, so the output is exactly
  fair. Unknown-law i.i.d. categorical sources are a special case. If n is
  the initial constant-run length, the stopping time satisfies
  tau<=n+2^nu_2(n)<=2n; unless n is a power of two it satisfies
  tau<=4n/3. This strengthens the adaptive part-(a) coin deadline and gives
  a lawful state-dependent alternative to a fixed response antipode; it does
  not supply exchangeability in LRC, GMC, or Jacobian carriers.
source: codex-2026-07-25-online-dyadic-contrast
related:
  - THM-2160-dyadic-checksum-extracts-a-fair-bit-under-the-critical-run-deadline
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
  - THM-2235-response-antipode-barriers-for-lrc-sheets-tournament-cycles-and-knot-kernels
script: 04-computation/categorical_coin_online_dyadic_contrast_thm2253.py
output: 05-knowledge/results/categorical_coin_online_dyadic_contrast_thm2253.out
script_sha256: 17f1d6653e09eab85f7cda370fe72abf1aeacbc8f9f87a48c20d47c511f988dd
output_sha256: 148b8a49b136ff77c491d33f17bca9c5f408318c902bce7628f36f0ba40bc8fc
hash_basis: working-tree bytes (LF)
---

# THM-2253 -- the first dyadic contrast is a fair tournament edge

THM-2225 answers the binary critical-run puzzle by a dyadic block
compression and, for the sharper deadline, a cyclic checksum. The compression
can be made genuinely online and its actual stopping time can be read from the
two-adic address of the first change. More generally, the ordered pair exposed
at the stopping node carries an intrinsic tournament.

## 1. The online rule

Let `A` be a finite alphabet and let

```text
X_1,X_2,... in A
```

have an exchangeable law: every finite coordinate permutation preserves
their distribution. Assume the source is nonconstant almost surely. An
i.i.d. categorical source with unknown law `pi` and `max_a pi(a)<1` is the
motivating special case.

For each aligned dyadic node `I`, fix any tournament `T_I` on `A`: for
distinct `a,b`, exactly one of `a->_I b` and `b->_I a` holds. The tournament
may depend on the node address, or on any observed context invariant under
swapping the two children of `I`.

An aligned dyadic interval is

```text
I=[2Lj+1,2L(j+1)],                 L a power of two. (1)
```

Its left and right children have length `L`. Call `I` a **contrast node**
when each child is constant and their two symbols are distinct.

Read the source from left to right. Stop at the first time `tau` which is
the right endpoint of a contrast node `I`. If its children carry `a,b`,
output heads when `a->_I b` and tails when `b->_I a`.

At a fixed right endpoint there is at most one contrast node. Aligned dyadic
intervals are laminar. If one contrast node were strictly contained in
another ending at the same time, the containing child would not be constant.
Thus the rule is unambiguous.

For a binary source the tournament is the single edge `0->1`.

## 2. Pathwise termination and its valuation deadline

Let `n` be the initial constant-run length:

```text
X_1=...=X_n=a,                 X_(n+1)=b!=a.          (2)
```

Put

```text
L=2^nu_2(n).                                           (3)
```

Because `n/L` is odd, the two length-`L` blocks on either side of the
boundary after `n` form the children of an aligned dyadic interval. Its
left child is constant `a`.

If the right child is constant, it is constant `b` and this interval is a
contrast node. If the right child is nonconstant, it contains a contrast
node of its own. Indeed, every nonconstant word of power-of-two length
contains such a node: at the root, either its two children are distinct
constants or at least one nonconstant child permits descent; the process
terminates at a differing adjacent pair.

In either case a contrast node has appeared by the end of the right block:

```text
tau<=n+2^nu_2(n)<=2n.                                 (4)
```

The constant continuation `a^n b^L` shows that (4) is sharp for this rule
at every `n`. If `n` is not a power of two, the odd integer `n/L` is at
least three, so

```text
tau<=n+L<=4n/3.                                       (5)
```

Thus the only critical runs on which this online rule can take the full
`2n` allowance are powers of two. THM-2225's checksum remains the separate
mechanism which obtains `max(2,2n-1)` on those dyadic boundaries.

For the sample `00001`, (4) gives the requested eight-flip ceiling. The
checksum rule is faster on this power-of-two boundary.

## 3. The terminal-cylinder involution

Let a terminal prefix end at time `tau`, and let `I` be its selected
contrast node with constant children `a^L,b^L`. Define `sigma` by swapping
those two children and leaving every other coordinate unchanged:

```text
... a^L b^L  |->  ... b^L a^L.                       (6)
```

The selected node remains a contrast node and its tournament orientation
reverses. No earlier node is created or destroyed:

- an earlier aligned interval disjoint from `I` is unchanged;
- an earlier aligned interval inside `I` lies wholly in one constant child,
  and after the swap it is still constant;
- an aligned ancestor has a nonconstant child containing `I` both before
  and after the swap.

The same observations handle other intervals ending at `tau`. Consequently
`sigma` maps a terminal prefix to a terminal prefix with the same stopping
time and selected node, and applying it twice restores the prefix.

The swap is a finite coordinate permutation, so exchangeability makes its
two terminal cylinders equiprobable. In the i.i.d. special case this is the
familiar composition identity

```text
P_pi[X_1...X_tau=w]=product_(a in A) pi(a)^N_a(w).   (7)
```

It is fixed-point-free and reverses heads and tails under the node-indexed
tournament. Hence it pairs all terminal cylinders and gives

```text
P_pi(heads)=P_pi(tails).                              (8)
```

By (4) the rule stops on every nonconstant sequence. These have total
probability one by assumption, so the two probabilities in (8) sum to one.
The extracted bit is exactly fair. In particular, this holds for every
unknown nondegenerate i.i.d. categorical law.

## 4. Why this is a genuine tournament and not a cosmetic one

The vertices are source symbols, the pairwise observable is the ordered
pair of distinct homogeneous child labels, and ties never occur at a
contrast node. Reversing the children reverses the node's tournament edge.
The swap preserves the full composition vector, terminal time, dyadic
address, and every allowed orientation context. It forgets the internal
coordinates only because each selected child is constant. Exchangeability,
not independence, is precisely the probability predicate it preserves.

Unlike THM-2225's cyclic checksum, this response is state-dependent. It is
not a homomorphism from one fixed source action to `F_2`, so it lawfully
evades THM-2235's odd-group/fixed-XOR obstruction.

That distinction identifies a possible but conditional LRC consumer.
Instead of asking a `13`-power sheet group for a nonexistent global
antipode, one may search a laminar carry/owner tree for the first pair of
sibling continuation states with an exact exchange involution. The required
lemma is strong: sibling exchange must preserve guard mass, labelled owner
incidence, and every future consumer. Bernoulli composition makes (7)
automatic; current LRC transfer states are not exchangeable, so this theorem
does not itself exclude any of the `166` scalar profiles. The likely lawful
replacement is a defect-bearing sibling swap, with the complementary owner
stratum paying the exact imbalance.

For GMC root packets, conjugate weights are likewise not equal scalar masses.
For the planar Jacobian response quotient, target shears preserve the present
response but need not be interchangeable for an unnamed continuation.
Thus the theorem contributes a new operation—**first laminar contrast and
local exchange**—together with the exact sidecar test that any transfer must
pass.

## 5. Exact referee

Run

```bash
python3 04-computation/categorical_coin_online_dyadic_contrast_thm2253.py
python3 -O 04-computation/categorical_coin_online_dyadic_contrast_thm2253.py
```

The companion exhausts binary terminal prefixes through length `16` and
ternary prefixes through length `10`, using a nontransitive three-cycle
tournament in the latter case. It checks node uniqueness, terminal-time
preservation, involutivity, output reversal, and exact bisection of every
observed composition shell. It also tests every continuation through the
deadline for binary `n<=16` and ternary `n<=8`, including the sharp constant
right-block witnesses and the nonpower `4n/3` inequality. The normal and
optimized transcripts are frozen byte for byte. QED.
