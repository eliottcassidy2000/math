# HYP-2459 - Parity projectors split midpoint odd channels from reversal even channels

**Status:** OPEN synthesis; projector lemma is exact, finite tournament audit verified for
labelled tournaments on `n=3,4,5`.
**Source:** codex-2026-06-13.
**Companions:** HYP-2458, HYP-2457, HYP-2456, HYP-2455, HYP-2444,
HYP-2443, HYP-2425, HYP-1308, HYP-1387, HYP-1477, HYP-1521, THM-466,
OCF/`Omega(T)`.
**Computation:** `04-computation/parity_projector_channel_atlas_codex.py`;
stored output `05-knowledge/results/parity_projector_channel_atlas_codex.out`.

This is the dual addendum to HYP-2458.  HYP-2458 isolates the midpoint
anti-symmetric side: odd Faulhaber moments survive.  This note adds the
tournament converse side: reversal-invariant scalars survive only on even
Walsh channels.  The two gates should be used together as a bookkeeping
principle for LRC14 and the other support-gate problems.

## Exact Projector Lemma

Let `J` be an involution on a finite state space and let

```text
P_+(f) = (f + f o J)/2,
P_-(f) = (f - f o J)/2.
```

Then `P_+` is the invariant projector and `P_-` is the anti-invariant
projector.  In particular, if a tournament on `n` labelled vertices is encoded
by edge signs `x_e in {+1,-1}`, the converse map is the global reversal

```text
J(x) = -x.
```

A Walsh monomial `x_S` picks up the factor `(-1)^|S|` under `J`.  Therefore:

```text
F(-x)= F(x)  => F has only even Walsh levels.
F(-x)=-F(x) => F has only odd Walsh levels.
```

This is the formal version of the prompt's split:

```text
midpoint scalar anti-gate     -> odd offset powers/moments survive
tournament reversal invariant -> even Walsh channels survive
```

## Scalar Midpoint Side

For paired offsets around a midpoint `c`,

```text
(c+z)^d + (c-z)^d keeps even powers of z,
(c+z)^d - (c-z)^d keeps odd powers of z.
```

The Faulhaber anchor is almost this pure anti-gate, except the interval balance
has one fixed central atom:

```text
D_p(c,n)
 = c^p + sum_{j=1}^n ((c-j)^p-(c+j)^p)
 = c^p - 2*sum_{r odd<=p} binom(p,r)c^(p-r)S_r(n).
```

So the surviving scalar data are:

```text
fixed-point atom: c^p
odd channels: S_1, S_3, S_5, ...
```

This refines HYP-2458: the odd moment inventory is not alone; it is attached
to a fixed midpoint atom.

## Tournament Reversal Side

The stored audit enumerates all labelled tournaments for `n=3,4,5` and computes
Walsh support under global converse.  It verifies:

```text
H(T)                  even Walsh only
c3(T)                 even Walsh only
writhe(T)             odd Walsh only
paths starting at 0   mixed
start0 + end0         even Walsh only
start0 - end0         odd Walsh only
raw edge flip delta   even Walsh only
oriented H-gradient   odd Walsh only
```

The important new refinement is the marked-viewpoint transport law:

```text
start0(T) = end0(T^op).
```

A rooted or observer-coupled statistic need not be even or odd by itself.
Reversal transports it to its dual perspective.  Only after pairing the two
perspectives do we obtain honest projectors:

```text
start0 + end0  -> even/invariant channel
start0 - end0  -> odd/anti-invariant channel
```

Similarly, the raw flip delta of `H` is still even, but the oriented derivative
obtained by multiplying by the marked edge sign is odd.  This matches the old
HYP-1521 slogan that the `H`-gradient is a Mobius-bundle section.

## LRC14 Transfer

For LRC, the immediate use is a clock taxonomy.

```text
even scalar clocks:
  reset period, denominator shell, observer distance, unmarked q-blocked total,
  complement/reversal-invariant tournament scalars such as H

odd marked clocks:
  owner support, carry residue, deletion derivative, pair sum/difference cut,
  oriented edge/runner pressure, writhe-like imbalance

transported clocks:
  source/sink, start/end, left/right boundary, observer-coupled perspectives

compatibility clocks:
  OCF alpha packets, Q27 owner/carry packets, code72 support-design packets,
  polynomial convolution factor grids
```

The proposed proof hygiene is:

1. quotient aggressively only on even scalar clocks;
2. split transported clocks into sum and difference before quotienting;
3. use odd marked clocks to prove pressure, descent, or owner-private openings;
4. attach compatibility packets before trusting an odd atom inventory.

This is directly aimed at the current LRC14 target.  A scalar row such as
"q is blocked" should not be treated as a complete state.  It should be split
into an even shell component and odd/marked fields: shell-27 class, divisor
fiber, owner support, carry residue, boundary atom, and deletion/opening target.

## Cross-Problem Reading

The same parity projector explains several recurring repo themes:

* **OCF and tournament enumeration.**  `H` is complement-invariant, so its
  scalar Walsh support is even.  Odd information reappears in gradients,
  rooted perspectives, deleted-card derivatives, and compatibility packets.

* **Unit distance.**  The metric equality graph is an even scalar channel.
  Unit-spine, orientation, nonunit tiling, and endpoint-compatible ears are
  marked support channels.

* **The `[72,36,16]` code.**  The Type II weight enumerator is an even scalar
  gate.  The open problem lives in the support/design/matroid lift.

* **Polynomial irreducibility.**  Sign-cube and residue scalars are quotient
  gates; convolution grids, Newton slopes, fixed-divisor rows, and
  factor-capture certificates are marked or compatibility channels.

* **Faulhaber/triangular towers.**  Odd moments survive the midpoint gate, but
  HYP-2458 remains essential because odd atom counts still need OCF-style
  packet compatibility.

## Tournament Analysis

Alternate vertex sets considered: runners, arcs, gaps, fixed circle sections,
section boundaries, denominator events, residues, cover arcs, Fourier modes,
edge signs, rooted perspectives, gradients, support packets, and proof
obligations.  The stored carrier tournament uses proof obligations as vertices.

Pairwise observable: lexicographic majority over LRC leverage, lift readiness,
parity clarity, exactness, algorithmic value, and cross-domain transfer.

Switch/gauge: higher retained-channel score wins; ties use the listed order as
the Hamiltonian path.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1]
hamiltonian_paths=1
leaders=[
  lrc14_q27_owner_carry_ledger,
  marked_section_parity_toggle,
  lrc_sign_cut_side_channel,
  ocf_compatibility_packets
]
```

This quotient preserves parity type, marked/transported viewpoint status, and
proof-carrier role.  It destroys raw time geometry, exact runner positions,
and unprojected circular endpoints.  The challenged assumption is that
"the tournament" must be a tournament on runners or arcs; here the useful
tournament is a ranker over parity-gated proof carriers.

## Open Problem

Build an LRC14 `Q27` ledger whose fields are explicitly typed as even scalar,
odd marked, transported, or compatibility-packet data.  Then test whether the
remaining primitive rows become pure after applying the projector protocol:

```text
even quotient first;
transported start/end or source/sink split second;
odd owner/carry/deletion pressure third;
OCF-style compatibility packet check last.
```

The target theorem would say that any row surviving the even quotient either
has a strict witness, descends to AP/`Vstar`/`2AP`, or exposes an odd marked
owner-private opening.
