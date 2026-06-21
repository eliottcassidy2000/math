---
id: HYP-2725
title: LRC14 generated-word compatibility hands off to relation-code packets
status: OPEN; exact bridge scout
source: codex-2026-06-21-S71
depends_on:
  - HYP-2724
  - HYP-2723
  - HYP-2722
  - HYP-2721
  - HYP-2719
  - HYP-2698
  - HYP-2702
related:
  - THM-558
  - THM-561
  - OPEN-Q-108
---

# HYP-2725 - Generated-Word / Relation-Code Bridge

## Claim

The HYP-2724 relation-code lens is not a replacement for HYP-2722 generated
miss-zeta word compatibility.  It is the next filter after compatibility has
reduced the problem to surviving low-support packets.

The proof stack should stay sequential:

```text
generated miss-zeta word compatibility
  -> low-support relation-code packet selection
  -> factorial q0 / Vitali atom boundary evaluation.
```

The tempting shortcut

```text
relation-code A3/dmin directly proves generated q0 compatibility
```

is too strong.  The exact bridge scout finds relation-code signal, but not a
single monotone scalar that replaces the generated-word death-chain and
context-merge lemmas.

## Exact Scout

Script:

```text
04-computation/lrc14_miss_zeta_relation_code_bridge_codex_s71.py
```

Stored output:

```text
05-knowledge/results/lrc14_miss_zeta_relation_code_bridge_codex_s71.out
```

The script reuses the HYP-2722 exact generated-word frontier:

```text
tests = 318
unique sparse-tail challenger shapes = 72
```

For every challenger shape, it aggregates the smallest generated-context
barriers:

```text
min q0,
min W1+W2 leakage,
min U4/q0,
min B2/q0,
min atom tax,
min distance to the cheap r=1 abstract LP direction.
```

It then attaches the KPS HYP-2724 relation-code observables for primitive
relations with `|coef| <= 2`:

```text
dmin, A2, A3, A4,
dA3 = A3(consecutive block) - A3(shape),
dA4 = A4(consecutive block) - A4(shape).
```

## Evidence

Globally, the relation-code features are real signal.  Across the `72` unique
frontier shapes:

```text
Pearson(A3, min_q0) = +0.830
Pearson(A3, min_U4) = +0.885
Pearson(A3, min_B2) = +0.904
Pearson(dmin, min_q0) = -0.778
Pearson(dmin, min_B2) = -0.800
```

But this global readout is size-confounded.  Within fixed shape size the signs
are mixed:

```text
size=4:
  Pearson(A3, min_q0)  = -0.492
  Pearson(A3, min_W12) = +0.568
  Pearson(A3, min_U4)  = +0.455
  Pearson(A3, min_d1)  = -0.945

size=5:
  Pearson(A3, min_q0)  = +0.039
  Pearson(A3, min_W12) = +0.160
  Pearson(A3, min_U4)  = +0.592
  Pearson(A3, min_d1)  = -0.006
```

The strongest warning is size `3`: it contains the worst generated
compatibility witnesses while the simple `|coef|<=2` relation spectrum is flat.
The smallest q0 and smallest U4 witness is still:

```text
size=3 shape=(0,1,3)
A2/A3/A4 = 0/0/0, dmin = 99
q0 = 20/16807
W1+W2 = 2199/614
U4/q0 = 2187/2005
B2/q0 = -533/280
```

So the generated death-chain/gap geometry is doing work that low-support
relation-code counts do not see.

The relation-code proxy is nevertheless useful on the finite packet side.  The
closest cheap `r=1` abstract LP direction is a size `5` low-support packet:

```text
shape=(0,1,2,3,5)
A2/A3/A4 = 1/10/31, dmin = 2
q0 = 199/24010
W1+W2 = 163599/61661
U4/q0 = 69930/61661
distance to cheap r=1 = 918112/308305
```

This is exactly the type of packet HYP-2724 can organize after HYP-2722 has
already shown it is not a silent generated q0-hiding move.

## Tournament Analysis

Vertices are sparse-tail challenger shapes, not runners or arcs.

Generated-risk edge:

```text
smaller (q0, W1+W2, U4/q0, distance to cheap r=1).
```

Relation-code proxy edge:

```text
smaller dmin, then larger A3, then larger A4.
```

Tie Hamiltonian path:

```text
lexicographic shape order.
```

Fingerprints from the scout:

```text
size=3: vertices=56, pairs=1540, relation-proxy flips=602
size=4: vertices=7,  pairs=21,   relation-proxy flips=6
size=5: vertices=8,  pairs=28,   relation-proxy flips=14
```

Every generated-risk tournament is transitive (`directed_3cycles=0`) because
the risk observable is a lexicographic scalar.  The high flip counts are the
point: the relation-code proxy and generated-word risk order are related but
not the same quotient.

## Proof Target

Use HYP-2725 as a guardrail for the LRC14 proof:

1. Prove the singleton generated-word exclusion using the HYP-2702 death-chain
   kernel.
2. Prove coherent context merging cannot erase the relevant q0, `W1+W2`, and
   `U4` witnesses.
3. Only then apply HYP-2724 to classify surviving low-support relation packets
   by `dmin`, `A3`, MDS/arc structure, and finite tournament-like types.
4. Evaluate the HYP-2721 q0/Vitali atom boundary after both compatibility and
   relation-support selection have been retained.

Rejected assumptions:

```text
A3 is a direct generated-compatibility certificate;
MDS/arc distance alone ranks the sparse-tail q0 bottlenecks;
the size-3 death-chain frontier can be seen by |coef|<=2 relation spectra;
full-row convex-order, product-word compatibility, and relation-code packets
are one scalar theorem.
```

The useful theorem target is narrower:

```text
Compatibility reduces surviving packets to a low-support relation ledger whose
leading finite cases are organized by A3/dmin/MDS.
```
