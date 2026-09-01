# THM-4311: endpoint-593 one-response size-preserving exchange

**FINITE-EXACT, relative to THM-4309 and THM-4310.  The fixed-pool result
below is proved by complete enumeration and exact certificate quotients.  It
does not give a physical entry or prove LRC(14).**

## Fixed objects

The inherited labelled pool has 30 speeds and hence
`binom(30,9)=14,307,150` labelled nine-bodies per row.  The carrier is
reconstructed from the canonical component ledgers of THM-4296, THM-4300,
THM-4302, and THM-4309; no serialized final carrier is trusted.  Its identity
is

```text
C_593^old: size 3,925, rank8 3,858, rank9 67,
FNV=6fbd0bffcf0ed78b, all 421 joint masks retained.
```

The exact THM-4310 endpoint-593 layer has 16 rows, FNV
`5424c07fa724011f`, SHA-256
`a4bf3d1e9aff29be4fb07d2644912cc098dbd996bc527b4f22dee980f36a6ce6`.
For preservation during exchange, use the complete carrier-closure target

```text
K_432 = THM-4309's 391 rows
        union THM-4310's complete 25-row endpoint-594 layer
        union the 16 endpoint-593 rows,
|K_432|=432, ordered row FNV=a7ed492c64d1c0d8.
```

The three endpoint-594 rows already typed before THM-4310 are deliberately
included. A preliminary 429-row calculation that omitted them was quarantined
during audit and is not part of this canonical packet.

## Exact hostile failure and response minimum

Independent O3 and O2 builds reconstruct the carrier, recompute activity, and
scan all

```text
16*binom(30,9)=228,914,400
```

row-body cases.  Their transcripts, pair ledgers, and failure ledgers are
byte-identical.  There is exactly one failure:

```text
(q,r,body)=(96,593,34087401),
failure-body FNV=22cd6ff93226979c,
ordered obligation FNV=643cc006d08e9a83,
pair-ledger FNV=ab23b8d0da917e96.
```

For an active rank-eight or rank-nine mask `m`, a response to this obligation
means `m & 34087401 = 0`.  A complete two-rank census gives

| rank | all masks | active at `(96,593)` | responders | responder FNV | least |
|---:|---:|---:|---:|---:|---:|
| 8 | 5,852,925 | 1,262,283 | 1,636 | `56f82f5dc11db83b` | `0134012c` |
| 9 | 14,307,150 | 7,218,133 | 16,209 | `0f615d806860553f` | `0036092c` |

Because the old carrier fails, zero additions cannot close the obligation.
Mask `0036092c` is active and disjoint, so one addition does.  The exact
rank-eight/rank-nine response minimum is therefore one.  A structurally
separate complement-subset enumerator starts from the body's 21-element
complement and reproduces both responder counts, FNVs, and least masks.

The selected response is

```text
A={0036092c}, rank9, singleton-list FNV=60873ef7a2b4ab90.
```

Among the old carrier's 3,925 masks, the number inactive on every row is zero
both on the prior 416-row carrier target and on all of `K_432`.  The exact
`3,925*432=1,695,600` sign census has FNV `c0a847088f7c904c` and no equality.
Thus strict-common-inactivity cannot supply a size-preserving deletion.

## Exact nonjoint one-deletion boundary after adjoining the response

Let `C^+=C_593^old union A`.  Then

```text
|C^+|=3,926, FNV=1dac1c2e57f3e682,
421 joint masks fixed, 3,505 nonjoint masks eligible for deletion.
```

On every row of `K_432`, the quotient enumerates every nine-body.  Bodies
already hit by an active joint mask are immutable while the joint deck is
held fixed.  For each remaining body it records a private obligation exactly
when there is one active disjoint nonjoint witness.  The complete result is

```text
body tests                         6,180,688,800
joint-exposed bodies                   1,509,470
private nonjoint obligations              1,520
private-obligation FNV          39fedd8ceb347304
protected nonjoint masks                     412
protected-mask FNV             7ee1d68a078d5b65
safe nonjoint single deletions              3,093
```

This is an exact one-deletion boundary: with the joint deck fixed, deleting a
nonjoint mask fails iff it is the unique nonjoint witness of some joint-exposed
body.  O2 and O3 reproduce the transcript, all 1,520 singleton rows, all 412
protected masks, and the selected deletion byte-for-byte.  The new response
`0036092c` is protected by exactly the original hostile obligation, an explicit
negative control.

The least safe old-carrier nonjoint deletion is

```text
D={0006e281}, rank8, singleton-list FNV=4c14214a64ec202c.
```

No simultaneous-deletion conclusion is inferred from this singleton quotient.

## Same-size exchange and complete raw closure

Define

```text
C_593^new=(C_593^old \ {0006e281}) union {0036092c}.
```

Then

```text
|C_593^new|=3,925, rank8=3,857, rank9=68,
FNV=c9e5faef52ca5707, all 421 joint masks retained.
```

The separate raw exchange program replays every one of

```text
432*binom(30,9)=6,180,688,800
```

cases at O3 and O2.  The outputs are byte-identical and give

```text
joint-exposed bodies=1,509,470,
nonjoint hit incidences=50,693,609,
failures=0,
pair-ledger FNV=48ad5090055eeeae.
```

Completeness is literal: for each of the 432 rows, Gosper enumeration visits
the full rank-nine universe from `(1<<9)-1` below `1<<30`, with the checked
count `14,307,150`.  Activity is recomputed exactly.  An active disjoint joint
mask closes a body; otherwise every active nonjoint mask is tested and failure
means zero disjoint hits.  Hence zero failures covers every labelled body on
the fixed 432-row target, not merely the bodies changed by the exchange.

## Typed consequence

Consuming the 16 endpoint-593 rows into THM-4310's audited partition gives

```text
typed union: 2,052 rows,
FNV=7c3f8bda9c37c5c2,
SHA256=ed15d387429df70be91a929e3f5821732cff05e1a0554ccc02c44c5ee1cf65c7;

residual: 20,595 rows,
FNV=0a9a532153a5e6dc,
SHA256=46339928143a508538adcfd7b24c3b0b208175246769f03a53a4acb62af22f64.
```

The residual maximum drops to endpoint 592 on 35 rows, FNV
`3eb23833c35b9266`, SHA-256
`74e0350a123ef06524e327a91d9cbea0c58eb0054fce176d426b2cdc1f9352dc`.

## Controls and scope firewall

- Baseline, singleton-quotient, and final-exchange O2/O3 artifacts agree
  byte-for-byte.
- Complement generation independently reproduces the complete response
  census for the unique failure.
- The original carrier failure, the response's unique private obligation, and
  zero strict-inactive capacity are hostile controls.
- `verify_endpoint593_packet.py` independently checks target seams, artifact
  equality, the hostile witness, the safe deletion, the zero-failure exchange,
  and the typed partition.
- The response minimum is only over rank-eight/rank-nine masks for this frozen
  failure.  The deletion boundary concerns one nonjoint deletion from `C^+`
  with all joint masks held fixed.  It is not a simultaneous-deletion theorem
  and does not prove a globally minimum carrier.
- Everything is finite-exact on the inherited labelled 30-speed pool and the
  fixed 432 rows.  No arbitrary-pair result, physical entry, terminating
  descent, or proof of LRC(14) follows.
