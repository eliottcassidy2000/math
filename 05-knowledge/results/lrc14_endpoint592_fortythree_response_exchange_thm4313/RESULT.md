# THM-4313: endpoint-592 forty-three-response size-preserving exchange

**FINITE-EXACT, relative to THM-4311. The fixed-pool result below is proved
by complete enumeration, explicit response certificates, and a direct
6.68-billion-case exchange replay. It does not give a physical entry or prove
LRC(14).**

## Fixed objects and hostile endpoint-592 boundary

Start from THM-4311's exchanged carrier

```text
C_592^old: size 3,925, rank8 3,857, rank9 68,
FNV=c9e5faef52ca5707, all 421 joint masks retained.
```

The inherited endpoint-592 layer has 35 rows, ordered FNV
`3eb23833c35b9266`. Independent O3 and O2 builds reconstruct the carrier and
scan all

```text
35*binom(30,9)=500,750,250
```

row-body cases. The outputs are byte-identical. Exactly 2,468 bodies fail on
seven rows:

| q | 96 | 100 | 105 | 192 | 210 | 256 | 416 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| failures at (q,592) | 13 | 3 | 48 | 1 | 288 | 2,101 | 14 |

The ordered obligation FNV is `2209b8d6760280cc`; the pair-ledger FNV is
`a3c8a93c0ea1ff45`. This is the hostile control closed below.

## Complete two-rank response census and multiplicity

For a frozen failure ((q,592,B)), a rank-eight or rank-nine mask responds
exactly when it is active at ((q,592)) and disjoint from (B). Complete
enumeration gives

| rank | responding masks | responder FNV | maximum total response |
|---:|---:|---:|---:|
| 8 | 188,462 | `bb1ff3319125eaae` | 409 |
| 9 | 2,205,776 | `ec506c87a74e64e1` | 536 |

Every obligation has a response. The exact per-obligation multiplicity ledger,
FNV `1a9370cc4b793710`, ranges from 423 to 34,654 total responders. Six
(q=256) bodies have no rank-eight response; this alone does not determine a
large response-cover obstruction.

An explicit pair-tagged packing supplies the rigorous full-universe lower
bound. It consists of

```text
(105,592,06d81088)
(256,592,0530e401) (256,592,06167400)
(256,592,07246088) (256,592,0f141208)
(256,592,13116408) (256,592,1430240b)
(256,592,16624401)
```

The certificate-only verifier scans all
(binom(30,8)+binom(30,9)=20,160,075) masks and finds zero co-responses
for all 28 pairs. Hence any rank-eight/rank-nine response cover of the frozen
2,468-obligation universe has at least eight masks. No maximum-packing or
exact-cover claim is made.

## A 43-mask cover and the retained-pool boundary

The deterministic retention rule keeps the best 20,000 masks by total
response, the best 3,000 per failed row, and a first-response sidecar for each
obligation. Their union has 37,497 masks. Optimization on this finite pool
produced an explicit cover

```text
A_592: 43 masks, rank8 2, rank9 41,
ordered FNV=ca3cb80f471f2e7e.
```

A direct arithmetic replay, separate from the optimizer's exported matrix,
recomputes all responses and obtains 7,842 incidences, all 2,468 obligations
covered, and zero missing. The compact retained-pool mask ledger has 37,497
distinct masks, ordered FNV `557039eeb8db4ed4`; direct set membership places
all 43 cover masks in that ledger.

An integer LP-dual certificate with denominator (10^9) has numerator
`35261737295`; every retained-pool column has score at most (10^9).
Consequently the optimum restricted to this particular 37,497-mask pool lies
in

```text
36 <= retained-pool optimum <= 43.
```

This lower bound does **not** extend to the complete responder universe. Exact
pricing of the same integer dual over all responding masks finds 20,986
violations: 508 rank-eight and 20,478 rank-nine omitted masks. The maximum
score is `3849498499` at mask `20188724`, and the violation-ledger FNV is
`9f5c2aabbea12034`. Thus the only proved full-universe interval here is

```text
8 <= complete response-cover minimum <= 43.
```

## Exact singleton boundary after adjoining the cover

Let (C^+=C_{592}^{old}cup A_{592}). Then

```text
|C^+|=3,968, FNV=b8221545b2d668d0.
```

The preservation target is the disjoint union

```text
K_467 = THM-4309's 391 rows
        union THM-4310's complete 25 endpoint-594 rows
        union THM-4311's 16 endpoint-593 rows
        union the 35 endpoint-592 rows,
|K_467|=467, ordered FNV=2d6aa988098aa5eb.
```

Adding masks cannot create a new singleton on an inherited row. The exact
quotient therefore filters THM-4311's complete 1,520-singleton ledger and
raw-scans all (35*binom(30,9)=500,750,250) new-row bodies. The cover resolves
852 inherited singletons and leaves 668. Together with the new layer:

```text
private obligations                       1,510
private-obligation FNV         4dfcadbd1d5c0c31
protected masks                              412
  protected old-carrier masks                369
  protected additions                         43
protected-mask FNV              6b6559e5bf1a529d
singleton-safe old nonjoint masks          3,135
safe-ledger FNV                 bd2d6c6eabf2203b
```

The deletion ledger takes the first 43 safe masks in the deterministic order
`(active-row count, rank, mask)`:

```text
D_592: 43 masks, ordered FNV=0dd14ef0fe3eec62.
```

The singleton quotient proves each member individually deletable from
(C^+). It makes no simultaneous-deletion inference.

## Size-preserving exchange and complete raw closure

Define

```text
C_592^new=(C_592^old \ D_592) union A_592.
```

Then

```text
|C_592^new|=3,925, rank8=3,818, rank9=107,
FNV=a0d08a38c10bdab7, all 421 joint masks retained.
```

A separate raw program reconstructs the carrier from canonical component
ledgers and replays every one of

```text
467*binom(30,9)=6,681,439,050
```

row-body cases. Fresh O3 and O2 builds give byte-identical transcripts,
467-row pair ledgers, and empty failure ledgers:

```text
joint-exposed bodies                 2,218,769
nonjoint hit incidences             74,178,388
failures                                     0
pair-ledger FNV          333e401fef6c7240
```

This direct simultaneous replay—not the singleton quotient—proves the
43-for-43 exchange.

## Typed consequence

Consuming the 35 endpoint-592 rows into THM-4311's exact partition gives

```text
typed union: 2,087 rows,
FNV=23e4136827b770a5,
SHA256=2d19070a08e6dd0b505070bea648cdcd7d6b8e3da1df18633ddc749fed8553ef;

residual: 20,560 rows,
FNV=8d797592a729e0b3,
SHA256=c23c20fb5552c12a9345ea34fba192e7817a3b332927ede639fefc81532a278a.
```

The residual maximum drops to endpoint 591 on 13 rows, FNV
`fc332c0697c671c7`, SHA-256
`b7c0cd2324bd4e7cb0819c30806faded369f5eba83f82e51d618be41e24bd211`.

## Controls and scope firewall

- Baseline, multiplicity, packing, direct-cover, singleton, and final-exchange
  O2/O3 artifacts match byte-for-byte.
- The packing lower bound uses the complete rank-eight/rank-nine universe.
  The stronger lower bound 36 is explicitly restricted to the retained pool;
  complete pricing supplies a hostile counter-certificate to globalization.
- All 43 additions are protected after adjoining the cover. The selected 43
  deletions come from old carrier masks only and retain all joint masks.
- Singleton safety was never promoted to simultaneous safety; the full
  6.68-billion-case replay is the simultaneous certificate.
- Everything is finite-exact on the inherited labelled 30-speed pool and the
  fixed 467 rows. There is no globally minimum carrier or response-cover
  theorem, arbitrary-pair theorem, physical entry, terminating descent, or
  proof of LRC(14).
