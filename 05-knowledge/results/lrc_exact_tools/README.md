# LRC Exact Tools State

Snapshot from `agent-v1410-codex-creative` on V1410-1.

## Verifier

Fingerprint:

```text
lrc_exact_tools.py lrc_exact_tools_v1 cc1ef7bad927e06117e0f04e
```

The fingerprint is the blake2b-12 hash of `lrc_exact_tools.py` bytes. Do not
edit that file without announcing a new fingerprint.

## Local Verification

Run from repository root:

```bash
PYTHONPATH=04-computation python3 04-computation/verify_certificate.py 05-knowledge/results/lrc_exact_tools/runs/certs/AP13_lrc_n14.json
PYTHONPATH=04-computation python3 04-computation/verify_certificate.py 05-knowledge/results/lrc_exact_tools/runs/certs/Vstar_lrc_n14.json
PYTHONPATH=04-computation python3 04-computation/verify_certificate.py 05-knowledge/results/lrc_exact_tools/runs/certs/B24_rank6_beta13_margin4_165.json
```

Observed results:

| cert | result | maximin | threshold | witnesses | candidates |
| --- | --- | --- | --- | --- | --- |
| `AP13_lrc_n14.json` | pass | `1/14` | `1/14` | 6 | 213 |
| `Vstar_lrc_n14.json` | pass | `1/14` | `1/14` | 6 | 437 |
| `B24_rank6_beta13_margin4_165.json` | pass | `1/11` | `1/15` | 14 | 491 |

## Minimal Sample

`runs/certs/minimal_contact_rank_sample.json` is a compact n=14 certificate
with the requested fields:

- `contact_graph`
- `rank_signature`
- `contact_matrix`
- `exact_candidate_format`

It uses speeds `[1,2,3,4,5,6,7,8,9,10,11,14,19,24]`, threshold `1/15`,
maximin `2/25`, witnesses `2/25,23/25`, and rank signature
`F2=1,beta1=1,contacts=2,edges=4`.

## Pilot Reproduction

Command shape:

```bash
PYTHONPATH=04-computation python3 04-computation/lrc_exact_tools.py pilot \
  --run-id 20260605Tpilot-v1410-B24-rankfields \
  --moving 14 \
  --n-total 15 \
  --delta 1/15 \
  --max-speed 24 \
  --shard-index 0 \
  --shards 256 \
  --time-budget-s 30
```

Shard params:

```json
{
  "moving_runners": 14,
  "n_total": 15,
  "delta": "1/15",
  "max_speed": 24,
  "shard_index": 0,
  "shards": 256,
  "time_budget_s": 30.0,
  "canonicalization": "sorted primitive gcd-normal form",
  "pruning": "stable blake2b shard filter; gcd primitive filter"
}
```

Pilot summary: processed 280 primitive rows in 30.03s at 9.32 rows/s; no
counterexamples; best margin `1/75`; best high-structure cert
`runs/certs/B24_rank6_beta13_margin4_165.json` has `F2=6`, `beta1=13`,
`D=[11,22]`, margin `4/165`.
