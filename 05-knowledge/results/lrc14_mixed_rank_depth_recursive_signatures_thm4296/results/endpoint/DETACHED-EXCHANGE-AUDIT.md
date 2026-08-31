# Incoming-signal probe of the THM-4287 endpoint-636 hostile layer

Status: `FINITE-EXACT SCRATCH RESULT`; not canon and not an LRC(14) proof.

## Result

Within the complete rank-eight universe, with activity typed separately at
`(100,636)` and `(256,636)`, the 64 and 37 frozen THM-4287 failure bodies have
exact append-only cover minima 8 and 6.  The combined 101-coordinate consumer
has exact minimum 14: no response mask can be shared so as to reduce the sum.
The 14 masks printed in `probe.out`, appended in that order to the maintained
9,006-mask carrier, give a 9,020-mask carrier with FNV
`aad14d1dcf2b3ebf`.  It covers all 101 hostile obligations.  Since the old
carrier already covers the other four endpoint-636 rows, this is a candidate
complete-layer carrier extension, subject to independent literal audit and
proof-graph integration by the root session.

The response quotient has 1,835 classes, ledger FNV `b6a561c395058684`, and
259 inclusion-maximal patterns.  The selected witness-response ledger has FNV
`2a4cd3ffe7eeca5f`.

## Detached audit of the zero-net exchange

The endpoint agent independently replaced 14 boundary-inactive nonjoint masks
by a different exact 14-mask cover, retaining size 9,006 and producing carrier
FNV `8062ce6d5728da1f`.  `detached_exchange_audit.cpp` wraps the maintained
source-independent literal-wall engine and verifies:

- all 14 deletions have strictly negative margin at all nine endpoint-637/636
  rows: 126 strict cells, zero equalities, ledger `640a94af980de74a`;
- the 14 additions' literal margins and response words cover all `64+37`
  hostile bodies, response-margin ledger `6c7b10e7ea30e546`;
- exhaustive replay of `9*C(30,9)=128,764,350` body instances has zero
  failures and zero activity equalities, boundary ledger `3bf60326e0da05c8`;
- every row's active-mask count and FNV agree with the endpoint agent's
  independently written scan.

O3 and O0 builds emit byte-identical `detached_exchange_audit.out`, SHA-256
`4c0dd784d05bdafed4581562d262f93d2adc1049e7f6debe0461a9dbc1f2892c`.

## Incoming-operation transfer

- Source: THM-4289's differentiated blowdown obstruction.
- Target/map: change endpoint `637 -> 636` at fixed first coordinate, and map
  each carrier mask to its exact normalized activity slack on the two layers.
- Preserved predicate: sign of exact wall activity at each fixed row.
- Destroyed data: body incidence and cross-row response; these must be retained
  as the 101-bit consumer sidecar.
- Hostile: for `(100,636)`, 2 of the 37 closest current-carrier responders stay
  inactive at endpoint 637; therefore layer difference alone is not a response
  congruence.  The full exact response quotient is necessary.

For first-order body deletion, price an atom by the number of active current
carrier masks which become disjoint after that atom is deleted.  Every hostile
body has positive price, but 8 bodies have tied maximizing atoms.  Hence this
operation is a useful sensitivity/branching heuristic, not an owner map.

## Reproduction

From the repository root on Windows with GCC available:

```powershell
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. `
  .scratch/incoming_lrc_signal_agent/probe.cpp `
  -o .scratch/incoming_lrc_signal_agent/probe.exe
.scratch/incoming_lrc_signal_agent/probe.exe `
  05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/inputs/reconstructed_final8951.txt `
  05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/inputs/additions45.txt `
  .scratch/incoming_lrc_signal_agent/suffix10.txt `
  05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/results/primary/failures.csv `
  unused
```

The final unused positional argument is retained only to keep this scratch
driver's invocation stable while it evolved.  Compare stdout byte-for-byte
with `probe.out` after normalizing platform newlines.

For the detached exchange replay:

```powershell
g++ -std=c++20 -O3 -DNDEBUG -pthread -I. `
  .scratch/incoming_lrc_signal_agent/detached_exchange_audit.cpp `
  -o .scratch/incoming_lrc_signal_agent/detached_exchange_audit.exe
.scratch/incoming_lrc_signal_agent/detached_exchange_audit.exe `
  05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/inputs/reconstructed_final8951.txt `
  05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/inputs/additions45.txt `
  05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/inputs/endpoint638_response_witness9.txt `
  05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/inputs/joint421_masks.txt `
  05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/inputs/current_residual22682.csv `
  05-knowledge/results/lrc14_repaired_carrier_endpoint637_descent_thm4287/results/primary/failures.csv `
  4
```

Scratch SHA-256 identities:

```text
probe.cpp                      9accb5e1f5d39c6d4ecfb81b3ec82d7b4eb0a77e64fbd6f4caf79e8e93c3c811
probe.out                      cd41f48d3f536c9d570cecbce7e8cd3200c9ec5312d091c9ed94d3e064287555
independent_witness_audit.py   fe9f55e90ca8313f4ffaedc7e105e7dfa86fd7da860d95a95050e365966b17ed
detached_exchange_audit.cpp    4c5176be7b7e56cdfc248fa8b32fb5362741c7b12efe64766eb86d1c2bb02b18
detached_exchange_audit.out    4c0dd784d05bdafed4581562d262f93d2adc1049e7f6debe0461a9dbc1f2892c
suffix10.txt                   38a40877b3618ed6967d387fd7af4dd6ec3281e30fabdb65d1f73cab1ecc5e69
```
