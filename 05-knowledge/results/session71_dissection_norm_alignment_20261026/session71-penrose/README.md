# Session71 Penrose and boundary-carrier lane

Read `REPORT.md` for theorem statements P1--P5, exact geometric witnesses,
primary literature and scope boundaries. Repository inputs are frozen at
`c2039c32f9726b09ed0e09d5f5c126fc72975897` in `SOURCE_PINS.json`.

## Reproduce

Python3, standard library only, from this directory:

```sh
python3 verify_packet.py
```

The verifier reruns `check_carriers.py` normally and with `-O`, requires
byte-identical output, checks the stored receipts, checks optimization-safe
injected failure exits, and verifies every local source-snapshot hash.
It does not access the network or mutate the repository.

Evidence scopes:

- **All-order proofs:** literal whole-side diagram/involution equivalence;
  relative-frame double cosets; native finite marked-tile conjunction guards;
  Euclidean development iff cycle consistency; abelian-field equal-area
  kite/dart obstruction.
- **Finite exact controls:** 6,574 involution-pair roundtrips through six
  side occurrences; `S4` matching double cosets; exact `Q(phi)` kite/dart
  placement; parity laws through seven independent gadgets; five-square
  nonoverlap hostile. These are not a Dudeney case enumeration.
- **Independent audit:** sibling `session71-matching-audit/check_penrose.py`
  independently uses `Q(sqrt(5))`, checks all angles and polygon simplicity,
  and checks 1,364 actual state/subset placements through five gadgets.
  It was replayed normally and with `-O` and matched its receipt exactly.
- **Injected failure:** the `--inject-fault` option is an explicit failing
  semantic gate demonstrating checks remain active under `-O`; it is not
  a claim of general artifact corruption coverage.

No source PDF is redistributed. `EXTERNAL_SOURCES.json` gives primary URL,
version/DOI, inspected location and hashes of retrieved bytes. P2 markings
were visually read from Penrose's printed page35, Fig9.

The parent owns publication. This lane made no shared checkout edits or
commits. Local named-tile attachment success is not infinite Penrose tiling
extendibility; geometric development is not global nonoverlap.
