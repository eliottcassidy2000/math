# Balanced fixed-path census: independent audit

Start with `AUDIT.md` for the proof audit, contribution boundary, prior-art
correction, and finite scopes. `receipt.json` seals the replays. `MANIFEST.json`
identifies every included executable and result; run `python3 verify_manifest.py`.

Independent methods are in `audit.py` and `parent_cards.py`. The producer's
immutable manifest-referenced actual copy is in `producer/`. No producer code
is imported to obtain the independent scalar counts. `compare.py` compares the
two already generated receipts and additionally verifies the divisor formula
against independently constructed rotation edge-orbits. `hostile.py` and
`hostile_packet.py` reject mathematical and receipt-tampering near misses.

Normal, optimized, and relocated results are retained. The optional extra
N=15 rotation fixed sets do not compute or claim a full n=14 census.
There is no new unmarked tournament class count or full metagraph reconstruction.
