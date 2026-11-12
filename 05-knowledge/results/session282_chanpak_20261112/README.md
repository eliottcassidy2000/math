# Session282 research packet

[Full synthesis](SYNTHESIS.md). Chan–Pak remains OPEN.

## Proof lanes

- [Renewable nine-bank, infinite eighteen-count intervals, exact receiver fibers](renewal/PROOF.md)
- [Global fixed-receiver witness compression](compression/PROOF.md)
- [Rational carries, local action, central169 switch and Fibonacci fibers](wildcard/PROOF.md)
- [Parametric reroot escape from all powers of13](wildcard/ROOT_RESET.md)
- [Central convolution and exact absolute target carry](central/PROOF.md)
- [Fixed-template target bounds and signed Jordan extraction](targeting/PROOF.md)
- [Receiver tuner and its scalar-neutrality price](targeting/RECEIVER_TUNER.md)

## Reproduce

From the repository root, using standard Python3 (no external packages):

```bash
P=05-knowledge/results/session282_chanpak_20261112
python3 "$P/renewal/verify.py"
python3 "$P/compression/verify.py"
python3 "$P/central/verify.py"
python3 "$P/wildcard/verify.py"
python3 "$P/wildcard/fiber_extension.py"
python3 "$P/wildcard/reroot_transfer.py"
python3 "$P/targeting/verify.py"
python3 "$P/targeting/verify_tuner.py"
```

Repeat with `python3 -O`; no theorem check relies on a removable assertion.
The renewal producer writes `RESULTS.json` next to itself. Other producers
print their result JSON; compare with their saved result. central/compression
reuse the explicitly credited Session275 carrier for some source operations;
independent literal count paths and the audit retain different recurrences.

[Independent review](audit/REVIEW.md) and source hashes bind the final restored
files. [Recovery record](RECOVERY.md) explains the environment rollback; fresh
replays, not lost old receipts, support this release. Search-derived witnesses
are finite controls; the all-order statements are proved separately.
