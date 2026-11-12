# Reproduction and validation scope

Python3 standard library only. Run these from the repository root; paths in
the scripts are derived from their own location. No external exact-count cache
is needed by either complete order64 producer/consumer.

```bash
P=05-knowledge/results/session279_interval_tightening_20261112
python3 "$P/enumeration/backbone.py"
python3 "$P/enumeration/second_jet.py"
python3 "$P/enumeration/integer_jet.py"
python3 "$P/enumeration/divisor_jet.py"
python3 "$P/marked/marked_graph.py"
python3 "$P/marked/marked_graph.py" --order64
python3 "$P/cross/ground_moments.py" --experiment
python3 "$P/cross/ground_moments.py" --consume "$P/cross/packet_n64_B3.json"
python3 "$P/cross/ground_moments.py" --consume "$P/cross/packet_n64_B11.json"
python3 "$P/physical/verify.py"
python3 "$P/audit/check_independent.py"
```

Most scripts print JSON. Timing fields vary; mathematical fields must agree
with their recorded results. The first and second jet test harnesses use full
partitions as an independent reference after computing their coefficient
outputs. The marked `--order64` producer does not. The ground-moment consumer
uses promised moments only; a separate producer currently builds those moments
by partitions. The independent audit uses normal Python with active assertions;
an optimized run is not a replacement for those checks.

The full repository documentation checker was attempted in this sparse
checkout. It requires omitted startup/hypothesis/canon files and therefore did
not provide a whole-repository pass. No existing truth-surface file was changed
except an appended method card in META-PATTERNS and the new S279 route. Local
packet links were separately resolved against the repository tree. This is a
research packet validation, not a re-audit of inherited mathematical claims or
a full repository CI claim.

The independent new-claim review and its precise controls are in
[audit/REVIEW.md](audit/REVIEW.md). Source and output hashes are retained there.
