"""Replay each scoped checker with and without optimization; freeze receipts."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys


root = Path(__file__).resolve().parent
checks = (
    ("independent_check.py", "independent_check.json"),
    ("arithmetic_check.py", "arithmetic_check.json"),
    ("square_counts.py", "square_counts.json"),
    ("reset_composition_check.py", "reset_composition_receipt.json"),
)
manifest = {}
for source, receipt in checks:
    outputs = []
    for flags in ([], ["-O"]):
        completed = subprocess.run([sys.executable, "-B", *flags, str(root/source)],
                                   check=True, capture_output=True, text=True, encoding="utf-8")
        outputs.append(completed.stdout.replace("\r\n", "\n"))
    if outputs[0] != outputs[1]:
        raise RuntimeError(f"normal/optimized mismatch: {source}")
    parsed = json.loads(outputs[0])
    # Canonicalize receipt bytes independently of host newline conversion.
    payload = json.dumps(parsed, indent=2, sort_keys=True) + "\n"
    (root/receipt).write_text(payload, encoding="utf-8", newline="\n")
    manifest[source] = {
        "source_sha256": hashlib.sha256((root/source).read_bytes()).hexdigest(),
        "receipt": receipt,
        "receipt_sha256": hashlib.sha256(payload.encode()).hexdigest(),
        "normal_optimized_stdout_identical": True,
    }
    print(f"PASS {source}: normal and optimized exact receipts agree")
(root/"validation_manifest.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
