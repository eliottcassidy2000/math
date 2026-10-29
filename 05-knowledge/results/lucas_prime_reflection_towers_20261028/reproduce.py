#!/usr/bin/env python3
"""Replay exact receipts without overwriting them."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

root=Path(__file__).resolve().parent
report=[]
for script,receipt in [('check.py','receipt.json'),('audit_alternate.py','audit_receipt.json')]:
    expected=(root/receipt).read_bytes()
    for optimized in [False,True]:
        command=[sys.executable]+(['-O'] if optimized else [])+[str(root/script)]
        actual=subprocess.check_output(command,cwd=root)
        if actual!=expected:
            raise RuntimeError(f'receipt mismatch: {script}, optimized={optimized}')
    report.append(dict(script=script,normal_and_optimized='PASS',
                       receipt_sha256=hashlib.sha256(expected).hexdigest()))
print(json.dumps(report,indent=2,sort_keys=True))
