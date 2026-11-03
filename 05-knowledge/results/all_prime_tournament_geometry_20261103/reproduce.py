#!/usr/bin/env python3
"""Normal/-O replay; receipts exclude wall-clock and platform noise."""
import json
from pathlib import Path
import subprocess
import sys

root=Path(__file__).resolve().parent
records=[]
for script,output in [('check.py','receipt.json'),('audit.py','audit.json')]:
    a=subprocess.check_output([sys.executable,'-B',str(root/script)],cwd=root)
    b=subprocess.check_output([sys.executable,'-B','-O',str(root/script)],cwd=root)
    if a!=b:
        raise RuntimeError(('optimized divergence',script))
    # check.py's output is the input receipt of the subsequent audit.
    (root/output).write_bytes(a)
    records.append({'script':script,'normal_optimized_identical':True,'valid_json':isinstance(json.loads(a),dict)})
print(json.dumps({'status':'PASS','runs':records},indent=2))
