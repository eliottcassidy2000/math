#!/usr/bin/env python3
"""Compare normal/-O mathematical receipts while excluding timing and mode."""
from pathlib import Path
import hashlib
import json

base = Path(__file__).resolve().parent

def normalize(name):
    record = json.loads((base/name).read_text())
    record.pop('optimized_python', None)
    for row in record.get('certificates', []):
        row.pop('seconds', None)
    return record

normal=normalize('exact_precision_results.json')
optimized=normalize('exact_precision_results_optimized.json')
if normal != optimized:
    raise RuntimeError('Normal and -O mathematical certificate receipts differ')
if json.loads((base/'hostile_results.json').read_text()) != json.loads((base/'hostile_results_optimized.json').read_text()):
    raise RuntimeError('Normal and -O hostile-test receipts differ')
blob=json.dumps(normal, sort_keys=True, separators=(',', ':')).encode()
receipt={
    'status':'PASS',
    'normal_and_optimized_mathematical_receipts_equal':True,
    'hostile_receipts_equal':True,
    'precision_inputs':[row['p'] for row in normal['certificates']],
    'hostile_test_count':len(normal['hostile_controls']),
    'normalized_sha256':hashlib.sha256(blob).hexdigest(),
}
(base/'normalized_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps(receipt, indent=2))
