"""Recompute the six receipts; explicit guards also execute under python -O."""
from pathlib import Path
import importlib
import json

HERE = Path(__file__).resolve().parent
total = 0
for name in ('charged', 'jets', 'memory', 'renewal', 'frobenius', 'rooks'):
    result = importlib.import_module(name).main(write=False)
    actual = json.loads(json.dumps(result))
    expected = json.loads((HERE / (name + '.json')).read_text(encoding='utf-8'))
    if actual != expected:
        raise RuntimeError('Receipt differs: ' + name)
    total += result['checks']
print(json.dumps({'status': 'PASS', 'receipts': 6, 'explicit_checks': total}))
