"""Reproduce all declared universes and compare the saved receipts exactly."""
from pathlib import Path
import json,importlib
HERE=Path(__file__).resolve().parent
for name in ('ports','paley','counting','arithmetic','crt','audit'):
 result=importlib.import_module(name).main(write=False)
 actual=json.loads(json.dumps(result));expected=json.loads((HERE/(name+'.json')).read_text())
 if actual!=expected:raise RuntimeError('receipt mismatch: '+name)
print('All six receipts agree exactly. Explicit guards remain active under -O.')
