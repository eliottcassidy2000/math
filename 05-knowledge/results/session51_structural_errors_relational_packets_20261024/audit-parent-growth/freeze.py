#!/usr/bin/env python3
"""Freeze only intentionally delivered sources, receipts, and dependency payloads."""
from pathlib import Path
import hashlib,json
B=Path(__file__).resolve().parent
names=['AUDIT.md','README.md','independent.cpp','incoming_trace.cpp','replay_independent.py','replay_producer.py','input_controls.py','mutations.py','reconcile_strata.py','verify_audit.py','replay.py','freeze.py','producer-pin.json','SOURCE_PINS.json','independent-replays.json','producer-replays.json','structural-replay.json','mutation-results.json','input-control-results.json','strata.json','strata-optimized.json','producer-verify-normal.json','producer-verify-optimized.json']
for label in ('O2','O3'):
 names += [f'independent-{label}.json',f'independent-{label}.out',f'incoming-{label}.json',f'incoming-{label}.out']
for label in ('O2','O3-ndebug-full'):
 names += [f'build-{label}.json',f'build-{label}.out',f'trace-{label}.json',f'trace-{label}.out']
for root,manifest in [('producer','MANIFEST.json'),('structural-audit','manifest.json')]:
 names.append(root+'/'+manifest);data=json.loads((B/root/manifest).read_text())['files'];names += [root+'/'+name for name in (data if isinstance(data,dict) else [item['path'] for item in data])]
names += ['source/session50-growth-THEOREM.md','source/session50-order-nine-THEOREM.md']
files={name:hashlib.sha256((B/name).read_bytes()).hexdigest() for name in sorted(set(names))}
(B/'manifest.json').write_text(json.dumps({'format':1,'files':files},indent=2)+'\n')
print(json.dumps({'manifest_sha256':hashlib.sha256((B/'manifest.json').read_bytes()).hexdigest(),'payload_files':len(files),'payload_bytes':sum((B/name).stat().st_size for name in files)},sort_keys=True))
