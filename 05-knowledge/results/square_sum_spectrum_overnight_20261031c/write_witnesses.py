#!/usr/bin/env python3
"""Reconstruct literal adjacency without importing the inverse or transfers."""
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031b'))
from audit_inverse_witnesses import adjacency

folder=HERE/'witnesses';folder.mkdir(exist_ok=True)
receipt=[]
for record in json.loads((HERE/'guarded_targets.json').read_text())['results']:
 if record['word'] is None:raise RuntimeError('missing positive witness')
 rows=adjacency(record['word']);N=len(rows);F=record['target']
 if N!=record['vertices']:raise RuntimeError('changed vertex clock')
 name=f'H{F}.txt'
 (folder/name).write_text('\n'.join([f'{N} {F}']+[''.join(str(r>>j&1) for j in range(N)) for r in rows])+'\n')
 receipt.append({'target':F,'vertices':N,'adjacency_file':name,'rows':rows})
print(json.dumps({'status':'reconstructed; require separate adjacency audit','witnesses':receipt},indent=2,sort_keys=True))
