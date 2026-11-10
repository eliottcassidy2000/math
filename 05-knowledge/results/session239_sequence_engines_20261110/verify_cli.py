#!/usr/bin/env python3
"""Validate convenience interface offsets and independently sourced OEIS prefixes."""
from generate import sequence
from pathlib import Path
import json
import subprocess
import sys

# Source pages checked 2026-11-10; integer sequences are factual data.
known={
 'A000568':(0,[1,1,1,2,4,12,56,456,6880,191536,9733056]),
 'A051337':(0,[1,1,0,1,1,6,35,353,6008,178133,9355949]),
 'A054946':(1,[1,0,2,24,544,22320,1677488,236522496,64026088576,33832910196480]),
 'A000088':(0,[1,1,2,4,11,34,156,1044,12346,274668,12005168]),
 'A000273':(0,[1,1,3,16,218,9608,1540944,882033440,1793359192848,13027956824399552,341260431952972580352]),
 'A001174':(1,[1,2,7,42,582,21480,2142288,575016219,415939243032,816007449011040]),
}
receipts=[]
for name,(offset,expect) in known.items():
    got=sequence(name,10)
    if got!=expect: raise ArithmeticError(('sequence or offset mismatch',name,got,expect))
    raw=subprocess.check_output([sys.executable,str(Path(__file__).with_name('generate.py')),'--sequence',name,'--through','10'])
    obj=json.loads(raw)
    if obj['offset']!=offset or [int(x) for x in obj['terms']]!=expect:
        raise ArithmeticError(('CLI mismatch',name))
    receipts.append({'sequence':name,'source':'https://oeis.org/'+name,'offset':offset,'through':10,'pass':True})
    z=sequence(name,0)
    if z != ([] if offset==1 else [1]):raise ArithmeticError(('zero horizon',name,z))
raw=subprocess.check_output([sys.executable,str(Path(__file__).with_name('generate.py')),'--through','20','--binary-precision','7'])
if list(map(int,json.loads(raw)['terms']))!=[a%128 for a in sequence('A000568',20)]:
    raise ArithmeticError('CLI residue mismatch')
large=json.loads(subprocess.check_output([sys.executable,str(Path(__file__).with_name('generate.py')),'--sequence','A054946','--through','200']))
if len(large['terms']) != 200 or len(large['terms'][-1]) <= 4300:
    raise ArithmeticError('large decimal serialization regression')
print(json.dumps({'large_decimal_serialization_through':200,'last_term_digits':len(large['terms'][-1]),'status':'PASS','sequence_checks':receipts,'empty_horizon_checks':6,'residue_prefix_through':20},indent=2))
