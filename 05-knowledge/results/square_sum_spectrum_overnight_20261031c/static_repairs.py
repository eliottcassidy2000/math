#!/usr/bin/env python3
"""Retain literal small repairs for the four extra restricted-grammar misses."""
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;OLD=HERE.parent/'square_sum_spectrum_overnight_20261030'
sys.path.insert(0,str(OLD))
from audit_feedback import decode,tournament,subset_dp,require
source=json.loads((OLD/'two_feedback_100000.json').read_text())
records=[];folder=HERE/'witnesses';folder.mkdir(exist_ok=True)
for target in (43,93,173,247):
 h,n,number=next(r for r in source['witnesses'] if r[0]==target)
 word=decode(n,number);rows=tournament(word);actual=subset_dp(rows)
 require(actual==target,'literal small repair')
 name=f'H{target}.txt';size=len(rows)
 (folder/name).write_text('\n'.join([f'{size} {target}']+[''.join(str(r>>j&1) for j in range(size)) for r in rows])+'\n')
 records.append({'target':target,'two_feedback_word':word,'vertices':size,'rows':rows,'literal_count':actual,'file':name})
print(json.dumps({'status':'passed','repairs':records,'scope':'Four actual small tournaments; combined with the audited finite arithmetic image only.'},indent=2,sort_keys=True))
