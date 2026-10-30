#!/usr/bin/env python3
"""Regression for the normalized parity repair of THM-722."""
from fractions import Fraction as F
from pathlib import Path
import json,runpy
root=Path(__file__).resolve().parents[3]
m=runpy.run_path(str(root/'04-computation/lrc14_leader_ledger_boxeph_S21.py'))
def need(b,msg):
    if not b:raise RuntimeError(msg)
positive=[]
for S in ([1,2],[2,4],[2,3],[4,6],[2,4,6],[3,6,9],[4,8,12],[2,5,8]):
    r=m['leader_ledger'](S)
    need(r['identity_ok'] and r['M']<F(1,2),'no-wrap ledger')
    positive.append(dict(speeds=S,M=str(r['M']),Hplus=r['Hplus']))
hostile=[]
for S in ([1,3],[2,6],[3,9],[6,18],[6,10]):
    try:m['leader_ledger'](S)
    except ValueError:hostile.append(S)
    else:raise RuntimeError('all-odd normalized core accepted')
# A small arithmetic 722 control; not used as a proof bridge.
need(722==2*19**2 and sum(d for d in range(1,723) if 722%d==0)==9*127,'aliquot identity')
need(pow(19,3,127)==1 and 19%127!=1,'exact order three')
out=dict(status='regression of corrected assumption, not new LRC coverage',
         accepted=positive,rejected=hostile,sigma722=1143,aliquot722=421)
Path(__file__).with_name('ledger_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print('ledger regression:',len(positive),'accepted,',len(hostile),'correctly rejected')
