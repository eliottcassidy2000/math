#!/usr/bin/env python3
"""Separate ROOT from unresolved mass in the actual native annular rows."""
import hashlib
import json
from pathlib import Path
import numpy as np
from collatz_compact_resets_20261015f import Graph, row_norms, need

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'05-knowledge/results'


def run(k):
    X,b=1<<k,1<<(k//2)
    g=Graph(8*X)
    cols,inv=np.unique(g.endpoint[1:],return_inverse=True)
    stages=[]
    frozen=json.loads((DATA/'collatz_survivor_packets_20261015e.json').read_text())
    expected=next(r for r in frozen['first_reset'] if r['k']==k)
    for j in range(1,k+9):
        g.extend(cols)
        g.close(cols,inv)
        cols,inv=np.unique(g.endpoint[1:],return_inverse=True)
        norms,rows=row_norms(k,g,inv,len(cols),matrix=True)
        need(norms==expected['stages'][j-1]['row_norms'],'fresh native norms equal frozen universe')
        root_pos=int(np.flatnonzero(cols==1)[0])
        a=[int(row[root_pos]) for row in rows]
        z=[u-abs(r) for u,r in zip(norms,a)]
        R=sum(r*r for r in a)
        Z=sum(v*v for v in z)
        C=2*sum(abs(r)*v for r,v in zip(a,z))
        need(R+Z+C==sum(u*u for u in norms),'full root/tail/cross energy invoice')
        checked=g.internal_clock_audit()
        g.replay(range(1,8*X+1,max(1,X//64)))
        stages.append(dict(stage=j,root_square_numerator=R,
                           unresolved_square_numerator=Z,cross_numerator=C,
                           denominator=b,root_signed_rows=a,unresolved_row_norms=z,
                           max_native_clock=int(g.clock.max()),
                           checked_internal_edges=checked))
        if R+Z+C <49*b*X*X:
            break
    return dict(k=k,source_count=8*X,stages=stages)


if __name__=='__main__':
    result=dict(status='FINITE-EXACT native k12,14,16; no universal inference',
                scales=[run(k) for k in [12,14,16]])
    out=DATA/'paired_root_profiles_20261015g.json'
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    for rec in result['scales']:
        s=rec['stages'][-1];X=1<<rec['k'];den=s['denominator']*X*X
        print(rec['k'],s['stage'],{v:s[v]/den for v in ['root_square_numerator','unresolved_square_numerator','cross_numerator']})
    print('sha256',hashlib.sha256(out.read_bytes()).hexdigest())
