#!/usr/bin/env python3
"""Point queries in the inherited seven-parent bank; retains physical words."""
import argparse
import json
from pathlib import Path
from source_inverse import Inverse, f

HERE=Path(__file__).resolve().parent
BANK=HERE.parent/'square_sum_spectrum_overnight_20261031c'


def profiles():
    rows=[]
    for name in ['D','DD']+[f'DT{k}D' for k in range(1,6)]:
        for line in (BANK/f'{name}_profiles.tsv').read_text().splitlines():
            rows.append((name,*map(int,line.split())))
    return rows


def block(name,mask):
    if name=='D':return [('D',mask)]
    k=0 if name=='DD' else int(name[2:-1])
    return [('D',mask&15)]+[('S',(mask>>(4+j))&1) for j in range(k)]+[('D',mask>>(4+k))]


def suffix(y):
    return [('S',1-((y>>j)&1)) for j in range(y.bit_length())]


def query(n,cap=56):
    solver=Inverse();seen=set();queries=0
    y=solver.solve(n,0,0,1)
    if y is not None:
        return {'target':n,'parent':'singleton','x':0,'y':y,'word':suffix(y),'queries':1,'nodes':solver.nodes}
    for x in range(cap+1):
        for name,h,e,p,q,s,t,F,mask in profiles():
            initial=h*f(x)+p*x+(h-s)*(x%2)+F-h
            if initial>n:continue
            Q=e*x+q;g=h-t
            key=(initial,Q,g,h)
            if key in seen:continue
            seen.add(key);queries+=1
            y=solver.solve(n-initial+h,Q,g,h)
            if y is None:continue
            if initial+Q*y+g*(y%2)+h*(f(y)-1)!=n:
                raise RuntimeError('bad scalar witness')
            word=([('S',int(b)) for b in bin(x)[2:]] if x else [])+block(name,mask)+suffix(y)
            return {'target':n,'parent':name,'mask':mask,'x':x,'y':y,'word':word,
                    'source':[initial,Q,g,h],'queries':queries,'nodes':solver.nodes}
    return {'target':n,'x':None,'cap':cap,'queries':queries,'nodes':solver.nodes}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('targets',type=int,nargs='+');parser.add_argument('--cap',type=int,default=56)
    args=parser.parse_args()
    if args.cap<0 or any(t<1 or t%2==0 for t in args.targets):raise ValueError('positive odd targets and nonnegative cap required')
    for t in args.targets:
        print(json.dumps(query(t,args.cap),sort_keys=True),flush=True)
