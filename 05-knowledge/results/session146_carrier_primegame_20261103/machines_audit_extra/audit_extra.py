#!/usr/bin/env python3
"""Independent transition/translation audit; reads ../machines/check.py."""
import hashlib, importlib.util, itertools, json
from pathlib import Path
p = Path(__file__).resolve().parent.parent / 'machines' / 'check.py'
spec = importlib.util.spec_from_file_location('checked_machine', p)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def need(condition, info):
    if not condition:
        raise RuntimeError(str(info))

traces=0
for n in range(1, 7):
    for t in range(4*n+5):
        for v in itertools.product([0,1], repeat=n):
            cells=[('_',None)]+[(x,None) for x in v]+[('_',None),('_',None)]
            cells[1]=(cells[1][0],())
            for k in range(n+1):
                heads=[(i,x[1]) for i,x in enumerate(cells) if x[1] is not None]
                need(heads==[(k+1,v[:k])], ('head',n,t,v,k,heads))
                shifted=[('_',None)]*3+cells+[('_',None)]*3
                a=m.ca_step(n,t,cells); b=m.ca_step(n,t,shifted)
                need(b[3:3+len(cells)]==a, ('translation',n,t,v,k))
                cells=a
            exp=next(((i,j) for i in range(n) for j in range(i,n)
                      if 2*i+2*j+6==t and v[i] and v[j]),None)
            marker=('R',) if exp is None else ('A',)+exp
            need(cells==[('_',None)]+[(x,None) for x in v]+[(marker,'H'),('_',None)],
                 ('output',n,t,v,cells,marker))
            traces+=1

n,t=2,8
states=[None,'H','E',(),(0,),(1,),(0,0),(0,1),(1,0),(1,1)]
syms=['_',0,1,('R',),('A',0,0),('A',0,1),('A',1,1),'#']
local=0
for ls,lq,s,q in itertools.product(syms,states,syms,states):
    cells=[('_',None),(ls,lq),(s,q),('_',None),('_',None)]
    a=m.ca_step(n,t,cells)
    b=m.ca_step(n,t,[('_',None)]*2+cells+[('_',None)]*2)
    need(a[1:4]==b[3:6], ('local_translation',cells))
    local+=1

def trial(x):
    return x>=2 and all(x%d for d in range(2,x) if d*d<=x)

stream=[]
for n in range(6,11):
    for name in ('owner','guy'):
        bits,steps=m.primegame_bits(n,name)
        need(bits==[trial(2*i+3) for i in range(n)], ('stream',n,name,bits))
        stream.append({'N':n,'list':name,'steps':steps})
print(json.dumps({'status':'PASS', 'checked_machine_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
  'single_head_traces_N_1_to_6':traces,'local_neighborhoods':local,
  'stream_checks':stream},indent=2,sort_keys=True))
