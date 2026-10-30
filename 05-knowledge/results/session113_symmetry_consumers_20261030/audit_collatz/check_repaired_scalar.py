#!/usr/bin/env python3
"""Small-precision scalar Fraction audit, independent from array checker."""
from fractions import Fraction as F
from pathlib import Path
import struct,json,subprocess,sys

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def height(z,E):
 if z==0:return E
 k=0
 while z%5==0:k+=1;z//=5
 return k
rs=[]
for E in [1,2,3]:
 q=5**E;D=8*q;lam=F(5,8)*(1-F(1,q));V=struct.unpack('<'+'q'*(5*q),Path(f'/tmp/repairedV{E}.bin').read_bytes());edges=0
 for u in range(1,5):
  for pending in [False,True]:
   if pending and u!=4:continue
   for z in range(q):
    if pending and z%5:continue
    idx=(4 if pending else u-1)*q+z
    need(V[idx]!=-(1<<60),'finite domain')
    for a in range(1,65):
     nu=3*u*pow(2**a,-1,5)%5;nz=(3*z+2**a-2)*pow(2**a,-1,q)%q;h=height(nz,E)
     npending=False;reward=F()
     if a!=1:
      if nu==1:reward=h
      elif nu==4 and h>0:npending=True
     elif pending:reward=h
     elif nu in [1,4] and h>0:
      reward=int(nu==1)+sum((F(1,2*5**(e-1)) for e in range(2,h+1)),F())
     dest=(4 if npending else nu-1)*q+nz;cost=D*(reward-lam*a)
     need(cost.denominator==1,'integer cost')
     need(V[dest]>=V[idx]+cost,('edge',E,u,pending,z,a))
     if z==0 and not pending and u in [1,4] and a>1:need(V[dest]>=cost,'entry')
     edges+=1
 rs.append(dict(E=E,scalar_edges=edges))
# Hostile certificate changes an ordinary valid state's value to positive.
original=bytearray(Path('/tmp/repairedV1.bin').read_bytes());original[:8]=struct.pack('<q',1)
bad=Path('/tmp/repaired_corrupted.bin');bad.write_bytes(original)
r=subprocess.run([sys.executable,'-O',str(Path(__file__).with_name('audit_repaired.py')),'1',str(bad)],capture_output=True,text=True)
need(r.returncode!=0 and 'nonpositive maximum' in r.stderr,'corruption was not rejected under -O')
report=dict(scalar_checks=rs,corrupted_vector_rejected_under_optimization=True)
Path(__file__).with_name('repaired_scalar_controls.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
