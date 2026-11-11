#!/usr/bin/env python3
"""Any finite number of synchronized fresh-return/birth steps; no regeneration claim."""
import json,argparse
from pathlib import Path
from ordinal_interfaces import profile,ordinal,selected_tail
A=(2,4,9,19,39,15);B=(2,4,57,3,11,27)
Q=(198,196,136,211,167,207,20,64)
def twin(T,v):
 ids=[i for i in range(len(T)) if i!=v]+[v,v];n=len(ids)
 return tuple(sum(1<<j for j in range(n) if i!=j and ((ids[i]==ids[j] and i<j) or (ids[i]!=ids[j] and T[ids[i]]>>ids[j]&1))) for i in range(n))
AT=twin(A,1);BT=twin(B,0)
def moment3(T):
 n=len(T);N=1<<n;ep=[[0]*n for _ in range(N)]
 for i in range(n):ep[1<<i][i]=1
 for mask in range(1,N):
  for i,z in enumerate(ep[mask]):
   if not z:continue
   bits=T[i]&~mask
   while bits:
    b=bits&-bits;bits-=b;ep[mask|b][b.bit_length()-1]+=z
 h=[sum(a) for a in ep];h[0]=1;ans=0;full=N-1
 for a in range(N):
  rem=full^a;b=rem
  while True:
   ans+=h[a]*h[b]*h[rem^b]
   if not b:break
   b=(b-1)&rem
 return h[-1],sum(h[a]*h[full^a] for a in range(N)),ans

def word(states,reverse=False):
 # 0 is an unused slot, 1/2 have consumed the A/B birth respectively.
 bank=(B if reverse else A,AT,BT);ls=list(states)
 if reverse:ls.reverse()
 T=()
 for state in ls:
  Z=bank[state]
  T=Z if not T else ordinal(T,0,Z,0)[0]
 return T

def substitute(left,right):
 modules=[(0,)]*8;modules[5]=left;modules[7]=right
 labels=[(i,j) for i,M in enumerate(modules) for j in range(len(M))]
 return tuple(sum(1<<b for b,(k,l) in enumerate(labels) if a!=b and ((i==k and modules[i][j]>>l&1) or (i!=k and Q[i]>>k&1))) for a,(i,j) in enumerate(labels))

def formula(states):
 h=q=r=1
 for state in states:
  x,y,z=((17,224,1449),(35,676,5649),(33,612,5163))[state];h*=x;q*=y;r*=z
 K=9*h*h-4*h*q+37*q*q
 U=-8*h*h+19*h*q+10*q*q-4*h*r+62*q*r
 E=7*h*q
 return dict(h=h,q=q,r=r,K=K,U=U,E=E,bare_apex=K+U-E,tail1=K+2*U)

def emit(states,t):
 if not states or any(s not in (0,1,2) for s in states):raise ValueError('a nonempty word over0,1,2 is required')
 T=substitute(word(states),word(states,True));return dict(states=states,profile=formula(states),source_rows=T,marked_root=0,tail_rows=selected_tail(T,1,t),tail_H=formula(states)['K']+2**t*formula(states)['U'])

def checks():
 local=[]
 for name,T,want in [('A',A,(17,224,1449)),('B',B,(17,224,1449)),('A_twin',AT,(35,676,5649)),('B_twin',BT,(33,612,5163))]:
  got=moment3(T)
  if got!=want:raise RuntimeError((name,got,want))
  local.append(dict(name=name,rows=T,moments=got))
 # Small quotient control, independently counting the inherited formula's origin.
 p=profile(Q,1)
 if (p['H'],p['U'],p['E'])!=(149,430,14):raise RuntimeError(p)
 out=dict(status='PASS',scope='Four local literal moment profiles and the inherited quotient initial profile; general Q formula is inherited, not recounted at large order',local_atoms=local,quotient=p,one_slot_states=[formula([s]) for s in range(3)])
 Path(__file__).with_name('fuelled_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--states');p.add_argument('--tail',type=int,default=1);a=p.parse_args()
 if a.states is None:checks()
 else:print(json.dumps(emit([int(x) for x in a.states.split(',')],a.tail)))
