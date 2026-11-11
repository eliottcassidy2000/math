#!/usr/bin/env python3
"""New two-step applications of the incoming sign-flip law; not its re-audit."""
import json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'constructor'))
from ordinal_interfaces import profile,selected_tail

def tt(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def step(B,O,a):
 if a not in (3,5):raise ValueError('a=3 or5 required')
 n=len(B);k=(a+1)//2;N=n+k
 out=list(B)+[x<<n for x in tt(k)]+[O]
 for i in range(n):
  out[i]|=((1<<k)-1)<<n
  if not O>>i&1:out[i]|=1<<N
 for i in range(n,N):out[i]|=1<<N
 return tuple(out),(1<<n)|(1<<(N-1))
def state(B,O):
 p=profile(B,O)
 return p['H']-p['W'],p['U'],p['H'],p['E']
def F(s,a):
 r,M,h,E=s
 return -r,a*(h+r)+(a*a+1)*M-E,r+(2*a-2)*M,a*M+r-h

def main():
 out=[]
 for a in (3,5):
  for b in (3,5):
   A,O=step((0,),1,a);B,P=step(A,O,b);ss=state(B,P);pred=F(F((1,1,1,1),a),b)
   if ss!=pred:raise RuntimeError('composite state mismatch')
   tailH=profile(selected_tail(B,P,1),0)['H']
   if tailH!=ss[0]+2*ss[1]:raise RuntimeError('tail mismatch')
   out.append({'word':[a,b],'rows':B,'marks':P,'state':ss,'order':len(B),'first_tail_H':tailH})
 receipt={'status':'PASS','input_status':'Derived from incoming33114ad gadget law; incoming independent audit pending','scope':'Four new composite initial cases and first positive tails; no inherited tower census','controls':out}
 HERE.joinpath('incoming_service_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
