from pathlib import Path
import sys,json,itertools
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'session275_exact_selection_20261111'/'switches'))
from carrier import profile,fuse
from escape_check import context,T41,T43,T45
from check import H

def common(a,b):
 n=len(a);inter=[x&y for x,y in zip(a,b)];dp=[0]*(1<<n);prev={}
 for i in range(n):dp[1<<i]|=1<<i
 for m in range(1,1<<n):
  v=dp[m]
  while v:
   z=v&-v;v^=z;i=z.bit_length()-1;rr=inter[i]&~m
   while rr:
    bit=rr&-rr;rr^=bit;j=bit.bit_length()-1
    if not dp[m|bit]&bit:prev[(m|bit,j)]=(m,i)
    dp[m|bit]|=bit
 if not dp[-1]:return None
 m=(1<<n)-1;i=(dp[m]&-dp[m]).bit_length()-1;out=[]
 while True:
  out.append(i)
  if (m,i) not in prev:break
  m,i=prev[m,i]
 return out[::-1]
for n in range(2,5):
 A=fuse(tuple(context(n,1)),tuple(T45));p=profile(A);target=57*(1<<n)-15;ed=list(itertools.combinations(range(len(A)),2));found=[]
 for es in itertools.chain(((e,) for e in ed),itertools.combinations(ed,2)):
  B=list(A)
  for i,j in es:B[i]^=1<<j;B[j]^=1<<i
  B=tuple(B)
  if H(B)!=p['H']:continue
  q=profile(B)
  if q['V']==target:
   frame=common(A,B)
   if frame:found.append(dict(edges=es,rows=B,profile=q,frame=frame))
 print(json.dumps(dict(n=n,H=p['H'],current_V=p['V'],target_V=target,repairs=found)),flush=True)
