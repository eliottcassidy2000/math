#!/usr/bin/env python3
"""Independent first-vertex subset recursion for final literal padded bank."""
import json
from functools import lru_cache
from pathlib import Path
from homometry_independent import must

def subset_counter(T):
 @lru_cache(None)
 def start(S,v):
  rest=S^(1<<v)
  if not rest:return 1
  ok=rest&T[v];result=0
  while ok:
   bit=ok&-ok;ok-=bit;result+=start(rest,bit.bit_length()-1)
  return result
 H=[1]
 for S in range(1,1<<len(T)):
  H.append(sum(start(S,v)for v in range(len(T))if S>>v&1))
 return H,start

def packet(T):
 n=len(T);full=(1<<n)-1;deletion=full^1
 reverse=[sum(1<<u for u in range(n)if T[u]>>v&1)for v in range(n)]
 H,start=subset_counter(T);_,end=subset_counter(reverse)
 U=V=d=e=Q=0
 for S in range(full+1):
  if S&1:
   sign=1 if S.bit_count()%2 else -1
   u=start(S,0)*H[full^S];v=end(S,0)*H[full^S]
   U+=u;V+=v;d+=sign*u;e+=sign*v
  else:Q+=H[S]*H[deletion^S]
 return [H[full],H[deletion],Q,U,V,d,e]

def main():
 p=Path(__file__).resolve().parent.parent
 j=json.loads((p/'additive/PADDED_BANK.json').read_text());out=[]
 for w in j['witnesses']:
  got=packet(w['rooted_rows']);must(got==w['packet'],'literal full padded packet');out.append(got)
 print(json.dumps(dict(status='PASS',literal_packets=out,method='Independent first-vertex subset recursion on supplied literal graph rows; no lane producer imports'),indent=2))
if __name__=='__main__':main()
