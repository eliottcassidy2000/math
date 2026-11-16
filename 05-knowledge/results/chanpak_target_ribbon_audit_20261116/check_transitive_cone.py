"""Exact BFS and conic separator; no LP solver is required."""
from collections import deque
import json
from pathlib import Path
from cap_check import transfer_system,need

def main():
 states,_,raw,ff=transfer_system(3)
 A=[raw[frozenset(i+1 for i in range(3) if m>>i&1)] for m in range(8)]
 d=len(states);initial=(1,)+(0,)*(d-1)
 # The historical separator uses the same reference state ordering.
 def step(v,m):
  z=[0]*d
  for (i,j),c in A[m].items(): z[i]+=c*v[j]
  return tuple(z)
 def H(v):return sum(v[i] for i in ff)
 stack=deque([initial]);seen={initial};tested=0;best=0
 while stack:
  y=stack.popleft()
  for m in range(8):
   v=step(y,m);h=H(v)
   if y!=initial:best=max(best,h);tested+=1
   if h==1 and v not in seen:seen.add(v);stack.append(v)
 seen.remove(initial);need(len(seen)==99 and tested==792 and best==39)
 C=json.loads(Path(__file__).with_name("TRANSITIVE_CONE_SEPARATOR.json").read_text())
 def dot(z):return sum(a*b for a,b in zip(C["separator"],z))
 need(all(dot(v)>=0 for v in seen))
 y=initial
 for m in C["word"]:y=step(y,m)
 need(H(y)==3 and dot(y)==C["value"]==-7)
 print("99 transitive rays, 792 exact extensions, maximum 39; reachable H=3 ray has separator −7.")
if __name__=='__main__':main()
