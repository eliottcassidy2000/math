#!/usr/bin/env python3
"""Audit the power sources using actual adjacency, without cover-weight code."""
import json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(p,msg):
 if not p:raise RuntimeError(msg)

def path_counts(rows):
 n=len(rows);size=1<<n;dp=[[0]*n for _ in range(size)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,size):
  for v in range(n):
   z=dp[mask][v]
   if not z:continue
   choices=rows[v] & (size-1-mask)
   while choices:
    bit=choices&-choices;choices-=bit;u=bit.bit_length()-1;dp[mask|bit][u]+=z
 return sum(dp[-1]),dp[-1]

def marked(parent,mask):
 rows=list(parent)+[mask];x=len(parent)
 for v in range(x):
  if not(mask>>v&1):rows[v]|=1<<x
 return rows,x

def append_zero(rows,x):
 rows=list(rows)+[1<<x];new=len(rows)-1
 for v in range(new):
  if v!=x:rows[v]|=1<<new
 return rows

def validate(rows):
 n=len(rows)
 need(all(0<=r<(1<<n) and not(r>>u&1) for u,r in enumerate(rows)),'row range/loop')
 need(all(((rows[u]>>v&1)+(rows[v]>>u&1))==1 for u in range(n) for v in range(u)),'complete orientation')

def main():
 input_file=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'power_sources.json'
 data=json.loads(input_file.read_text());checks=0;special=[];fvs=None
 for rec in data['hits']:
  parent,mask=rec['rows'],rec['mask'];validate(parent)
  h,ends=path_counts(parent);g=sum(ends[v] for v in range(len(parent)) if mask>>v&1)
  rows,x=marked(parent,mask);F,_=path_counts(rows)
  tail=append_zero(rows,x);H1,_=path_counts(tail);q=H1-F-g
  need([F,q,g,h]==rec['state'],'independent actual source state')
  need(q>0 and q&(q-1)==0 and q-F-g==rec['deficit'],'dyadic deficit')
  checks+=1
  if rec['deficit']==9:
   for r in range(1,9):
    if r>1:tail=append_zero(tail,x)
    actual,_=path_counts(tail);expected=q*(1<<r)-9
    need(actual==expected,'deficit-nine actual adjacency family')
    special.append({'r':r,'vertices':len(tail),'actual_H':actual,'exponent':q.bit_length()-1+r})
   # Source minimum feedback-vertex deletion is exhaustively checked.
   import itertools
   def acyclic(keep):
    return all(not(parent[a]>>b&1 and parent[b]>>c&1 and parent[c]>>a&1) and not(parent[b]>>a&1 and parent[a]>>c&1 and parent[c]>>b&1) for a,b,c in itertools.combinations(keep,3))
   fvs=next(k for k in range(len(parent)+1) if any(acyclic([v for v in range(len(parent)) if v not in D]) for D in itertools.combinations(range(len(parent)),k)))
   need(fvs==3,'deficit-nine parent is outside two-feedback-vertex carriers')
 need(len(special)==(8 if any(r['deficit']==9 for r in data['hits']) else 0),'deficit-nine audit scope')
 print(json.dumps({'status':'passed','input_file':input_file.name,'sources_checked':checks,'method':'Literal endpoint subset DP on the parent, marked extension and first zero append; no cover-weight or transfer implementation imported.','deficit9_actual_adjacency':special,'deficit9_parent_minimum_feedback_vertices':fvs,'scope':'Every accepted finite source yields its stated all-length pure-tail family by the proved recurrence; search completeness is reported only by the input census; sampled witnesses make no absence claim.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
