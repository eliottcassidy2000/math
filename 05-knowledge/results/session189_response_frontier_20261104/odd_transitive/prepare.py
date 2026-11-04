#!/usr/bin/env python3
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def tt(n):return [sum(1<<j for j in range(i+1,n))for i in range(n)]
def reverse(r):return [((1<<len(r))-1)^s^(1<<i)for i,s in enumerate(r)]
def graph(m,external,into):
 out=tt(m)+[s<<m for s in external]
 for j in range(len(external)):
  for i in range(m):
   if into>>j&1:out[m+j]|=1<<i
   else:out[i]|=1<<(m+j)
 return out

def main():
 cases=[];graphs=[]
 def add(m,external,into,label):
  source=graph(m,external,into);target=list(source);mask=(1<<m)-1
  for i in range(m):target[i]^=mask^(1<<i)
  ids=[]
  for kind,rows in [('source',source),('local_converse',target)]:
   ids.append(len(graphs));graphs.append({'id':len(graphs),'case':len(cases),'kind':kind,'rows':rows})
  cases.append({'id':len(cases),'m':m,'exterior_rows':external,'into_mask':into,'label':label,'graph_ids':ids})
  return cases[-1]
 for m in range(3,8):
  for edge in range(2):
   external=[2,0]if edge else[0,1]
   for into in range(4):add(m,external,into,f'all_rooted_two_exterior_{edge}_{into}')
 for m in [5,7]:
  case=add(m,tt(6),56,'strong_C3_TT_m_TT3_TT3')
  case['witness_groups']=[list(range(m,m+3)),list(range(m+3,m+6))]
 # A fully covered root XOR in a fixed strong ambient, with TT2 witness children.
 m=5;core_n=9;source=tt(m)+[0]*6
 X=[5,6];Y=[7,8];outer=[9,10]
 for i in range(m):
  for x in X:source[i]|=1<<x
  for y in Y:source[y]|=1<<i
 source[5]|=1<<6;source[7]|=1<<8
 for x in X:
  for y in Y:source[x]|=1<<y
 for i in range(core_n):source[i]|=1<<9;source[10]|=1<<i
 source[9]|=1<<10
 exterior=[row>>m for row in source[m:]];into=sum(1<<j for j,row in enumerate(source[m:])if row&1)
 case=add(m,exterior,into,'fully_covered_root_in_strong_fixed_ambient')
 if graphs[case['graph_ids'][0]]['rows']!=source:raise RuntimeError('ambient source mismatch')
 case['witness_groups']=[X,Y];case['ambient_names']=outer
 # Every witness case also retains quotient-only converse on the same names.
 for case in cases:
  if 'witness_groups' not in case:continue
  m=case['m'];source=graphs[case['graph_ids'][0]]['rows'];target=list(source)
  groups=[list(range(m))]+case['witness_groups']
  for a in range(3):
   for b in range(a+1,3):
    for u in groups[a]:
     for v in groups[b]:target[u]^=1<<v;target[v]^=1<<u
  idx=len(graphs);graphs.append({'id':idx,'case':case['id'],'kind':'quotient_only','rows':target});case['graph_ids'].append(idx)

 (HERE/'controls.json').write_text(json.dumps({'cases':cases,'graphs':graphs},indent=2)+'\n')
 with (HERE/'graphs.txt').open('w')as f:
  f.write(str(len(graphs))+'\n')
  for z in graphs:f.write(str(len(z['rows']))+' '+' '.join(map(str,z['rows']))+'\n')
 print(json.dumps({'cases':len(cases),'graphs':len(graphs)}))
if __name__=='__main__':main()
