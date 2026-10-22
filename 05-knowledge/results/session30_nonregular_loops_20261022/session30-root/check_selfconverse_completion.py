"""Independent converse testing on supplied complete census, no producer imports."""
import json,itertools,sys
from pathlib import Path

def matrix(n,bits):
 s=[[0]*n for _ in range(n)]
 for b,(i,j) in zip(bits,itertools.combinations(range(n),2)):
  s[i][j]=1 if b=='1' else -1;s[j][i]=-s[i][j]
 return s

def anti_exists(s,root=None):
 n=len(s);scores=list(map(sum,s));mapping=[]
 def search(used):
  i=len(mapping)
  if i==n:return True
  for j in range(n):
   if used>>j&1 or scores[j]!=-scores[i]:continue
   if root is not None and ((i==root)!=(j==root)):continue
   if any(s[j][mapping[k]]!=-s[i][k] for k in range(i)):continue
   mapping.append(j)
   if search(used|1<<j):return True
   mapping.pop()
  return False
 return search(0)

base=Path(__file__).resolve().parent.parent
bank=json.loads((base/'session30-almost-regular/census8_bank.json').read_text())
result={}
for key,n in [('regular9',9),('almostregular8',8)]:
 rec=[]
 for row in bank[key]:
  s=matrix(n,row['bits']);sc=anti_exists(s)
  if n%2:
   fixed=[v for v in range(n) if anti_exists(s,v)]
   if bool(fixed)!=sc:raise RuntimeError('fixed root existence')
  else:
   d=list(map(sum,s));completion=[r+[-d[i]] for i,r in enumerate(s)]+[d+[0]]
   if anti_exists(completion,n)!=sc:raise RuntimeError('completion mismatch')
   fixed=None
  rec.append({'bits':row['bits'],'self_converse':sc,'admissible_fixed_vertices':fixed})
 result[key]={'count':len(rec),'self_converse':sum(r['self_converse'] for r in rec),'records':rec}
if result['regular9']['self_converse']!=result['almostregular8']['self_converse']:raise RuntimeError('class count mismatch')
print(json.dumps({'status':'PASS','scope':'independent converse tests on producer mass-certified complete census','results':result},indent=2))
