"""Exact one-transitive-block endpoint run compiler, exploratory search."""
from pathlib import Path
import json
from hashlib import sha256
from math import comb
HERE=Path(__file__).resolve().parent

def add(a,b):
 for i,v in enumerate(b): a[i]+=v

def run_profile(rows,e,roots):
 fixed=[v for v in range(len(rows)) if v!=e];m=len(fixed);lim=m+2;full=(1<<m)-1
 def paths(start=None):
  dp=[[[0]*lim for _ in range(m+1)] for _ in range(full+1)]
  if start is None:dp[0][m][1]=1
  for i,v in enumerate(fixed):
   if start is None or v==start: dp[1<<i][i][0]=1
  for mask in range(full+1):
   for j in range(lim):
    for last in range(m+1):
     val=dp[mask][last][j]
     if not val:continue
     old=e if last==m else fixed[last]
     if last!=m and j+1<lim and rows[old]>>e&1:dp[mask][m][j+1]+=val
     for b,w in enumerate(fixed):
      if not(mask>>b&1) and rows[old]>>w&1:dp[mask|1<<b][b][j]+=val
  return dp
 dp=paths();p=[]
 for mask in range(full+1):
  a=[0]*lim
  for last in range(m+1):add(a,dp[mask][last])
  p.append(a)
 p[0][0]+=1
 out={}
 for a in roots:
  ia=fixed.index(a);F=[0]*(2*lim-1)
  for mask in range(full+1):
   if not mask>>ia&1:continue
   for j,x in enumerate(dp[mask][ia]):
    if x:
     for k,y in enumerate(p[full^mask]):
      F[j+k]+=x*y
  ds=paths(a);S=[0]*(2*lim-1)
  for last in range(m+1):
   for j,x in enumerate(ds[full][last]):S[j]+=x
  out[a]=(S[1:],F[1:])
 return p[full][1:],out

def evalrun(c,k):
 # coefficient j multiplies onto maps from k labelled TT vertices to j ordered runs.
 return sum(v*sum((-1)**(j-i)*__import__('math').comb(j,i)*i**k for i in range(j+1)) for j,v in enumerate(c,1))

def literal(rows):
 n=len(rows);f=(1<<n)-1;dp=[[0]*n for _ in range(f+1)];ds=[[0]*n for _ in range(f+1)];h=[0]*(f+1);h[0]=1
 for mask in range(1,f+1):
  for a in range(n):
   if mask>>a&1:
    q=mask^(1<<a);dp[mask][a]=1 if not q else sum(dp[q][b] for b in range(n) if rows[b]>>a&1);ds[mask][a]=1 if not q else sum(ds[q][b] for b in range(n) if rows[a]>>b&1)
  h[mask]=sum(dp[mask])
 F=[sum(dp[mask][a]*h[f^mask] for mask in range(1,f+1)) for a in range(n)]
 return h[f],ds[f],F

def inflate(rows,e,k):
 owner=[i for i in range(len(rows)) for _ in range(k if i==e else 1)]
 rr=tuple(sum((int(x<y) if owner[x]==owner[y] else((rows[owner[x]]>>owner[y])&1))<<y for y in range(len(owner)) if x!=y)for x in range(len(owner)))
 return rr,{i:owner.index(i)for i in range(len(rows))}

def main():
 src=HERE.parent/'returns/port_collisions.json'
 data=json.loads(src.read_text())
 hits=[];controls=0;tested=0
 for row in data['hits']:
  rows=row['parent'];a,b=row['ports']
  for e in range(len(rows)):
   if e in (a,b):continue
   H,out=run_profile(rows,e,(a,b));Sa,Fa=out[a];Sb,Fb=out[b]
   if tested<20:
    for k in(1,2,3):
     rr,idx=inflate(rows,e,k);h,s,f=literal(rr)
     if (evalrun(H,k),evalrun(Sa,k),evalrun(Fa,k),evalrun(Sb,k),evalrun(Fb,k))!=(h,s[idx[a]],f[idx[a]],s[idx[b]],f[idx[b]]):raise RuntimeError((rows,e,k,'literal'))
     controls+=1
   tested+=1
   D=[x-y-z+w for x,y,z,w in zip(Fa,Sa,Fb,Sb)]
   if not any(D):
    hits.append(dict(parent=rows,parent_H=row['parent_H'],ports=[a,b],inflate=e,strong_parent=row['strong_parent'],H_runs=H,Sa=Sa,Sb=Sb,Fa=Fa,Fb=Fb,samples=[dict(k=k,H=evalrun(H,k)+evalrun(Fa,k)-evalrun(Sa,k),defects=[evalrun(Sa,k),evalrun(Sb,k)])for k in(1,2,3,4,8)]))
 report=dict(status='FINITE-EXACT finite catalogue; all-module identities follow from the run-cover bijection',source_bank_sha256=sha256(src.read_bytes()).hexdigest(),universe='Every unequal-start equal-entry-score pair in the complete ordinary parent bank through order 8; each non-port vertex chosen as the one substituted nonempty module',tested_charts=tested,literal_controls=controls,neutral_infinite_charts=len(hits),strong_neutral_charts=sum(x['strong_parent'] for x in hits),hits=hits)
 (HERE/'inflation_probe.json').write_text(json.dumps(report,indent=2)+'\n')
 def modes(c):
  d={}
  for j,v in enumerate(c,1):
   for i in range(1,j+1):d[i]=d.get(i,0)+v*(-1)**(j-i)*comb(j,i)
  return {i:v for i,v in d.items()if v}
 ready=[]
 for q in hits:
  if not q['strong_parent']:continue
  if not any(c[0]in(1,2)and not any(c[1:])for c in(q['Sa'],q['Sb'])):continue
  cc=[(q['H_runs'][i]if i<len(q['H_runs'])else 0)+a-b for i,(a,b)in enumerate(zip(q['Fa'],q['Sa']))]
  ready.append({**{key:q[key]for key in('parent','ports','inflate')},'H':modes(cc),'defects':[modes(q['Sa']),modes(q['Sb'])]})
 (HERE/'ready_families.json').write_text(json.dumps(ready,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if k!='hits'}))

if __name__=='__main__':main()
