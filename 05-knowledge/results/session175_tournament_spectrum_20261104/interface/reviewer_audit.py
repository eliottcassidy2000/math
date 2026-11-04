#!/usr/bin/env python3
"""Independent completeness/maximality/strong/witness certificate audit."""
import json,itertools,hashlib
from pathlib import Path
P=Path(__file__).resolve().parent

def need(b,m):
 if not b:raise ValueError(m)
def matrix(f,k):
 a=[[False]*7 for _ in range(7)]
 for off,typ in [(0,f[0]),(3,f[1])]:
  for i in range(3):
   for j in range(i+1,3):
    forward=typ=='T' or (i,j)!=(0,2)
    a[off+i][off+j]=forward;a[off+j][off+i]=not forward
 for i in range(3):
  a[6][i]=True;a[i+3][6]=True
  for j in range(3):a[i][j+3]=bool(k>>(3*i+j)&1);a[j+3][i]=not a[i][j+3]
 return a
def isstrong(a):
 r=[row[:] for row in a]
 for k in range(7):
  for i in range(7):
   for j in range(7):r[i][j]|=r[i][k] and r[k][j]
 return all(all(row) for row in r)
def literal(a):return sum(all(a[p[i]][p[i+1]] for i in range(6)) for p in itertools.permutations(range(7)))
def components(h,allowed):
 unseen=set(allowed);out=[]
 while unseen:
  x=unseen.pop();stack=[x];c={x}
  while stack:
   x=stack.pop()
   for b in range(9):
    y=x^(1<<b)
    if y in unseen and abs(h[x]-h[y])<=2:unseen.remove(y);c.add(y);stack.append(y)
  out.append(frozenset(c))
 return set(out)
def prime(n):return n>=2 and all(n%d for d in range(2,n))
def main():
 j=json.loads((P/'certificate.json').read_text());need(j==json.loads((P/'certificate.optimized.json').read_text()),'optimized equality')
 details=[];matrices={};lane_values=set()
 for f in j['families']:
  h=[];strong=[]
  for k in range(512):
   a=matrix(f['name'],k);matrices[f['name'],k]=a
   need(all(not a[i][i] for i in range(7)) and all(a[i][l]!=a[l][i] for i in range(7) for l in range(i)),'tournament')
   h.append(literal(a));strong.append(isstrong(a))
   if strong[-1] and (prime(h[-1]) or (h[-1]%7==0 and prime(h[-1]//7) and h[-1]//7!=3)):lane_values.add(h[-1])
  need(h==f['H'],'all literal values');need(sum(strong)==f['strong_states'],'strong census')
  need(sorted(set(h))==f['distinct_H'],'distinct H')
  for key,allowed in [('components',range(512)),('strong_components',[k for k in range(512) if strong[k]])]:
   got=components(h,allowed);listed=[frozenset(c['states']) for c in f[key]]
   need(len(listed)==len(set(listed)) and set(listed)==got,'maximal exact component partition')
   for c in f[key]:
    vals=sorted({h[k] for k in c['states']})
    need(c['size']==len(c['states'])==len(set(c['states'])),'size')
    need(vals==c['values']==list(range(c['min'],c['max']+1,2)),'exact interval metadata')
  details.append({'family':f['name'],'literal_states':512,'strong':sum(strong),'full_components':len(f['components']),'strong_components':len(f['strong_components'])})
 need({int(x) for x in j['witnesses']}==lane_values,'all lane values witnessed')
 for key,w in j['witnesses'].items():
  a=matrices[w['family'],w['interface']];rows=[sum(int(v)<<i for i,v in enumerate(row)) for row in a]
  need(rows==w['rows'],'witness rows');need(isstrong(a) and w['strong'] is True,'witness strong')
  h=literal(a);need(h==int(key),'witness literal')
  need(prime(h) if w['lane']=='prime' else h%7==0 and prime(h//7) and h//7!=3,'witness lane')
 print(json.dumps({'status':'PASS','families':details,'witnesses':len(lane_values),'certificate_sha256':hashlib.sha256((P/'certificate.json').read_bytes()).hexdigest()},sort_keys=True,indent=2))
if __name__=='__main__':main()
