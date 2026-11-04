#!/usr/bin/env python3
import itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent

def profile(r):
 n=len(r); q=[0]*n
 for p in itertools.permutations(range(n)):
  for cuts in range(1<<(n-1)):
   if all((cuts>>i&1) or (r[p[i]]>>p[i+1]&1) for i in range(n-1)):q[cuts.bit_count()]+=1
 return q

def rows(code,n=4):
 r=[0]*n;b=0
 for i in range(n):
  for j in range(i+1,n):
   if code>>b&1:r[i]|=1<<j
   else:r[j]|=1<<i
   b+=1
 return r

def modules(r):
 n=len(r);return [m for m in range(1,1<<n) if m.bit_count()>1 and all((r[v]&m) in (0,m) for v in range(n) if not(m>>v&1))]

def xor(r,fam):
 r=r[:];n=len(r)
 for i in range(n):
  for j in range(i+1,n):
   if sum((m>>i&1) and (m>>j&1) for m in fam)%2:r[i]^=1<<j;r[j]^=1<<i
 return r

def main():
 records=[]
 for code in range(64):
  r=rows(code);q=profile(r);ms=modules(r);targets=[]
  for f in range(1<<len(ms)):
   fam=[m for j,m in enumerate(ms) if f>>j&1];t=xor(r,fam);targets.append(dict(family=fam,rows=t,profile=profile(t)))
  records.append(dict(code=code,rows=r,degrees=sorted(x.bit_count() for x in r),profile=q,nontransitive=q[0]>1,modules=ms,targets=targets))
 (P/'four_child_records.json').write_text(json.dumps(records,indent=2)+'\n')
 print('nontransitive',sum(r['nontransitive'] for r in records),'targets',sum(len(r['targets']) for r in records if r['nontransitive']))
 if not all(t['profile']==r['profile'] for r in records if r['nontransitive'] for t in r['targets']):raise RuntimeError('profile changed')
if __name__=='__main__':main()
