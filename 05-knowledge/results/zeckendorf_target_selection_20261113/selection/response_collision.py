from pathlib import Path
from math import gcd
from itertools import permutations
import subprocess,json
REF='ed79a99437c3bf19509b2f9b7b696ab9349032d3'
BANK='05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/extended_census.json'
def moments(T):
 n=len(T);full=(1<<n)-1;E=[[0]*n for _ in range(full+1)];H=[1]+[0]*full
 for v in range(n):E[1<<v][v]=1
 for m in range(1,full+1):
  for v,z in enumerate(E[m]):
   if not z:continue
   bits=T[v]&~m
   while bits:
    b=bits&-bits;bits-=b;E[m|b][b.bit_length()-1]+=z
  H[m]=sum(E[m])
 return H,E

def profiles(T):
 n=len(T);full=(1<<n)-1;H,E=moments(T);_,S=moments(tuple(full^(1<<i)^row for i,row in enumerate(T)))
 return [(H[-1],S[-1][r],H[full^(1<<r)],sum(E[m][r]*H[full^m]for m in range(1,full+1)))for r in range(n)]
def main():
 bank=json.loads(subprocess.check_output(['git','show',REF+':'+BANK]))['records'];seen={};hits=[];tested=0
 for rec in bank:
  if not rec['strong']:continue
  T=tuple(rec['adjacency_outmasks']);Ts=[T]
  if rec['ordinary_multiplicity']==2:
   full=(1<<len(T))-1;Ts.append(tuple(full^(1<<i)^row for i,row in enumerate(T)))
  for T in Ts:
   for r,z in enumerate(profiles(T)):
    tested+=1
    if gcd(gcd(z[0],z[1]),gcd(z[2],z[3]))!=1:continue
    for S,q in seen.get(z,[]):
     if len(S)!=len(T):hits.append(dict(state=z,left=[S,q],right=[T,r]))
    seen.setdefault(z,[]).append((T,r))
 for hit in hits:
  S,r=hit['left'];T,q=hit['right'];best=None
  for omit in range(len(T)):
   if omit==q:continue
   for mapping in permutations([v for v in range(len(T))if v!=omit]):
    if mapping[r]!=q:continue
    if any(((S[r]>>i)&1)!=((T[q]>>mapping[i])&1)for i in range(len(S))if i!=r):continue
    changes=[(i,j)for i in range(len(S))for j in range(i+1,len(S))if((S[i]>>j)&1)!=((T[mapping[i]]>>mapping[j])&1)]
    if best is None or len(changes)<best[0]:best=(len(changes),omit,mapping,changes)
  hit['star_preserving_minimum_embedding']=best
 print('tested',tested,'cross-order primitive hits',len(hits))
 print(hits[:3])
 Path(__file__).with_name('response_collision.json').write_text(json.dumps(dict(ref=REF,scope='Inherited ordinary strong classes through7, primitive rooted4-state collisions across orders',tested=tested,hits=hits),indent=2)+'\n')
if __name__=='__main__':main()
