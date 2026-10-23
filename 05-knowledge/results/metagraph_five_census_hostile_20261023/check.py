"""Independent Python finite audit of the explicit order-23 hostile pair."""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
from hashlib import sha256
import argparse,json
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(ok,why):
 global CHECKS
 CHECKS+=1
 if not ok:raise ValueError(why)

def tables(k,root=False):
 pairs=list(combinations(range(k),2));out=[]
 perms=[(0,)+p for p in permutations(range(1,k))] if root else list(permutations(range(k)))
 for code in range(1<<len(pairs)):
  a=[[0]*k for _ in range(k)]
  for bit,(i,j) in enumerate(pairs):a[i][j]=(code>>bit)&1;a[j][i]=1-a[i][j]
  out.append(min(sum(a[p[i]][p[j]]<<bit for bit,(i,j) in enumerate(pairs)) for p in perms))
 return out

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
 p=23;cases=[];alltab={k:tables(k) for k in range(1,6)};roottab={k:tables(k,True) for k in range(1,6)}
 need(len(set(alltab[5]))==12 and len(set(roottab[5]))==48,'five-type universes')
 for mask in (382046,431774):
  S=[i for i in range(1,p) if mask>>i&1]
  need(len(S)==11 and all(((mask>>i)&1)+((mask>>(p-i))&1)==1 for i in range(1,p)),'Cayley tournament')
  a=[[int(i!=j and ((mask>>((j-i)%p))&1)) for j in range(p)] for i in range(p)]
  need(all(a[i][j]==a[(i+1)%p][(j+1)%p] for i in range(p) for j in range(p)),'translation automorphism')
  need(all(a[i][(i+1)%p] for i in range(p)),'strong Hamiltonian cycle')
  need(all(sum(row)==11 for row in a),'regular degree11')
  glob={};root={}
  for k in range(1,6):
   hist=Counter();rhist=Counter()
   for subset in combinations(range(p),k):
    code=sum(a[i][j]<<bit for bit,(i,j) in enumerate(combinations(subset,2)))
    hist[alltab[k][code]]+=1
    if subset[0]==0:rhist[roottab[k][code]]+=1
   glob[k]=dict(sorted(hist.items()));root[k]=dict(sorted(rhist.items()))
  tt6=0;root_tt6=0
  for subset in combinations(range(p),6):
   degrees=sorted(sum(a[i][j] for j in subset) for i in subset)
   if degrees==list(range(6)):
    tt6+=1
    if subset[0]==0:root_tt6+=1
  # Prime-order circulant isomorphism classes have canonical multiplier representative.
  orbit=sorted({sum(1<<((u*s)%p) for s in S) for u in range(1,p)})
  need(mask==min(orbit),'canonical multiplier representative')
  cases.append({'mask':mask,'difference_set':S,'multiplier_orbit':orbit,
    'global_censuses':glob,'root0_censuses':root,'transitive_six_sets':tt6,'root0_in_transitive_six':root_tt6})
 need(cases[0]['global_censuses']==cases[1]['global_censuses'],'every induced census throughfive equal')
 need(cases[0]['root0_censuses']==cases[1]['root0_censuses'],'every rooted census throughfive equal')
 need(set(cases[0]['multiplier_orbit']).isdisjoint(cases[1]['multiplier_orbit']),'multiplier orbits disjoint')
 need([c['transitive_six_sets'] for c in cases]==[1173,1265],'six-motif separates')
 need([c['root0_in_transitive_six'] for c in cases]==[306,330],'rooted six incidence separates')
 # An independent source-extension double count: a k-set is in choose(23-k,5-k) five-sets.
 from math import comb
 for k in range(1,6):
  need(sum(cases[0]['global_censuses'][k].values())==comb(23,k),'global census mass')
  need(sum(cases[0]['root0_censuses'][k].values())==comb(22,k-1),'rooted census mass')
 result={'status':'FINITE-EXACT explicit strongly connected order23 pair; all-size padding theorem in README',
  'orientation':'i -> j iff (j-i) mod23 belongs to difference_set','cases':cases,'checks':CHECKS,
  'source_sha256_lf':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()}
 data=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode()
 if args.write:(HERE/'receipt.json').write_bytes(data)
 else:need((HERE/'receipt.json').read_bytes().replace(b'\r\n',b'\n')==data,'frozen receipt')
 print(json.dumps({'checks':CHECKS,'rooted5_nonzero_types':len(cases[0]['root0_censuses'][5]),'transitive6':[1173,1265],'receipt_sha256':sha256(data).hexdigest()}))
if __name__=='__main__':main()
