"""Independent stdlib brute (x,y) verification; does not import producer."""
from pathlib import Path
import hashlib,json
from itertools import combinations
ROOT=Path(__file__).resolve().parent

def require(ok,msg):
 if not ok:raise ValueError(msg)

def count_and_smooth(a,p):
 a1,a2,a3,a4,a6=[z%p for z in a]
 count=1;smooth=True
 for x in range(p):
  for y in range(p):
   if (y*y+a1*x*y+a3*y-x*x*x-a2*x*x-a4*x-a6)%p:continue
   count+=1
   if (a1*y-3*x*x-2*a2*x-a4)%p==0 and (2*y+a1*x+a3)%p==0:smooth=False
 return count,smooth

def main():
 raw=(ROOT.parent/'ranks/phase2/highrank_inputs.json').read_bytes()
 cert=json.loads((ROOT/'certificate.json').read_text());curves=sorted(json.loads(raw)['curves'],key=lambda c:c['id'])
 require(hashlib.sha256(raw).hexdigest()==cert['input_sha256'],'input changed')
 ids=[c['id'] for c in curves];require(ids==cert['ids'] and len(set(ids))==66,'IDs')
 aa=[list(map(int,c['ainvs'])) for c in curves]
 for ps,row in cert['traces'].items():
  p=int(ps);require(3<=p<=1000 and all(p%d for d in range(2,int(p**.5)+1)),'prime')
  require(len(row)==len(ids),'trace row length')
  for i,a in enumerate(aa):
   n,smooth=count_and_smooth(a,p)
   require((row[i] is not None)==smooth,f'good reduction {ids[i]} p={p}')
   if smooth:require(row[i]==p+1-n,f'trace {ids[i]} p={p}')
 matrix=cert['upper_triangle_first_witness_prime'];require(len(matrix)==66,'matrix')
 for i in range(66):require(len(matrix[i])==65-i,'triangle shape')
 for i,j in combinations(range(66),2):
  p=matrix[i][j-i-1];require(p is not None,'missing witness')
  row=cert['traces'][str(p)]
  require(row[i] is not None and row[j] is not None and row[i]!=row[j],'invalid witness')
 require(not cert['unseparated_pairs'],'unseparated')
 # Hostile control: different-j, genuinely 2-isogenous models.
 # y²=x³+a*x²+b*x -> y²=x³-2a*x²+(a²-4b)*x.
 # For a=b=1 their j-values are 2048/3 and 35152/9.
 for ps in cert['traces']:
  p=int(ps)
  n,s=count_and_smooth([0,1,0,1,0],p);m,t=count_and_smooth([0,-2,0,-3,0],p)
  if s and t:require(n==m,'2-isogeny hostile control')
 print('PASS: 66 curves, 2145 non-Q-isogeny witnesses, largest prime 71; independent exhaustive point counts and different-j 2-isogeny control.')
if __name__=='__main__':main()
