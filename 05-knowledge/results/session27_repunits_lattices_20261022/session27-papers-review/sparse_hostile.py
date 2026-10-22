from pathlib import Path
import itertools,json
G=[[int(x) for x in row.split(',')] for row in open(str(Path(__file__).with_name('downloaded_G96_reduced.txt'))).read().strip().strip('[]').split(';')];n=len(G);best=(min(G[i][i] for i in range(n)),None)
for r in [2,3,4]:
 for I in itertools.combinations(range(n),r):
  for ss in itertools.product([-1,1],repeat=r-1):
   S=(1,)+ss;q=sum(G[i][i] for i in I)+2*sum(S[a]*S[b]*G[I[a]][I[b]] for a in range(r) for b in range(a))
   if q<best[0]:best=(q,(I,S));print('best',best,flush=True)
 print('finished support',r,'best',best,flush=True)
