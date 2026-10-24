#!/usr/bin/env python3
"""Independent exact labeled tournament degree-sequence dynamic programming."""
import functools,itertools,json,math,argparse
@functools.lru_cache(None)
def count(degrees):
 n=len(degrees)
 if not n:return 1
 if min(degrees)<0 or max(degrees)>=n or sum(degrees)!=n*(n-1)//2:return 0
 # Landau inequalities are only a safe rejection filter, not a count oracle.
 if any(sum(degrees[:k])<k*(k-1)//2 for k in range(1,n+1)):return 0
 d=degrees[0];groups=[(x,len(tuple(g))) for x,g in itertools.groupby(degrees[1:])]
 total=0
 # Choose root outneighbors in each equal-degree group. A vertex beating the
 # deleted root loses one unit of required degree; root outneighbors do not.
 for wins in itertools.product(*(range(c+1) for x,c in groups)):
  if sum(wins)!=d:continue
  rem=[];mult=1
  for (x,c),w in zip(groups,wins):
   rem.extend([x]*w);rem.extend([x-1]*(c-w));mult*=math.comb(c,w)
  total+=mult*count(tuple(sorted(rem)))
 return total

def run():
 out=[]
 for n in (3,5,7,9,11):
  k=(n-1)//2
  regular=count((k,)*n)
  particular=count((k-1,)*k+(k,)*k)
  almost=math.comb(2*k,k)*particular
  if regular!=almost:raise RuntimeError('regular/almost-regular labeled bijection')
  out.append({'regular_order':n,'labelled_regular_count':regular,'almost_regular_order':2*k,'fixed_assignment_almost_regular_count':particular,'labelled_almost_regular_count':almost})
 parser=argparse.ArgumentParser();parser.add_argument('--output',default='degree-mass.json');args=parser.parse_args()
 result={'status':'PASS','method':'independent residual-degree dynamic programming; equal-degree subset multiplicities','results':out,'cache':str(count.cache_info())}
 with open(args.output,'w') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
 print(json.dumps(result,sort_keys=True))
if __name__=='__main__':run()
