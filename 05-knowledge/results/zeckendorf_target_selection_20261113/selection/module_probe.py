from response_collision import profiles
from itertools import permutations
import json
from pathlib import Path
S=(46,12,56,48,35,2);mapping=(6,0,5,4,3,2,1);oldT=(16,1,3,71,14,31,55)
T=tuple(sum(((oldT[mapping[i]]>>mapping[j])&1)<<j for j in range(7))for i in range(7));ROOT=4

def tt(n):return tuple(((1<<n)-1)^((1<<(i+1))-1)for i in range(n))
def inflate(rows,v,module):
 blocks=[module if i==v else (0,) for i in range(len(rows))];starts=[];z=0
 for B in blocks:starts.append(z);z+=len(B)
 out=[]
 for i,B in enumerate(blocks):
  ext=sum((((1<<len(blocks[j]))-1)<<starts[j])for j in range(len(rows))if rows[i]>>j&1)
  out.extend((r<<starts[i])|ext for r in B)
 return tuple(out),starts[ROOT]
def main():
 print('normalizedtarget',T)
 records=[]
 for v in range(6):
  if v==ROOT:continue
  # At most min(in,out)+2 runs in an ordered2cover.
  bound=max(min(len(R)-1-R[v].bit_count(),R[v].bit_count())+2 for R in (S,T))
  tests=[]
  for m in range(1,bound+1):
   A,r=inflate(S,v,tt(m));B,q=inflate(T,v,tt(m));x=profiles(A)[r];y=profiles(B)[q]
   tests.append(dict(m=m,left=x,right=y,equal=x==y))
  print(v,bound,[(z['m'],z['equal'])for z in tests])
  records.append(dict(vertex=v,run_bound=bound,tests=tests,universal=all(z['equal']for z in tests)))
 Path(__file__).with_name('module_probe.json').write_text(json.dumps(dict(source=S,target=T,root=ROOT,records=records),indent=2)+'\n')
if __name__=='__main__':main()
