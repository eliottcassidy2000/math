import json,itertools
from pathlib import Path
D=json.load(open(Path(__file__).resolve().parent/'RESULT.json'))['dictionary'];out=[]
for h,(a,b,n,code) in D.items():
 edges=list(itertools.combinations(range(n),2)); E=set((i,j) if (code>>k)&1 else (j,i) for k,(i,j) in enumerate(edges))
 count=sum(all((p[i],p[i+1]) in E for i in range(n-1)) for p in itertools.permutations(range(n)))
 assert count==int(h);assert 11*n>=16*a+24*b;out.append([int(h),n,count])
assert 375**3*405**2==15**11
print(json.dumps({'all_permutation_counts':out,'integer_dual_checks':'pass','40_vertex_identity':'pass'}))
