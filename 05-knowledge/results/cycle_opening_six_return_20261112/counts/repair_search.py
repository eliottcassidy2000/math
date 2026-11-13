import sys,itertools,time,json
from pathlib import Path
sys.path.insert(0,'05-knowledge/results/cycle_opening_six_return_20261112/counts')
from inflation_probe import literal,inflate,run_profile,evalrun
P=(78,12,24,16,35,79,30)
def birth(rows,a):
 n=len(rows);return tuple(r|((1<<n)if i==a else 0)for i,r in enumerate(rows))+(((1<<n)-1)^(1<<a),)
X=birth(birth(P,4),7);edges=list(itertools.combinations(range(len(X)),2));found=[];checked=0;t=time.time()
for cost in range(3):
 for changes in itertools.combinations(edges,cost):
  q=list(X)
  for u,v in changes:q[u]^=1<<v;q[v]^=1<<u
  h,s,f=literal(q);checked+=1
  for a in range(len(q)):
   if a==1 or s[a]!=1 or h+f[a]-s[a]!=293:continue
   H,out=run_profile(q,1,[a]);S,F=out[a]
   coeff=[(H[i]if i<len(H)else 0)+F[i]-S[i]for i in range(len(F))]
   good=coeff[:2]==[293,135] and not any(coeff[2:])and S[0]==1 and not any(S[1:])
   found.append({'cost':cost,'changes':changes,'entry':a,'parent':q,'H_runs':coeff,'S_runs':S,'uniform':good})
 print('cost',cost,'checked',checked,'found',len(found),'seconds',round(time.time()-t,2),flush=True)
out={'checked':checked,'target':'135*2^k+23','hits':found}
Path('05-knowledge/results/cycle_opening_six_return_20261112/counts/repair_search.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(found))
