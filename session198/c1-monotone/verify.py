"""Reproduce exact finite coverage using only independently replayable sources."""
from pathlib import Path
import hashlib,json
from core import require
from source_oracle import source_bank,build_source
import bilinear,c3
from monotone import inventory,readout,ONE
P=Path(__file__).resolve().parent
GAPS=[7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821,1067,1339,1861,2647]

def replay(recipe):
 D=E=U=d=1
 for kind,L in recipe:
  if kind=='O':E=0;U*=2;continue
  trail=[];k=L
  while k>1:trail.append(k);k=(k+1)//2
  g=Q=B=1
  for k in reversed(trail):
   if k&1:ell=(k+1)//2;g+=ell-ell%2;Q+=B
   B*=2
  D,E,U,d=g*D-(L&1)*E+L*U,D-(E if L==1 else 0),Q*D+B*U,D
 return D,E,U,d

def hc(T):
 n=len(T);N=1<<n;A=[[0]*n for _ in range(N)];C=[[0]*n for _ in range(N)]
 for i in range(n):A[1<<i][i]=1
 C[1][0]=1
 for S in range(1,N):
  for u in range(n):
   if not S>>u&1:continue
   outs=T[u]&~S
   while outs:
    bit=outs&-outs;outs-=bit;v=bit.bit_length()-1
    A[S|bit][v]+=A[S][u];C[S|bit][v]+=C[S][u]
 return sum(A[-1]),sum(C[-1][u]for u in range(n)if T[u]&1)

def main():
 # Hostile/boundary probes do not rely on positive-range saturation.
 for N in range(1,30):
  support,seen,count=inventory(N,True)
  reference=source_bank(N)
  expected={s for s in reference if readout(s)<=N}
  require(seen==expected,'small target bank boundary')
  require(set(support)=={readout(s)for s in expected if s!=ONE},'small target image boundary')
 support,seen,count=inventory(5001,True)
 reference=source_bank(5001);expected={s for s in reference if readout(s)<=5001}
 require(seen==expected,'complete source set comparison at5001')
 small_count=count;old_count=len(reference)
 del support,seen,reference,expected
 N=50001;support,_,count=inventory(N)
 missing=[n for n in range(3,N+1,2)if n not in support]
 require(missing==GAPS,'exact C1 missing list')
 digest=hashlib.sha256()
 for H,r in sorted(support.items()):
  s=replay(r['recipe'])
  require(s==r['source'] and readout(s)==H and H%2==1,'independent positive source replay')
  require(3**len(r['recipe']) <= (2**len(r['recipe']))*H,'short physical recipe growth')
  require(all(L is None or L<=s[2] for _,L in r['recipe']),'bounded binary gadget labels')
  digest.update((json.dumps([H,s,r['recipe']],separators=(',',':'))+'\n').encode())
 repairs=json.loads((P/'REPAIR_WITNESSES.json').read_text());positive_repairs=[]
 for text,r in sorted(repairs.items(),key=lambda q:int(q[0])):
  H=int(text);record=r['record'];C=r['cycles'];v=replay(record['left_recipe']);w=replay(record['right_recipe'])
  require(v==tuple(record['left']) and w==tuple(record['right']),'repair source replay')
  module=bilinear if C==2 else c3
  require(module.value(record['kind'],v,w)==H,'repair scalar')
  T,root=module.realize(record);require(hc(T)==(H,C),'literal repair H and C')
  positive_repairs.append(dict(H=H,C=C,order=len(T),rows=T,root=root))
 require(set(map(int,repairs))==set(missing)-{7,21},'all exceptional repairs supplied')
 s2,_,_=bilinear.support(67);s3,_,_=c3.support(67)
 require(all(n not in s2 and n not in s3 for n in (7,21)),'complete small low-cycle exclusions')
 require(67 not in s2 and 67 in s3,'C3-only repair control')
 result=dict(status='PASS',bound=N,C1_values=len(support),C1_missing=missing,source_states=count,comparison_bound=5001,comparison_monotone_source_states=small_count,comparison_U_bounded_source_states=old_count,independently_replayed_positive_recipes=len(support),short_recipe_growth_checks=len(support),positive_certificate_sha256=digest.hexdigest(),C2_repairs=sum(x['C']==2 for x in positive_repairs),C3_repairs=sum(x['C']==3 for x in positive_repairs),maximum_repair_graph_order=max(x['order']for x in positive_repairs),C123_missing=[7,21])
 encoded=json.dumps(result,indent=2,sort_keys=True)+'\n';(P/'VALIDATION.json').write_text(encoded)
 (P/'REPAIR_GRAPHS.json').write_text(json.dumps(positive_repairs,indent=2,sort_keys=True)+'\n')
 print(encoded,end='')
if __name__=='__main__':main()
