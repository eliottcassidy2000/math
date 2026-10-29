"""Same-agent alternate implementations; explicit guards survive python -O."""
from pathlib import Path
from itertools import permutations,product
from collections import defaultdict,Counter
from fractions import Fraction as Q
from math import gcd
import json
import ports,paley,crt
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(v,m):
 global CHECKS
 CHECKS+=1
 if not v:raise RuntimeError(m)
def half_paths(T,length):
 n=len(T);layer={1<<v:[int(i==v) for i in range(n)] for v in range(n)}
 for size in range(1,length):
  nxt={}
  for mask,row in layer.items():
   for last,count in enumerate(row):
    if not count:continue
    options=T[last]&~mask
    while options:
     bit=options&-options;options-=bit;j=bit.bit_length()-1;newmask=mask|bit
     if newmask not in nxt:nxt[newmask]=[0]*n
     nxt[newmask][j]+=count
  layer=nxt
 return layer
def meet_hamilton(T):
 n=len(T)
 if n<2:return 1
 k=n//2;forward=half_paths(T,k)
 backwards=tuple(sum(1<<i for i in range(n) if T[i]>>j&1) for j in range(n))
 reverse=half_paths(backwards,n-k);full=(1<<n)-1;answer=0
 for mask,row in forward.items():
  other=reverse.get(full^mask)
  if other is None:continue
  for a,c in enumerate(row):
   if c:answer+=c*sum(d for b,d in enumerate(other) if T[a]>>b&1)
 return answer
def convolution(a,b):
 out=defaultdict(Q)
 for x,p in a.items():
  for y,q in b.items():out[x+y]+=p*q
 return out
def main(write=True):
 global CHECKS
 CHECKS=0
 # Global converse preserves the C3 parent; actual complement need not do so.
 S4=(2,12,9,1);P4=(3,0,1,2);C4=ports.core.complement(S4,P4)
 need(ports.core.is_strong(S4) and all(S4[a]>>b&1 for a,b in zip(P4,P4[1:])),'strong four-vertex native control')
 child4=tuple(row&7 for row in C4[:3]);op=ports.core.converse(S4)
 need(not ports.core.is_strong(child4) and ports.core.is_strong(tuple(row&7 for row in op[:3])),'native complement differs from global converse')
 left=ports.frames.p.controls.build_hostile(False);right=ports.frames.p.controls.build_hostile(True)
 def cycles(T):
  return [(a,b,c) for a in range(12) for b in range(a+1,12) for c in range(b+1,12)
          if ((T[a]>>b)&1)==((T[b]>>c)&1)==((T[c]>>a)&1)]
 need(cycles(left)==cycles(right),'literal 81/27 labelled cyclic triples')
 for T in (ports.tournament(1),ports.endpoint_tournament(5),paley.atom(7)):
  brute=sum(all(T[a]>>b&1 for a,b in zip(P,P[1:])) for P in permutations(range(len(T))))
  need(meet_hamilton(T)==brute,'middle-split paths versus literal paths')
 data=json.loads((HERE/'ports.json').read_text());hs=[]
 for row in data['cuts']['smaller_context']:
  S=row['source'];P=row['path'];n=len(S);arcs=set(zip(P,P[1:]));C=[]
  for a in range(n):
   C.append(sum(1<<b for b in range(n) if a!=b and (((S[a]>>b)&1) if ((a,b) in arcs or (b,a) in arcs) else ((S[b]>>a)&1))))
  h=meet_hamilton(tuple(C));need(h==row['H'],'literal complement and independent middle split');hs.append(h)
 p=7;A=paley.affine(p);H=paley.squares(p);tt=paley.neg(p,2);t=paley.neg(p,1)
 chars=Counter();fixed=Counter();total=0
 for root_a in H:
  for zero_a in H:
   for gs in product(A,repeat=3):
    child=[None]*p;child[0]=tuple(zero_a*x%p for x in range(p))
    for i,g in zip(H,gs):child[i]=g;child[-i%p]=paley.compose(t,paley.compose(g,t))
    c=paley.lift(p,root_a,0,child)
    need(paley.compose(c,tt)==paley.compose(tt,c),'literal centralizer permutation')
    chars[tuple(paley.characters(c,p,2))]+=1;fixed[sum(c[i]==i for i in range(49))]+=1;total+=1
 need(total==83349 and len(chars)==9 and set(chars.values())=={9261},'exact multiplier fibres')
 need(Q(sum(x*c for x,c in fixed.items()),total)==5,'Paley7 depth-two first moment')
 need(Q(sum(x*x*c for x,c in fixed.items()),total)==49,'Paley7 depth-two second moment')
 for r in json.loads((HERE/'crt.json').read_text())['rows']:
  n=r['period'];orb=sum(crt.lucas(gcd(n,j)) for j in range(n))//n
  need(orb==r['chart_count'],'independent divisor Burnside')
 local=Counter()
 for x in product(range(8),repeat=2):
  y=crt.golden_mul(x,(0,1),8);period=1
  while y!=x:y=crt.golden_mul(y,(0,1),8);period+=1
  local[period]+=1
 need(local=={1:1,3:3,6:12,12:48},'ring8 phase strata')
 X={1:Q(1)};Y={1:Q(1)};spines=[]
 for k in range(4):
  biased={x:x*q for x,q in X.items() if x};mean=sum(x*q for x,q in Y.items())
  vy=sum((x-mean)**2*q for x,q in Y.items());vb=sum((x-mean)**2*q for x,q in biased.items())
  need(sum(biased.values())==1 and sum(x*q for x,q in biased.items())==mean==2*k+1,'same spine means')
  need((vy,vb)==(4*k*(k-1),2*k*(k-1)),'shared-copy variance differs')
  spines.append({'depth':k,'mean':str(mean),'mirrored_variance':str(vy),'ordinary_variance':str(vb)})
  Y=convolution(Y,{2*x:q for x,q in X.items()})
  q=convolution(convolution(X,X),X);X={x:p/3 for x,p in q.items()};X[0]=X.get(0,0)+Q(2,3)
 out={'status':'PASS; same-agent alternate implementation, not an independent researcher audit','checks':CHECKS,
      'native_vs_converse':{'source':S4,'path':P4,'native_partner':C4,'native_parent':child4,'labelled_cyclic_triples_81_27':len(cycles(left))},
      'sixteen_vertex_H':hs,'Paley7_depth2':{'order':total,'character_fibres':[[list(k),v] for k,v in sorted(chars.items())],
                                            'fixed_counts':sorted(fixed.items())},
      'ring8_periods':sorted(local.items()),'spines':spines}
 if write:(HERE/'audit.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'audit_checks':CHECKS,'H':hs,'Paley7_centralizer':total}))
 return out
if __name__=='__main__':main()
