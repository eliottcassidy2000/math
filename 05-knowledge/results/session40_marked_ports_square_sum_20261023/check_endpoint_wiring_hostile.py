import itertools,json
checks=0

def data(n,T):
 es=list(itertools.combinations(range(n),2)); ix={e:i for i,e in enumerate(es)}
 def arc(a,b): return bool((T>>ix[min(a,b),max(a,b)])&1)==(a<b)
 ps=[p for p in itertools.permutations(range(n)) if all(arc(a,b) for a,b in zip(p,p[1:]))]
 return es,ix,arc,ps

for n in range(1,5):
 for T in range(1<<(n*(n-1)//2)):
  es,ix,arc,ps=data(n,T)
  ends=[(p[0],p[-1]) for p in ps]
  if len(ends)!=len(set(ends)):raise ValueError('small endpoint collision')
  checks+=1
es,ix,arc,ps=data(5,4)
if len(ps)!=5:raise ValueError('source count')
rows=[]
for p,expected,H in [((4,0,3,2,1),871,15),((4,2,0,3,1),733,11)]:
 if p not in ps:raise ValueError('invalid source path')
 pe={frozenset(e) for e in zip(p,p[1:])}
 U=sum(1<<i for i,(a,b) in enumerate(es) if (arc(a,b) if frozenset((a,b)) in pe else not arc(a,b)))
 if U!=expected:raise ValueError('complement target')
 _,_,ua,ups=data(5,U)
 if len(ups)!=H:raise ValueError('target count')
 V=((1<<10)-1)^U
 if len(data(5,V)[3])!=H:raise ValueError('converse count')
 scores=sorted(sum(ua(a,b) for b in range(5) if b!=a) for a in range(5))
 if scores!=[1,2,2,2,3]:raise ValueError('target scores')
 rows.append(dict(path=p,target=U,H=H,scores=scores))
 checks+=5
print(json.dumps(dict(status='PASS',checks=checks,source=4,source_H=5,targets=rows,small_order_universe='all labelled tournaments n1..4'),indent=2))
