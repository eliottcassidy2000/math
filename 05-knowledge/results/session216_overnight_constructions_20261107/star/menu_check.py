import json,itertools,pathlib
# Reuse only the independently implemented literal subset-DP functions, not a recurrence counter.
ns={};exec((pathlib.Path(__file__).with_name('check.py')).read_text().split('checks=0')[0],ns)
tables,whole=ns['tables'],ns['whole']
def req(b,m):
 if not b:raise RuntimeError(m)
def state(a,S):
 h,s,e=tables(a);N=(1<<len(a))-1;n=len(a)
 alpha=sum(h[N^q]*sum(s[q][v] for v in range(n) if S>>v&1) for q in range(1,N+1))
 tau=sum(sum(e[p][v] for v in range(n) if S>>v&1)*sum(s[N^p][v] for v in range(n) if S>>v&1) for p in range(1,N))
 en=sum(e[N][v] for v in range(n) if S>>v&1)
 return (h[N],alpha,tau,en)
seeds=[([0,5,9,19,39,15,63],105,150),([0,1,3,23,71,31,175,63],77,152),([0,5,17,7,43,15,63],101,154),([0,1,11,35,15,87,31],109,156)]
checks=0;seed_receipts=[]
for a,S,tau0 in seeds:
 req(state(a,S)==(9,245,tau0,9),'seed state')
 core=whole(a,[],S,0,1);n=len(core)
 cycle=next((list((0,)+p) for p in itertools.permutations(range(1,n)) if all(core[x]>>y&1 for x,y in zip((0,)+p,p+(0,)))),None)
 req(cycle is not None,'initial strong cycle')
 seed_receipts.append({'rows':a,'cut':S,'tau':tau0,'core_cycle':cycle})
 for m in range(3):
  for bits in itertools.product([0,1],repeat=m):
   b=a.copy();cut=S;st=(9,245,tau0,9)
   for eps in bits:
    n=len(b);b=[x|(1<<n) for x in b]+[0];cut|=eps<<n
    h,al,ta,en=st;st=(h,2*al+eps*h,ta+eps*(al+en),eps*h)
   req(state(b,cut)==st,'recurrence')
   for k in [1,2]:
    actual=tables(whole(b,[],cut,0,k))[0][-1]
    req(actual==st[1]*(1<<k)+9-st[2],'closed count');checks+=1
print(json.dumps({'status':'PASS','literal_closed_checks':checks,'seeds':seed_receipts},indent=2))
