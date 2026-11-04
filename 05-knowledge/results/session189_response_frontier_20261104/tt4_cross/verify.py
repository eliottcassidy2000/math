import collections,hashlib,itertools,json,pathlib
from exterior import capacities,exterior,mul,plus
from prepare import ORDERS,COMBOS,CONSTANTS,kernels
P=pathlib.Path(__file__).resolve().parent
def load(n):return json.loads((P/n).read_text())
def need(b,msg):
    if not b:raise RuntimeError(msg)
def poly(p,tau):
    out={};power={0:1}
    for z in p:out=plus(out,power,z);power=mul(power,tau)
    return out
def deriv(p):return[(i+1)*p[i+1]for i in range(len(p)-1)]
inner=load('inner_kernels.json');lookup={tuple(x['order']):{(v['p'],v['q'],v['r']):(v['L'],v['R'])for v in x['coefficients']}for x in inner}
for order,items in lookup.items():
    for (p,q,r),(L,R)in items.items():need((L,R)==kernels(order,p,q,r),'complete literal kernel')
certificates=[]
for target,combo,const in zip(ORDERS[1:],COMBOS,CONSTANTS):
    B=tuple(reversed(target));lp=[0]*3;rp=[0]*3;terms=[]
    for p,q,r,w in combo:
        L,R=lookup[ORDERS[0]][p,q,r];l,rr=lookup[B][p,q,r];dp=[x-y for x,y in zip(L,l)];dq=[x-y for x,y in zip(R,rr)]
        lp=[x+w*y for x,y in zip(lp,dp)];rp=[x+w*y for x,y in zip(rp,dq)];terms.append(dict(p=p,q=q,r=r,weight=w,P=dp,Q=dq))
    need(lp==rp==[const,0,0],'constant symmetric kernel certificate')
    certificates.append(dict(target=target,local_replacement=B,terms=terms,constant=const))
bank=load('bank.json');controls=load('controls.json');values=list(map(int,(P/'independent_values.txt').read_text().split()));degree=collections.Counter();pos=0
for b in bank:
    masks=b['masks'];L=len(masks);index={x:i for i,x in enumerate(masks)};need(values[pos:pos+L*L]==b['values'],'independent complete physical counts');pos+=L*L
    for i,a in enumerate(masks):
        for j,c in enumerate(masks):
            x=sum((-1)**((a^u).bit_count()+(c^v).bit_count())*b['values'][index[u]*L+index[v]]for u in masks if u&a==u for v in masks if v&c==v)
            need(x==b['coefficients'][i*L+j],'full named Mobius inversion');degree[a.bit_count()+c.bit_count()]+=1
need(pos==len(values),'complete data size')
def coef(g,p,q,r,i):
    b=bank[g];index={x:j for j,x in enumerate(b['masks'])};return b['coefficients'][index[(1<<p)|(1<<q)]*len(index)+index[(1<<r)|(1<<i)]]
individual=[];capacities_saved=[];outcomes=[]
for caseid,case in enumerate(controls['cases']):
    src=bank[case['source']]['rows'];n=len(src)-4;R=[x>>4 for x in src[4:]];into=sum((x&1)<<i for i,x in enumerate(src[4:]));alpha,beta,tau,S,E,ell,rho=exterior(R,into);full=(1<<n)-1
    V=case['virtual'][:]
    if case['ambient']:
        N=len(V);V=[x|(1<<N)for x in V]+[1<<(N+1),(1<<N)-1]
    C=capacities(V);G=[i-3 for i in case['group']];sym=sum(C[0][i]+C[i][0]for i in G);need(sym>0,'strict symmetric capacity')
    capacities_saved.append(dict(case=caseid,rows=V,matrix=C,group=G,symmetric_capacity=sym))
    for h,(target,combo,const)in enumerate(zip(ORDERS[1:],COMBOS,CONSTANTS)):
        local=case['local_replacements'][h];B=tuple(reversed(target))
        for p,q,r in lookup[ORDERS[0]]:
            L,RR=lookup[ORDERS[0]][p,q,r];l,rr=lookup[B][p,q,r];dp=[x-y for x,y in zip(L,l)];dq=[x-y for x,y in zip(RR,rr)]
            for i in range(n):
                f=plus(mul(mul(alpha,S[i]),poly(dp,tau)),mul(mul(beta,E[i]),poly(dq,tau)))
                f=plus(f,mul(mul(alpha,beta),plus(mul(ell[i],poly(deriv(dp),tau)),mul(rho[i],poly(deriv(dq),tau)))))
                got=coef(case['source'],p,q,r,i+4)-coef(local,p,q,r,i+4);need(got==f.get(full,0),'general two-kernel replacement identity')
                individual.append([caseid,h,p,q,r,i+4,got])
        actual=sum(w*sum(coef(case['source'],p,q,r,i)-coef(case['targets'][h],p,q,r,i)for i in case['group'])for p,q,r,w in combo)
        loc=sum(w*sum(coef(case['source'],p,q,r,i)-coef(local,p,q,r,i)for i in case['group'])for p,q,r,w in combo)
        need(actual==loc==const*sym,'quotient restoration, fixed ambient and symmetric constant')
        outcomes.append(dict(case=caseid,target=h,contrast=actual,symmetric_capacity=sym,constant=const))
need(load('individual_quartics.json')==individual,'all individual marked contrast data')
need(load('capacity_controls.json')==capacities_saved,'complete capacity data')
physical=load('physical_queries.json');actual_values=list(map(int,(P/'physical_values.txt').read_text().split()));need(len(actual_values)==len(physical['queries']),'physical query size')
for q,v in zip(physical['queries'],actual_values):
    b=bank[q['graph']];L=len(b['masks']);expect=sum(b['coefficients'][i*L+j]for i,a in enumerate(b['masks'])if a&q['a']==a for j,c in enumerate(b['masks'])if c&q['c']==c);need(v==expect,'actual nonuniform query count')
for c in physical['cases']:
    got=sum(t['weight']*sum(q['sign']*(actual_values[q['source']]-actual_values[q['target']])for q in t['queries'])for t in c['terms'])
    expected=next(x['contrast']for x in outcomes if x['case']==c['case']and x['target']==c['target']);need(got==expected,'actual physical finite difference')
print(json.dumps(dict(status='PASS',kernel_certificates=certificates,full_named_graphs=len(bank),full_named_coefficients=pos,degree_counts=dict(sorted(degree.items())),individual_replacement_quartics=len(individual),physical_queries=len(actual_values),outcomes=outcomes,scope='Crossed TT4 target branch, profile-preserving other children and every fixed ambient. General all-transitive exclusion requires separately reviewed incoming common-end theorem.'),indent=2,sort_keys=True))
