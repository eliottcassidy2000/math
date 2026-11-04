import collections,hashlib,itertools,json,pathlib
from quotient_counts import converse,endpoint,repeated_endpoint
P=pathlib.Path(__file__).resolve().parent
def load(n):return json.loads((P/n).read_text())
def need(b,msg):
    if not b:raise RuntimeError(msg)
def module(r,m):return all((r[i]&m)in(0,m)for i in range(len(r))if not(m>>i&1))
def valid(r):
    for i,x in enumerate(r):
        need(not(x>>i&1),'loop')
        for j in range(i):need(((x>>j)&1)+((r[j]>>i)&1)==1,'not tournament')
def inflate(r,cap):
    off=[];n=0
    for k in cap:off.append(n);n+=k
    out=[0]*n
    for i,k in enumerate(cap):
        for u in range(k):
            out[off[i]+u]=sum(1<<(off[i]+v)for v in range(u+1,k))
            for j,l in enumerate(cap):
                if r[i]>>j&1:out[off[i]+u]|=((1<<l)-1)<<off[j]
    return out
selected=load('selected.json');weighted=load('selected.out');runword=load('runword.out');need(len(selected)==len(weighted)==len(runword)==4,'four certificates')
certificates=[];run_skews=[]
for q,w,r in zip(selected,weighted,runword):
    Q=q['rows'];valid(Q);need(len(Q)==9,'quotient order')
    need(not any(module(Q,m)for m in range(1,511)if 1<m.bit_count()<9),'quotient prime')
    need(any(x.bit_count()!=4 for x in Q),'noncirculant degree witness')
    E=endpoint(Q);S=endpoint(converse(Q));need(E==q['E']and S==q['S'],'endpoints')
    need(not any((a+b)%2 for a,b in zip(S,E)),'p zero')
    skew=[[0]*9 for _ in range(9)]
    for j in range(9):
        a=repeated_endpoint(converse(Q),j);b=repeated_endpoint(Q,j)
        for i in range(9):skew[i][j]=a[i]-b[i]
    sig=[(S[i]-E[i]+2*sum(skew[i]))%4 for i in range(9)]
    need(sig==q['sigma']==[0]*9,'complete sigma zero')
    run_skews.append(dict(node=q['node'],start_run_skews=skew))
    witnesses=[]
    for i,j in itertools.permutations(range(1,9),2):
        if Q[0]>>i&1 and Q[i]>>j&1 and Q[j]&1 and all(((Q[u]>>i)&1)==((Q[u]>>j)&1)for u in range(1,9)if u not in(i,j)):witnesses.append([i,j])
    need(witnesses==q['old_ordered_witnesses_at_z0']==[],'ordered witness absent')
    need(q['contracted_rows']==inflate(Q,[1]+[2]*8),'contracted physical rows')
    g=q['certificate_group'];k=w['skews'][g]
    need(len(k)==16 and len(w['skews'])==9,'complete weighted capacities')
    need(r['target']==g and r['skew']==k+[0],'independent run-word derivative')
    need(k[0]>0 and all(x>=0 for x in k),'coefficientwise positivity')
    certificates.append(dict(node=q['node'],group=g,k_coefficients=k,strict_constant=k[0]))
need(load('quotient_run_skews.json')==run_skews,'saved full run skew matrix')
controls=load('compression_controls.json');vals=list(map(int,(P/'compression_output.txt').read_text().split()));need(len(controls)==len(vals)==100,'compression controls')
for c,value in zip(controls,vals):
    need(c['physical_rows']==inflate(c['rows'],c['capacities']),'compression physical rows')
    need(sum(endpoint(c['physical_rows']))==c['H']==value,'independent literal physical Hamilton count')
physical=load('physical_queries.json');values=list(map(int,(P/'physical_values.txt').read_text().split()));need(len(values)==len(physical['queries'])==208,'physical query count')
for q in physical['queries']:
    valid(q['physical_rows']);need(q['physical_rows']==inflate(q['compressed_rows'],q['capacities']),'physical-port expansion')
    n=len(q['physical_rows'])-2
    need(q['physical_rows'][n]>>(n+1)&1,'port order')
    need(sum(((q['physical_rows'][v]>>n)&1)<<v for v in range(n))==q['a_mask'],'physical a mask')
    need(sum(((q['physical_rows'][v]>>(n+1))&1)<<v for v in range(n))==q['c_mask'],'physical c mask')
idx={q['node']:i for i,q in enumerate(selected)};actual=[]
for c in physical['cases']:
    got=sum(t['weight']*sum(q['sign']*(values[q['source']]-values[q['target']])for q in t['queries'])for t in c['terms'])
    k=weighted[idx[c['node']]]['skews'][c['group']]
    want=2*(3*k[0]+k[1]if c['ambient']else k[0]);need(got==want and got>0,'physical grouped mixed binomial')
    actual.append(dict(node=c['node'],m=c['m'],ambient=c['ambient'],contrast=got))
bank=load('small_bank.json');ind=list(map(int,(P/'small_values.txt').read_text().split()));pos=0;degrees=collections.Counter()
for b in bank:
    masks=b['masks'];L=len(masks);need(ind[pos:pos+L*L]==b['values'],'independent small full response counts');pos+=L*L;index={v:i for i,v in enumerate(masks)}
    for i,a in enumerate(masks):
        for j,c in enumerate(masks):
            z=sum((-1)**((a^u).bit_count()+(c^v).bit_count())*b['values'][index[u]*L+index[v]]for u in masks if u&a==u for v in masks if v&c==v)
            need(z==b['coefficients'][i*L+j],'full named coefficient inversion');degrees[a.bit_count()+c.bit_count()]+=1
need(pos==len(ind),'full small bank count')
small=load('small_controls.json');sw=load('small_weighted.out');small_actual=[]
def coef(g,a,c):
    b=bank[g];ix={v:i for i,v in enumerate(b['masks'])};return b['coefficients'][ix[a]*len(ix)+ix[c]]
for c,k in zip(small['cases'],sw):
    s,t,l=c['graph_ids'];m=c['m'];terms=[(m-1,1)]+([(2,3)]if m==5 else [])
    for g,names in enumerate(c['groups'],1):
        ds=sum(w*sum(coef(s,3,(1<<r)|(1<<i))-coef(t,3,(1<<r)|(1<<i))for i in names)for r,w in terms)
        dl=sum(w*sum(coef(s,3,(1<<r)|(1<<i))-coef(l,3,(1<<r)|(1<<i))for i in names)for r,w in terms)
        kp=k['skews'][g];pred=2*(3*kp[0]+kp[1]if c['ambient']else kp[0]);need(ds==dl==pred,'small full-bank restoration and ambient transfer')
        small_actual.append(dict(graph=s,group=g,contrast=ds))
print(json.dumps(dict(status='PASS',certificates=certificates,physical_actual_queries=len(values),physical_controls=actual,independent_compression_controls=100,small_complete_graphs=len(bank),small_complete_named_coefficients=pos,small_degree_counts=dict(sorted(degrees.items())),small_mixed_controls=small_actual,scope='General sign-definite excursion-polynomial sufficient criterion; four explicit prime-nine odd-transitive-root families with eight TT2 children, every fixed ambient. No unrestricted prime-nine closure or response collision.'),indent=2,sort_keys=True))
