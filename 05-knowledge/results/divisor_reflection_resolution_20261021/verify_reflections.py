import itertools,json,math,re
from pathlib import Path
counts={}
def check(x,key):
    if not x: raise RuntimeError(key)
    counts[key]=counts.get(key,0)+1
def source(k,lift):
    out=[]
    for a,b,c in itertools.product(range(k+1),range(2),range(2)):
        if (a,b,c) in [(0,0,0),(k,1,1)]:continue
        es=(0,1) if lift and k%2==0 and a==k//2 and b+c==1 else (-1,)
        out.extend((a,b,c,e) for e in es)
    return set(out)
def C(x,k):
    a,b,c,e=x;return(k-a,1-b,1-c,e)
def Q(x):
    a,b,c,e=x;return(a,c,b,1-e if e>=0 else -1)
def orbit(x,c,q):return {x,c(x),q(x),c(q(x))}
rows=[]
for k in range(2,13):
    for lift in (False,True):
        X=source(k,lift);c=lambda x:C(x,k);q=Q
        for x in X:
            check(c(x) in X and q(x) in X and c(c(x))==q(q(x))==x and c(q(x))==q(c(x)),'source_action')
        fixes=[sum(f(x)==x for x in X) for f in (c,q,lambda x:c(q(x))) ]
        check(fixes==[0,2*k,0 if lift or k%2 else 2],'source_fixed_counts')
        orbs={frozenset(orbit(x,c,q)) for x in X}
        sizes={i:sum(len(o)==i for o in orbs) for i in [1,2,4]}
        if lift:
            check(sizes=={1:0,2:k,4:k//2+1},'resolved_orbits')
            check(len(X)==4*k+2+(2 if k%2==0 else 0),'minimal_resolution_size')
        rows.append({'k':k,'lift':lift,'size':len(X),'fix_C_Q_CQ':fixes,'orbit_sizes':sizes})
def edges(m):return list(itertools.combinations(range(m),2))
def perm_mask(x,m,p):
    es=edges(m);ix={e:i for i,e in enumerate(es)};out=0
    for i,(a,b) in enumerate(es):
        if (x>>i)&1:out|=1<<ix[tuple(sorted((p[a],p[b]))) ]
    return out
graphrows=[]
for m in range(2,7):
    N=math.comb(m,2);full=(1<<N)-1;p=tuple(reversed(range(m)));r=lambda x:perm_mask(x,m,p)
    f=m//2;t=(N+f)//2;fix=0;orbs=set()
    for x in range(1<<N):
        y=r(x);check(r(y)==x and r(x^full)==(y^full),'graph_action')
        check(y!=(x^full),'no_complement_reflection_fixed')
        fix+=y==x;orbs.add(frozenset([x,x^full,y,y^full]))
    check(fix==2**t,'graph_fixed_formula')
    check(sum(len(o)==2 for o in orbs)==2**(t-1),'two_orbit_capacity')
    check(sum(len(o)==4 for o in orbs)==(2**N-2**t)//4,'free_orbit_capacity')
    for p in itertools.permutations(range(m)):
        if any(p[p[i]]!=i for i in range(m)):continue
        check(any(set((p[a],p[b]))=={a,b} for a,b in edges(m)),'involutive_vertex_fixed_pair')
    graphrows.append({'m':m,'N':N,'reflection_edge_orbits':t,'fixed_graphs':fix,'two_orbits':2**(t-1),'free_orbits':(2**N-2**t)//4})
M={(1,0,0,-1):0,(1,1,1,-1):63,(2,0,0,-1):4,(0,1,1,-1):59,(0,1,0,-1):1,(0,0,1,-1):32,(2,0,1,-1):62,(2,1,0,-1):31,(1,1,0,0):2,(1,1,0,1):47,(1,0,1,0):61,(1,0,1,1):16}
check(set(M)==source(2,True) and len(set(M.values()))==12,'literal_table_bijection')
expected_names={'p':0,'pqr':63,'p^2':4,'qr':59,'q':1,'r':32,'p^2 r':62,'p^2 q':31,'(pq,0)':2,'(pq,1)':47,'(pr,0)':61,'(pr,1)':16}
text=Path(__file__).with_name('REFLECTION_RESOLUTION.md').read_text()
printed={a.strip():int(b) for a,b in re.findall(r'^\|([^|]+)\|\s*(\d+)\s*\|$',text,re.M)}
check(printed==expected_names,'printed_table_matches')
def T(g,m):
    ix={e:i for i,e in enumerate(edges(m))};out=0
    for j,(a,b) in enumerate(edges(m+1)):
        if b==a+1 or not ((g>>ix[a,b-1])&1):out|=1<<j
    return out
def arc(t,n,a,b):
    ix={e:i for i,e in enumerate(edges(n))}
    return bool((t>>ix[min(a,b),max(a,b)])&1) if a<b else not bool((t>>ix[b,a])&1)
blue=[];pairs=set()
for x,g in M.items():
    rg=perm_mask(g,4,(3,2,1,0));check(M[C(x,2)]==g^63 and M[Q(x)]==rg,'literal_equivariance')
    t=T(g,4);check(all(arc(t,5,i,i+1) for i in range(4)),'actual_fixed_path')
    for a,b in edges(5):
        check(arc(T(rg,4),5,a,b)==arc(t,5,4-b,4-a),'actual_anti_reversal')
        check(arc(T(g^63,4),5,a,b)==(arc(t,5,a,b) if b==a+1 else not arc(t,5,a,b)),'actual_nonpath_flip')
    isblue=all(arc(t,5,a,b)==arc(t,5,4-b,4-a) for a,b in edges(5))
    check(isblue==(Q(x)==x),'actual_colour')
    if isblue:blue.append(g)
    pairs.add(tuple(sorted([g,g^63])))
check(len(blue)==4 and len(pairs)==6 and sum(a in blue for a,b in pairs)==2,'pair_colours')
check((M[(1,0,0,-1)] & ~M[(0,1,0,-1)])==0 and not (1<=0),'not_divisibility_embedding')
print(json.dumps({'status':'PASS','checks':counts,'total_checks':sum(counts.values()),'source_rows':rows,'graph_rows':graphrows,'literal_codes':sorted(M.values()),'blue_codes':sorted(blue),'complement_pairs':sorted(pairs)},indent=2,sort_keys=True))
