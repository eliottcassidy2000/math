"""A000568 by parent-only attachment states and root-switch rows.

No child canonicalization or supplied sequence values enter generation.
Exact directed colour refinement is a filter, followed by exact isomorphism.
Run from anywhere; max order defaults to 8. Stdlib only.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
import gzip
from math import factorial, gcd
from pathlib import Path
from time import monotonic

CHECKS = 0
def need(b, msg):
    global CHECKS
    CHECKS += 1
    if not b: raise RuntimeError(msg)

def colours(a):
    n = len(a)
    if not n: return (),()
    c = [x.bit_count() for x in a]
    for _ in range(n):
        keys = [(c[i], tuple(sorted(c[j] for j in range(n) if a[i]>>j&1))) for i in range(n)]
        ids = {x:i for i,x in enumerate(sorted(set(keys)))}
        d = [ids[x] for x in keys]
        if len(set(c))==len(set(d)):
            return tuple(d), tuple(sorted(keys))
        c = d
    raise RuntimeError('refinement failed to stabilize')

def isomorphisms(a,b,all_maps=False):
    ca,ka = colours(a); cb,kb = colours(b)
    if ka!=kb: return []
    n=len(a); mapping=[-1]*n; answers=[]
    candidates=[[j for j in range(n) if ca[i]==cb[j]] for i in range(n)]
    def visit(todo, used):
        if not todo:
            answers.append(tuple(mapping))
            return not all_maps
        choices=[]
        for i in todo:
            js=[j for j in candidates[i] if not (used>>j&1) and all(((a[i]>>u)&1)==((b[j]>>mapping[u])&1) for u in range(n) if mapping[u]>=0)]
            if not js: return False
            choices.append((len(js),i,js))
        _,i,js=min(choices)
        for j in js:
            mapping[i]=j
            if visit([u for u in todo if u!=i],used|(1<<j)):
                mapping[i]=-1
                return True
            mapping[i]=-1
        return False
    visit(list(range(n)),0)
    return answers

def move_mask(mask,p):
    return sum(1<<p[i] for i in range(len(p)) if mask>>i&1)

def extend(a,s):
    n=len(a)
    return tuple(row|((1-(s>>i&1))<<n) for i,row in enumerate(a))+(s,)

def delete(a,v):
    rest=[i for i in range(len(a)) if i!=v]
    return tuple(sum(((a[i]>>j)&1)<<k for k,j in enumerate(rest)) for i in rest),sum(((a[v]>>j)&1)<<k for k,j in enumerate(rest))

def hamilton(a):
    n=len(a)
    if not n: return 1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,1<<n):
        for i in range(n):
            if dp[s][i]:
                choices=a[i]&~s
                while choices:
                    bit=choices&-choices; choices^=bit
                    dp[s|bit][bit.bit_length()-1]+=dp[s][i]
    return sum(dp[-1])

def is_prime(a):
    n=len(a)
    if n<3: return False
    for s in range(1,(1<<n)-1):
        if s.bit_count()<2: continue
        if all((a[v]&s) in (0,s) for v in range(n) if not (s>>v&1)):
            return False
    return True

def odd_parts(n,least=1):
    if not n:
        yield ()
    for j in range(least,n+1,2):
        for tail in odd_parts(n-j,j): yield (j,)+tail

def burnside(n):
    ans=Fraction(0)
    for parts in odd_parts(n):
        multiplicities=Counter(parts)
        z=1
        for k,m in multiplicities.items(): z*=k**m*factorial(m)
        exponent=sum((k-1)//2 for k in parts)+sum(gcd(a,b) for a,b in combinations(parts,2))
        ans+=Fraction(1<<exponent,z)
    need(ans.denominator==1,'Burnside integrality')
    return int(ans)

def generate(parents,N):
    start=monotonic(); n=N-1
    auts=[isomorphisms(a,a,True) for a in parents]
    need(all(len(g)%2==1 for g in auts),'tournament odd automorphisms')
    canonmasks=[tuple(min(move_mask(s,g) for g in group) for s in range(1<<n)) for group in auts]
    buckets=defaultdict(list)
    for k,a in enumerate(parents): buckets[colours(a)[1]].append(k)
    @lru_cache(None)
    def classify(a):
        for k in buckets[colours(a)[1]]:
            maps=isomorphisms(a,parents[k])
            if maps:
                p=maps[0]
                need(all(((a[i]>>j)&1)==((parents[k][p[i]]>>p[j])&1) for i in range(n) for j in range(n)),'invalid parent transport')
                return k,p
        raise RuntimeError('missing deletion parent')
    states=[(k,s) for k,masks in enumerate(canonmasks) for s in sorted(set(masks))]
    stateid={s:i for i,s in enumerate(states)}
    rows=[]; representatives={}; root_auts={}; trace=0
    for idx,(k,s) in enumerate(states):
        child=extend(parents[k],s); row=Counter()
        for v in range(N):
            deletion,mask=delete(child,v)
            parent,p=classify(deletion)
            key=(parent,canonmasks[parent][move_mask(mask,p)])
            row[stateid[key]]+=1
        signature=tuple(sorted(row.items()))
        need(sum(row.values())==N and all(w%2 for w in row.values()),'row mass / odd vertex orbit')
        trace+=row[idx]; rows.append(signature)
        representatives.setdefault(signature,child)
        if signature not in root_auts:
            root_auts[signature]=row[idx]*sum(move_mask(s,g)==s for g in auts[k])
        if N>=9 and (idx+1)%100000==0:
            print(json.dumps({'order':N,'processed_pointed_states':idx+1,'total':len(states),'seconds':round(monotonic()-start,1)}),flush=True)
    for i,row in enumerate(rows):
        need(all(rows[j]==row for j,w in row),'root-switch block law')
        need(any(j==i for j,w in row),'missing own root')
    children=list(representatives.values())
    need(trace==N*len(children),'trace identity')
    need(len(children)==burnside(N),'independent Burnside count')
    # New parents are chosen before this audit. No full-child canon is used.
    hist=Counter(); aut_hist=Counter(root_auts.values()); prime_count=None
    if N<=8:
        prime_count=0
        for signature,a in representatives.items():
            aut=len(isomorphisms(a,a,True)); paths=hamilton(a)
            need(aut==root_auts[signature],'root stabilizer transport')
            need(paths%aut==0,'path action is free')
            hist[paths//aut]+=1
            prime_count+=is_prime(a)
        need(sum(k*v for k,v in hist.items())==1<<((N-1)*(N-2)//2),'fixed-path cube mass')
    result={'order':N,'pointed_states':len(states),'trace':trace,'classes':len(children),'prime_classes':prime_count,'rank_over_F2':len(children),'rank_Q_squared_over_F2':len(children)*(N%2),'tc_histogram':dict(sorted(hist.items())),'automorphism_order_histogram':dict(sorted(aut_hist.items())),'parent_cache_size':classify.cache_info().currsize,'seconds':round(monotonic()-start,3)}
    print(json.dumps({k:v for k,v in result.items() if 'histogram' not in k}),flush=True)
    return children,result

def hostile():
    def block(kind):
        a=[2,4,1]
        return tuple(x|(8 if kind=='A' else 0) for x in a)+(0 if kind=='A' else 7,)
    def compose(word):
        a=[]
        for k,kind in enumerate(word):
            for row in block(kind): a.append((row<<(4*k))|(15<<(4*((k+1)%3))))
        return tuple(a)
    def cyclic(a):
        return [c for c in combinations(range(12),3) if all(sum(a[v]>>u&1 for u in c)==1 for v in c)]
    u,v=compose('AAA'),compose('BAA')
    need(cyclic(u)==cyclic(v),'labelled triples hostile')
    au,av=len(isomorphisms(u,u,True)),len(isomorphisms(v,v,True))
    hu,hv=hamilton(u),hamilton(v)
    need((au,av,hu,hv)==(81,27,462915,462915),'hostile invariant')
    need(all(u[i]>>(i+1)&1 for i in range(11)) and not all(v[i]>>(i+1)&1 for i in range(11)),'named Hamilton guard')
    return {'rows_AAA':u,'rows_BAA':v,'labelled_cyclic_triples':len(cyclic(u)),'automorphism_orders':[au,av],'Hamilton_paths':[hu,hv],'identity_order_legal':[True,False]}

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--max-order',type=int,default=8); args=parser.parse_args()
    need(1<=args.max_order<=9,'explicit finite universe 1..9')
    parents=[()]; results=[]
    for N in range(1,args.max_order+1):
        parents,r=generate(parents,N); results.append(r)
    out={'status':'FINITE-EXACT','generation':'parent attachment orbits and identical root-switch rows; no child canon/count input','results':results,'hostile_81_27':hostile(),'checks':CHECKS}
    Path(__file__).with_name('root_switch.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    bank=(json.dumps({'order':args.max_order,'rows':parents},separators=(',',':'))+'\n').encode()
    Path(__file__).with_name('tournament_bank.json.gz').write_bytes(gzip.compress(bank,mtime=0))
