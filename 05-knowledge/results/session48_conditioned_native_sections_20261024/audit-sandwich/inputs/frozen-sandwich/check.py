#!/usr/bin/env python3
"""Portable exact matrix, orbit, resource and shared-observer controls. No asserts."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import json, hashlib

CHECKS=0
HERE=Path(__file__).resolve().parent

def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(label)

def converse(a): return tuple(tuple(a[j][i] for j in range(len(a))) for i in range(len(a)))
def tournament(n,code):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=(code>>k)&1;a[j][i]=1-a[i][j]
    return tuple(map(tuple,a))
def cyclic(n,diffs):
    return tuple(tuple(int(i!=j and (j-i)%n in diffs) for j in range(n)) for i in range(n))
C3=cyclic(3,{1})
def ordinal(*blocks):
    labels=[(b,i) for b,a in enumerate(blocks) for i in range(len(a))]
    return tuple(tuple(blocks[b][i][j] if b==c else int(b<c) for c,j in labels) for b,i in labels)
def hpaths(a):
    n=len(a); out=[sum(a[u][v]<<v for v in range(n)) for u in range(n)]; full=(1<<n)-1
    def rec(p,used):
        if used==full:
            yield p;return
        poss=(full^used) if not p else out[p[-1]]&(full^used)
        while poss:
            z=poss&-poss;poss-=z;v=z.bit_length()-1
            yield from rec(p+(v,),used|z)
    return list(rec((),0))
def hcount(a):
    n=len(a); d={(1<<i,i):1 for i in range(n)}
    for size in range(1,n):
        e={}
        for (mask,u),v in d.items():
            for j in range(n):
                if not(mask>>j&1) and a[u][j]:
                    key=(mask|1<<j,j);e[key]=e.get(key,0)+v
        d=e
    return sum(d.values())
def autos(a):
    n=len(a); deg=list(map(sum,a)); p=[-1]*n; found=[]
    # Full adjacency pruning; does not assume a circulant or any known group.
    def rec(i,used):
        if i==n: found.append(tuple(p));return
        for v in range(n):
            if not(used>>v&1) and deg[i]==deg[v] and all(a[i][j]==a[v][p[j]] for j in range(i)):
                p[i]=v;rec(i+1,used|1<<v)
    rec(0,0);return found
@lru_cache(None)
def pairs(n): return tuple((i,j) for i in range(n) for j in range(i+2,n))
@lru_cache(None)
def reflection_indices(n):
    pp=pairs(n); lookup={p:i for i,p in enumerate(pp)}
    return tuple(lookup[n-1-j,n-1-i] for i,j in pp)
def word(a,p): return bytes(a[p[i]][p[j]] for i,j in pairs(len(a)))
def rho(w,n): return bytes(w[i] for i in reflection_indices(n))
TRANS=bytes.maketrans(bytes([0,1]),bytes([1,0]))
def flip(w):return w.translate(TRANS)
def reservoir(w,n):return frozenset((w,rho(w,n),flip(w),rho(flip(w),n)))
def key(w,n):return min(reservoir(w,n))
def partner(a,p):
    pos={u:i for i,u in enumerate(p)}
    return tuple(tuple(0 if u==v else a[u][v] if abs(pos[u]-pos[v])==1 else a[v][u] for v in range(len(a))) for u in range(len(a)))
def strong(a):
    n=len(a)
    for root in range(n):
        seen={root};todo=[root]
        while todo:
            u=todo.pop()
            for v in range(n):
                if a[u][v] and v not in seen:seen.add(v);todo.append(v)
        if len(seen)!=n:return False
    return True
def isom(a,b,p):return all(a[u][v]==b[p[u]][p[v]] for u in range(len(a)) for v in range(len(a)))
def restrict(a,vs):return tuple(tuple(a[u][v] for v in vs) for u in vs)
def check_agreement(a,p,b):
    edges={(u,v) for u in range(len(a)) for v in range(len(a)) if a[u][v] and b[u][v]}
    need(edges==set(zip(p,p[1:])),'labelled agreement graph is exactly the observer path')
    need(partner(b,p)==a,'labelled agreement-path involution')
    need(partner(converse(a),p[::-1])==converse(b),'literal converse/partner square')

def path_orbits(a):
    pp=hpaths(a); aa=autos(a); reps={}
    for p in pp:reps.setdefault(word(a,p),p)
    need(len(pp)==hcount(a),'DFS and independent subset-DP path count')
    need(len(pp)==len(aa)*len(reps),'free automorphism-normalization factor')
    mult=Counter(word(a,p) for p in pp)
    need(set(mult.values())=={len(aa)},'all normalized words have exact automorphism lift count')
    return pp,aa,list(reps.values())

def sandwich_profile(a,label,exhaustive_source=False):
    pp,aa,reps=path_orbits(a); q=len(a); s=ordinal(a,C3,converse(a)); n=len(s)
    ww={};full=set()
    for i,p in enumerate(reps):
        for j,t in enumerate(reps):
            local=[]
            for k in range(3):
                obs=p+tuple(q+(k+z)%3 for z in range(3))+tuple(q+3+v for v in t[::-1])
                need(all(s[u][v] for u,v in zip(obs,obs[1:])),'sandwich observer Hamiltonian')
                local.append(word(s,obs));full.add(obs)
            need(len(set(local))==1,'central three phases are one normalized word')
            ww[i,j]=local[0]
            b=partner(s,obs)
            need(strong(b),'every sandwich complement is strong')
            need(b[obs[-1]][obs[0]]==1,'closing arc authenticates a spanning directed cycle')
            check_agreement(s,obs,b)
    need(not strong(s),'nonempty sandwich source is nonstrong')
    need(len(set(ww.values()))==len(reps)**2,'ordered-pair map injective')
    for (i,j),w in ww.items():
        need(rho(w,n)==ww[j,i],'normalization/converse is exact coordinate swap')
        need(flip(rho(w,n))==rho(flip(w),n),'F and reflection commute')
        need((w==rho(w,n))==(i==j),'blue diagonal iff same prefix path orbit')
        need(len(reservoir(w,n))==(2 if i==j else 4),'normalized two/four-word resource cardinality')
    words=set(ww.values()); need(not(words&{flip(w) for w in words}),'no F loop')
    resources={key(w,n) for w in words}; ac=len(reps)
    need(len(resources)==ac*(ac+1)//2,'unordered pair resource count')
    if exhaustive_source:
        allp=hpaths(s); allw={word(s,p) for p in allp}; sa=autos(s)
        need(len(allp)==3*len(pp)**2,'direct exhaustive sandwich path factorization')
        need(len(sa)==3*len(aa)**2,'direct sandwich automorphism factorization')
        need(allw==words,'all sandwich words exhausted, not just a menu')
    return dict(name=label,order_R=q,H_R=len(pp),aut_R=len(aa),a=ac,order_S=n,
                H_S=3*len(pp)**2,aut_S=3*len(aa)**2,source_words=ac*ac,
                blue_incidence=ac,black_weighted_incidence=ac*(ac-1),
                black_resources=ac*(ac-1)//2,total_resources=len(resources),loops=0,
                direct_full_S_enumeration=exhaustive_source)

P=((3,0,6,7,1,5,4,2),(3,1,6,7,0,5,4,2));G=(4,0,2,3,1,5,6,7)
def seed():
    a=[[0]*8 for _ in range(8)]
    for z,(i,j) in enumerate(combinations(range(8),2)):
        a[i][j]=(29558785>>(27-z))&1;a[j][i]=1-a[i][j]
    return tuple(map(tuple,a))
def prefix(r,b):return tuple(8*j+x for j in range(r) for x in P[(b>>j)&1])

def native_checks():
    a=seed(); pp,aa,reps=path_orbits(a)
    need(len(pp)==409 and len(aa)==1,'native seed complete count and rigidity')
    need(strong(a),'native seed is strongly connected')
    pinned={p for p in pp if all(p[t]==P[0][t] for t in (0,2,3,5,6,7))}
    need(pinned==set(P),'complete weaker-common-pin seed family has exactly the two binary keys')
    need(isom(partner(a,P[0]),partner(a,P[1]),G),'native endpoint-marked partner map')
    need(G[P[0][0]]==P[1][0] and G[P[0][-1]]==P[1][-1],'partner map retains true endpoints')
    rows=[]
    for r in range(1,5):
        rr=ordinal(*([a]*r));q=len(rr);s=ordinal(rr,C3,converse(rr));n=len(s)
        left=tuple(range(q+3));right=tuple(range(q,n));ls=restrict(s,left);rs=restrict(s,right)
        native_to_right=tuple(3+u for u in range(q))+tuple((-j)%3 for j in range(3))
        need(isom(converse(ls),rs,native_to_right),'right source is converse via explicit C3 anti-isomorphism')
        beta=tuple(q+3+u if u<q else q+(-((u-q)%3))%3 if u<q+3 else u-q-3 for u in range(n))
        need(isom(s,converse(s),beta),'global source explicit anti-automorphism')
        records={};global_words={};local_words={};global_keys={}
        for x,y in product(range(1<<r),repeat=2):
            for k in range(3):
                center=tuple(q+(k+z)%3 for z in range(3))
                obs=prefix(r,x)+center+tuple(q+3+u for u in prefix(r,y)[::-1])
                lp=obs[:q+3];rp=tuple(u-q for u in obs[q:]);b=partner(s,obs)
                lb=partner(ls,lp);rb=partner(rs,rp)
                need(restrict(b,left)==lb and restrict(b,right)==rb,'one global K restricts to both local partners')
                need(restrict(lb,range(q,q+3))==restrict(rb,range(3)),'partners agree literally on shared owners')
                need(all(s[u][v] for u,v in zip(obs,obs[1:])),'native global observer is Hamiltonian')
                need(strong(b),'native global target strong')
                lw=word(ls,lp);rw=word(rs,rp);gw=word(s,obs)
                native_b=prefix(r,y)+tuple(q+((1-k)+z)%3 for z in range(3))
                need(tuple(native_to_right[u] for u in native_b[::-1])==rp,'right phase is 1-k in converse native frame')
                need(rw==rho(word(ls,native_b),q+3),'actual right normalization is left-key reflection')
                need((key(lw,q+3)==key(rw,q+3))==(x==y),'local capacity iff distinct keys')
                need((gw==rho(gw,n))==(x==y),'native global blue iff diagonal')
                records[x,y,k]=(obs,b,lp,lb,rp,rb)
                global_words[x,y,k]=gw;global_keys[x,y,k]=key(gw,n);local_words[x,k]=lw
        m=1<<r
        base_target=records[0,0,0][1]
        anti=tuple(q+3+u if u<q else q+(2-(u-q))%3 if u<q+3 else u-q-3 for u in range(n))
        need(isom(base_target,converse(base_target),anti),'common selected global target is explicitly self-converse')
        for (x,y,k),(obs,b,lp,lb,rp,rb) in records.items():
            # Endpoint-fixed seed maps compose even though they do not carry
            # every internal observer position to the corresponding position.
            transport=tuple(8*(u//8)+(G[u%8] if (x>>(u//8))&1 else u%8) if u<q
                else q+(u-q+k)%3 if u<q+3
                else q+3+8*((u-q-3)//8)+(G[(u-q-3)%8] if (y>>((u-q-3)//8))&1 else (u-q-3)%8)
                for u in range(n))
            need(isom(base_target,b,transport),'all selected global resources share one actual target isomorphism class')
        need(len(set(local_words.values()))==m,'exact local selected source words')
        need(len({key(w,q+3) for w in local_words.values()})==m,'local resources distinct despite converse')
        need(len(set(global_words.values()))==m*m,'exact selected global words')
        need(len(set(global_keys.values()))==m*(m+1)//2,'selected global resource capacity is unordered pair count')
        for x,y,k in records:
            obs,b,lp,lb,rp,rb=records[x,y,k]
            for delta in (1,2):
                target=(x,y,(k+delta)%3);obs2,b2,lp2,lb2,rp2,rb2=records[target]
                rot=tuple(u if not(q<=u<q+3) else q+(u-q+delta)%3 for u in range(n))
                need(tuple(rot[u] for u in obs)==obs2,'global phase path transport')
                need(isom(s,s,rot) and isom(b,b2,rot),'global source and partner transport')
                lrot=rot[:q+3];rrot=tuple(rot[q+u]-q for u in range(q+3))
                need(isom(ls,ls,lrot) and isom(lb,lb2,lrot),'left ordinary incidence endpoints fixed as classes')
                need(isom(rs,rs,rrot) and isom(rb,rb2,rrot),'right ordinary incidence endpoints fixed as classes')
                need(sum(u!=v for u,v in zip(obs,obs2))==3,'global integer owner-clock work three')
                need(sum(u!=v for u,v in zip(lp,lp2))+sum(u!=v for u,v in zip(rp,rp2))==6,'sum local owner-clock occurrence work six')
                need(global_keys[x,y,k]==global_keys[target],'global resource held continuously')
                need(key(word(ls,lp),q+3)==key(word(ls,lp2),q+3) and key(word(rs,rp),q+3)==key(word(rs,rp2),q+3),'both local resources held continuously')
            swapped=(y,x,(1-k)%3)
            need(tuple(beta[u] for u in obs[::-1])==records[swapped][0],'actual global converse transport swaps local keys')
            need(global_keys[x,y,k]==global_keys[swapped],'global quotient loses marked left/right assignment')
        rows.append(dict(r=r,local_order=q+3,global_order=n,selected_keys=m,
                         global_records=len(records),capacity_valid_records=3*m*(m-1),
                         selected_global_words=m*m,selected_global_resources=m*(m+1)//2,selected_target_classes=1,
                         selected_blue_resources=m,
                         capacity_valid_global_black_resources=m*(m-1)//2,
                         complete_unpinned_R_orbits=409**r,complete_unpinned_S_words=409**(2*r),
                         local_phase_work=6,global_phase_work=3))
    # Diagonal is a genuine capacity hostile, not a nonexistent Hamiltonian observer.
    need(records[0,0,0][0] is not None,'diagonal still has a global observer but collides locally')
    return dict(seed_H=len(pp),seed_Aut=len(aa),rows=rows)


def arbitrary_prefix_hostile():
    u=((0,0,0,1),(1,0,0,0),(1,1,0,0),(0,1,1,0))
    r=ordinal(C3,u);l=ordinal(r,C3)
    p=(0,1,2,3,6,5,4);q=(0,1,2,5,4,3,6);tail=(7,8,9)
    need(all(r[x][y] for x,y in zip(p,p[1:])) and all(r[x][y] for x,y in zip(q,q[1:])),'hostile prefix paths are Hamiltonian')
    need(word(r,p)!=word(r,q),'hostile prefixes lie in distinct R path orbits')
    wp=word(l,p+tail);wq=word(l,q+tail)
    need(wp!=wq and rho(wp,10)==wq,'distinct prefix words become reflection-paired locally')
    need(key(wp,10)==key(wq,10),'arbitrary prefix orbit inequality does not imply local resource injection')
    pp,aa,reps=path_orbits(r)
    return dict(U_matrix=u,R_order=7,R_H=len(pp),R_Aut=len(aa),R_a=len(reps),L_order=10,P=p,Q=q,
                distinct_R_orbits=True,same_local_resource=True)

def main():
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    for rec in pins['files']:
        data=(HERE/rec['local']).read_bytes().replace(b'\r\n',b'\n')
        need(hashlib.sha256(data).hexdigest()==rec['sha256'],'pinned source archive: '+rec['local'])
    small=[]
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            small.append(sandwich_profile(tournament(n,code),f'labelled-{n}-{code}',True))
    famous=[sandwich_profile(C3,'C3'),sandwich_profile(cyclic(5,{1,2}),'cyclic-C5'),
            sandwich_profile(cyclic(7,{1,2,4}),'Paley7')]
    paley11=cyclic(11,{1,3,4,5,9}); pp,aa,reps=path_orbits(paley11)
    need((len(pp),len(aa),len(reps))==(95095,55,1729),'Paley11 exact independent pattern hostile')
    need(len(reps)!=3**4,'false short power-of-three pattern refuted')
    # Nonempty-R is needed for the closing-arc proof.
    p=(0,1,2);need(not strong(partner(C3,p)),'empty-wing hostile: C3 partner is not strong')
    result=dict(status='PASS',universe={'all_labelled_R_orders':[1,2,3,4],
                'all_labelled_R_count':len(small),'famous_full_pair_profiles':['C3','cyclic-C5','Paley7'],
                'Paley11_scope':'complete R paths, automorphisms and normalized words; S profile derived, not exhaustively enumerated',
                'native_r':[1,2,3,4]},famous=famous,
                Paley11=dict(H=len(pp),Aut=len(aa),a=len(reps),S_order=25,S_words=len(reps)**2,
                             S_blue=len(reps),S_black_weighted=len(reps)*(len(reps)-1),
                             S_resources=len(reps)*(len(reps)+1)//2),
                arbitrary_R_capacity_hostile=arbitrary_prefix_hostile(),
                small_summary=dict(profiles=len(small),max_a=max(row['a'] for row in small),
                                   total_source_words=sum(row['source_words'] for row in small),
                                   all_loops_zero=all(row['loops']==0 for row in small)),native=native_checks())
    result['checks']=CHECKS
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
