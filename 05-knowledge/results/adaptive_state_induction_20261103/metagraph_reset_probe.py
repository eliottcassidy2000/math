#!/usr/bin/env python3
"""Exact named native resets and small complete ordinary carriers.
No repository imports. Bit1 means earlier labelled vertex beats later.
"""
from itertools import combinations,permutations
from collections import Counter,deque
from functools import lru_cache
from math import comb
import json


def need(p,m):
    if not p:raise RuntimeError(m)


@lru_cache(None)
def pairs(n):return tuple(combinations(range(n),2))


def decode(n,c):
    T=[0]*n
    for k,(i,j) in enumerate(pairs(n)):
        u,v=(i,j) if c>>k&1 else (j,i);T[u]|=1<<v
    return tuple(T)


def encode(T,p=None):
    if p is None:p=tuple(range(len(T)))
    return sum(((T[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(pairs(len(T))))


def converse(T):return tuple(((1<<len(T))-1)^r^(1<<i) for i,r in enumerate(T))


def relabel(T,p):return tuple(sum(((T[u]>>v)&1)<<j for j,v in enumerate(p)) for u in p)


def paths(T):
    return tuple(p for p in permutations(range(len(T))) if all(T[u]>>v&1 for u,v in zip(p,p[1:])))


@lru_cache(None)
def H(T):
    n=len(T);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for S in range(1,1<<n):
        for v,c in enumerate(D[S]):
            if c:
                A=T[v]&~S
                while A:
                    b=A&-A;A-=b;D[S|b][b.bit_length()-1]+=c
    return sum(D[-1]) if n else 1


def native(T,p):
    need(sorted(p)==list(range(len(T))) and all(T[u]>>v&1 for u,v in zip(p,p[1:])),'named native guard')
    edges={frozenset((a,b)) for a,b in zip(p,p[1:])}
    R=list(T)
    for i,j in pairs(len(T)):
        if frozenset((i,j)) not in edges:R[i]^=1<<j;R[j]^=1<<i
    return tuple(R)


def ear(T,S):
    """New vertex beats old set S, all old names retained."""
    n=len(T)
    return tuple(r|(0 if S>>i&1 else 1<<n) for i,r in enumerate(T))+(S,)


def birth(T,v):return ear(T,((1<<len(T))-1)^(1<<v))


def induced(T,V):return tuple(sum(((T[u]>>v)&1)<<j for j,v in enumerate(V)) for u in V)


def root_profile(T,v):
    n=len(T);ps=paths(T)
    F=sum(all(T[p[i]]>>p[i+1]&1 for i in range(n-1) if i!=p.index(v)) for p in permutations(range(n)))
    return H(T),sum(p[0]==v for p in ps),H(induced(T,tuple(i for i in range(n) if i!=v))),F


def cover_profile(T):
    n=len(T);out=[0]*n
    for p in permutations(range(n)):
        bad=sum(not(T[u]>>v&1) for u,v in zip(p,p[1:]))
        for k in range(bad+1,n+1):out[k-1]+=comb(n-1-bad,k-1-bad)
    return tuple(out)


def degree_canon(T):
    cells={}
    for i,row in enumerate(T):cells.setdefault(row.bit_count(),[]).append(i)
    from itertools import product
    best=min((encode(T,sum(blocks,())),sum(blocks,())) for blocks in product(*(permutations(cells[d]) for d in sorted(cells))))
    return decode(len(T),best[0]),best[1]


def same_H_component(T):
    h=H(T);seed,_=degree_canon(T);todo=deque([seed]);seen={seed};rows=[]
    while todo:
        S=todo.popleft();out=Counter()
        for p in paths(S):
            U=native(S,p)
            if H(U)!=h:continue
            D,_=degree_canon(U);out[D]+=1
            if D not in seen:seen.add(D);todo.append(D)
        a=attachment(S)
        rows.append({'rows':S,'Aut':sum(encode(S,p)==encode(S) for p in permutations(range(len(S)))),
                     'proper_modules':[M for M in modules(S) if M!=(1<<len(S))-1],
                     'outgoing_path_masses':[{'rows':D,'W':count} for D,count in sorted(out.items())],
                     'ternary_hits':{str(x):[c for c,z in a.items() if z==x] for x in (3*h,3*h+2,3*h+4)}})
    return rows


def two_step_returns(T):
    h=H(T);mids=set();ends={};repairs=[]
    for P in paths(T):
        U=native(T,P);C,pi=degree_canon(U)
        if C in mids:continue
        mids.add(C)
        for Qref in paths(C):
            Q=tuple(pi[v] for v in Qref);V=native(U,Q)
            if H(V)!=h:continue
            D,_=degree_canon(V)
            if D in ends:continue
            a=attachment(V);hits={str(z):[c for c,v in a.items() if v==z] for z in (3*h,3*h+2,3*h+4)}
            ends[D]=hits
            if any(hits.values()):
                for z,cuts in hits.items():
                    if cuts:need(len(paths(ear(V,cuts[0])))==int(z),'literal eight-vertex repair')
                repairs.append({'first_path':P,'second_path':Q,'middle_rows':U,'middle_H':H(U),'final_rows':V,'hits':hits})
    return {'first_step_ordinary_classes':len(mids),'return_ordinary_classes':len(ends),
            'endpoint_classes':[{'rows':D,'hits':hits} for D,hits in sorted(ends.items())],'repairs':repairs}


def attachment(T):return {m:H(ear(T,m)) for m in range(1<<len(T))}


def modules(T):
    n=len(T)
    return [M for M in range(1,1<<n) if M.bit_count()>1 and
        all(((T[v]&M)==0 or (T[v]&M)==M) for v in range(n) if not M>>v&1)]


def reverse_module(T,M):
    R=list(T)
    for i,j in pairs(len(T)):
        if M>>i&1 and M>>j&1:R[i]^=1<<j;R[j]^=1<<i
    return tuple(R)


def ordinary_bank(n):
    bank={};orbit={}
    for c in range(1<<len(pairs(n))):
        if c in orbit:continue
        T=decode(n,c);bank[c]=T
        for p in permutations(range(n)):
            # Current labels i map to canonical representative labels p[i].
            orbit.setdefault(encode(T,p),(c,p))
    need(len(bank)=={3:2,4:4,5:12,6:56}[n],'complete ordinary universe')
    return bank,orbit


def class_path(adj,a,b,predicate=lambda c:True):
    q=deque([a]);parent={a:None}
    while q:
        c=q.popleft()
        if c==b:
            out=[]
            while parent[c] is not None:
                prev,p=parent[c];out.append((prev,c,p));c=prev
            return list(reversed(out))
        for z,p in adj[c].items():
            if z not in parent and predicate(z):parent[z]=(c,p);q.append(z)
    return None


def literal_route(T,route,orbit):
    out=[]
    for a,b,p in route:
        c,transport=orbit[encode(T)];need(c==a,'route source class')
        inv={v:i for i,v in enumerate(transport)}
        named=tuple(inv[v] for v in p);U=native(T,named)
        need(orbit[encode(U)][0]==b,'route target class')
        out.append({'path':named,'source_rows':T,'target_rows':U,'H':[H(T),H(U)]});T=U
    return T,out


def main():
    result={'status':'FINITE-EXACT; declared complete small carriers and selected named sources'}
    layers=[];last=None
    for n in range(3,7):
        bank,orbit=ordinary_bank(n);hc={c:H(T) for c,T in bank.items()};adj={c:{} for c in bank}
        for c,T in bank.items():
            ps=paths(T);need(len(ps)==hc[c],'literal Hamilton count')
            for p in ps:adj[c].setdefault(orbit[encode(native(T,p))][0],p)
        fibres={h:sorted(c for c in bank if hc[c]==h) for h in sorted(set(hc.values()))}
        comps={}
        for h,verts in fibres.items():
            unseen=set(verts);parts=[]
            while unseen:
                seed=min(unseen);q=[seed];seen={seed};unseen.remove(seed)
                while q:
                    c=q.pop()
                    for z in adj[c]:
                        if z in unseen and hc[z]==h:unseen.remove(z);seen.add(z);q.append(z)
                parts.append(sorted(seen))
            comps[h]=parts
        layers.append({'n':n,'ordinary_classes':len(bank),'H_fibres':fibres,'strict_H_native_components':comps,
            'strict_H_support_edges':[[c,z] for c in adj for z in adj[c] if c<=z and hc[c]==hc[z]]})
        if n==6:last=bank,orbit,hc,adj
    result['complete_layers']=layers
    bank,orbit,hc,adj=last
    A=(40,17,3,6,13,30);B=(34,24,3,5,13,30)
    ca=orbit[encode(A)][0];cb=orbit[encode(B)][0]
    route=class_path(adj,ca,cb);need(route is not None,'ordinary native route found')
    final,steps=literal_route(A,route,orbit)
    result['H37_hostile']={'A_rows':A,'A_root':2,'B_rows':B,'B_root':3,'ordinary_codes':[ca,cb],
        'proper_modules_A':[M for M in modules(A) if M!=(1<<6)-1],
        'proper_modules_B':[M for M in modules(B) if M!=(1<<6)-1],
        'strict_H_route':class_path(adj,ca,cb,lambda c:hc[c]==37),
        'shortest_ordinary_native_route':steps,'final_to_B_isomorphisms':[p for p in permutations(range(6)) if relabel(B,p)==final],
        'A_next_root_H':[H(birth(A,v)) for v in range(6)],'B_next_root_H':[H(birth(B,v)) for v in range(6)]}
    h37=[c for c in bank if hc[c]==37];rooted={}
    need(all(sum(encode(bank[c],p)==c for p in permutations(range(6)))==1 for c in h37),'rigid H37 root gauges')
    for c in h37:
        for v in range(6):
            rooted[c,v]={}
            for p in paths(bank[c]):
                U=native(bank[c],p);z,pi=orbit[encode(U)]
                if hc[z]==37:rooted[c,v].setdefault((z,pi[v]),p)
    src=(ca,orbit[encode(A)][1][2]);dst=(cb,orbit[encode(B)][1][3])
    rp=class_path(rooted,src,dst);need(rp is not None and len(rp)==4,'shortest rooted neutral route')
    T,steps=literal_route(A,[(a[0],b[0],p) for a,b,p in rp],orbit)
    need(any(p[2]==3 and relabel(B,p)==T for p in permutations(range(6))),'physical root2 arrives at B root3')
    all_graphs=[A]+[tuple(s['target_rows']) for s in steps]
    need(all(any(S[p[-1]]>>p[0]&1 for p in paths(S)) for S in all_graphs),'all neutral-reset sources strong')
    need(any(bank[144][p[-1]]>>p[0]&1 for p in paths(bank[144])) and orbit[encode(converse(bank[144]))][0]==144,
         'isolated H37 source strong and self-converse')
    profiles=[cover_profile(S) for S in all_graphs]
    need(len(set(profiles))==1 and profiles[0]==(37,302,1044,1896,1800,720),'universal homogeneous reset certificate')
    U=birth(birth(T,2),4);need(H(U)==307 and len(paths(U))==307,'literal reset repair307')
    before=[H(birth(birth(A,2),w)) for w in range(7)];need(307 not in before,'original same-root source misses307')
    result['H37_rooted_neutral_reset']={'steps':steps,'root':2,'root_profiles':[root_profile(S,2) for S in all_graphs],
        'full_cover_profile':profiles[0],'final_to_B_maps':[p for p in permutations(range(6)) if relabel(B,p)==T],
        'after_reset_birth_roots':[2,4],'after_reset_rows':U,'after_reset_H':307,'original_two_birth_H':before,
        'isolated_H37_class_rows':bank[144],'isolated_cover_profile':cover_profile(bank[144])}
    # A genuine module-reset positive/hostile bank, complete through order six.
    controls=Counter();witness=None;profile_witness=None
    for c,T in bank.items():
        a=attachment(T)
        for M in modules(T):
            if M==(1<<6)-1:continue
            U=reverse_module(T,M);need(H(U)==H(T),'module reversal preserves H')
            b=attachment(U);controls['proper_module_resets']+=1
            for cut in a:
                if cut&M in (0,M):need(a[cut]==b[cut],'uniform cut survives module reversal');controls['uniform_cuts']+=1
            for v in range(6):
                if not M>>v&1:need(H(birth(T,v))==H(birth(U,v)),'outside-root covariance')
            if set(a.values())!=set(b.values()) and witness is None:
                witness={'source_rows':T,'module_mask':M,'reset_rows':U,'H':H(T),
                    'source_ear_support':sorted(set(a.values())),'reset_ear_support':sorted(set(b.values())),
                    'new_target_cut':next((s for s,h in b.items() if h not in a.values()),None)}
            if orbit[encode(U)][0]!=c and [H(birth(T,v)) for v in range(6)]!=[H(birth(U,v)) for v in range(6)] and profile_witness is None:
                profile_witness={'source_rows':T,'module_mask':M,'reset_rows':U,'H':H(T),
                    'before':[H(birth(T,v)) for v in range(6)],'after':[H(birth(U,v)) for v in range(6)]}
    result['module_resets']={'controls':dict(controls),'ear_support_witness':witness,'named_adaptive_response_witness':profile_witness}
    h25=[]
    for c,T in bank.items():
        if H(T)!=25:continue
        a=attachment(T);hits={str(h):[s for s,z in a.items() if z==h] for h in (75,77,79)}
        need(all(hits.values()),'order-six H25 arbitrary-ear branches')
        for h,cuts in hits.items():need(len(paths(ear(T,cuts[0])))==int(h),'literal H25 branch')
        h25.append({'ordinary_code':c,'rows':T,'hits':hits})
    result['complete_order6_H25_arbitrary_ears']=h25
    # Fixed incoming H93/H95 representatives, no all-order-seven census.
    seed=ear(decode(5,759),2)
    selected=[]
    for m in (50,24):
        T=ear(seed,m);h=H(T);targets={3*h,3*h+2,3*h+4};a=attachment(T)
        resets=[]
        for p in paths(T):
            U=native(T,p)
            if H(U)!=h:continue
            b=attachment(U)
            resets.append({'path':p,'rows':U,'ternary_hits':{str(x):[s for s,z in b.items() if z==x] for x in sorted(targets)}})
        selected.append({'source_rows':T,'H':h,'source_hits':{str(x):[s for s,z in a.items() if z==x] for x in sorted(targets)},
                         'same_H_native_resets':resets,'proper_modules':[M for M in modules(T) if M!=(1<<len(T))-1],
                         'closed_same_H_component':same_H_component(T),'two_step_endpoint_returns':two_step_returns(T)})
    result['selected_incoming']=selected
    return result


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
