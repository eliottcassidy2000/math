#!/usr/bin/env python3
"""Independent edge-set audit, no imports from the theorem's implementation.

Small paths use complete permutations; large path counts use bitmask DP.
Automorphisms use simultaneous candidate propagation and MRV branching.
Native controls enumerate the complete common-pin permutation universe.
All checks survive python -O.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
import json

TESTS = 0

def require(condition, label):
    global TESTS
    TESTS += 1
    if not condition:
        raise RuntimeError(label)

def cycle(n, differences):
    return n, frozenset((u,v) for u in range(n) for v in range(n)
                        if (v-u) % n in differences)

def from_code(n, c):
    return n, frozenset((u,v) if (c>>k)&1 else (v,u)
                        for k,(u,v) in enumerate(combinations(range(n),2)))

def opposite(t):
    return t[0], frozenset((v,u) for u,v in t[1])

def stack(*tt):
    edges = set()
    start = 0
    for n,e in tt:
        edges.update((start+u,start+v) for u,v in e)
        edges.update((u,v) for u in range(start) for v in range(start,start+n))
        start += n
    return start, frozenset(edges)

def induced(t, labels):
    return len(labels), frozenset((i,j) for i,u in enumerate(labels)
                                  for j,v in enumerate(labels) if (u,v) in t[1])

def path_ok(t, p):
    return len(p)==t[0] and len(set(p))==t[0] and all(z in t[1] for z in zip(p,p[1:]))

def brute_paths(t):
    return [p for p in permutations(range(t[0])) if all(z in t[1] for z in zip(p,p[1:]))]

def walk_paths(t):
    n,e=t
    todo=[((),frozenset(range(n)))]; out=[]
    while todo:
        p,free=todo.pop()
        if not free:
            out.append(p)
        else:
            todo.extend((p+(u,),free-{u}) for u in free if not p or (p[-1],u) in e)
    return out

def dp_count(t):
    n,e=t
    incoming=[sum(1<<u for u in range(n) if (u,v) in e) for v in range(n)]
    @lru_cache(None)
    def endings(mask,last):
        if mask==1<<last:
            return 1
        prev=mask^(1<<last)
        choices=prev&incoming[last]
        return sum(endings(prev,u) for u in range(n) if choices>>u&1)
    return sum(endings((1<<n)-1,u) for u in range(n))

def is_maps(t, u):
    """Complete edge-isomorphism search by propagated candidate domains."""
    n,e=t; m,f=u
    if n!=m:
        return []
    de=[sum((a,b) in e for b in range(n)) for a in range(n)]
    df=[sum((a,b) in f for b in range(n)) for a in range(n)]
    initial={a:{b for b in range(n) if de[a]==df[b]} for a in range(n)}
    ans=[]
    def rec(mapping, candidates):
        if not candidates:
            ans.append(tuple(mapping[a] for a in range(n)))
            return
        a=min(candidates, key=lambda x:(len(candidates[x]),x))
        for b in sorted(candidates[a]):
            nxt={x:{y for y in ys if y!=b and (((a,x) in e)==((b,y) in f))}
                 for x,ys in candidates.items() if x!=a}
            if all(nxt.values()):
                rec(dict(mapping,**{})|{a:b},nxt)
    rec({},initial)
    return ans

def image(t, p):
    return t[0], frozenset((p[u],p[v]) for u,v in t[1])

def K(t,p):
    keep=frozenset(zip(p,p[1:]))
    return t[0], frozenset((u,v) if (u,v) in keep else (v,u) for u,v in t[1])

@lru_cache(None)
def cells(n):
    # Deliberately use gap-major, rather than the worker's row-major word order.
    return tuple((i,i+d) for d in range(2,n) for i in range(n-d))

def word(t,p):
    return tuple(int((p[i],p[j]) in t[1]) for i,j in cells(t[0]))

def reflect(w,n):
    values=dict(zip(cells(n),w))
    return tuple(values[n-1-j,n-1-i] for i,j in cells(n))

def flip(w):
    return tuple(1-z for z in w)

def resource(w,n):
    return frozenset((w,flip(w),reflect(w,n),flip(reflect(w,n))))

def reach_strong(t):
    n,e=t
    reachable=[{u}|{v for v in range(n) if (u,v) in e} for u in range(n)]
    for k in range(n):
        for u in range(n):
            if k in reachable[u]:
                reachable[u]|=reachable[k]
    return all(len(z)==n for z in reachable)

C3=cycle(3,{1})
P0=(3,0,6,7,1,5,4,2)
P1=(3,1,6,7,0,5,4,2)
G=(4,0,2,3,1,5,6,7)
SEED_ROWS=("00001110","10000110","11000010","11100000",
           "01110000","00111000","00011101","11111100")
SEED=(8,frozenset((u,v) for u,row in enumerate(SEED_ROWS) for v,z in enumerate(row) if z=='1'))

def audit_small():
    data=[]
    for q in range(1,5):
        for c in range(1<<(q*(q-1)//2)):
            R=from_code(q,c)
            pp=brute_paths(R); aa=is_maps(R,R)
            reps={word(R,p):p for p in pp}
            a=len(reps)
            require(len(pp)==dp_count(R)==len(aa)*a,'small free path action')
            S=stack(R,C3,opposite(R)); n=S[0]
            allp=walk_paths(S); autS=is_maps(S,S)
            require(len(allp)==3*len(pp)**2,'whole-S literal path factorization')
            require(len(autS)==3*len(aa)**2,'whole-S full automorphism factorization')
            W={word(S,p) for p in allp}; pairwords={}
            for x,p in enumerate(reps.values()):
                for y,t in enumerate(reps.values()):
                    obs=p+tuple(q+i for i in range(3))+tuple(q+3+u for u in t[::-1])
                    pairwords[x,y]=word(S,obs)
            require(W==set(pairwords.values()) and len(W)==a*a,'square bijection, not just product')
            for (x,y),w in pairwords.items():
                require(reflect(w,n)==pairwords[y,x],'reflection exactly swaps the coordinates')
                require((reflect(w,n)==w)==(x==y),'blue iff diagonal')
                require(resource(w,n)&W=={pairwords[x,y],pairwords[y,x]},'exact resource intersection')
            for p in allp:
                target=K(S,p)
                require(target[1]&S[1]==frozenset(zip(p,p[1:])),'complete agreement graph')
                require(K(target,p)==S,'literal K involution')
                require(K(opposite(S),p[::-1])==opposite(target),'literal converse square')
                require(reach_strong(target),'every target is strong')
                require((p[-1],p[0]) in target[1],'retained path closes globally')
            require(not reach_strong(S),'nonempty source nonstrong')
            require(not(W&{flip(w) for w in W}),'no ordinary or merged loop')
            require(sum(w==reflect(w,n) for w in W)==a,'blue count')
            require(len({resource(w,n) for w in W})==a*(a+1)//2,'all resource count')
            data.append((q,c,a,len(W)))
    return {'all_labelled_R_orders':[1,2,3,4],'count':len(data),
            'total_source_words':sum(z[3] for z in data),'max_a':max(z[2] for z in data)}

def audit_famous():
    results=[]
    for n,d,h,g,a in [(3,{1},3,3,1),(5,{1,2},15,5,3),
                       (7,{1,2,4},189,21,9),(11,{1,3,4,5,9},95095,55,1729)]:
        T=cycle(n,d); count=dp_count(T); maps=is_maps(T,T)
        require((count,len(maps),count//len(maps))==(h,g,a),'independent famous H/Aut control')
        if n<=7:
            pp=brute_paths(T)
            require(len(pp)==h and len({word(T,p) for p in pp})==a,'full famous words by permutations')
        else:
            pp=walk_paths(T)
            words=Counter(word(T,p) for p in pp)
            require(len(pp)==h and len(words)==a and set(words.values())=={g},
                    'Paley11 complete normalized words and exact lift multiplicity')
        results.append({'R_order':n,'H':count,'Aut':len(maps),'a':a,
                        'S_order':2*n+3,'S_words':a*a,'blue':a,'black_weighted':a*(a-1),
                        'resources':a*(a+1)//2})
    require(not reach_strong(K(C3,(0,1,2))),'empty-R strong-target hostile')
    return results

def pinned_prefixes(r):
    # Enumerate all owner assignments under the six fixed clocks in each block.
    result=[]
    for bits in product((0,1),repeat=r):
        result.append(tuple(8*j+v for j,b in enumerate(bits) for v in (P0,P1)[b]))
    return result

def clock_work(p,q):
    first={u:i for i,u in enumerate(p)}
    second={u:i for i,u in enumerate(q)}
    return sum(first[u]!=second[u] for u in first)

def audit_native():
    encoded=0
    for u,v in combinations(range(8),2):
        encoded=2*encoded+int((u,v) in SEED[1])
    require(encoded==29558785,'literal seed is authenticated canonical bit code')
    pp=brute_paths(SEED)
    require(len(pp)==dp_count(SEED)==409,'seed H independently brute forced')
    aa=[p for p in permutations(range(8)) if image(SEED,p)==SEED]
    require(aa==[tuple(range(8))],'seed rigidity by all 8! permutations')
    require(reach_strong(SEED),'seed strongly connected')
    fixed=(0,2,3,5,6,7)
    require({p for p in pp if all(p[i]==P0[i] for i in fixed)}=={P0,P1},'complete pinned seed paths')
    require(image(K(SEED,P0),G)==K(SEED,P1),'literal seed partner map')
    require(G[3]==3 and G[2]==2,'both actual junction endpoints fixed')
    rows=[]
    for r in range(1,5):
        R=stack(*([SEED]*r)); q=R[0]; S=stack(R,C3,opposite(R)); n=S[0]
        L=induced(S,tuple(range(q+3))); D=induced(S,tuple(range(q,n)))
        right_map=tuple(range(3,q+3))+tuple((-i)%3 for i in range(3))
        require(image(opposite(L),right_map)==D,'right local source actual anti-isomorphism')
        beta=tuple(q+3+u for u in range(q))+tuple(q+(-i)%3 for i in range(3))+tuple(range(q))
        require(image(S,beta)==opposite(S),'whole source actual anti-automorphism')
        prefixes=pinned_prefixes(r); m=len(prefixes); allrecords={}; GW={}; LW={}; RW={}
        base=prefixes[0]+tuple(range(q,q+3))+tuple(q+3+v for v in prefixes[0][::-1])
        base_target=K(S,base)
        # Every unpinned assignment is considered, including all six C3 orders.
        accepted=0
        for x,p in enumerate(prefixes):
            for y,t in enumerate(prefixes):
                for middle in permutations(range(q,q+3)):
                    obs=p+middle+tuple(q+3+v for v in t[::-1])
                    if not path_ok(S,obs):
                        require((middle[0],middle[1]) not in S[1] or (middle[1],middle[2]) not in S[1],
                                'only cycle direction rejects a common-pin assignment')
                        continue
                    accepted+=1
                    k=middle[0]-q; target=K(S,obs)
                    lp=obs[:q+3]; rp=tuple(u-q for u in obs[q:])
                    lk=K(L,lp); dk=K(D,rp)
                    require(induced(target,tuple(range(q+3)))==lk,'left complete adjacency flags')
                    require(induced(target,tuple(range(q,n)))==dk,'right complete adjacency flags')
                    require(induced(lk,tuple(range(q,q+3)))==induced(dk,(0,1,2)),'shared partner literally equal')
                    native_path=t+tuple(q+((1-k+j)%3) for j in range(3))
                    require(tuple(right_map[v] for v in native_path[::-1])==rp,'right phase exactly 1-k')
                    lw=word(L,lp); rw=word(D,rp); gw=word(S,obs)
                    require(rw==reflect(word(L,native_path),q+3),'right native word reflected without F')
                    require((resource(lw,q+3)==resource(rw,q+3))==(x==y),'capacity-one iff off-diagonal')
                    require((gw==reflect(gw,n))==(x==y),'global blue iff capacity-invalid diagonal')
                    require(reach_strong(target),'native target strong')
                    # Read the choice from actual owner order, independently of key-bit order.
                    maps=[]
                    for side,path in enumerate((p,t)):
                        side_map=[]
                        for j in range(r):
                            block=tuple(v-8*j for v in path[8*j:8*j+8])
                            gm=G if block==P1 else tuple(range(8))
                            side_map.extend(8*j+gm[v] for v in range(8))
                        maps.append(tuple(side_map))
                    transport=maps[0]+tuple(q+(j+k)%3 for j in range(3))+tuple(q+3+v for v in maps[1])
                    require(image(base_target,transport)==target,'actual common-target blockwise transport')
                    allrecords[x,y,k]=(obs,lp,rp,target,lk,dk)
                    GW[x,y,k]=gw; LW[x,k]=lw; RW[y,k]=rw
        require(accepted==3*m*m,'complete common-pin global family')
        require(len(set(LW.values()))==m,'selected local source-word injection')
        require(len({resource(w,q+3) for w in LW.values()})==m,'selected local resource injection')
        require(len(set(GW.values()))==m*m,'selected global word square')
        require(len({resource(w,n) for w in GW.values()})==m*(m+1)//2,'selected global resource triangle')
        diagonal_anti=tuple(q+3+u for u in range(q))+tuple(q+(2-i)%3 for i in range(3))+tuple(range(q))
        require(image(base_target,diagonal_anti)==opposite(base_target),
                'common target self-converse through explicit phase-corrected diagonal map')
        for (x,y,k),(obs,lp,rp,target,lk,dk) in allrecords.items():
            swapped=allrecords[y,x,(1-k)%3][0]
            require(tuple(beta[v] for v in obs[::-1])==swapped,'literal reflection swaps marked key order')
            require(reflect(GW[x,y,k],n)==GW[y,x,(1-k)%3],'global quotient records precisely the swap')
            for delta in (1,2):
                nxt=(x,y,(k+delta)%3); o2,l2,r2,t2,lk2,dk2=allrecords[nxt]
                rot=tuple(range(q))+tuple(q+(j+delta)%3 for j in range(3))+tuple(range(q+3,n))
                require(tuple(rot[v] for v in obs)==o2,'literal phase observer motion')
                require(image(S,rot)==S and image(target,rot)==t2,'global source and target phase transport')
                lrot=rot[:q+3]; rrot=tuple(rot[q+j]-q for j in range(q+3))
                require(image(L,lrot)==L and image(lk,lrot)==lk2,'left incidence constant by marked map')
                require(image(D,rrot)==D and image(dk,rrot)==dk2,'right incidence constant by marked map')
                require(clock_work(obs,o2)==3,'physical owner-clock work three')
                require(clock_work(lp,l2)+clock_work(rp,r2)==6,'local occurrence owner-clock work six')
                require(resource(GW[x,y,k],n)==resource(GW[nxt],n),'global resource continuously retained')
                require(resource(LW[x,k],q+3)==resource(LW[x,nxt[2]],q+3)
                        and resource(RW[y,k],q+3)==resource(RW[y,nxt[2]],q+3),'two local resources continuously retained')
        rows.append({'r':r,'common_pin_records':accepted,'keys':m,'capacity_valid_records':3*m*(m-1),
                     'source_words':m*m,'blue_resources':m,'black_resources':m*(m-1)//2,
                     'target_classes':1,'local_work':6,'global_work':3})
    return {'seed_H':409,'seed_Aut':1,'rows':rows}

def audit_hostile_injection():
    # Search rather than assume the worker's four-owner code or path pair.
    for c in range(64):
        U=from_code(4,c)
        if not reach_strong(U) or not is_maps(U,opposite(U)):
            continue
        R=stack(C3,U); L=stack(R,C3)
        pp=[(0,1,2)+tuple(3+v for v in p) for p in brute_paths(U)]
        seen={}
        for p in pp:
            w=word(L,p+(7,8,9)); rr=resource(w,10)
            if rr in seen and seen[rr][1]!=word(R,p):
                p0,w0=seen[rr]
                require(word(R,p0)!=word(R,p),'generic prefix path orbits are distinct')
                require(resource(word(L,p0+(7,8,9)),10)==rr,'distinct prefix orbits collide as resources')
                require(reflect(word(L,p0+(7,8,9)),10)==w,'collision is literal reflection')
                return {'R_description':'C3 -> U_strong4','U_code_lsb':c,'R_paths':[list(p0),list(p)],
                        'L_order':10,'prefix_words_distinct':True,'local_resources_equal':True}
            seen[rr]=(p,word(R,p))
    raise RuntimeError('no generic-injection hostile found')

def main():
    small=audit_small()
    famous=audit_famous()
    native=audit_native()
    hostile=audit_hostile_injection()
    print(json.dumps({'status':'PASS','checks':TESTS,'small':small,'famous':famous,
                      'native':native,'generic_injection_hostile':hostile},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
