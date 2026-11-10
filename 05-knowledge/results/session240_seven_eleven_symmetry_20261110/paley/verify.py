"""Dependency-free exhaustive incidence automorphisms and skew block labellings."""
from collections import Counter
import json

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def blocks(p,S): return [frozenset((i+s)%p for s in S) for i in range(p)]
def automorphisms(B):
    n=len(B); ans=[]
    # Cells contain source/target blocks with identical membership signatures.
    def dfs(perm,left,cells):
        d=len(perm)
        if d==n:
            beta=[None]*n
            for a,b in cells:
                check(len(a)==len(b)==1, 'verify.py validation originally at line 14')
                beta[a[0]]=b[0]
            ans.append((tuple(perm),tuple(beta)));return
        source=[([a for a in x if d not in B[a]],[a for a in x if d in B[a]],y) for x,y in cells]
        for v in left:
            nxt=[]
            for a0,a1,y in source:
                b0=[b for b in y if v not in B[b]]; b1=[b for b in y if v in B[b]]
                if len(a0)!=len(b0):break
                if a0:nxt.append((a0,b0))
                if a1:nxt.append((a1,b1))
            else: dfs(perm+[v],left-{v},nxt)
    dfs([],set(range(n)),[(list(range(n)),list(range(n)))])
    return ans

def skew_labels(B):
    n=len(B); ans=[]
    def dfs(f,left):
        i=len(f)
        if i==n: ans.append(tuple(f));return
        for b in left:
            if i not in B[b] and all((j in B[b]) != (i in B[f[j]]) for j in range(i)):
                dfs(f+[b],left-{b})
    dfs([],set(range(n)));return ans

def orbit_lengths(perms,n):
    unseen=set(range(n)); sizes=[]
    while unseen:
        x=min(unseen); orb={g[x] for g in perms};sizes.append(len(orb));unseen-=orb
    return sorted(sizes)

def path_endpoints(B):
    n=len(B); dp={(1,0):1}
    for size in range(1,n):
        nxt={}
        for (mask,last),count in dp.items():
            for v in B[last]:
                if not (mask>>v)&1:
                    key=(mask|(1<<v),v);nxt[key]=nxt.get(key,0)+count
        dp=nxt
    return [dp.get(((1<<n)-1,v),0) for v in range(n)]

def run(p,S,tag):
    B=blocks(p,S); G=automorphisms(B); K=[a for a,b in G if a==b]
    skew=skew_labels(B)
    # Transport canonical matching i -> B_i by (alpha,beta).
    transported=set()
    for a,b in G:
        f=[None]*p
        for i in range(p):f[a[i]]=b[i]
        transported.add(tuple(f))
    check(transported <= set(skew), 'verify.py validation originally at line 65')
    check(len(transported)*len(K)==len(G), 'verify.py validation originally at line 66')
    # Independent direct check on completed permutations.
    for a,b in G:
        check(all(frozenset(a[x] for x in B[i])==B[b[i]] for i in range(p)), 'verify.py validation originally at line 69')
    check(all(all(((j-i)%p in S)==((a[j]-a[i])%p in S) for i in range(p) for j in range(p)) for a in K), 'verify.py validation originally at line 70')
    result=dict(tag=tag,p=p,S=sorted(S),design_aut=len(G),tournament_aut=len(K),canonical_matching_orbit=len(transported),all_skew_matchings=len(skew),all_skew_in_one_orbit=set(skew)==transported,design_point_stabilizer_orbits=orbit_lengths([a for a,b in G if a[0]==0],p),tournament_point_stabilizer_orbits=orbit_lengths([a for a in K if a[0]==0],p),design_pair_stabilizer=sum(a[0]==0 and a[1]==1 for a,b in G),tournament_pair_stabilizer=sum(a[0]==0 and a[1]==1 for a in K),skew_matchings=[list(f) for f in skew] if p<=11 else 'omitted')
    if p<=11:
        matchings=sorted(skew); idx={f:i for i,f in enumerate(matchings)}
        frame_perms=[]
        for a,b in G:
            action=[]
            for f in matchings:
                new=[None]*p
                for i in range(p):new[a[i]]=b[f[i]]
                action.append(idx[tuple(new)])
            frame_perms.append(action)
        zero=idx[tuple(range(p))]
        result['matching_stabilizer_orbits']=orbit_lengths([a for a in frame_perms if a[zero]==zero],len(skew))
        result['matching_action_kernel']=sum(a==list(range(len(skew))) for a in frame_perms)
        ends=path_endpoints(B)
        result['paths_from_0_by_endpoint']=ends
        result['H']=p*sum(ends)
        result['H_tournament_orbits']=result['H']//len(K)
        result['H_mod_design_aut']=result['H']%len(G)
    print(json.dumps(result,sort_keys=True),flush=True)
    return result
if __name__=='__main__':
    for p in [3,7,11,19]:run(p,{x*x%p for x in range(1,p)},'Paley')
    for p in [7,11]:run(p,set(range(1,(p+1)//2)),'rotational')
    print('ALL CHECKS PASSED')
