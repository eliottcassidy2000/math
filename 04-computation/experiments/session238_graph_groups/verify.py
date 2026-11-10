#!/usr/bin/env python3
"""Exact finite controls, no assertions (also runs under -O)."""
import itertools, json, sys

def check(ok, msg):
    if not ok: raise RuntimeError(msg)

def test(mods, differences, inject=False):
    group=list(itertools.product(*(range(n) for n in mods)))
    zero=group[0]; m=len(group)
    def add(a,b): return tuple((x+y)%n for x,y,n in zip(a,b,mods))
    def neg(a): return tuple((-x)%n for x,n in zip(a,mods))
    edges={frozenset((v,add(v,d))) for v in group for d in differences}
    edges=sorted(edges,key=lambda e: sorted(e))
    vertices=[('s',tuple(sorted(e))) for e in edges]
    degree=len(differences); priv=5+degree%2
    vertices += [('p',v,i) for v in group for i in range(priv)]
    index={v:i for i,v in enumerate(vertices)}; N=len(vertices)
    def move(v,a):
        if v[0]=='p': return ('p',add(v[1],a),v[2])
        return ('s',tuple(sorted(add(x,a) for x in v[1])))
    permutations=[[index[move(v,a)] for v in vertices] for a in group]
    check(all(all(p[i]!=i for i in range(N)) for p in permutations[1:]),'free action')
    cycles=[]; triangles=[]; arcs=set()
    def arc(x,y):
        check((y,x) not in arcs, 'prescription conflict'); arcs.add((x,y))
    baseedges=[e for e in edges if zero in e]
    base=[index[('s',tuple(sorted(e)))] for e in baseedges]
    base += [index[('p',zero,i)] for i in range(priv)]
    for p in permutations:
        c=[p[x] for x in base]; cycles.append(c)
        for x,y in zip(c,c[1:]+c[:1]): arc(x,y)
        tri=[p[index[('p',zero,i)]] for i in range(3)]
        arc(tri[2],tri[0]); triangles.append(tri)
    forced=len(arcs)
    check(forced==2*len(edges)+6*m+m*(degree%2),'forced count')
    expected=(N*(N-1)//2-forced)//m
    freeorbits=0
    for x in range(N):
        for y in range(x+1,N):
            if (x,y) in arcs or (y,x) in arcs: continue
            orbit={(p[x],p[y]) for p in permutations}
            check(len(orbit)==m and not any((b,a) in orbit for a,b in orbit),'pair orbit')
            for a,b in orbit: arc(a,b)
            freeorbits+=1
    check(freeorbits==expected,'completion exponent')
    if inject:
        a,b=next(iter(arcs)); arcs.remove((a,b)); arcs.add((b,a))
    for p in permutations:
        check(all((p[a],p[b]) in arcs for a,b in arcs),'group invariance')
    for c in cycles+triangles:
        check(all((a,b) in arcs for a,b in zip(c,c[1:]+c[:1])),'directed cycle')
    for i,j in itertools.combinations(range(m),2):
        check(bool(set(cycles[i])&set(cycles[j]))==(frozenset((group[i],group[j])) in edges),'selected conflict')
    for i,j in itertools.product(range(m),repeat=2):
        check(bool(set(cycles[i])&set(triangles[j]))==(i==j),'corona attachment')
        if i!=j: check(not set(triangles[i])&set(triangles[j]),'private triangle')
    for i,c in enumerate(cycles):
        a=index[('p',group[i],3)]
        check(any((a,b) in arcs and (b,d) in arcs and (d,a) in arcs for b,d in itertools.permutations(set(c)-{a},2)),'extra triangle')
    result={'group':mods,'degree':degree,'m':m,'edges':len(edges),'N':N,'completion_exponent':expected}
    if m<=3 and not edges:
        # Independent subset DP, not the cycle formula.
        incoming=[sum(1<<a for a,b in arcs if b==v) for v in range(N)]
        dp=[{} for _ in range(1<<N)]
        for v in range(N): dp[1<<v][v]=1
        for mask in range(1,1<<N):
            if mask&(mask-1)==0: continue
            for v in range(N):
                if mask>>v&1:
                    prev=mask^(1<<v)
                    dp[mask][v]=sum(cnt for u,cnt in dp[prev].items() if incoming[v]>>u&1)
        H=sum(dp[-1].values()); lower=5**m+2*m
        check(H% m==0 and H>=lower and H//m not in (7,21),'H bound/group quotient')
        result.update(H=H,lower_bound=lower)
    return result

def main():
    rows=[]
    for n in (1,3,5,7):
        for mask in range(1<<((n-1)//2)):
            ds={(i,) for j in range((n-1)//2) if mask>>j&1 for i in (j+1,n-j-1)}
            rows.append(test((n,),ds, '--inject-failure' in sys.argv and n==3 and mask==0))
    nonzero=set(itertools.product(range(3),repeat=2))-{(0,0)}
    rows.append(test((3,3),set())); rows.append(test((3,3),nonzero))
    # Hostile even action: its involution swaps a pair, making invariant orientation impossible.
    check({(0,1),(1,0)}=={((0+a)%2,(1+a)%2) for a in range(2)},'even hostile')
    print(json.dumps({'status':'PASS','cases':len(rows),'even_hostile':'reversed pair detected','rows':rows},indent=2))
if __name__=='__main__': main()
