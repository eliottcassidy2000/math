from itertools import combinations, permutations
from collections import Counter

def verify_words(words):
    vertices=set(words[0]); pairs=set(combinations(sorted(vertices),2))
    initial={}; used=Counter()
    for i,p in enumerate(words):
        assert set(p)==vertices and len(p)==len(vertices)
        for a,b in zip(p,p[1:]):
            e=tuple(sorted((a,b)))
            direction=(a,b) if (i-used[e])%2==0 else (b,a)
            assert e not in initial or initial[e]==direction
            initial[e]=direction;used[e]+=1
    assert set(initial)==pairs and all(v%2==1 for v in used.values())
    t=set(initial.values());start=set(t);states=[]
    for p in words:
        assert t not in states
        states.append(set(t))
        keep=set(zip(p,p[1:]));assert keep<=t
        t=keep|{(b,a) for a,b in t-keep}
    assert t==start
    return used,states

root=[(3,4,1,2,5),(4,5,1,2,3),(5,3,1,2,4)]
used,states=verify_words(root)
assert [e for e,c in used.items() if c==3]==[(1,2)]
print('ROOT_N5 valid directed triangle, triple edge (1,2), states', [sorted(s) for s in states])
for n in range(6,103,4):
    k=n//2
    words=[[i]+[z%n for d in range(1,k) for z in (i-d,i+d)]+[(i-k)%n] for i in range(k)]
    used,states=verify_words(words)
    assert set(used.values())=={1}
print('WALECKI independent directed-set validation n=6,10,...,102 PASSED')
# Complete labelled n5 graph; attach path edge sets to every edge independently.
n=5;pairs=list(combinations(range(n),2));words=list(permutations(range(n)))
adj={}; paths={}
for bits in range(1<<len(pairs)):
    t={e if bits>>i&1 else e[::-1] for i,e in enumerate(pairs)}
    adj[bits]=set()
    for p in words:
        path=set(zip(p,p[1:]))
        if not path<=t:continue
        new=path|{(b,a) for a,b in t-path}
        other=sum(1<<i for i,e in enumerate(pairs) if e in new)
        assert other!=bits and other not in adj[bits]
        adj[bits].add(other);paths[bits,other]=p
triangles=0
for a in adj:
    for b in adj[a]:
        if b<=a:continue
        for c in adj[a]&adj[b]:
            if c<=b:continue
            ws=[paths[a,b],paths[b,c],paths[c,a]]
            used,states=verify_words(ws)
            assert Counter(used.values())=={1:9,3:1}
            endpoints=Counter(v for p in ws for v in (p[0],p[-1]))
            assert sorted(endpoints.values())==[2,2,2]
            triple=next(e for e,k in used.items() if k==3)
            assert not any(v in endpoints for v in triple)
            triangles+=1
print('N5 exhaustive labelled graph:',len(adj),'vertices',sum(map(len,adj.values()))//2,'edges',triangles,'triangles; equality anatomy passed in every triangle')
print('ALL INDEPENDENT CHECKS PASSED')
ps9=[(5,0,8,6,1,2,3,7,4),(4,2,3,6,0,1,5,7,8),(5,4,8,3,0,1,7,2,6),(6,5,2,8,1,3,4,0,7),(7,6,4,1,0,2,3,5,8)]
used,states=verify_words(ps9)
assert {e for e,c in used.items() if c==3}=={(0,1),(2,3)}
pairs=list(combinations(range(9),2))
bits=[sum(1<<i for i,e in enumerate(pairs) if e in s) for s in states]
assert bits==[5927068074,45360910277,59874146331,9281284065,58367449686]
print('N9 independent directed-set certificate passed; odd girth exactly 5; states',bits)
