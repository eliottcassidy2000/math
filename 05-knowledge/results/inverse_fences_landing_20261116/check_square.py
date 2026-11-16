from math import gcd

def graph(m,s):
    return {x:{t-x for t in (s,m,m+s) if 1<=t-x<m and 2*x!=t}
            for x in range(1,m)}

def components(adj):
    missing=set(adj)
    out=[]
    while missing:
        root=next(iter(missing))
        stack=[root]
        seen={root}
        while stack:
            for v in adj[stack.pop()]-seen:
                seen.add(v)
                stack.append(v)
        missing-=seen
        out.append((len(seen),all(len(adj[x])==2 for x in seen)))
    return sorted(out)

def walk(adj,start,banned=None):
    ans=[start]
    prev=banned
    while nxt:=adj[ans[-1]]-{prev}:
        assert len(nxt)==1
        prev=ans[-1]
        ans.append(nxt.pop())
        if ans[-1]==banned:break
    return ans

def repair(m,s,u):
    d=gcd(s,m)
    assert d%2==1 and gcd(d,u)==1 and ((d-1)//2)*u<s
    adj=graph(m,s)
    ans=walk(adj,s)[::-1]
    for i in range(1,(d-1)//2+1):
        x,y=m-s+i*u,s-i*u
        ans.extend(walk(adj,x,y))
    assert len(ans)==m-1 and set(ans)==set(range(1,m))
    sums=[a+b for a,b in zip(ans,ans[1:])]
    assert set(sums)<={s,m,m+s,m+u}
    assert sums.count(m+u)==(d-1)//2
    return ans

for m in range(2,81):
    for s in range(1,m):
        d=gcd(m,s); L=m//d
        expected=[(L-1,False)]+[(L,False)]*(not d%2)+[(2*L,True)]*((d-1)//2)
        assert components(graph(m,s))==sorted(expected)
        if d%2:
            for u in range(1, min(11,2*(s-1)//max(1,d-1)+1)):
                if gcd(u,d)==1 and ((d-1)//2)*u<s and u!=s:
                    repair(m,s,u)
example=repair(2025,576,91)
print("Component formula checked for all 3160 pairs through m=80; patch checked at admissible small u.")
print("Square-only example: size", len(example), "endpoints", example[0],example[-1])
print("New seams", [(a,b) for a,b in zip(example,example[1:]) if a+b==2116])
