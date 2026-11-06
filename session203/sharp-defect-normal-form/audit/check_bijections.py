"""Independent literal path/cover sets for endpoint equality and the one-edge repair."""
from itertools import permutations
import json

def need(b,s):
    if not b:raise RuntimeError(s)
def good(g,p):return all(g[x]>>y&1 for x,y in zip(p,p[1:]))
def paths(g,V):return {p for p in permutations(V) if good(g,p)}
def covers(g,V,k=2):
    ans=set()
    for p in permutations(V):
        def splits(start,left,bins):
            if left==1:
                q=p[start:]
                if good(g,q):ans.add(bins+(q,))
                return
            for stop in range(start,len(p)+1):
                q=p[start:stop]
                if good(g,q):splits(stop,left-1,bins+(q,))
        splits(0,k,())
    return ans
roots=bad=repairs=lossH=gainH=lossC=gainC=0
for n in range(1,6):
 for code in range(1<<(n*(n-1)//2)):
    g=[0]*n;b=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>b&1:g[i]|=1<<j
            else:g[j]|=1<<i
            b+=1
    V=tuple(range(n));H=paths(g,V);C=covers(g,V)
    for v in V:
        L={(q,p) for q,p in C if q and q[-1]==v}
        image={(p[:p.index(v)+1],p[p.index(v)+1:]) for p in H}
        complement=L-image;need(complement=={(q,p) for q,p in L if p and g[p[0]]>>v&1},'bad cover set')
        A=tuple(x for x in V if x!=v);incoming={x for x in A if g[x]>>v&1}
        expect={((),p) for p in paths(g,A) if p and p[0] in incoming}|{(q,p) for q,p in covers(g,A) if q and p and q[-1] in incoming and p[0] in incoming}
        need({(q[:-1],p) for q,p in complement}==expect,'bad deletion exact bijection')
        roots+=1;bad+=len(complement)
    if n>4:continue
    for x in V:
        if any(p[0]==x or p[-1]==x for p in H):continue
        # a=n, v=n+1, literal unflipped and flipped graphs.
        a=n;v=n+1
        base=[row|(1<<v) for row in g]+[(1<<n)-1,1<<a]
        flip=base.copy();flip[a]^=1<<x;flip[x]|=1<<a
        VV=V+(a,v);old=paths(base,VV);new=paths(flip,VV)
        lost=old-new;gained=new-old
        R={(q,p) for q,p in C if q and q[0]==x}
        L={(q,p) for q,p in C if q and q[-1]==x}
        need(lost=={p+(v,a)+q for q,p in R},'lost Hamilton paths')
        need(gained=={q+(a,)+p+(v,) for q,p in L},'gained Hamilton paths')
        oldC=covers(base,V+(a,));newC=covers(flip,V+(a,));lostC=oldC-newC;newlyC=newC-oldC
        image=set()
        for q,p in R:image|={((a,)+q,p),(p,(a,)+q)}
        need(lostC==image,'lost two-cover bin labels')
        E3={(q,p,r) for q,p,r in covers(g,V,3) if q and q[-1]==x}
        image=set()
        for q,p,r in E3:image|={(q+(a,)+p,r),(r,q+(a,)+p)}
        need(newlyC==image and len(image)==2*len(E3),'gained two-cover deletion inverse')
        repairs+=1;lossH+=len(lost);gainH+=len(gained);lossC+=len(lostC);gainC+=len(newlyC)
print(json.dumps(dict(status='PASS',marked_inputs=roots,bad_cover_bijections=bad,physical_flip_cases=repairs,lost_Hamilton_paths=lossH,gained_Hamilton_paths=gainH,lost_two_covers=lossC,gained_two_covers=gainC),sort_keys=True,indent=2))
