#!/usr/bin/env python3
from verify_extension import *
from functools import lru_cache

@lru_cache(None)
def endpoints(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,c in enumerate(d[mask]):
            if not c:continue
            todo=a[i]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;d[mask|bit][bit.bit_length()-1]+=c
    return d
@lru_cache(None)
def all_bins(a):
    n=len(a);h=[sum(v) for v in endpoints(a)];h[0]=1;c=[[0]*(n+1) for _ in h];c[0][0]=1
    for mask in range(1,len(h)):
        bit=mask&-mask;s=mask
        while s:
            if s&bit:
                for k in range(1,n+1):c[mask][k]+=h[s]*c[mask^s][k-1]
            s=(s-1)&mask
    falls=[[1]]
    for k in range(n):falls.append(mul(falls[-1],[-k,1]))
    return [add(*[scale(falls[k],v) for k,v in enumerate(row)]) for row in c]

def marked(a,anchors,roles):
    n=len(a);reverse=tuple(sum(1<<j for j in range(n) if a[j]>>i&1) for i in range(n))
    ep=endpoints(a);sp=endpoints(reverse);bins=all_bins(a);amask=sum(1<<v for v in anchors)
    @lru_cache(None)
    def rec(j,mask):
        if j==len(anchors):return bins[mask]
        root=anchors[j];tab=sp if roles[j]=='S' else ep;out=[0];s=mask
        while s:
            if s>>root&1 and s&amask==1<<root and tab[s][root]:out=add(out,scale(rec(j+1,mask^s),tab[s][root]))
            s=(s-1)&mask
        return out
    return rec(0,(1<<n)-1)

def substitute(q,r,a):
    outside=[i for i in range(len(q)) if i!=r];N=len(outside);o=[0]*(N+len(a))
    for i,u in enumerate(outside):
        for j,v in enumerate(outside):
            if q[u]>>v&1:o[i]|=1<<j
        for j in range(len(a)):
            if q[u]>>r&1:o[i]|=1<<(N+j)
            else:o[N+j]|=1<<i
    for j,v in enumerate(a):o[N+j]|=v<<N
    validate(o);return tuple(o)

def test_module(m,roles):
    e=roles.count('E');s=roles.count('S');ei=0;si=e+m;anchors=[]
    for role in roles:
        if role=='E':anchors.append(ei);ei+=1
        else:anchors.append(si);si+=1
    return transitive(m+e+s),tuple(anchors)

def positive_extended(p):
    if p[0]<0:return False
    odd=[v if j%2 else 0 for j,v in enumerate(p)];even=[v if j and not j%2 else 0 for j,v in enumerate(p)]
    for f,par in ((trim(odd),1),(trim(even),0)):
        if f!=[0] and any(v<0 for v in coordinates(f,par)):return False
    return True

def markedD(f,s,e):return add(D(f),scale(mul(mon(1),f),s+e),scale(shift(f,1),s+e),scale(f,s*e))

def run():
    inputs=[]
    for n in range(1,4):
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask)
            for k in range(1,min(2,n)+1):
                for anchors in permutations(range(n),k):
                    for roles in product('SE',repeat=k):inputs.append((a,anchors,roles))
    contexts=[(tournament(n,m),r) for n in range(1,4) for m in range(1<<(n*(n-1)//2)) for r in range(n)]
    checks=0;virtual=0
    for q,r in contexts:
        cols={}
        for roles in [('S',),('E',),('S','S'),('S','E'),('E','S'),('E','E')]:
            images=[]
            for m in range(5):
                a,anchors=test_module(m,roles);require(marked(a,anchors,roles)==mon(m),'test-module identification')
                out=substitute(q,r,a);images.append(marked(out,tuple(len(q)-1+v for v in anchors),roles))
            cols[roles]=images
            for par,k in product((0,1),range(2)):
                b=basis(k,par);result=add(*[scale(images[m],v) for m,v in enumerate(b) if v])
                require(positive_extended(result),('virtual positivity',q,r,roles,par,k));virtual+=1
        for a,anchors,roles in inputs:
            f=marked(a,anchors,roles);out=substitute(q,r,a);actual=marked(out,tuple(len(q)-1+v for v in anchors),roles)
            expected=add(*[scale(cols[roles][m],v) for m,v in enumerate(f) if v])
            require(actual==expected,('full marked linearity',q,r,a,anchors,roles));require(positive_extended(f),'input shifted positivity');checks+=1
    concrete=0
    for n in range(1,4):
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask);w=wrapper(a)
            for k in range(1,n+1):
                for anchors in permutations(range(n),k):
                    for roles in product('SE',repeat=k):
                        f=marked(a,anchors,roles);s=roles.count('S');e=k-s
                        require(marked(w,anchors,roles)==markedD(f,s,e),'closed W marked transfer');concrete+=1
            f=bin_polynomial(a);na=n;nb=n+1
            # wrapper() uses T -> vertex n -> vertex n+1 -> T, so a=n+1, b=n.
            expected_a=add(mul(mon(1),f),mul([1,1],add(shift(f,1),scale(f,-1))))
            expected_b=add(mul(mon(1),f),shift(f,1))
            require(marked(w,(nb,),('S',))==expected_a,'new apex a start')
            require(marked(w,(na,),('S',))==expected_b,'new apex b start')
            require(marked(w,(nb,),('E',))==expected_b,'new apex a end')
            require(marked(w,(na,),('E',))==expected_a,'new apex b end')
    for par,k in product((0,1),range(31)):
        b=basis(k,par);require(positive_extended(add(shift(b,1),scale(b,-1))),'positive forward difference')
    tri=(2,4,1);f=marked(tri,(0,),('S',));require(f==[1,1,1],'C3 reserved shift hostile')
    print(json.dumps({'input_marked_states':len(inputs),'rooted_contexts':len(contexts),'physical_marked_context_checks':checks,'formal_basis_context_checks':virtual,'explicit_W_role_checks':concrete,'forward_difference_basis_checks':62,'C3_shifted_root_start':f},sort_keys=True,indent=2))
if __name__=='__main__':run()
