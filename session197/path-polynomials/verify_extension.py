#!/usr/bin/env python3
"""Independent cyclic-skeleton/context controls plus exact F/G physical checks."""
from verify import *
from functools import lru_cache

def delta(p):return add(shift(p,1),scale(shift(p,-1),-1))
def lap(p):return add(shift(p,1),shift(p,-1),scale(p,-2))
def F(p):return add(mul(mon(4),p),scale(mul(mon(1),add(mul([1,1,1],shift(p,1)),scale(mul([1,-1,1],shift(p,-1)),-1))),2))
def G(p):return mul(mon(1),add(mul([0,0,-2,0,1],p),mul([1,-3,4,-2],shift(p,-1)),mul([1,3,4,2],shift(p,1))))
def L(p):
    xp=mul(mon(2),p)
    return add(D(xp),mul(mon(2),D(p)),scale(mul(mon(4),p),-1),scale(xp,-2))

def clone(a,r,m):
    vs=[(i,0) for i in range(len(a)) if i!=r]+[(r,j) for j in range(m)];o=[0]*len(vs)
    for i,j in combinations(range(len(vs)),2):
        v,k=vs[i];w,l=vs[j];u,t=(i,j) if (k<l if v==w else a[v]>>w&1) else (j,i);o[u]|=1<<t
    return tuple(o)
def port_wrapper(a,ports,s):
    n=len(a);m=len(s);w=n+m;o=list(a)+[v<<n for v in s]+[0]
    for i in range(n):
        if i in ports:o[i]|=1<<w
        else:o[w]|=1<<i
    for j in range(n,n+m):o[j]|=(1<<n)-1;o[w]|=1<<j
    validate(o);return tuple(o)

def lmul(p,q):
    o={}
    for i,a in p.items():
        for j,b in q.items():o[i+j]=o.get(i+j,0)+a*b
    return {i:c for i,c in o.items() if c}
def kernel(r,eps):
    out={}
    for j in range(r+1):
        c=(-1)**(r-j)*comb(r,j);out[j]=out.get(j,0)+c;out[-j]=out.get(-j,0)+(-1)**eps*c
    return {i:c for i,c in out.items() if c}

def skeleton_operator(a,root):
    """Enumerate cyclic skeletons literally; canonical least outside vertex."""
    outside=[i for i in range(len(a)) if i!=root];N=len(outside);bybit=[[] for _ in outside]
    for size in range(1,N+1):
        for subset in combinations(range(N),size):
            sm=sum(1<<i for i in subset)
            for tail in permutations(subset[1:]):
                word=(subset[0],)+tail
                for marks in range(1<<size):
                    r=marks.bit_count()
                    if not r and (size<3 or not size%2):continue
                    ok=True
                    for j in range(size):
                        u=outside[word[j]];v=outside[word[(j+1)%size]]
                        if marks>>j&1:
                            if not(a[u]>>root&1 and a[root]>>v&1):ok=False;break
                        elif not a[u]>>v&1:ok=False;break
                    if not ok:continue
                    ker=kernel(r,(1-size)%2) if r else {0:2}
                    bybit[subset[0]].append((sm,ker))
    @lru_cache(None)
    def pack(mask):
        if not mask:return {(0,0):1}
        bit=mask&-mask;i=bit.bit_length()-1;out={}
        for (power,j),c in pack(mask^bit).items():out[(power+1,j)]=out.get((power+1,j),0)+c
        for sm,ker in bybit[i]:
            if sm&mask!=sm:continue
            for (power,j),c in pack(mask^sm).items():
                for s,b in ker.items():out[(power+1,j+s)]=out.get((power+1,j+s),0)+c*b
        return {key:c for key,c in out.items() if c}
    return pack((1<<N)-1)
def apply_skeleton(op,f):
    out=[0]
    for (power,j),c in op.items():out=add(out,scale(mul(mon(power),shift(f,j)),c))
    return out

def run():
    report={};lower=0
    require(delta([1])==lap([1])==[0],'constant boundary')
    for k in range(31):
        q=basis(k,1);r=basis(k,0)
        require(lap(q)==(scale(basis(k-1,1),2*k*(2*k+1)) if k else [0]),'lap q')
        require(lap(r)==(scale(basis(k-1,0),(2*k+2)*(2*k+1)) if k else [2]),'lap r')
        require(delta(q)==(scale(basis(k-1,0),2*(2*k+1)) if k else [2]),'delta q')
        require(delta(r)==add(scale(q,4*(k+1)),scale(basis(k-1,1),2*k*(k+1)*(2*k+1)) if k else [0]),'delta r');lower+=4
        for p in (q,r):
            require(F(p)==add(scale(mul(mon(1),D(mul(mon(1),p))),2),scale(D(p),2),scale(mul(mon(4),p),-1),scale(mul(mon(2),p),-2)),'F identity')
            require(G(p)==L(mul(mon(1),p)),'G identity')
            require(F(p)==add(mul(mon(4),p),scale(mul(mon(2),p),4),scale(mul(mon(2),lap(p)),2),scale(mul(add(mon(3),mon(1)),delta(p)),2)),'F manifest positivity')
            require(G(p)==add(mul(mon(5),p),scale(mul(mon(3),p),6),scale(mul(mon(1),p),2),mul(add(mon(1),scale(mon(3),4)),lap(p)),mul(add(scale(mon(2),3),scale(mon(4),2)),delta(p))),'G manifest positivity')
        for par in (0,1):
            b=basis(k,par);cs=coordinates(F(b),par);require(all(c>=0 for c in cs),'F positivity')
            cs=coordinates(G(b),1-par);require(all(c>=0 for c in cs),'G positivity')
            top=2*k*k+14*k+(11 if par else 15)
            middle=k**4+20*k**3+(40 if par else 52)*k*k+(32 if par else 52)*k+(9 if par else 19)
            down=2*k*(k if par else k+1)*(2*k+1)*(2*k*k+2*k+1)
            rhs=add(basis(k+2,par),scale(basis(k+1,par),top),scale(b,middle),scale(basis(k-1,par),down) if k else [0])
            require(L(b)==rhs,'L bands')
    report['exact_lowering_identities']=lower
    # Check parity kernels directly against explicit assignments of small clone sets.
    kernelchecks=0
    for m,r,eps,d in product(range(7),range(1,4),(0,1),range(4)):
        literal=0
        for word in product(range(d+r),repeat=m):
            used={v for v in word if v>=d}
            if len(used)==r and sum(v>=d for v in word)%2==eps:literal+=1
        computed=sum(c*(d+j)**m for j,c in kernel(r,eps).items())
        require(computed==2*literal,('parity kernel',m,r,eps,d));kernelchecks+=1
    report['literal_parity_kernel_checks']=kernelchecks
    graphchecks=0;basischecks=0;contextcount=0
    for n in range(1,5):
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask)
            for root in range(n):
                op=skeleton_operator(a,root);physical={}
                for m in range(5):
                    physical[m]=bin_polynomial(clone(a,root,m)) if n-1+m else [1]
                    require(apply_skeleton(op,mon(m))==physical[m],('skeleton graph',a,root,m));graphchecks+=1
                for par,k in product((0,1),range(2)):
                    b=basis(k,par);out=apply_skeleton(op,b);co=coordinates(out,(n-1+par)%2)
                    require(all(v>=0 for v in co),('skeleton basis',a,root,par,k));basischecks+=1
                contextcount+=1
    report.update(rooted_contexts=contextcount,skeleton_physical_profiles=graphchecks,skeleton_basis_tests=basischecks)
    cores=[((4,1,2),(0,1),F,13),((0,1,3,7),(0,2),G,19)]
    witnesses=[];portchecks=0
    for core,ports,operator,target in cores:
        witness=port_wrapper(core,ports,(0,));require(h_subsets(witness)[-1]==target,'prime witness')
        witnesses.append({'core':core,'ports':ports,'rows':witness,'H':target})
        for n in range(1,5):
            for mask in range(1<<(n*(n-1)//2)):
                a=tournament(n,mask);require(bin_polynomial(port_wrapper(core,ports,a))==operator(bin_polynomial(a)),('physical port transfer',core,ports,a));portchecks+=1
    report.update(prime_witnesses=witnesses,physical_port_profiles=portchecks)
    print(json.dumps(report,sort_keys=True,indent=2))
if __name__=='__main__':run()
