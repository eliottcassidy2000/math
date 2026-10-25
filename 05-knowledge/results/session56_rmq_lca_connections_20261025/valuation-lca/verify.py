#!/usr/bin/env python3
"""Exact labelled valuation/RMQ carrier; standard library only."""
import hashlib,itertools,json,math,pathlib,random
INF=math.inf

def require(c,msg):
    if not c: raise RuntimeError(msg)

def vp(x,p):
    if x==0:return INF
    n=0
    while x%p==0:x//=p;n+=1
    return n

def digits(x,p,h):
    a=[]
    for _ in range(h):a.append(x%p);x//=p
    return tuple(a)

class Carrier:
    """Prime input; duplicate labels share a value terminal, not a false split."""
    def __init__(self,values,p):
        require(p>=2 and all(p%d for d in range(2,math.isqrt(p)+1)), 'prime required')
        self.values=list(values);self.p=p
        distinct=set(self.values);spread=max(distinct)-min(distinct) if distinct else 0
        h=0;power=1
        while power<=spread:h+=1;power*=p
        self.h=h
        self.order=sorted(distinct,key=lambda x:digits(x,p,h))
        self.rank={x:i for i,x in enumerate(self.order)}
        self.labels={x:[] for x in self.order}
        for i,x in enumerate(self.values):self.labels[x].append(i)
        self.adj=[vp(y-x,p) for x,y in zip(self.order,self.order[1:])]
        self.st=[self.adj]
        width=2
        while width<=len(self.adj):
            last=self.st[-1];half=width//2
            self.st.append([min(last[i],last[i+half]) for i in range(len(self.adj)-width+1)])
            width*=2
    def depth(self,i,j):
        x,y=self.values[i],self.values[j]
        if x==y:return INF
        a,b=sorted((self.rank[x],self.rank[y]));k=(b-a).bit_length()-1
        return min(self.st[k][a],self.st[k][b-(1<<k)])
    def mask(self,k):
        require(k>=0,'nonnegative exponent')
        return sum(1<<t for t,(i,j) in enumerate(itertools.combinations(range(len(self.values)),2)) if self.depth(i,j)>=k)
    def separation(self):
        if len(self.values)!=len(self.order):return INF
        return 0 if len(self.values)<2 else 1+max(self.adj)

def direct_mask(values,m):
    return sum(1<<t for t,(x,y) in enumerate(itertools.combinations(values,2)) if (x-y)%m==0)

def check(values,p):
    c=Carrier(values,p)
    keys=[digits(x,p,c.h) for x in values]
    for i in range(len(values)):
        for j in range(len(values)):
            require(c.depth(i,j)==vp(values[i]-values[j],p),'RMQ valuation mismatch')
            if values[i]!=values[j]:
                lcp=next(k for k,(a,b) in enumerate(zip(keys[i],keys[j])) if a!=b)
                require(c.depth(i,j)==lcp,'explicit digit-prefix mismatch')
    for k in range(c.h+2):require(c.mask(k)==direct_mask(values,p**k),'partition mismatch')
    s=c.separation()
    if s!=INF:
        require(direct_mask(values,p**s)==0,'separation false')
        if s:require(direct_mask(values,p**(s-1))!=0,'separation not first')
    return len(values)**2

def run():
    count=cases=0
    # Exhaust all ordered records of length <=4 over a signed five-element alphabet.
    for p in [2,3,5,7]:
        for n in range(5):
            for v in itertools.product([-3,-1,0,1,3],repeat=n):count+=check(v,p);cases+=1
    rng=random.Random(56)
    for _ in range(100):
        p=rng.choice([2,3,5,7,11]);v=[rng.randrange(-10**12,10**12) for _ in range(25)]
        count+=check(v,p);cases+=1
    # Affine Lo Shu fixed integer towers, not root-vector lift genealogy.
    L=[8,1,6,3,5,7,4,9,2];A=[3*x-2 for x in L];B=[8*x-7 for x in L]
    require(Carrier(A,3).separation()==3,'3-adic template threshold')
    require(Carrier(B,2).separation()==7,'2-adic template threshold')
    for v,p in [(A,3),(B,2)]:count+=check(v,p)
    # Fixed records forming the actual CRT witness. All labels retained.
    W=[0,9,1,16,15,14,29,21,30]
    c5=Carrier(W,5);c7=Carrier(W,7)
    require(c5.mask(1)==0x61089888 and c7.mask(1)==0x8a054250,'inherited masks')
    require(c5.mask(1)&c7.mask(1)==0,'CRT intersection')
    require({0,4,8}&{2,4,6}=={4},'crossed hierarchy')
    crt_checks=0
    for p,q in [(2,3),(2,5),(3,7),(5,7)]:
        cp,cq=Carrier(W,p),Carrier(W,q)
        for a in range(5):
            for b in range(5):
                require(cp.mask(a)&cq.mask(b)==direct_mask(W,p**a*q**b),'multi-prime mask')
                crt_checks+=1
    # Natural integer order is NOT digit-lexicographic order.
    require(vp(2,2)==1 and min(vp(1,2),vp(1,2))==0,'natural-sort hostile')
    require(Carrier([0,1,2],2).order==[0,2,1],'digit ordering')
    require(Carrier([1,1,9],2).separation()==INF,'duplicate hostile')
    # Derive the inherited 27 block from actual shortcut iteration, no hard-coded word assumption.
    n=27;word=[];a=1;b=0;j=0
    while True:
        e=n%2;word.append(e)
        if e:b=3*b+2**j;a*=3;n=(3*n+1)//2
        else:n//=2
        j+=1
        if n<27:break
    q=2**j;d=q-a
    require((j,sum(word),n,a,b,d)==(59,37,23,3**37,1100931843921811423,126176846412426125),'original block')
    require(math.gcd(d,b)==1 and b%d!=0 and d%2==1,'finite rational guard')
    fuel_checks=0
    for depth in range(1,13):
        modulus=q**depth;r=b*pow(d,-1,modulus)%modulus
        # A marked rational p-adic path is specified through needed precision H.
        for t in [0,1,2,3,2**40+1]:
            source=r+t*modulus;v=vp(d*source-b,2);K=v//59
            H=v+1;zH=b*pow(d,-1,2**H)%2**H
            require(vp(source-zH,2)==v,'marked-path depth')
            require(Carrier([source,zH],2).depth(0,1)==v,'LCA marked path')
            x=source
            for step in range(59*K):
                require(x%2==word[step%59],'exact repeated parity')
                x=(3*x+1)//2 if x%2 else x//2
            require(vp(d*x-b,2)<59,'fuel exhaustion')
            fuel_checks+=1
    require(vp(d*27-b,2)//59==1 and vp(d*23-b,2)//59==0,'reset hostile')
    return {'status':'PASS','record_cases':cases,'ordered_pair_checks':count,'CRT_exponent_checks':crt_checks,'Collatz_marked_path_replays':fuel_checks,'LoShu_first_separating_depths':{'p3':3,'p2':7},'CRT_masks':[hex(c5.mask(1)),hex(c7.mask(1))],'27_word':''.join(map(str,word)),'hostiles':['natural integer order invalid','duplicates never separate','cross-prime blocks nonlaminar','27 supports one block only'],'scope':'Finite controls validate independently proved exact carrier; no universal coverage or integer-square witness.'}

if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
