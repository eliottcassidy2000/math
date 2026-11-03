#!/usr/bin/env python3
"""Standalone exact controls; every validation remains active under python -O."""
import itertools, math, json

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

B = ((16,1,3,7,14), (20,1,2,7,14))
def edge(g,i,j): return (g[i]>>j)&1

def tournament(g):
    return all(not edge(g,i,i) for i in range(len(g))) and all(edge(g,i,j)+edge(g,j,i)==1 for i in range(len(g)) for j in range(i))

def ham(g):
    n=len(g); dp={}
    for j in range(n): dp[(1<<j,j)]=1
    for mask in range(1,1<<n):
        for j in range(n):
            count=dp.get((mask,j),0)
            if count:
                rem=g[j]&~mask
                while rem:
                    bit=rem&-rem; t=bit.bit_length()-1; rem-=bit
                    key=(mask|bit,t); dp[key]=dp.get(key,0)+count
    return sum(dp.get(((1<<n)-1,j),0) for j in range(n))

def cover(g):
    n=len(g); ordered=[0]*(n+1)
    for p in itertools.permutations(range(n)):
        for cuts in range(1<<(n-1)):
            if all((cuts>>j)&1 or edge(g,p[j],p[j+1]) for j in range(n-1)):
                ordered[1+cuts.bit_count()]+=1
    require(all(ordered[j]%math.factorial(j)==0 for j in range(1,n+1)), 'cover division')
    return [ordered[j]//math.factorial(j) for j in range(1,n+1)]

def read_seed(g):
    return tuple(g[i]|(0 if i==0 else 1<<5) for i in range(5))+(1,)

def ordinal(blocks):
    total=sum(map(len,blocks)); rows=[]; offset=0
    for block in blocks:
        after=((1<<total)-1)^((1<<(offset+len(block)))-1)
        rows.extend((row<<offset)|after for row in block)
        offset+=len(block)
    return tuple(rows)

for g in B:
    require(tournament(g),'seed tournament')
    require(cover(g)==[9,33,31,10,1], 'cover profile')
    for u in range(8):
        val=sum(cover(g)[j-1]*math.prod(range(u-j+1,u+1)) for j in range(1,6))
        require(val==u**5+6*u**3+2*u,'slot profile')
    for i in range(5):
        seen=1<<i
        while True:
            nxt=seen
            for j in range(5):
                if seen>>j&1: nxt|=g[j]
            if nxt==seen: break
            seen=nxt
        require(seen==31,'strongness')
require(not any(all(edge(B[0],i,j)==edge(B[1],p[i],p[j]) for i in range(5) for j in range(5)) for p in itertools.permutations(range(5))), 'nonisomorphism')
require([ham(read_seed(g)) for g in B]==[17,19], 'read heads')
for word in itertools.product(range(2),repeat=2):
    require(ham(ordinal([B[b] for b in word]))==81,'literal tape')
    for j in range(2):
        blocks=[read_seed(B[b]) if i==j else B[b] for i,b in enumerate(word)]
        require(ham(ordinal(blocks))==(17+2*word[j])*9,'literal addressed context')
    require(ham(ordinal([read_seed(B[b]) for b in word]))==17**(2-sum(word))*19**sum(word),'literal simultaneous read')

def step(n): return (3*n+1)//2 if n&1 else n//2

def data(n,k):
    word=[]; vals=[]
    for _ in range(k):
        word.append(n&1); n=step(n); vals.append(n)
    return tuple(word), vals

def cocycles(word):
    a=s=0; out=[]
    for j,b in enumerate(word):
        a=3**b*a+b*2**j; s+=b
        out.append((s,a))
    return out

checked=0; exceptions=0; cylinder_totals={}
for k in range(1,13):
    words=set(); affine_maps=set(); offsets_by_weight={}; survivors=0
    for r in range(2**k):
        w,_=data(r,k); words.add(w); cs=cocycles(w)
        s,a=cs[-1]
        affine_maps.add((s,a))
        offsets_by_weight.setdefault(s,set()).add(a)
        require((-a*pow(3**s,-1,2**k))%(2**k)==r,'inverse residue')
        thresholds=[a//(2**j-3**s) for j,(s,a) in enumerate(cs,1) if 2**j>3**s]
        h=min(thresholds) if thresholds else None
        if h is None: survivors+=1
        rpos=r if r else 2**k
        badcount=max(0,(h-rpos)//(2**k)+1) if h is not None else None
        if badcount: exceptions+=badcount
        for q in range(6):
            n=rpos+q*2**k; ww, vals=data(n,k)
            require(ww==w,'residue itinerary')
            for j,((s,a),v) in enumerate(zip(cs,vals),1):
                require(2**j*v==3**s*n+a,'affine identity')
                require((v<n)==((2**j-3**s)*n>a),'endpoint threshold')
            require(any(v<n for v in vals)==(h is not None and n>h),'cylinder threshold')
            if h is not None:
                require((q<badcount)==all(v>=n for v in vals),'exception count')
            checked+=1
        n=rpos
        w1,_=data(n,k+1); w2,_=data(n+2**k,k+1)
        require(w1[:k]==w2[:k] and w1[k]!=w2[k], 'precision loss')
    require(len(words)==2**k,'parity bijection')
    require(len(affine_maps)==2**k,'affine interface injectivity')
    require(all(len(offsets_by_weight[s])==math.comb(k,s) for s in range(k+1)),'fixed-weight offset injectivity')
    cylinder_totals[k]=survivors

v={0:1}; selected={}
for j in range(1,257):
    v={s:v.get(s,0)+v.get(s-1,0) for s in range(j+1) if 3**s>2**j}
    num=sum(v.values())
    if j<=12: require(num==cylinder_totals[j],'survival DP')
    if j in (1,2,4,8,12,16,32,64,128,256):
        selected[str(j)]={'survivors':num,'denominator':2**j,'density':num/2**j}
for k in range(1,65):
    for q in (1,2,7):
        n=q*2**k-1; w,vals=data(n,k)
        require(all(w),'all-one cylinder')
        require(all(t>n for t in vals),'no bounded descent')
        require(vals==[3**j*q*2**(k-j)-1 for j in range(1,k+1)],'all-one formula')
for k in range(1,9):
    N=2**k-1
    counts={17**(N-z)*19**z for z in range(2**k)}
    require(len(counts)==2**k,'rich output injectivity')
print(json.dumps({'status':'PASS','seed_profiles':[9,33,31,10,1], 'read_counts':[17,19], 'positive_orbit_representatives_checked':checked,'cylinders_checked':sum(2**k for k in range(1,13)), 'finite_exceptions_summed_across_depths':exceptions,'survival_counts':selected},sort_keys=True,indent=2))
