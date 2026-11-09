#!/usr/bin/env python3
"""Independent subset-path / Newton-profile audit. No producer imports or inputs."""
import itertools, math, json, hashlib
from functools import lru_cache

def require(x,s):
    if not x: raise ValueError(s)
def graphs(n):
    edges=list(itertools.combinations(range(n),2))
    for bits in range(1<<len(edges)):
        rows=[0]*n
        for k,(i,j) in enumerate(edges):
            a,b=(i,j) if (bits>>k)&1 else (j,i);rows[a]|=1<<b
        yield tuple(rows)
def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def substitute(template,parts):
    offsets=[];n=0
    for p in parts:offsets.append(n);n+=len(p)
    rows=[]
    for i,p in enumerate(parts):
        outside=sum(((1<<len(q))-1)<<offsets[j] for j,q in enumerate(parts) if template[i]>>j&1)
        rows.extend((r<<offsets[i])|outside for r in p)
    return tuple(rows)
def path_counts(rows):
    n=len(rows);d=[[0]*n for _ in range(1<<n)]
    for v in range(n):d[1<<v][v]=1
    for m in range(1,1<<n):
        for v in range(n):
            if m>>v&1:
                prev=m^(1<<v)
                if prev:d[m][v]=sum(d[prev][u] for u in range(n) if prev>>u&1 and rows[u]>>v&1)
    return (1,)+tuple(sum(row) for row in d[1:])
def profile(rows):
    h=path_counts(rows)
    @lru_cache(None)
    def partitions(m):
        if m==0:return (1,)
        least=m&-m;out=[0]*(m.bit_count()+1);s=m
        while s:
            if s&least:
                for k,c in enumerate(partitions(m^s)):out[k+1]+=h[s]*c
            s=(s-1)&m
        return tuple(out)
    return partitions(len(h)-1)
def q(pc,t):
    f=1;out=pc[0]
    for k in range(1,len(pc)):
        f*=t-k+1;out+=pc[k]*f
    return out
def from_values(vals):
    out=[];k=0
    while vals:
        require(vals[0]%math.factorial(k)==0,'Newton integral')
        out.append(vals[0]//math.factorial(k));vals=[y-x for x,y in zip(vals,vals[1:])];k+=1
    return tuple(out)
@lru_cache(None)
def ordered(a,b):return from_values([q(a,j)*q(b,j) for j in range(len(a)+len(b)-1)])
@lru_cache(None)
def collar_pair(a,b):
    return from_values([j*(q(a,j)*q(b,j)+(q(a,j+1)-q(a,j))*(q(b,j+1)-q(b,j))+(q(a,j)-q(a,j-1))*(q(b,j)-q(b,j-1))) for j in range(len(a)+len(b))])
def c3_h(a,b,c):
    fs=[tuple(math.factorial(k)*v for k,v in enumerate(p)) for p in (a,b,c)]
    total=0
    for i,j,k in itertools.product(*(range(1,len(f)) for f in fs)):
        spread=max(i,j,k)-min(i,j,k)
        if spread<=1:total+=(3 if spread==0 else 1)*fs[0][i]*fs[1][j]*fs[2][k]
    return total

def compile_old(seed,limit=4095):
    N=1+(limit-1).bit_length();bank=[{} for _ in range(N+1)]
    def put(n,p,d):
        if p[1]<=limit:bank[n][p]=min(d,bank[n].get(p,10**9))
    for p in seed:put(len(p)-1,p,0)
    for n in range(2,N+1):
        for i in range(1,n):
            for a,da in list(bank[i].items()):
                for b,db in list(bank[n-i].items()):
                    if a[1]*b[1]<=limit:put(n,ordered(a,b),max(da,db))
        for i in range(1,n-1):
            for a,da in bank[i].items():
                for b,db in bank[n-i-1].items():
                    h=2*a[1]*b[1]+(q(a,2)-a[1])*(q(b,2)-b[1])
                    if h<=limit:put(n,collar_pair(a,b),1+max(da,db))
    generators={p[1]:0 for p in seed}
    # Every bounded-order expression is valid as an exterior factor; seeds too.
    for row in bank:
        for p,d in row.items():generators[p[1]]=min(d,generators.get(p[1],10**9))
    depth={1:0}
    for h in range(2,limit+1):
        candidates=[generators[h]] if h in generators else []
        for g,d in generators.items():
            if 1<g<h and h%g==0 and h//g in depth:candidates.append(max(d,depth[h//g]))
        if candidates:depth[h]=min(candidates)
    return bank,depth

def main():
    seeds={};small=[];count=0
    for n in range(1,6):
        for r in graphs(n):
            p=profile(r);seeds.setdefault(p,r);count+=1
            if n<=3:small.append((r,p))
    bank,depth=compile_old(seeds)
    require(len(depth)==577 and 137 not in depth,'prior grammar exact membership')
    require({d:sum(v==d for v in depth.values()) for d in set(depth.values())}=={0:79,1:390,2:107,3:1},'prior depths')
    controls=0
    for (ra,a),(rb,b),(rc,c) in itertools.product(small,repeat=3):
        r=substitute((2,4,1),(ra,rb,rc));h=path_counts(r)[-1]
        require(h==c3_h(a,b,c),'C3 literal versus kernel')
        require(h>=sum(2**len(p) for p in (ra,rb,rc))-3,'strong cycle split bound')
        controls+=1
    family=[]
    for r,p in seeds.items():
        a=p;pc=r
        out=substitute((2,4,1),(a,transitive(2),transitive(2)))
        expected=4*q(pc,3)+5*q(pc,2)-11*pc[1]
        require(path_counts(out)[-1]==expected==c3_h(pc,profile(transitive(2)),profile(transitive(2))),'family identity')
        family.append(expected)
    witness=substitute((2,4,1),(transitive(2),transitive(2),transitive(3)))
    require(path_counts(witness)[-1]==137,'137 witness')
    primitive=set();witnesses={}
    for a,b,c in itertools.combinations_with_replacement([p for row in bank[1:6] for p in row],3):
        h=c3_h(a,b,c)
        if h<=4095 and h not in depth:
            primitive.add(h);witnesses.setdefault(h,(seeds[a],seeds[b],seeds[c]))
    # Independent literal realization of every escape, not just its scalar profile.
    for h,parts in witnesses.items():
        require(path_counts(substitute((2,4,1),parts))[-1]==h,'every escape graph witness')
    generic=0;template_count=0
    for k in (3,4):
        for template in graphs(k):
            cycle=next((p for p in itertools.permutations(range(k)) if all(template[p[i]]>>p[(i+1)%k]&1 for i in range(k))),None)
            if cycle is None:continue
            template_count+=1
            for sizes in itertools.product((1,2),repeat=k):
                h=path_counts(substitute(template,tuple(transitive(n) for n in sizes)))[-1]
                n=sum(sizes);a,b=divmod(n,k);balanced=(k-b)*2**a+b*2**(a+1)-k
                require(h>=sum(2**s for s in sizes)-k>=balanced,'generic cycle bound')
                generic+=1
    sharp=0
    for sizes in itertools.product(range(1,9),repeat=3):
        pcs=[from_values([j**n for j in range(n+1)]) for n in sizes]
        h=c3_h(*pcs);lower=2**(sum(sizes)-2)+1
        require(h>=lower,'sharp C3 bound')
        if sorted(sizes)[:2]==[1,1]:require(h==lower,'sharp equality')
        sharp+=1
    encoded=json.dumps(sorted(depth.items()),separators=(',',':')).encode()
    result={'status':'PASS','seed_graphs':count,'seed_profiles':len(seeds),'old_profile_counts':[len(row) for row in bank[1:]],'old_scalar_count':len(depth),'old_scalar_depth_sha256':hashlib.sha256(encoded).hexdigest(),'literal_C3_controls':controls,'literal_generic_template_controls':generic,'cycle_certified_templates':template_count,'sharp_C3_size_controls':sharp,'literal_escape_controls':len(witnesses),'family_controls':len(family),'family_H_values':sorted(set(family)),'escape_witness_rows':{str(h):substitute((2,4,1),parts) for h,parts in sorted(witnesses.items())},'137_rows':witness,'primitive_escape_values':sorted(primitive)}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
