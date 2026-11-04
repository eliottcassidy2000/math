"""Independent literal-query audit of exposed-anchor laminar module reversals.
No producer imports and no polynomial response banks are used.
"""
from itertools import combinations,product
from functools import lru_cache
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent.parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

@lru_cache(None)
def rows(n,code):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        u,v=(i,j) if code>>k&1 else (j,i);a[u]|=1<<v
    return tuple(a)

def hp(a,start=None):
    n=len(a);full=(1<<n)-1
    if not n:return 1
    @lru_cache(None)
    def tail(mask,v):
        if mask==full:return 1
        avail=a[v]&~mask;s=0
        while avail:
            bit=avail&-avail;avail-=bit;s+=tail(mask|bit,bit.bit_length()-1)
        return s
    return sum(tail(1<<u,u) for u in range(n)) if start is None else tail(1<<start,start)

@lru_cache(None)
def value(n,code,x,c):
    a=list(rows(n,code))+[1<<(n+1),0]
    for u in range(n):
        for mask,port in [(x,n),(c,n+1)]:
            if mask>>u&1:a[u]|=1<<port
            else:a[port]|=1<<u
    return hp(tuple(a))

@lru_cache(None)
def modules(n,code):
    a=rows(n,code);out=[]
    for m in range(1,1<<n):
        if m.bit_count()<2:continue
        if all(not(m>>u&1) and a[u]&m in (0,m) or m>>u&1 for u in range(n)):out.append(m)
    return tuple(out)

def flip(n,m):
    return sum(1<<i for i,(a,b) in enumerate(combinations(range(n),2)) if m>>a&1 and m>>b&1)

@lru_cache(None)
def anchor_parity(n,code,B,v):
    a=rows(n,code);vs=[u for u in range(n) if B>>u&1];core=tuple(sum(1<<j for j,w in enumerate(vs) if a[u]>>w&1) for u in vs)
    rev=tuple(sum(1<<j for j in range(len(vs)) if core[j]>>i&1) for i in range(len(vs)))
    k=vs.index(v);return (hp(core,k)+hp(rev,k))%2

def battery(n,code,B,v):
    a=rows(n,code);vs=[u for u in range(n) if B>>u&1]
    outside=sum(1<<u for u in range(n) if not(B>>u&1) and a[u]>>vs[0]&1)
    d=(a[v]&B).bit_count()
    qs=[(outside|(1<<u)|((1<<v) if a[v]>>u&1 else 0),outside|(1<<v)) for u in vs if u!=v]
    if d%2:qs.append((outside|(1<<v),outside|(1<<v)))
    if (len(vs)-1+d)%2:qs.append((outside,outside|(1<<v)))
    return qs,(outside|(1<<v),outside),d

def masks(m):
    return [s for s in range(m+1) if s&m==s]

def main():
    expected=set();counts=[];families=covered=0
    for n in range(1,6):
        nf=na=nc=0
        for code in range(1<<(n*(n-1)//2)):
            ms=modules(n,code)
            for choice in range(1,1<<len(ms)):
                selected=tuple(m for k,m in enumerate(ms) if choice>>k&1)
                if any(a&b and a&b not in (a,b) for a,b in combinations(selected,2)):continue
                exposed=[]
                for B in selected:
                    if any(B!=a and B&a==B for a in selected):continue
                    union=0
                    for a in selected:
                        if a!=B and a&B==a:union|=a
                    exposed.extend((B,v) for v in range(n) if (B&~union)>>v&1)
                if exposed:
                    nf+=1
                    for B,v in exposed:expected.add((n,code,selected,B,v));na+=1
                else:nc+=1
        families+=nf;covered+=nc;counts.append({'n':n,'families':nf,'anchors':na,'covered':nc})
    data=json.loads((ROOT/'simultaneous_response/anchored.json').read_text());seen=set();literal=0
    for item in data['records']:
        n,code=item['n'],item['source_code'];ms=tuple(item['module_masks']);B,v=item['chosen_outer_module'],item['anchor'];key=n,code,ms,B,v
        need(key in expected and key not in seen,('unexpected/duplicate anchor',key));seen.add(key)
        target=code
        for m in ms:target^=flip(n,m)
        need(target==item['target_code'],('edge XOR',key))
        qs,ep,d=battery(n,code,B,v);parity=anchor_parity(n,code,B,v)
        need([list(q) for q in qs]==item['query_masks'] and list(ep)==item['endpoint_query'],('physical masks',key))
        need(d==item['outdegree_in_outer'] and parity==item['endpoint_parity'],('local data',key))
        whole=code^flip(n,B);sums=[sum(value(n,k,*q) for q in qs) for k in (code,whole,target)]
        points=[value(n,k,*ep) for k in (code,whole,target)]
        need(sums==[item['source_query_sum'],item['single_outer_query_sum'],item['all_reversals_query_sum']],('literal sums',key))
        need(points==item['endpoint_counts_source_single_all'],('literal endpoint',key))
        need(sums[1]==sums[2] and points[1]==points[2],('exact invisibility',key))
        if B.bit_count()>=3:need((sums[0]-sums[2])%4==2*(1-parity),('battery residue',key))
        else:need(parity==1,('n2 endpoint boundary',key))
        if parity:need((points[0]-points[2])%4==2,('endpoint residue',key))
        literal+=3*(len(qs)+1)
    need(seen==expected,'anchor completeness')
    avg=0
    for n in range(1,5):
        full=(1<<n)-1
        for code in range(1<<(n*(n-1)//2)):
            for A in modules(n,code):
                ext=full^A;vs=[u for u in range(n) if A>>u&1]
                for x,c,cy in product(masks(ext),masks(ext),(0,A)):
                    sums=[sum(value(n,k,x|(1<<u),c|cy) for u in vs) for k in (code,code^flip(n,A))]
                    need(sums[0]==sums[1],('averaging',n,code,A,x,c,cy));avg+=1
    source={f:hashlib.sha256((ROOT/'simultaneous_response'/f).read_bytes()).hexdigest() for f in ['ANCHORED.md','anchored.py','search.py','anchored.json']}
    return {'status':'PASS','method':'Independent original-module/family enumeration and literal subset Hamilton-path counting; no producer imports or polynomial banks','exposed_families':families,'anchor_records':len(seen),'covered_families_not_promoted':covered,'literal_query_occurrences':literal,'distinct_literal_queries':value.cache_info().currsize,'exact_averaging_checks':avg,'by_order':counts,'source_sha256':source}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
