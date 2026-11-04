#!/usr/bin/env python3
"""Standalone exact controls; no producer imports or third-party dependencies."""
import itertools, json, math

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def graph(n,k):
    a=[[0]*n for _ in range(n)]
    for bit,(i,j) in enumerate(itertools.combinations(range(n),2)):
        a[i][j]=(k>>bit)&1;a[j][i]=1-a[i][j]
    return a

def count(a):
    n=len(a); dp={}
    for v in range(n):dp[(1<<v,v)]=1
    for s in range(1,1<<n):
        for v in range(n):
            c=dp.get((s,v),0)
            for w in range(n):
                if c and not (s>>w)&1 and a[v][w]:
                    key=(s|(1<<w),w);dp[key]=dp.get(key,0)+c
    return sum(dp.get(((1<<n)-1,v),0) for v in range(n))

def components(a):
    n=len(a);r=[[bool(i==j or a[i][j]) for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):r[i][j]=r[i][j] or (r[i][k] and r[k][j])
    left=set(range(n));out=[]
    while left:
        v=min(left);s=sorted(w for w in left if r[v][w] and r[w][v]);left.difference_update(s);out.append(s)
    return out

def induce(a,vs):return [[a[v][w] for w in vs] for v in vs]
def blow(a,rs):
    labels=[(i,j) for i,r in enumerate(rs) for j in range(r)]
    return [[int(x<y) if i==j else a[i][j] for j,y in labels] for i,x in labels]
def prime(x):return x>=2 and all(x%d for d in range(2,math.isqrt(x)+1))
def reverse(a):return list(map(list,zip(*a)))
def canon(a):
    n=len(a)
    return min(tuple(b[p[i]][p[j]] for i in range(n) for j in range(i+1,n)) for b in (a,reverse(a)) for p in itertools.permutations(range(n)))
def lower(n):
    m=n//2;base=2*m;b=math.comb(2*m,m);t=0;x=1
    while x<b:x*=base;t+=1
    return max(0,t-2)

def main():
    rows=brute=cores=0; witnessed=set();bounds=0
    for n in range(1,5):
        for k in range(1<<(n*(n-1)//2)):
            a=graph(n,k);cs=components(a);non=[s for s in cs if len(s)>1]
            for rs in itertools.product((1,2),repeat=n):
                g=blow(a,rs);h=count(g);rows+=1
                factors=[count(blow(induce(a,s),[rs[i] for i in s])) for s in cs]
                need(h==math.prod(factors),'SCC factorization')
                if len(non)==0:need(h==1,'transitive padding')
                elif len(non)==1:
                    s=non[0];core=induce(a,s);hc=count(blow(core,[rs[i] for i in s]));need(h==hc,'exact core reduction');cores+=1
                    if prime(h):witnessed.add(h)
                    L=(h-1).bit_length()-1
                    need(all(rs[i]<=L for i in s),'triangle size cap')
                else:need(not prime(h),'multiple strong components cannot yield prime')
                if len(g)<=5:
                    hp=sum(all(g[p[i]][p[i+1]] for i in range(len(g)-1)) for p in itertools.permutations(range(len(g))))
                    need(h==hp,'independent literal path count');brute+=1
            need(canon(a)==canon(reverse(a)),'reverse canonical type')
    # Genuine hostile examples beyond the small quotient universe: C3 + C3,
    # with arbitrary source/sink padding. Extracting either core alone would be unsound.
    c=graph(3,2)
    need(len(components(c))==1,'chosen triangle')
    hostile=[]
    for padding in range(3):
        n=6+padding;a=[[0]*n for _ in range(n)]
        blocks=[list(range(padding)),list(range(padding,padding+3)),list(range(padding+3,n))]
        for b,vs in enumerate(blocks):
            for i,v in enumerate(vs):
                for j,w in enumerate(vs):a[v][w]=int(i<j) if b==0 else c[i][j]
                for ws in blocks[b+1:]:
                    for w in ws:a[v][w]=1
        h=count(a);need(h==9 and not prime(h),'two-core hostile');hostile.append(h)
    for n in range(3,1001):
        actual=sum(prime(p) and p!=7 for p in range(3,n+1));need(lower(n)<=actual,'elementary exact prime bound');bounds+=1
    table=[]
    for n in (10,100,1000,10000,100000):
        L=(n-1).bit_length()-1;b=lower(n)
        for K in (1,10,100):
            d=0
            while K*L**d<b:d+=1
            table.append({'N':n,'cores':K,'elementary_required_prime_values':b,'L':L,'capacity_necessary_order':d})
    print(json.dumps({'status':'PASS','inflated_quotients':rows,'literal_path_crosschecks':brute,'single_core_cases':cores,'finite_prime_values':sorted(witnessed),'two_core_hostile_counts':hostile,'elementary_bound_cutoffs':bounds,'capacity_table':table},indent=2,sort_keys=True))
if __name__=='__main__':main()
