"""Independent exact checks. Standard library only; run python3 -B verify.py."""
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import comb, factorial, gcd
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def require(ok, message):
    if not ok: raise RuntimeError(message)

def decode(n, code):
    rows=[0]*n; b=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>b&1: rows[i]|=1<<j
            else: rows[j]|=1<<i
            b+=1
    return rows

def induced(rows, keep):
    return [sum(((rows[i]>>j)&1)<<k for k,j in enumerate(keep)) for i in keep]

def tables(rows):
    n=len(rows); full=(1<<n)-1
    end=[[0]*n for _ in range(1<<n)]
    start=[[0]*n for _ in range(1<<n)]
    h=[0]*(1<<n);h[0]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not mask>>v&1:continue
            rest=mask^(1<<v)
            if not rest:end[mask][v]=start[mask][v]=1
            else:
                end[mask][v]=sum(end[rest][u] for u in range(n) if rest>>u&1 and rows[u]>>v&1)
                start[mask][v]=sum(start[rest][u] for u in range(n) if rest>>u&1 and rows[v]>>u&1)
        h[mask]=sum(end[mask])
    return end,start,h

def signature(rows):
    end,start,h=tables(rows); n=len(rows);full=(1<<n)-1
    U=[sum(start[s][v]*h[full^s] for s in range(1,full+1)) for v in range(n)]
    W=[[0]*n for _ in range(n)]
    for s in range(1,full):
        for v in range(n):
            if not end[s][v]:continue
            for w in range(n):
                q=end[s][v]*start[full^s][w]
                W[v][w]+=q;W[w][v]+=q
    return h[-1],U,W

def H(rows):return tables(rows)[2][-1]

def path_dfs(rows):
    n=len(rows);full=(1<<n)-1
    def walk(mask,v):
        if mask==full:return 1
        q=rows[v]&~mask;ans=0
        while q:
            bit=q&-q;q-=bit;ans+=walk(mask|bit,bit.bit_length()-1)
        return ans
    return sum(walk(1<<v,v) for v in range(n)) if n else 1

def closure(rows,O,t):
    # B -> TT_t -> root; root -> O and B\O -> root.
    n=len(rows);r=n+t;out=list(rows)+[0]*(t+1)
    for i in range(n):
        for z in range(n,n+t):out[i]|=1<<z
        if O>>i&1:out[r]|=1<<i
        else:out[i]|=1<<r
    for z in range(n,n+t):
        for w in range(z+1,r+1):out[z]|=1<<w
    return out

def automorphisms(rows):
    # Degree-partitioned backtracking, with adjacency tests at each assignment.
    n=len(rows);degs=[r.bit_count() for r in rows]
    candidates=[[j for j in range(n) if degs[i]==degs[j]] for i in range(n)]
    order=sorted(range(n),key=lambda i:len(candidates[i]));p=[-1]*n;ans=[]
    def rec(k,used):
        if k==n:ans.append(tuple(p));return
        i=order[k]
        for j in candidates[i]:
            if used>>j&1:continue
            if any(((rows[i]>>a)&1)!=((rows[j]>>p[a])&1) for a in order[:k]):continue
            p[i]=j;rec(k+1,used|1<<j);p[i]=-1
    rec(0,0);return ans

def verify_seeds():
    found=[json.loads(s) for s in (HERE/'search.jsonl').read_text().splitlines() if s.strip()]
    result=[]
    for item in found:
        rows=item['rows'];O=item['O'];h,U,W=signature(rows)
        ids=[v for v in range(len(rows)) if O>>v&1]
        T=sum(W[v][w] for v in ids for w in ids if v<w);M=sum(U[v] for v in ids)
        require(h==item['h'] and M==item['M'] and h-T==1,'search seed failed')
        group=automorphisms(rows)
        stab=[p for p in group if {p[i] for i in ids}==set(ids)]
        require(len(stab)==1,'unit-offset stabilizer')
        rays=[]
        for t in [1,2]:
            graph=closure(rows,O,t);count=H(graph)
            require(count==1+(1<<t)*M,'literal ray failed')
            require(count==path_dfs(graph),'independent DFS failed')
            rays.append(count)
        power=(M&-M).bit_length()-1
        require(M>>power==item['c'],'odd slope failed')
        result.append(dict(item,T=T,U=U,selected=ids,aut_base=len(group),stabilizer=len(stab),s0=power+1,ray_counts=rays))
    return result

def inherited_seeds():
    wanted={1,3,5,7,13,19,27,1039};result=[]
    bank=ROOT/'05-knowledge/results/kps_s21_two_adic_lines_20261110/line_seeds.txt'
    for line in bank.read_text().splitlines():
        if line.startswith('#'):continue
        c,s0,kind,n,code,v,*old=line.split();c=int(c)
        if c not in wanted:continue
        require(kind=='tail','expected tail certificate')
        n=int(n);v=int(v);A=decode(n,int(code));keep=[i for i in range(n) if i!=v]
        B=induced(A,keep);O=sum(((A[v]>>i)&1)<<j for j,i in enumerate(keep))
        h,U,W=signature(B);ids=[i for i in range(n-1) if O>>i&1]
        T=sum(W[i][j] for i in ids for j in ids if i<j);M=sum(U[i] for i in ids)
        require(h-T==1 and M==c*(1<<(int(s0)-1)),'inherited certificate failed')
        require(H(closure(B,O,1))==1+2*M,'inherited literal closure')
        result.append(dict(c=c,s0=int(s0),base_rows=B,O=O,h=h,T=T,M=M,root_order=n))
    return result

def small_symmetry():
    checked=units=nontrivial=0
    # Every labelled tournament through order 5, every selected subset.
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            rows=decode(n,code);h,U,W=signature(rows);group=automorphisms(rows)
            for O in range(1<<n):
                ids=[i for i in range(n) if O>>i&1]
                T=sum(W[i][j] for i in ids for j in ids if i<j);M=sum(U[i] for i in ids)
                g=sum({p[i] for i in ids}==set(ids) for p in group)
                require(h%g==T%g==M%g==(h-T)%g==0,'divisibility failed')
                if h-T==1:units+=1;require(g==1,'rigidity failed')
                nontrivial+=g>1;checked+=1
    return dict(coloured_tournaments=checked,unit_sets=units,nontrivial_stabilizers=nontrivial)

def extension_check():
    cases=0;polynomial_cases=0;collision=None
    for n,k in [(1,1),(2,1),(3,1),(3,2),(4,1)]:
        for code in range(1<<(n*(n-1)//2)):
            base=decode(n,code);N=n+k
            variable=[(i,j) for i in range(N) for j in range(i+1,N) if j>=n]
            edge_index={e:i for i,e in enumerate(variable)};r=len(variable)
            cylinders=[];cover_orders=Counter()
            for p in permutations(range(N)):
                mask=value=0;valid=True
                for a,b in zip(p,p[1:]):
                    if a<n and b<n:
                        if not base[a]>>b&1:valid=False;break
                    else:
                        bit=1<<edge_index[min(a,b),max(a,b)];mask|=bit
                        if a<b:value|=bit
                if valid:
                    cylinders.append((mask,value))
                    runs=sum(p[i]<n and (i==0 or p[i-1]>=n) for i in range(N))
                    cover_orders[runs]+=1
            counts=[]
            for bits in range(1<<r):
                rows=base+[0]*k
                for i,(a,b) in enumerate(variable):
                    if bits>>i&1:rows[a]|=1<<b
                    else:rows[b]|=1<<a
                count=sum(bits&mask==value for mask,value in cylinders)
                require(count==H(rows),'extension cylinder formula');counts.append(count);polynomial_cases+=1
            mean=Fraction(sum(counts),1<<r)
            expected=sum(Fraction(num,1<<(j+k-1)) for j,num in cover_orders.items())
            require(mean==expected,'extension mean law')
            cases+=1
            # Transitive3 extension hostiles: finite positive mean does not fill every odd total.
            if n==3 and k==1 and H(base)==1 and collision is None:
                collision=dict(n=n,k=k,base=base,mean=str(mean),histogram=dict(sorted(Counter(counts).items())))
    return dict(base_families=cases,extensions=polynomial_cases,hostile=collision)

def prime_geometry():
    result=[]
    for p in [3,5,7,13,19,1039]:
        require(all(p%d for d in range(2,int(p**.5)+1)),'not prime')
        v=2%p;order=1
        while v!=1:v=v*2%p;order+=1
        result.append(dict(p=p,ord2=order,doubling_invariant_sets=(1<<((p-1)//(2*order))) if order%2 else 0,square_affine_order=p*(p-1)//2 if p%4==3 else None))
    squares={x*x%7 for x in range(1,7)}
    P7=[sum(1<<j for j in range(7) if (j-i)%7 in squares) for i in range(7)]
    h,U,W=signature(P7)
    require(h==189,'P7 path count');require(len(set(U))==1,'P7 vertex response')
    w=W[0][1];require(all(W[i][j]==w for i in range(7) for j in range(i+1,7)),'P7 pair homogeneity')
    require(len(automorphisms(P7))==21,'P7 automorphisms')
    possible=[m for m in range(8) if comb(m,2)*w==h-1]
    require(not possible,'expected P7 no seed')
    return dict(primes=result,P7=dict(H=h,U=U[0],W=w,translation_path_orbits=h//7,affine_path_orbits=h//21,unit_set_sizes=possible))

def main():
    report=dict(new_seeds=verify_seeds(),inherited_seeds=inherited_seeds(),small_symmetry=small_symmetry(),extension=extension_check(),prime_geometry=prime_geometry())
    (HERE/'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['inherited_seeds']},indent=2))

if __name__=='__main__':main()
