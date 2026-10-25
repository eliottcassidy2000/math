"""Independent audit; does not import producer code. Python standard library."""
from math import comb
from collections import Counter
from pathlib import Path
import hashlib,json

def require(x,msg):
    if not x: raise RuntimeError(msg)

def convolution(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out

def shifted_from_formula(p,k,s):
    # Construct shifted Burnside factors directly; not shift an orbit polynomial.
    n=p**k; N=s*n; f=[0]+[(-1)**(p-j)*comb(p,j) for j in range(1,p+1)]
    power=[1]
    for _ in range(s*n//p): power=convolution(power,f)
    numerator=[(n-1)*v for v in power]; numerator[N]+=1
    require(all(v%n==0 for v in numerator),'integer division')
    return [v//n for v in numerator]

def order(a,q=0):
    return next(i for i,v in enumerate(a) if (v%q if q else v))

def audit_arithmetic():
    rows=[]
    # Independent finite box; includes p | s and p^2 | s, k=1 and k>=2.
    for p,k in [(3,1),(3,2),(3,3),(5,1),(5,2),(7,1),(7,2),(11,1),(13,1),(17,1),(19,1)]:
        for s in sorted(set([1,2,3,p,p*p])):
            if s*p**k>800: continue
            a=shifted_from_formula(p,k,s); n=p**k; N=n*s
            tests=[]
            for q in [0,2,3,5,7,11,13,17,19,23,29,31,37,41,43]:
                expected=s*n//p if q==0 else (N if (n-1)%q==0 else ((N if s%p==0 else N-p+1) if q==p else s*n//p))
                require(order(a,q)==expected,f'depth {(p,k,s,q)}')
                tests.append([q,expected])
            h=[0]+[(-1)**(p-j)*comb(p,j)//p for j in range(1,p)]
            exact=[0]*(N+1);exact[N]=1
            for j,v in enumerate(h): exact[N-p+j]-=s*v
            require([v%p for v in a]==[v%p for v in exact],'whole polynomial characteristic p')
            if s==1: require(a[n//p]==(n-1)*p**(n//p-k),'rational leading coefficient')
            rows.append({'p':p,'k':k,'s':s,'points':N,'tests':tests})
    return rows

def adjacency():
    S={(0,1),(0,2),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2),(2,3),(2,4),(4,0)}
    xy=[divmod(i,5) for i in range(25)]
    return [[int(((b[0]-a[0])%5,(b[1]-a[1])%5) in S) for b in xy] for a in xy]

def homogeneous(M):
    return [(i,j) for i in range(len(M)) for j in range(i+1,len(M)) if all(M[k][i]==M[k][j] for k in range(len(M)) if k not in (i,j))]

def refine_partition(M):
    # Refine blocks by incoming-neighbor counts; no use of producer colors/ranks.
    partition=[(0,),tuple(range(1,len(M)))]; trace=[partition]
    while True:
        blocks=[]
        for block in partition:
            groups={}
            for v in block:
                key=tuple(sum(M[w][v] for w in B) for B in partition)
                groups.setdefault(key,[]).append(v)
            blocks += [tuple(groups[key]) for key in sorted(groups)]
        if len(blocks)==len(partition): return trace
        partition=blocks; trace.append(partition)

def is_auto(M,g):
    return all(M[i][j]==M[g[i]][g[j]] for i in range(len(M)) for j in range(len(M)))

def cycle_sizes(g):
    unseen=set(range(len(g))); sizes=[]
    while unseen:
        u=next(iter(unseen)); v=u; c=0
        while v in unseen:
            unseen.remove(v);c+=1;v=g[v]
        require(v==u,'permutation');sizes.append(c)
    return tuple(sorted(sizes))

def spectrum_from_cycles(perms):
    N=len(perms[0]); total=[0]*(N+1)
    for g in perms:
        f=[1]
        for m in cycle_sizes(g): f=convolution(f,[1]+[0]*(m-1)+[1])
        total=[a+b for a,b in zip(total,f)]
    require(all(v%len(perms)==0 for v in total),'Burnside integrality')
    c=[v//len(perms) for v in total]
    shifted=[sum(c[r]*comb(r,j)*(-1)**(r-j) for r in range(j,N+1)) for j in range(N+1)]
    return c,shifted

def audit_tournament():
    M=adjacency()
    require(all(M[i][i]==0 for i in range(25)),'no loops')
    require(all(M[i][j]+M[j][i]==1 for i in range(25) for j in range(i)),'tournament')
    require(set(map(sum,M))=={12},'regular')
    require(homogeneous(M)==[],'pair free')
    trace=refine_partition(M)
    require([len(P) for P in trace]==[2,3,15,21,25],'independent root refinement')
    permutations=[]
    for a in range(5):
        for b in range(5):
            g=[5*((i//5+a)%5)+(i%5+b)%5 for i in range(25)]
            require(is_auto(M,g),'translation preserves all arcs');permutations.append(g)
    require({g[0] for g in permutations}==set(range(25)),'transitive translations')
    T=[[M[u//2][v//2] if u//2!=v//2 else int(u%2<v%2) for v in range(50)] for u in range(50)]
    pairs=homogeneous(T)
    require(pairs==[(2*i,2*i+1) for i in range(25)],'exact pairs')
    require([sum(T[i]) for i in range(50)]==[25,24]*25,'degree layers canonically identify pair positions')
    # First degree layer induces R; pair uniqueness fixes the rest of an automorphism.
    require([[T[2*i][2*j] for j in range(25)] for i in range(25)]==M,'degree layer induces core')
    for g in permutations:
        require(is_auto(T,[2*g[i//2]+i%2 for i in range(50)]),'lifted automorphism')
    poly,shifted=spectrum_from_cycles(permutations)
    require(shifted==shifted_from_formula(5,2,1),'generic cycle census vs shifted formula')
    require(sum(poly)==1342208,'orbit count')
    require([order(shifted,q) for q in [0,5,2,3,7]]==[5,21,25,25,5],'required spectrum')
    return {'core_size':25,'tournament_size':50,'root_partition_class_counts':[len(P) for P in trace],
            'root_partition_trace':trace,'core_pair_count':0,'inflated_pair_count':len(pairs),
            'full_automorphism_count':25,'gap_orbits':1,'degrees':dict(Counter(map(sum,T))),
            'cycle_types':{str(c):v for c,v in Counter(cycle_sizes(g) for g in permutations).items()},
            'integer_orbit_coefficients':poly,'shifted_coefficients':shifted,
            'depths':{str(q):order(shifted,q) for q in [0,2,3,5,7]},
            'adjacency_sha256':hashlib.sha256(bytes(sum(M,[]))).hexdigest()}

def general_odd_controls():
    groups={}
    for n in [1,3,5,7,9,15,21,25]: groups['cyclic'+str(n)]=[[(i+a)%n for i in range(n)] for a in range(n)]
    groups['nonabelian_order21_degree7']=[[(a*i+b)%7 for i in range(7)] for a in [1,2,4] for b in range(7)]
    rows=[]
    for name,G in groups.items():
        c,a=spectrum_from_cycles(G);N=len(c)-1;d=min(len(cycle_sizes(g)) for g in G)
        require(order(a)==d,'minimum cycle formula')
        require([v%2 for v in c]==[comb(N,j)%2 for j in range(N+1)],'odd orbit mod2 law')
        for q in [2,3,5,7,11,13,17,19,23,29,31,37,41,43]:
            require((order(a,q)>d)==(a[d]%q==0),'prime support criterion')
            if q!=2: require(order(a,q)%2==N%2,'reciprocity parity')
        rows.append({'name':name,'points':N,'order':len(G),'depth_Q':d,'kappa':a[d], 'depth_2':order(a,2),'depth_5':order(a,5)})
    return rows

def main():
    arithmetic=audit_arithmetic();witness=audit_tournament();odd=general_odd_controls()
    even,even_shifted=spectrum_from_cycles([[0,1],[1,0]])
    require(even==[1,1,1] and order(even_shifted,2)==0,'even transposition hostile')
    independent=[1]
    for _ in range(3): independent=convolution(independent,shifted_from_formula(3,1,1))
    diagonal=shifted_from_formula(3,1,3)
    require(order(independent,3)==3 and order(diagonal,3)==9,'independent versus synchronized action')
    payload={'verdict':'PASS','arithmetic_case_count':len(arithmetic),'arithmetic_cases':arithmetic,'witness':witness,'general_odd_controls':odd,
             'even_transposition_polynomial':even,'even_transposition_depth_F2':order(even_shifted,2),
             'independent_versus_diagonal_three_C3_char3_depths':[order(independent,3),order(diagonal,3)]}
    print(json.dumps(payload,indent=2,sort_keys=True))
if __name__=='__main__': main()
