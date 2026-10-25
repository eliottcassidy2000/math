"""Independent exact Session63 audit; no producer imports; Python stdlib only."""
from itertools import permutations, product
from math import comb, gcd
import json

checks = 0
def check(test, label):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(label)

def walsh(r):
    return [[(-1)**((i & j).bit_count()) for j in range(1 << r)] for i in range(1 << r)]

def transform(H, v):
    return [sum(a*b for a,b in zip(row,v)) for row in H]

def tensor(A,B):
    return [[x*y for x in ar for y in br] for ar in A for br in B]

def multiply(A,B):
    return [[sum(x*y for x,y in zip(ar,bc)) for bc in zip(*B)] for ar in A]

def main():
    # Integral jets from Burnside coefficients, not producer output.
    a=[]
    for i in range(26):
        numerator=comb(25,i)+(24*comb(5,i//5) if i%5==0 else 0)
        check(numerator%25==0,'integer orbit coefficient')
        a.append(numerator//25)
    b=[sum((-1)**(i-j)*comb(i,j)*a[i] for i in range(j,26)) for j in range(26)]
    delta=[]; d=0
    for x in b:
        d=gcd(d,x); delta.append(d)
    check(delta==[0]*5+[3000]*8+[600]*4+[120]*4+[24]*4+[1], 'jet gcd chain')
    # All signed involutions of ranks 1..3, all coefficient rings Z/m.
    signed_cases=0
    for rank in range(1,4):
        for p in permutations(range(rank)):
            for s in product((-1,1),repeat=rank):
                if any(p[p[i]]!=i or s[i]*s[p[i]]!=1 for i in range(rank)):
                    continue
                seen=set(); good=bad=0
                for i in range(rank):
                    if i in seen: continue
                    seen.update((i,p[i]))
                    if p[i]==i and s[i]==-1: bad+=1
                    else: good+=1
                for mod in (2,3,4,6):
                    count=0
                    for v in product(range(mod),repeat=rank):
                        count+=all(v[p[i]]==s[i]*v[i]%mod for i in range(rank))
                    check(count==mod**good*gcd(2,mod)**bad,'signed orbit basechange')
                    signed_cases+=1
    # All graphs and tournaments n=2..4, endpoint-reversal stabilizer signs.
    graph_results=[]
    for n in range(2,5):
        edges=list((i,j) for i in range(n) for j in range(i+1,n)); edge_index={e:i for i,e in enumerate(edges)}
        perms=list(permutations(range(n))); maps=[]
        for p in perms:
            mp=[]; rev=0
            for k,(i,j) in enumerate(edges):
                mp.append(edge_index[tuple(sorted((p[i],p[j])))])
                if p[i]>p[j]: rev|=1<<k
            maps.append((mp,rev))
        def relabel(v, mp):
            return sum(1<<mp[k] for k in range(len(edges)) if v>>k&1)
        unseen=set(range(1<<len(edges))); good=bad=0
        while unseen:
            v=min(unseen); orbit=set(); odd=False
            for mp,rev in maps:
                w=relabel(v,mp); orbit.add(w)
                if w==v and (v&rev).bit_count()%2: odd=True
            unseen-=orbit
            if odd: bad+=1
            else: good+=1
        unseen=set(range(1<<len(edges))); tournaments=0
        while unseen:
            v=min(unseen)
            orbit={relabel(v^rev,mp) for mp,rev in maps}
            unseen-=orbit; tournaments+=1
        check(good==tournaments,'even graphs equal tournament types')
        graph_results.append({'n':n,'trivial_sign_orbits':good,'nontrivial_sign_orbits':bad,'point_lattice_invariant_rank':tournaments,'Fourier_mod2_invariant_rank':good+bad})
    # Walsh inversion, positivity criterion, and normalized convolution.
    fourier_cases=0
    snf_controls=[]
    for r in range(0,5):
        H=walsh(r); N=1<<r
        U=V=[[1]]
        for _ in range(r):
            U=tensor(U,[[1,0],[1,-1]])
            V=tensor(V,[[1,-1],[0,1]])
        D=multiply(multiply(U,H),V)
        check(D==[[2**i.bit_count() if i==j else 0 for j in range(N)] for i in range(N)],'explicit unimodular Walsh diagonalization')
        factors=sorted(D[i][i] for i in range(N))
        check(factors==sum(([2**j]*comb(r,j) for j in range(r+1)),[]),'Walsh Smith multiplicities')
        snf_controls.append({'rank':r,'factors':factors})
        for j in range(N):
            unit=[int(k==j) for k in range(N)]
            check(transform(H,transform(H,unit))==[N*x for x in unit],'Walsh orthogonality')
        for trial in range(8):
            f=[(j*j+trial*j+trial)%5 for j in range(N)]
            g=[(3*j+2*trial+1)%4 for j in range(N)]
            A=transform(H,f); B=transform(H,g)
            check(transform(H,A)==[N*x for x in f], 'histogram integrality and positivity')
            convolution=[sum(A[i]*B[i^j] for i in range(N)) for j in range(N)]
            check(convolution==[N*x for x in transform(H,[x*y for x,y in zip(f,g)])],'normalized convolution')
            # C_{ij}=A_{i xor j}; character eigenvalue is N f_s.
            for s in range(N):
                v=H[s]
                check([sum(A[i^j]*v[j] for j in range(N)) for i in range(N)]==[N*f[s]*z for z in v], 'PSD eigenvalue identity')
            fourier_cases+=1
        if r:
            f=[1]*N; h=[1+H[1][j] for j in range(N)]
            A=transform(H,f); B=transform(H,h)
            check(all(A[j]==B[j] for j in range(N) if j!=1) and A[1]!=B[1], 'missing mode not context complete')
    # PSD does not imply integral histogram; lattice does not imply positivity.
    H=walsh(1)
    check(transform(H,[1,0])==[1,1], 'PSD fractional histogram hostile')
    check(transform(H,transform(H,[-1,2]))==[-2,4], 'integral negative histogram hostile')
    # Legendre branch-free twists by literal y enumeration.
    legendre=[]
    for p in (5,7):
        nonsquare=next(d for d in range(1,p) if pow(d,(p-1)//2,p)==p-1)
        for lam in (2,3):
            counts=[]
            for d in (1,nonsquare):
                count=sum(y*y%p==d*x*(x-1)*(x-lam)%p for x in range(p) if x not in (0,1,lam) for y in range(p))
                counts.append(count)
            expected={(5,2):[4,0],(5,3):[0,4],(7,2):[4,4],(7,3):[0,8]}[p,lam]
            check(counts==expected,'Legendre branch-free twist counts')
            legendre.append({'p':p,'lambda':lam,'counts':counts,'mean':sum(counts)//2,'variance':(counts[0]-counts[1])**2//4})
    return {'status':'PASS','checks':checks,'signed_cases':signed_cases,'fourier_cases':fourier_cases,'Smith_controls':snf_controls,'jet_coefficients':b,'gcd_chain':delta,'graph_controls':graph_results,'Legendre_controls':legendre}

if __name__=='__main__':
    print(json.dumps(main(),sort_keys=True,indent=2))
