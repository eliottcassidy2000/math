#!/usr/bin/env python3
"""Independent exact controls; no repository imports or optimization-sensitive gates."""
from fractions import Fraction as F
from itertools import product, combinations, permutations
from math import gcd, lcm, prod, factorial
import hashlib, json, sys

COUNTS={}
def need(ok,label):
    if not ok: raise RuntimeError(label)
def sig(ds,x): return tuple((x[j]-x[i])%gcd(ds[i],ds[j]) for i,j in combinations(range(len(ds)),2))
def canonical(ds,x):
    return min(tuple((a+t)%d for a,d in zip(x,ds)) for t in range(lcm(*ds)))
def triangular(ds,s):
    edge=dict(zip(combinations(range(len(ds)),2),s))
    return all((edge[i,j]+edge[j,k]-edge[i,k])%gcd(ds[i],ds[j],ds[k])==0 for i,j,k in combinations(range(len(ds)),3))
def crt_pair(a,m,b,n):
    need(gcd(m,n)==1,'coprime CRT input')
    return (a+m*((b-a)*pow(m,-1,n)%n))%(m*n),m*n

def reconstruct(ds,s):
    e=dict(zip(combinations(range(len(ds)),2),s))
    def diff(i,j):
        return 0 if i==j else (e[i,j] if i<j else -e[j,i])
    rem=lcm(*ds); primes=[]; p=2
    while p*p<=rem:
        if rem%p==0:
            primes.append(p)
            while rem%p==0: rem//=p
        p+=1
    if rem>1: primes.append(rem)
    xs=[0]*len(ds); mods=[1]*len(ds)
    for p in primes:
        pp=[]
        for d in ds:
            z=1
            while d%p==0: z*=p; d//=p
            pp.append(z)
        anchor=max(range(len(ds)),key=lambda i:pp[i])
        for i,n in enumerate(pp):
            if n>1: xs[i],mods[i]=crt_pair(xs[i],mods[i],diff(anchor,i)%n,n)
    need(tuple(mods)==tuple(ds),'CRT reconstruction moduli')
    need(sig(ds,xs)==s,'triangle-compatible signature reconstructs')
    return tuple(xs)

def phase_tests():
    profiles=[ds for r in range(1,5) for ds in product(range(1,5),repeat=r)]
    profiles += [(4,6,8),(3,5,15),(4,6,10,15),(6,10,15),(3,9,15,45),(15,21,35),(3,9,9)]
    states=signatures=assignments=edge_graphs=0
    for ds in profiles:
        bysig={}; counts={}; L=lcm(*ds)
        # Each first encountered state independently generates its whole translation orbit.
        orbit_seen={}
        for x in product(*(range(d) for d in ds)):
            states+=1; s=sig(ds,x)
            if s not in bysig:
                bysig[s]=x
                orb={tuple((a+t)%d for a,d in zip(x,ds)) for t in range(L)}
                need(len(orb)==L,'cyclic orbit length')
                for y in orb:
                    need(y not in orbit_seen,'distinct signatures imply disjoint orbit')
                    orbit_seen[y]=s
            need(orbit_seen.get(x)==s,'pair phases determine orbit')
            need(triangular(ds,s),'native triangle compatibility')
            counts[s]=counts.get(s,0)+1
        need(len(bysig)==prod(ds)//L,'all-arity phase count')
        need(all(n==L for n in counts.values()),'all phase fibres exact orbit')
        for s in bysig: reconstruct(ds,s)
        signatures+=len(bysig)
        edge_pairs=list(combinations(range(len(ds)),2))
        # Independent signature collision test for EVERY edge subgraph.
        for mask in range(1<<len(edge_pairs)):
            E=[ij for k,ij in enumerate(edge_pairs) if mask>>k&1]
            projected={tuple(s[k] for k in range(len(edge_pairs)) if mask>>k&1) for s in bysig}
            connected=True
            for p in range(2,L+1):
                if any(p%q==0 for q in range(2,__import__('math').isqrt(p)+1)):continue
                pa=p
                while pa<=L:
                    I={i for i,d in enumerate(ds) if d%pa==0}
                    if I:
                        seen={min(I)}
                        while True:
                            nxt=seen|{j for i,j in E if i in seen and j in I}|{i for i,j in E if j in seen and i in I}
                            if nxt==seen:break
                            seen=nxt
                        if seen!=I:connected=False
                    pa*=p
            need((len(projected)==len(bysig))==connected,'prime-layer graph completeness iff')
            edge_graphs+=1
        ranges=[range(gcd(ds[i],ds[j])) for i,j in combinations(range(len(ds)),2)]
        if prod(map(len,ranges))<=50000:
            for s in product(*ranges):
                assignments+=1
                need(triangular(ds,s)==(s in bysig),'triangle image is exact')
    # A complete graph is required: only anchor pair phases fail for (2,4,4).
    ds=(2,4,4); x=(0,0,0); y=(0,0,2)
    need(sig(ds,x)[:2]==sig(ds,y)[:2] and canonical(ds,x)!=canonical(ds,y),'star-only phase hostile')
    # Multidimensional abelian action: pairwise orbit data do not suffice.
    G=list(product(range(3),repeat=2))
    def orbit(x): return {((x[0]+u)%3,(x[1]+v)%3,(x[2]+u+v)%3) for u,v in G}
    x=(0,0,0); y=(0,0,1)
    need(y not in orbit(x),'F3^2 triple invariant hostile')
    for i,j in combinations(range(3),2):
        need({(z[i],z[j]) for z in orbit(x)}==set(product(range(3),repeat=2)),'every pair orbit full')
    # A nonabelian boundary as an independent control, NOT a tournament automorphism claim.
    A4=[p for p in permutations(range(4)) if sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))%2==0]
    x=(0,1,2); y=(0,2,1)
    need(y not in {tuple(p[a] for a in x) for p in A4},'A4 triple orientation')
    for i,j in combinations(range(3),2): need(any((p[x[i]],p[x[j]])==(y[i],y[j]) for p in A4),'A4 pair orbit')
    COUNTS['cyclic_phase']={'profiles':len(profiles),'states':states,'signatures':signatures,'candidate_edge_assignments':assignments,'edge_subgraphs':edge_graphs,'all_pass':True}

def field_mul(a,b):
    out=0
    while b:
        if b&1: out^=a
        b>>=1; a<<=1
        if a&8: a^=11
    return out

def fpow(a,n):
    out=1
    while n:
        if n&1:out=field_mul(out,a)
        a=field_mul(a,a);n//=2
    return out

def rank_mod(rows,p):
    a=[list(x) for x in rows]; r=0
    for j in range(len(a[0]) if a else 0):
        k=next((k for k in range(r,len(a)) if a[k][j]%p),None)
        if k is None: continue
        a[r],a[k]=a[k],a[r]; inv=pow(a[r][j]%p,-1,p); a[r]=[(z*inv)%p for z in a[r]]
        for k in range(len(a)):
            if k!=r:
                v=a[k][j]%p;a[k]=[(z-v*w)%p for z,w in zip(a[k],a[r])]
        r+=1
    return r

def finite_field_tests():
    basis=[3,5,7]
    need(rank_mod([[(b>>j)&1 for b in basis] for j in range(3)],2)==3,'F8 normal basis')
    for w in product(range(2),repeat=3):
        x=0
        for a,b in zip(w,basis):
            if a:x^=b
        shifted=(w[-1],)+w[:-1]; y=0
        for a,b in zip(shifted,basis):
            if a:y^=b
        need(field_mul(x,x)==y,'rotation Frobenius intertwiner')
    # Tensor F8^tensor3 -> F8^9: evaluate independent embeddings.
    triples=list(product(range(3),repeat=3)); phases=[(0,i,j) for i in range(3) for j in range(3)]
    def evaluate(t,phase):
        out=1
        for exponent,r in zip(t,phase):out=field_mul(out,fpow(2,exponent*(2**r)))
        return out
    cols=[[evaluate(t,p) for p in phases] for t in triples]
    matrix=[[(col[k]>>bit)&1 for col in cols] for k in range(9) for bit in range(3)]
    need(rank_mod(matrix,2)==27,'tensor evaluation isomorphism rank27')
    tests=0
    for i,t in enumerate(triples):
        for j,u in enumerate(triples):
            coeffs=[fpow(2,a+b) for a,b in zip(t,u)]
            terms=[v for v in triples if all((c>>a)&1 for c,a in zip(coeffs,v))]
            for k,p in enumerate(phases):
                right=0
                for v in terms:right^=evaluate(v,p)
                need(field_mul(cols[i][k],cols[j][k])==right,'tensor evaluation multiplicativity')
            tests+=1
    COUNTS['finite_field']={'normal_basis_words':8,'tensor_dimension':27,'tensor_factors':9,'basis_multiplication_pairs':tests}

def seed(bit):
    a=[[False]*4 for _ in range(4)]
    for i,j in [(0,1),(0,2),(0,3),(1,2),(2,3),(3,1)]:a[i][j]=True
    return a if bit==0 else [list(x) for x in zip(*a)]
def join(parts):
    n=sum(map(len,parts)); a=[[False]*n for _ in range(n)]; off=0
    for part in parts:
        m=len(part)
        for i in range(m):
            for j in range(m):a[off+i][off+j]=part[i][j]
            for j in range(off+m,n):a[off+i][j]=True
        off+=m
    return a

def scc(a):
    n=len(a); reach=[sum((1<<j) for j in range(n) if a[i][j]) | (1<<i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    remaining=set(range(n)); comps=[]
    while remaining:
        i=min(remaining); c={j for j in remaining if reach[i]>>j&1 and reach[j]>>i&1};comps.append(c);remaining-=c
    saved=comps[:]
    comps.sort(key=lambda c:-sum(a[next(iter(c))][next(iter(d))] for d in saved if d!=c))
    return comps

def carousel(rows):
    r=len(rows);m=len(rows[0]); parts=[join([seed(rows[j][i]) for j in range(r)]) for i in range(m)];b=4*r;n=m*b
    a=[[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            bi,li=divmod(i,b);bj,lj=divmod(j,b)
            a[i][j]=parts[bi][li][lj] if bi==bj else 1<=((bj-bi)%m)<=(m-1)//2
    return a

def module_closure(a,u,v):
    n=len(a); mask={u,v}; changed=True
    while changed:
        changed=False
        for x in range(n):
            if x in mask:continue
            values={a[x][y] for y in mask}
            if len(values)>1:mask.add(x);changed=True
    return mask

def aut_count(a):
    n=len(a); deg=[sum(row) for row in a];can=[[j for j in range(n) if deg[j]==deg[i]] for i in range(n)]
    order=sorted(range(n),key=lambda i:len(can[i])); assignment={}; used=set()
    def visit(k):
        if k==n:return 1
        i=order[k];total=0
        for j in can[i]:
            if j in used or any(a[i][u]!=a[j][v] for u,v in assignment.items()):continue
            assignment[i]=j;used.add(j);total+=visit(k+1);used.remove(j);del assignment[i]
        return total
    return visit(0)

def tournament_tests():
    colors=0
    for r in range(1,6):
        signatures={}
        for e in product(range(2),repeat=r):
            a=join([seed(b) for b in e]); seq=tuple(map(len,scc(a)))
            need(seq==tuple(z for bit in e for z in ((1,3) if bit==0 else (3,1))),'SCC colour decoder')
            need(seq not in signatures,'distinct bit colours');signatures[seq]=e;colors+=1
    need(aut_count(seed(0))==aut_count(seed(1))==3,'seed automorphism counts')
    matrices=closures=autos=0
    examples=[]
    for r,m in [(1,3),(2,3),(1,5)]:
        for flat in product(range(2),repeat=r*m):
            rows=tuple(tuple(flat[j*m:(j+1)*m]) for j in range(r));a=carousel(rows);b=4*r;n=len(a);matrices+=1
            # Minimal-module closure is an adjacency-only recovery of the outer blocks.
            for i,j in combinations(range(n),2):
                c=module_closure(a,i,j);need((len(c)<n)==(i//b==j//b),'outer block recovery');closures+=1
            periods=[next(d for d in range(1,m+1) if m%d==0 and all(row[i]==row[(i+d)%m] for i in range(m))) for row in rows]
            expected=3**(r*m)*(m//lcm(*periods))
            rotations=sum(all(row[i]==row[(i+t)%m] for row in rows for i in range(m)) for t in range(m))
            need(expected==3**(r*m)*rotations,'joint period automorphism formula')
            if r==1 and m==3 or r==2 and m==3 and flat in [(0,)*6,(1,0,0,0,1,0),(1,0,0,1,0,0),(1,0,1,0,1,0)]:
                direct=aut_count(a);need(direct==expected,'independent adjacency aut enumeration');autos+=1
                examples.append({'rows':rows,'order':n,'aut':direct})
    COUNTS['native_tournament']={'bit_colors':colors,'outer_tournaments':matrices,'pair_module_closures':closures,'direct_automorphism_censuses':autos,'examples':examples}

def affine(bits,base_odd=False):
    a=1;c=0;b=1
    for bit in bits:
        if base_odd:
            c=3*c+b;a*=3;b*=2**bit
        else:
            c=3*c+bit*b;a*=3;b*=2
    return a,c,b

def v2(n):
    need(n!=0,'finite valuation input');n=abs(n);return (n&-n).bit_length()-1

def arithmetic_tests():
    words=[w for n in range(5) for w in product(range(2),repeat=n)];tests=0
    for u in words:
        au,cu,bu=affine(u)
        for v in words:
            av,cv,bv=affine(v)
            need(affine(u+v)==(av*au,av*cu+bu*cv,bv*bu),'carry cocycle');tests+=1
    u=(1,7);v=(7,1);a,c,b=affine(u,True);A,C,B=affine(v,True)
    need((a,b)==(A,B) and c!=C,'trace blind to carry')
    need(F(a*483+c,b)==F(A*469+C,B)==17,'actual odd Collatz merge')
    # Corrected complete matching branch.
    word=(1,0,0,0);a,c,b=affine(word)
    need(a*5+c==432 and a*21+c==1728,'rank full-match witness')
    need(len(word)-v2(a*5+c)==0 and len(word)-v2(a*21+c)==-2,'residual valuation sidecar')
    z=7;q=1;prefix=[]
    for i in range(110):
        if i==86:z86=z;q86=q
        bit=z%2;need(not(q==1 and bit==0),'root7 legal prefix');prefix.append(bit);z=(3*z+bit)//2;q+=3*bit-1
    need(z86==10952205219944621 and q86==47,'root7 macro86')
    P=(1,0,0);Q=(1,0,1,0,0,0);a,c,b=affine(P);A,C,B=affine(Q);D=a-b;E=A-B;delta=C*D-c*E
    y=(a*z86+c)//b
    need(delta==684 and v2(D*z86+c)==5 and v2(D*y+c)==2 and v2(E*y+C)==9,'seven bit refill')
    need(b*D*(E*y+C)==a*E*(D*z86+c)+b*delta,'switch cocycle')
    COUNTS['arithmetic']={'binary_word_pairs':tests,'legal_root7_macros':110,'refill':{'old_fuel':5,'after_P':2,'after_switch':9,'delta':delta},'Collatz_merge':17}

def partial(word):
    # A word of reflections x->L-x: exact open-interval domain.
    e=1;c=F(0);lo=None;hi=None
    for L in word:
        bounds=sorted([F(-c,e),F(L-c,e)])
        lo=bounds[0] if lo is None else max(lo,bounds[0]);hi=bounds[1] if hi is None else min(hi,bounds[1])
        e=-e;c=L-c
    return e,c,lo,hi

def geometry_tests():
    nwords=probes=0
    for n in range(7):
        for w in product(range(1,5),repeat=n):
            e,c,lo,hi=partial(w);nwords+=1
            for x in (F(i,2) for i in range(-2,11)):
                actual=x;legal=True
                for L in w:
                    if not (0<actual<L):legal=False;break
                    actual=L-actual
                predicted=(lo is None or lo<x) and (hi is None or x<hi)
                need(legal==predicted,'partial affine exact domain')
                if legal:need(actual==e*x+c,'partial affine exact value')
                probes+=1
    need(partial((2,2))[:2]==partial((3,3))[:2] and partial((2,2))[2:]!=partial((3,3))[2:],'equal formal identity unequal guard')
    def positive_scale_comp(a,upper,b,other_upper):
        return a*b,min(F(upper),F(other_upper,a))
    gf=positive_scale_comp(2,1,3,2);fg=positive_scale_comp(3,2,2,1)
    need(gf==(6,F(1)) and fg==(6,F(1,3)),'same formal centre unequal partial composition domains')
    coords=[(i,j) for i in range(-5,6) for j in range(-5,6)];pairs=0
    def phi(u):m,n=u;return(m*m-n*n,2*m*n,m*m+n*n)
    def lor(u,v):return u[2]*v[2]-u[0]*v[0]-u[1]*v[1]
    for u in coords:
        for v in coords:
            need(lor(phi(u),phi(v))==2*(u[0]*v[1]-u[1]*v[0])**2,'Gaussian Lorentz determinant');pairs+=1
    need(phi((8,1))[2]==phi((7,4))[2]==65,'norm65 equal norm')
    need(lor(phi((8,1)),phi((1,0)))==2 and lor(phi((7,4)),phi((1,0)))==32,'norm65 distinct adjacency')
    origins=[(0,0),(1,0),(1,1),(0,1),(0,0)]
    need(all(abs(x-X)+abs(y-Y)==1 for (x,y),(X,Y) in zip(origins,origins[1:])),'local square seams')
    need(len(set(origins))==4,'flat tree nonoverlap hostile')
    # Named involutions: product alone loses component connectivity.
    need(tuple((0,1)[(0,1)[i]] for i in range(2))==tuple((1,0)[(1,0)[i]] for i in range(2)),'involution products identical')
    def matching_components(a,b):
        edges=[((0,min(i,a[i])),(1,min(i,b[i]))) for i in range(2)]
        unseen={v for edge in edges for v in edge}; count=0
        while unseen:
            seen={min(unseen)}
            while True:
                nxt=seen|{v for u,v in edges if u in seen}|{u for u,v in edges if v in seen}
                if nxt==seen:break
                seen=nxt
            unseen-=seen;count+=1
        return count
    need(matching_components((0,1),(0,1))==2 and matching_components((1,0),(1,0))==1,'matching product connectivity hostile')
    COUNTS['geometry']={'formal_center_partial_domain_hostile':True,'matching_product_hostile':{'identity_pair_components':2,'swap_pair_components':1},'split_words':nwords,'rational_probes':probes,'Gaussian_vector_pairs':pairs,'flat_square_union_area':4,'piece_area_sum':5}

def det_bareiss(matrix):
    a=[row[:] for row in matrix];n=len(a);prev=1;sign=1
    for k in range(n-1):
        row=next((i for i in range(k,n) if a[i][k]),None)
        if row is None:return 0
        if row!=k:a[k],a[row]=a[row],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                val=a[i][j]*pivot-a[i][k]*a[k][j];need(val%prev==0,'Bareiss integral division');a[i][j]=val//prev
            a[i][k]=0
        prev=pivot
    return sign*a[-1][-1]

def moment_tests():
    for d in range(11):
        H=[[factorial(i+j) for j in range(d+1)] for i in range(d+1)]
        need(det_bareiss(H)==prod(factorial(j)**2 for j in range(d+1)),'factorial Hankel determinant')
    primes=[2,3,5,7,11,13]
    for p in primes:
        H=[[factorial(i+j)%p for j in range(p)] for i in range(p)]
        Hp=[[factorial(i+j)%p for j in range(p+1)] for i in range(p+1)]
        need(rank_mod(H,p)==p and rank_mod(Hp,p)==p,'char-p sharp carrier')
    need(factorial(1)-factorial(14)==-87178291199,'cyclic13 alias')
    need(6*factorial(1)-factorial(3)==0,'nonideal scalar kernel first moment')
    need(36*factorial(2)-12*factorial(4)+factorial(6)==504,'kernel not multiplicatively stable')
    # Exact law-level hostile, not a deterministic orbit obstruction.
    for k in range(2,8):
        laws=[[b for b in product(range(2),repeat=k) if sum(b)%2==e] for e in range(2)]
        for mask in range(2**k-1):
            vals=[F(sum(all(b[j] for j in range(k) if mask>>j&1) for b in law),len(law)) for law in laws]
            need(vals[0]==vals[1],'proper subset moments agree')
        vals=[F(sum(all(b) for b in law),len(law)) for law in laws]
        need(abs(vals[0]-vals[1])==F(1,2**(k-1)),'full conjunction separates law')
    COUNTS['moments']={'factorial_hankel_degrees':list(range(11)),'characteristics':primes,'cyclic_alias':-87178291199,'parity_law_dimensions':list(range(2,8))}

def verify_sources():
    from pathlib import Path
    root=Path(__file__).resolve().parent
    path=root/'SOURCE_PINS.json'
    need(path.exists(),'required source pin ledger exists')
    data=json.loads(path.read_text())
    for row in data['files']:
        blob=(root/row['snapshot']).read_bytes()
        need(hashlib.sha256(blob).hexdigest()==row['sha256'],'source pin '+row['path'])
        need(hashlib.sha1(b'blob '+str(len(blob)).encode()+b'\0'+blob).hexdigest()==row['git_blob'],'Git blob '+row['path'])
    COUNTS['source_files_verified']=len(data['files'])

if __name__=='__main__':
    if '--inject-failure' in sys.argv:need(False,'injected explicit failure')
    verify_sources();phase_tests();finite_field_tests();tournament_tests();arithmetic_tests();geometry_tests();moment_tests()
    print(json.dumps({'status':'PASS','checks':COUNTS},sort_keys=True,indent=2))
