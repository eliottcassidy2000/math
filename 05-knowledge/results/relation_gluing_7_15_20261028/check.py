"""Exact finite controls for relation gluing, arithmetic projective geometry, and fusion."""
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd, isqrt, prod
import json
import sys
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

sys.stdout.reconfigure(encoding="utf-8",newline="\n")

def need(b,msg):
    if not b:raise RuntimeError(msg)

def factor(n):
    f={};p=2
    while p*p<=n:
        while n%p==0:f[p]=f.get(p,0)+1;n//=p
        p+=1
    if n>1:f[n]=1
    return f

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def sigma(n):return prod(sum(p**j for j in range(e+1)) for p,e in factor(n).items())
def sf(n):return prod(p for p,e in factor(n).items() if e%2)
def is_square(n):return isqrt(n)**2==n

def span(columns):
    s={0}
    for c in columns:s|={x^c for x in list(s)}
    return s

def brank(columns):
    basis={}
    for c in columns:
        while c:
            i=c.bit_length()-1
            if i not in basis:basis[i]=c;break
            c^=basis[i]
    return len(basis)

def rank_mod(rows,p):
    a=[[x%p for x in row] for row in rows]
    if not a:return 0
    r=0
    for j in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][j]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        u=pow(a[r][j],-1,p);a[r]=[x*u%p for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][j]:
                v=a[i][j];a[i]=[(x-v*y)%p for x,y in zip(a[i],a[r])]
        r+=1
    return r

def integer_profile(rows):
    M=Matrix(rows);D=smith_normal_form(M,domain=ZZ)
    ds=[abs(int(D[i,i])) for i in range(min(D.rows,D.cols)) if D[i,i]]
    for p in (2,3,5):
        null_mod=M.cols-rank_mod(rows,p)
        need(null_mod==M.cols-len(ds)+sum(d%p==0 for d in ds),"coefficient-change exact-sequence dimensions")
    return {"smith_nonzero":ds,"rational_nullity":M.cols-len(ds),
            "mod2_nullity":M.cols-rank_mod(rows,2)}

def subset_syndromes(columns):
    out=[0]*(1<<len(columns))
    for m in range(1,len(out)):
        b=m&-m;out[m]=out[m^b]^columns[b.bit_length()-1]
    return out

def projective_bank(r):
    coordinates=[2,3,5,7][:r];rows=[]
    for mask in range(1,1<<r):
        d=prod(p for i,p in enumerate(coordinates) if mask>>i&1)
        for a in range(1,10001):
            p=d*a*a-1
            if prime(p):break
        else:raise RuntimeError("bounded realization failed")
        need(sf(p+1)==d,"exact prime squareclass")
        rows.append({"vector":mask,"squareclass":d,"square_factor":a,"prime":p})
    need(len({row["prime"] for row in rows})==len(rows),"distinct realizing primes")
    return rows

def code_receipt(rows,r):
    n=len(rows);synd=subset_syndromes([row["vector"] for row in rows])
    words=[];circuits=[];weights=Counter()
    totals=[1]*(1<<n)
    for m in range(1<<n):
        if m:
            b=m&-m;totals[m]=totals[m^b]*(rows[b.bit_length()-1]["prime"]+1)
        need(is_square(totals[m])==(synd[m]==0),"all arithmetic totals versus binary syndrome")
        if synd[m]==0:
            words.append(m);weights[m.bit_count()]+=1
            if m and m.bit_count()<=r+1:
                sub=(m-1)&m
                while sub and synd[sub]:sub=(sub-1)&m
                if not sub:circuits.append(m)
    # Independent character-product calculation in the polynomial ring Z[z].
    coefficients=[0]*(n+1)
    for y in range(1<<r):
        poly=[1]
        for row in rows:
            sign=(-1)**((y&row["vector"]).bit_count()%2)
            poly=[(poly[j] if j<len(poly) else 0)+(sign*poly[j-1] if j else 0)
                  for j in range(len(poly)+1)]
        coefficients=[a+b for a,b in zip(coefficients,poly)]
    need(all(c%(1<<r)==0 for c in coefficients),"Fourier projection divisibility")
    need({i:c//(1<<r) for i,c in enumerate(coefficients) if c}==dict(weights),"Fourier weight enumerator")
    need(len(words)==1<<(n-r),"complete Hamming kernel size")
    return {"length":n,"rank":r,"kernel_size":len(words),"weight_enumerator":dict(sorted(weights.items())),
            "circuit_weights":dict(sorted(Counter(m.bit_count() for m in circuits).items()))},synd

def duad(v):
    mask=0
    for i in range(4):
        if v>>i&1:mask^=1|(1<<(i+1))
    if mask.bit_count()==4:mask^=63
    need(mask.bit_count()==2,"unique duad representative")
    return tuple(i for i in range(6) if mask>>i&1)

def line_spreads(lines):
    point_lines={v:[i for i,L in enumerate(lines) if v in L] for v in range(1,16)}
    answers=set()
    def visit(unused,chosen):
        if not unused:answers.add(tuple(sorted(chosen)));return
        v=min(unused)
        for i in point_lines[v]:
            if set(lines[i])<=unused:visit(unused-set(lines[i]),chosen+[i])
    visit(set(range(1,16)),[])
    return sorted(answers)

def geometry(rows):
    lines=sorted({tuple(sorted((a,b,a^b))) for a,b in combinations(range(1,16),2)})
    isotropic=[L for L in lines if not (set(duad(L[0]))&set(duad(L[1])))]
    anisotropic=[L for L in lines if L not in isotropic]
    need(len(lines)==35 and len(isotropic)==15 and len(anisotropic)==20,"35=15+20")
    for L in isotropic:need(sorted(x for v in L for x in duad(v))==list(range(6)),"syntheme line")
    for L in anisotropic:
        need(len(set(x for v in L for x in duad(v)))==3,"triangle line")
    iso_spreads=line_spreads(isotropic)
    need(len(iso_spreads)==6 and all(len(s)==5 for s in iso_spreads),"six pentads")
    need(all(len(set(a)&set(b))==1 for a,b in combinations(iso_spreads,2)),"pentads pairwise meet once")
    # A literal four-colour certificate for same-leaf hierarchy covering.
    adjacency=[{j for j,K in enumerate(isotropic) if i!=j and set(L)&set(K)}
               for i,L in enumerate(isotropic)]
    colours={}
    def colour():
        if len(colours)==15:return True
        v=max((v for v in range(15) if v not in colours),
              key=lambda v:(len({colours[w] for w in adjacency[v] if w in colours}),len(adjacency[v])))
        taken={colours[w] for w in adjacency[v] if w in colours}
        for c in range(4):
            if c not in taken:
                colours[v]=c
                if colour():return True
                del colours[v]
        return False
    need(colour(),"four-tree doily certificate")
    all_spreads=line_spreads(lines)
    need(len(all_spreads)==56,"all projective line spreads")
    spread_masks=[sum(1<<i for i in s) for s in all_spreads]
    def packing(unused,chosen):
        if not unused:return chosen
        b=unused&-unused
        for i,s in enumerate(spread_masks):
            if s&b and s&unused==s:
                found=packing(unused^s,chosen+[i])
                if found is not None:return found
        return None
    parallelism=packing((1<<35)-1,[])
    need(parallelism is not None and len(parallelism)==7,"seven spreads partition all35 lines")
    fano_lines=[L for L in lines if max(L)<=7]
    need(len(fano_lines)==7 and all(set(a)&set(b) for a,b in combinations(fano_lines,2)),
         "Fano line hierarchy cover needs seven trees")
    need(sum(L in isotropic for L in fano_lines)==3,"tetrahedral3+4 split")
    def named(L):return [rows[v-1]["prime"] for v in L]
    return {"duad_prime_labels":[{"duad":duad(row["vector"]),"prime":row["prime"],"vector":row["vector"]} for row in rows],
            "isotropic_prime_triples":[named(L) for L in isotropic],
            "fano_isotropic_triples":[named(L) for L in fano_lines if L in isotropic],
            "six_pentads":[[named(isotropic[i]) for i in s] for s in iso_spreads],
            "four_hierarchy_classes":[[list(isotropic[i]) for i,c in sorted(colours.items()) if c==k] for k in range(4)],
            "projective_spreads":len(all_spreads),
            "seven_spread_parallelism":[[list(lines[j]) for j in all_spreads[i]] for i in parallelism]}

def rank_table(columns):
    return [brank(columns[i] for i in range(len(columns)) if m>>i&1) for m in range(1<<len(columns))]

def tree_receipt(r):
    columns=list(range(1,1<<r));n=len(columns);full=(1<<n)-1;ranks=rank_table(columns)
    cut=[ranks[m]+ranks[full^m]-r for m in range(full+1)]
    def construct(indices):
        if len(indices)==1:return indices[0]
        # Split first by highest coordinate with variation; zeros/ones give linear/affine sections.
        for bit in reversed(range(r)):
            a=[v for v in indices if not (v>>bit&1)];b=[v for v in indices if v>>bit&1]
            if a and b:return [construct(a),construct(b)]
        raise RuntimeError("distinct projective points")
    tree=construct(columns);node_receipts=[]
    def visit(t):
        if isinstance(t,int):return 1<<(t-1),cut[1<<(t-1)],0
        a,wa,ga=visit(t[0]);b,wb,gb=visit(t[1]);s=a|b
        gluing=ranks[a]+ranks[b]-ranks[s]
        node_receipts.append({"subset":s,"gluing_dimension":gluing,"separator_dimension":cut[s]})
        return s,max(wa,wb,cut[s]),ga+gb+gluing
    _,width,total=visit(tree)
    need(width==r-1 and total==n-r,"projective tree construction and LCA nullity telescoping")
    need(all(cut[s]>=r-1 for s in range(full+1) if (n+2)//3<=s.bit_count()<=2*n//3),
         "all balanced cuts meet the universal lower bound")
    result={"rank":r,"leaves":n,"optimal_separator_dimension":width,"number_of_separator_values":1<<width,
            "kernel_dimension":total,"tree":tree,"node_quotients":node_receipts}
    if r==3:
        @lru_cache(None)
        def counts(s):
            if s&(s-1)==0:return Counter({cut[s]:1})
            first=s&-s;out=Counter();a=(s-1)&s
            while a:
                b=s^a
                if b and a&first:
                    for u,cu in counts(a).items():
                        for v,cv in counts(b).items():out[max(cut[s],u,v)]+=cu*cv
                a=(a-1)&s
            return out
        dist=counts(full)
        need(sum(dist.values())==10395 and min(dist)==2,"all10395 Fano leaf trees")
        result["all_rooted_binary_tree_widths"]=dict(sorted(dist.items()))
    return result

def gluing_controls():
    cuts=0
    for columns in product(range(8),repeat=4):
        global_synd=subset_syndromes(columns);kernel=[m for m,s in enumerate(global_synd) if s==0]
        for a in range(16):
            b=15^a
            left=[i for i in range(4) if a>>i&1];right=[i for i in range(4) if b>>i&1]
            sa=span(columns[i] for i in left);sb=span(columns[i] for i in right)
            image=Counter(global_synd[m&a] for m in kernel)
            ka=sum(global_synd[m]==0 for m in range(16) if m&a==m)
            kb=sum(global_synd[m]==0 for m in range(16) if m&b==m)
            need(set(image)==sa&sb and set(image.values())=={ka*kb},"exact gluing image and every fibre")
            cuts+=1
    return cuts

def fusion_controls():
    cases=0
    for p in (2,3,5,7):
        def f(e):return sigma(p**e)
        def delta(a,b):return sf(f(a+b)*f(a)*f(b))
        for a,b,c in product(range(6),repeat=3):
            need(sf(delta(a,b)*delta(a+b,c))==sf(delta(b,c)*delta(a,b+c)),"fusion cocycle identity")
            cases+=1
    for a,b in product(range(1,81),repeat=2):
        fa,fb=factor(a),factor(b)
        defect=1
        for p in fa.keys()&fb.keys():
            defect*=sigma(p**(fa[p]+fb[p]))*sigma(p**fa[p])*sigma(p**fb[p])
        need(sf(sigma(a*b))==sf(sigma(a)*sigma(b)*defect),"all-product repaired squareclass law")
    need(sf(sigma(8)*sigma(2)*sigma(4))==35,"2/4 fusion defect")
    need(sf(sigma(520))==35,"formal square receives missing35")
    # Full prime-power choice presentation, including exponents0 and repeated-prime fusion.
    choices=0
    for exps in product(range(6),repeat=4):
        n=prod(p**e for p,e in zip((2,3,5,7),exps))
        local=prod(sigma(p**e) for p,e in zip((2,3,5,7),exps))
        need(local==sigma(n),"canonical prime-power choice presentation")
        tau=prod(e+1 for e in exps);omega=sum(e>0 for e in exps)
        if n>1:
            F=tau-2;U=omega-int(omega==1 and max(exps)==1)
            S=2**omega-1-int(max(exps)<=1)
            shape=sorted(e for e in exps if e)
            need((F==U+S)==(shape in ([1],[3],[1,1,2])),"FUS evaluation in full presentation")
        choices+=1
    return {"local_cocycle_cases":cases,"all_product_cases":80**2,"prime_power_choices":choices}

def main():
    rows3=projective_bank(3);rows4=projective_bank(4)
    c3,_=code_receipt(rows3,3);c4,_=code_receipt(rows4,4)
    generators=(7,25,42,75)
    words=[]
    for v in range(16):
        c=0
        for i,g in enumerate(generators):
            if v>>i&1:c^=g
        need(subset_syndromes(list(range(1,8)))[c]==0,"relation-as-object kernel word")
        words.append(c)
    need(len(set(words))==16,"chosen Hamming kernel basis")
    sources=[prod(row["prime"] for i,row in enumerate(rows3) if c>>i&1) for c in words]
    for a,b in product(range(16),repeat=2):
        need(sources[a^b]==sources[a]*sources[b]//gcd(sources[a],sources[b])**2,
             "old arithmetic relations form the new point group")
    for a,b,c in combinations(range(1,16),3):
        need((a^b^c==0)==is_square(prod(rows4[v-1]["prime"]+1 for v in (a,b,c))),
             "relations become arithmetic points")
    need(c3["circuit_weights"]=={3:7,4:7},"Fano14 circuits")
    need(c4["circuit_weights"]=={3:35,4:105,5:168},"rank-four308 circuits")
    coordinate_primes=[2,3,5]
    A=[[factor(row["prime"]+1).get(p,0) for row in rows3] for p in coordinate_primes]
    subset_profiles={}
    for s in range(1,128):
        subset_profiles[s]=integer_profile([[row[i] for i in range(7) if s>>i&1] for row in A])
    byprime={row["prime"]:i for i,row in enumerate(rows3)}
    types={}
    for primes in ((2,5,7),(7,19,89)):
        s=sum(1<<byprime[p] for p in primes);types[str(primes)]=subset_profiles[s]
    need(types["(2, 5, 7)"]["smith_nonzero"]==[1,1],"lifted integer relation")
    need(types["(7, 19, 89)"]["smith_nonzero"]==[1,1,6],"torsion-born parity relation")
    need(subset_profiles[127]["smith_nonzero"]==[1,1,1],"full arithmetic lattice is saturated")
    smith={s:tuple(v["smith_nonzero"]) for s,v in subset_profiles.items()};smith[0]=()
    stabilizer=[];linear_permutations=0
    for basis in product(range(1,8),repeat=3):
        if brank(basis)!=3:continue
        linear_permutations+=1
        perm={}
        for v in range(1,8):
            image=0
            for i,c in enumerate(basis):
                if v>>i&1:image^=c
            perm[v]=image
        if all(smith[s]==smith[sum(1<<(perm[v]-1) for v in range(1,8) if s>>(v-1)&1)]
               for s in range(128)):
            stabilizer.append([perm[v] for v in range(1,8)])
    need(linear_permutations==168 and stabilizer==[[1,2,3,4,5,6,7],[1,2,3,5,4,7,6]],
         "integral subset modules reduce symmetry168 to2")
    G=Matrix([[1,0,-3],[0,1,-2],[0,0,-1]])
    need(G*G==Matrix.eye(3) and G.det()==-1,"integral reflection")
    perm=stabilizer[1]
    for i,j in enumerate(perm):
        need(G*Matrix(A)[:,i]==(1 if i<3 else -1)*Matrix(A)[:,j-1],
             "signed integral lift of the surviving involution")
    triangle=integer_profile([[1,1,0],[0,1,1],[1,0,1]])
    need(triangle["smith_nonzero"]==[1,1,2],"odd-cycle torsion")
    for entries in product((-1,0,1),repeat=6):
        integer_profile([list(entries[:3]),list(entries[3:])])
    # Named Eisenstein/Fano chart, inherited marking, independently checked against actual sigma.
    residues=[0,1,6,2,5,3,4];marking=[1,2,4,5,6,3,7]
    inv={x:i for i,x in enumerate(residues)};norm7=[]
    for t in range(7):
        vertices=[inv[t],inv[(t+1)%7],inv[(t+3)%7]]
        primes=[rows3[marking[v]-1]["prime"] for v in vertices]
        need(is_square(prod(p+1 for p in primes)),"norm-seven face to arithmetic square")
        norm7.append(sorted(primes))
    print(json.dumps({"prime_bank7":rows3,"prime_bank15":rows4,"code7":c3,"code15":c4,
        "geometry":geometry(rows4),"tree7":tree_receipt(3),"tree15":tree_receipt(4),
        "binary_gluing_matrix_cut_cases":gluing_controls(),
        "coefficient_change":{"integer_2_by_3_matrices":729,"arithmetic_subsets":127,
                              "line_types":types,"triangle":triangle,
                              "binary_automorphisms":linear_permutations,
                              "integer_profile_stabilizer":stabilizer,
                              "integral_reflection":[[int(x) for x in G.row(i)] for i in range(3)]},
        "fusion":fusion_controls(),"norm7_face_arithmetic_triples":norm7,
        "relation_to_object":{"basis_support_masks":generators,
            "group_law_pairs":256,"three_point_tests":455,
            "points":[{"vector":v,"old_support_mask":words[v],"old_square_total_source":sources[v],
                       "new_prime":rows4[v-1]["prime"]} for v in range(1,16)]}},
        indent=2,sort_keys=True))

if __name__=="__main__":main()
