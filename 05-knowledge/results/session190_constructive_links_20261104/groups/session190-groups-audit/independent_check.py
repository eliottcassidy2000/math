#!/usr/bin/env python3
"""Independent audit: directed arc sets, DFS paths, binary matrix enumeration.
Does not import the producer, its verifier, or its numerical result files.
"""
from itertools import combinations,permutations,product
from collections import Counter
from pathlib import Path
import json,math

def ck(b,m):
    if not b:raise RuntimeError(m)

def trans(arcs,p):return frozenset((p[u],p[v]) for u,v in arcs)
def converse(arcs):return frozenset((v,u) for u,v in arcs)
def comp(p,q):return tuple(p[q[i]] for i in range(len(p)))
def inv(p):return tuple(p.index(i) for i in range(len(p)))
def closure(gens,identity,op):
    seen={identity};todo=[identity]
    while todo:
        a=todo.pop()
        for b in gens:
            c=op(a,b)
            if c not in seen:seen.add(c);todo.append(c)
    return seen

def dfs_paths(n,A):
    found=[]
    def go(p,left):
        if not left:found.append(p);return
        for v in sorted(left):
            if not p or (p[-1],v) in A:go(p+(v,),left-{v})
    go((),set(range(n)));return found

def components(items,gens,action):
    remaining=set(items);out=[]
    while remaining:
        seed=min(remaining);seen={seed};todo=[seed]
        while todo:
            a=todo.pop()
            for g in gens:
                b=action(g,a)
                if b not in seen:seen.add(b);todo.append(b)
        ck(seen<=remaining,'component partition');remaining-=seen;out.append(seen)
    return out

def signed_edges():
    rows=[]
    for n in range(3,6):
        off=[(i,j) for i,j in combinations(range(n),2) if j!=i+1];P={(i,i+1) for i in range(n-1)}
        r=tuple(reversed(range(n)));perms=list(permutations(range(n)));identity=tuple(range(n));signed=Counter()
        bitpairs=set();signed_classes=set();fixed=0
        for word in product([0,1],repeat=len(off)):
            A=frozenset(P|{(u,v) if b else(v,u) for (u,v),b in zip(off,word)})
            B=frozenset(P|{(v,u) if b else(u,v) for (u,v),b in zip(off,word)})
            ck(A&B==P,'common path');ordinary=[];negative=[]
            for p in perms:
                ap,bp=trans(A,p),trans(B,p)
                if {ap,bp}=={A,B}:ordinary.append(p)
                if {converse(ap),converse(bp)}=={A,B}:negative.append(p);ck(converse(ap)==A,'negative swaps')
            ck(ordinary==[identity],'ordinary edge stabilizer')
            ck(negative in ([],[r]),'signed edge stabilizer')
            ck(bool(negative)==(trans(A,r)==converse(A)),'reflection criterion')
            signed[len(negative)+1]+=1
            complement=tuple(1-x for x in word);key=min(word,complement);bitpairs.add(key)
            lookup={e:b for e,b in zip(off,word)}
            reflected=tuple(lookup[(n-1-v,n-1-u)] for u,v in off)
            rkey=min(reflected,tuple(1-x for x in reflected));signed_classes.add(min(key,rkey))
            fixed+=reflected==word
        q=len(off);f=(n-1)//2;c=(q+f)//2
        ck(len(bitpairs)==2**(q-1),'edge orbits');ck(len(signed_classes)==(2**(q-1)+2**(c-1))//2,'signed edge orbits')
        rows.append({'n':n,'normalized_tournaments':2**q,'all_permutations_per_pair':math.factorial(n),'signed_stabilizers_on_words':dict(signed),'ordinary_edge_orbits':len(bitpairs),'signed_edge_orbits':len(signed_classes),'reflection_fixed_words':fixed})
    return rows

def paley():
    A=frozenset((u,v) for u in range(7) for v in range(7) if (v-u)%7 in(1,2,4));reverse=converse(A);ps=list(permutations(range(7)))
    aut=[p for p in ps if trans(A,p)==A];anti=[p for p in ps if trans(A,p)==reverse]
    ck(len(aut)==21 and len(anti)==21,'Paley full groups')
    ck(set(aut)=={tuple((a*x+b)%7 for x in range(7)) for a in(1,2,4) for b in range(7)},'Paley affine equality')
    t=tuple((x+1)%7 for x in range(7));d=tuple(2*x%7 for x in range(7));minus=tuple(-x%7 for x in range(7))
    paths=dfs_paths(7,A);ck(len(paths)==189,'Paley path DFS');E=[[0]*7 for _ in range(7)]
    for p in paths:E[p[0]][p[-1]]+=1
    ck(all(E[u][v]==(0 if u==v else 1 if(u,v)in A else 8) for u in range(7) for v in range(7)),'endpoint 1/8')
    orbits=components(paths,[t,d],lambda g,p:tuple(g[x] for x in p))
    signed_orbits=components(paths,[(t,0),(d,0),(minus,1)],lambda g,p:tuple(g[0][x] for x in(p[::-1] if g[1] else p)))
    ck(sorted(map(len,orbits))==[21]*9,'ordinary path orbits');ck(sorted(map(len,signed_orbits))==[21]*3+[42]*3,'signed path orbits')
    anti_fixed={}
    for g in anti:
        if comp(g,g)==tuple(range(7)):
            anti_fixed[str(g)]=sum(tuple(g[x] for x in p[::-1])==p for p in paths)
    ck(len(anti_fixed)==7 and set(anti_fixed.values())=={9},'anti involution fixed paths')
    S=[[0 if u==v else 1 if (u,v)in A else-1 for v in range(7)] for u in range(7)]
    ck(all(sum(S[u][w]*S[w][v] for w in range(7))==1-7*(u==v) for u in range(7) for v in range(7)),'Paley skew square')
    neigh=[]
    for p in paths:
        arcs=set(zip(p,p[1:]));B=frozenset(arcs|{(v,u) for u,v in A-arcs});neigh.append(B)
        ck(all(trans(B,g)!=B for g in aut if g!=tuple(range(7))),'Paley edge symmetry')
    ck(len(set(neigh))==189,'distinct labelled neighbors')
    return {'aut':len(aut),'signed_aut':len(aut)+len(anti),'H':len(paths),'endpoint_matrix':E,'ordinary_path_orbits':sorted(map(len,orbits)),'signed_path_orbits':sorted(map(len,signed_orbits)),'anti_involution_fixed_counts':sorted(anti_fixed.values()),'distinct_neighbors':len(set(neigh))}

def bitrank(cols):
    basis={}
    for x in cols:
        while x:
            i=x.bit_length()-1
            if i in basis:x^=basis[i]
            else:basis[i]=x;break
    return len(basis)
def mapply(M,x):
    r=0
    for i,c in enumerate(M):
        if x>>i&1:r^=c
    return r
def mmul(A,B):return tuple(mapply(A,b) for b in B)
def binary():
    I=(1,2,4);GL=[M for M in product(range(8),repeat=3) if bitrank(M)==3];ck(len(GL)==168,'GL3 count')
    C=(2,4,3);F=(1,4,6);C7=closure([C],I,mmul);G=closure([C,F],I,mmul)
    ck(len(C7)==7 and len(G)==21,'binary generated groups')
    inverses={M:next(N for N in GL if mmul(M,N)==I) for M in GL}
    ck(mmul(mmul(F,C),inverses[F])==mmul(C,C),'Frobenius conjugacy')
    normalizer={M for M in GL if {mmul(mmul(M,N),inverses[M]) for N in C7}==C7}
    ck(normalizer==G,'binary Singer normalizer')
    # Seven-dimensional permutation modules, independent product factorization and subspace check.
    t=tuple((i+1)%7 for i in range(7));d=tuple(2*i%7 for i in range(7));minus=tuple(-i%7 for i in range(7))
    def a(v,p):return sum(((v>>i)&1)<<p[i] for i in range(7))
    def poly(v,exps):
        z=0;w=v
        for k in range(max(exps)+1):
            if k in exps:z^=w
            w=a(w,t)
        return z
    kernels=[{v for v in range(128) if poly(v,e)==0} for e in[(0,1),(0,1,3),(0,2,3)]]
    ck([len(k) for k in kernels]==[2,8,8],'binary module dimensions')
    ck(len({x^y^z for x in kernels[0] for y in kernels[1] for z in kernels[2]})==128,'binary direct sum')
    ck(all(a(v,d) in k for k in kernels for v in k),'multiplier invariance')
    ck({a(v,minus) for v in kernels[1]}==kernels[2],'negative swaps summands')
    return {'GL3_order':len(GL),'C7_order':len(C7),'F21_order':len(G),'normalizer_order':len(normalizer),'C_matrix_columns':C,'F_matrix_columns':F,'module_kernel_sizes':list(map(len,kernels))}

def modular_nonsplit():
    # An invariant linear functional of a cyclic permutation has constant coefficients.
    # Explicitly enumerate coefficients in dimensions 2,3; use row reduction for 7.
    rows=[]
    for p in(2,3,7):
        M=[]
        for i in range(p):
            row=[0]*(p+1);row[i]=1;row[(i+1)%p]=-1;M.append([x%p for x in row])
        M.append([1]*p+[1]);r=0
        for c in range(p):
            k=next((k for k in range(r,len(M)) if M[k][c]),None)
            if k is None:continue
            M[r],M[k]=M[k],M[r];invv=pow(M[r][c],-1,p);M[r]=[v*invv%p for v in M[r]]
            for k in range(len(M)):
                if k!=r:
                    q=M[k][c];M[k]=[(x-q*y)%p for x,y in zip(M[k],M[r])]
            r+=1
        inconsistent=any(not any(row[:p]) and row[p] for row in M)
        ck(inconsistent,'nonsplit constant line check');rows.append({'characteristic':p,'cyclic_permutation_dimension':p,'invariant_retraction_exists':False})
    return rows

def su3_symbolic():
    # A monomial matrix sends e_i to zeta^(e_i) e_(p_i).
    def op(A,B):
        p,a=A;q,b=B
        return(tuple(p[q[i]] for i in range(3)),tuple((b[i]+a[q[i]])%7 for i in range(3)))
    I=((0,1,2),(0,0,0));D=((0,1,2),(1,2,4));P=((2,0,1),(0,0,0));G=closure([D,P],I,op)
    ck(len(G)==21,'SU3 symbolic group order')
    ck(op(op(P,D),op(P,P))==op(D,D),'P D P^-1 = D^2')
    ck(all(sum(e)%7==0 and sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))%2==0 for p,e in G),'determinants SU3')
    orbits=[sorted({k,2*k%7,4*k%7}) for k in(1,3)];ck(set(orbits[0]).isdisjoint(orbits[1]),'real conjugate orbits')
    return {'symbolic_group_order':len(G),'nontrivial_C7_character_orbits':orbits,'minimum_real_dimension_from_conjugate_orbits':6}

if __name__=='__main__':
    result={'status':'PASS','signed_edges':signed_edges(),'Paley7':paley(),'binary':binary(),'modular_nonsplit':modular_nonsplit(),'SU3':su3_symbolic()}
    Path(__file__).with_name('independent_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
