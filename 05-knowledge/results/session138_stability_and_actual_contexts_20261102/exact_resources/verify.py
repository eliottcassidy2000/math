#!/usr/bin/env python3
"""Standalone exact verifier for the one-edge arithmetic reflector; stdlib only."""
import hashlib,itertools,json,math
from collections import Counter
from functools import lru_cache

QMIN=(16,17,3,7,12)
QPLUS=(16,25,3,5,12)

def need(x,msg):
    if not x:raise RuntimeError(msg)

def tournament(n,code):
    A=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            u,v=(i,j) if code>>k&1 else (j,i)
            A[u]|=1<<v;k+=1
    return tuple(A)

def chain(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def reverse(A):
    m=(1<<len(A))-1
    return tuple(m^a^(1<<i) for i,a in enumerate(A))

def flip(A,i,j):
    B=list(A);B[i]^=1<<j;B[j]^=1<<i;return tuple(B)

def relabel(A,p):
    return tuple(sum(1<<j for j in range(len(p)) if A[p[i]]>>p[j]&1) for i in range(len(p)))

def substitute(Q,U):
    m=len(U);B=list(U)+[0]*(len(Q)-1)
    labs=[list(range(m))]+[[m+i-1] for i in range(1,len(Q))]
    for i in range(len(Q)):
        for j in range(len(Q)):
            if Q[i]>>j&1:
                for u in labs[i]:
                    for v in labs[j]:B[u]|=1<<v
    return tuple(B)

def join(A,B):
    n=len(A);m=len(B);tail=((1<<m)-1)<<n
    return tuple(a|tail for a in A)+tuple(b<<n for b in B)

def strong(A):
    if len(A)<2:return False
    for start in range(len(A)):
        seen=1<<start;todo=seen
        while todo:
            b=todo&-todo;todo-=b;v=b.bit_length()-1
            new=A[v]&~seen;seen|=new;todo|=new
        if seen!=(1<<len(A))-1:return False
    return True

@lru_cache(None)
def htable(A):
    n=len(A);dp=[[0]*n for _ in range(1<<n)];H=[0]*(1<<n);H[0]=1
    for s in range(1,1<<n):
        for v in range(n):
            if s>>v&1:
                t=s^(1<<v)
                dp[s][v]=1 if not t else sum(dp[t][u] for u in range(n) if t>>u&1 and A[u]>>v&1)
        H[s]=sum(dp[s])
    return tuple(H)

def H(A):return htable(A)[-1]

def literal_H(A):
    return sum(all(A[a]>>b&1 for a,b in zip(w,w[1:])) for w in itertools.permutations(range(len(A))))

def profile(A):
    n=len(A);h=htable(A);P=[[0]*(n+1) for _ in range(1<<n)];P[0][0]=1
    for s in range(1,1<<n):
        bit=s&-s;t=s
        while t:
            if t&bit:
                for d in range(1,n+1):P[s][d]+=h[t]*P[s^t][d-1]
            t=(t-1)&s
    return tuple(P[-1])

def E(A):
    h=htable(A);mask=(1<<len(A))-1
    return sum(h[s]*h[mask^s] for s in range(mask+1))

def trim(a):
    a=list(a)
    while a and not a[-1]:a.pop()
    return tuple(a)

def add(*terms):
    a=[0]*max((len(t) for t in terms),default=0)
    for t in terms:
        for i,c in enumerate(t):a[i]+=c
    return trim(a)

def scale(a,c):return trim([x*c for x in a])
def times(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)

def shift(a,s):
    out=[0]*len(a)
    for j,c in enumerate(a):
        for i in range(j+1):out[i]+=c*math.comb(j,i)*s**(j-i)
    return trim(out)

def cover_poly(p):
    out=();fall=(1,)
    for k,c in enumerate(p):
        out=add(out,scale(fall,c));fall=times(fall,(-k,1))
    return out

def pred_minus(F):
    return add(times((0,0,-2,0,1),F),times((0,1,3,2),shift(F,1)),times((0,-1,3,-2),shift(F,-1)))

def pred_plus(F):
    return add(times((0,0,0,0,1),F),times((0,2,2,2),shift(F,1)),times((0,-2,2,-2),shift(F,-1)))

def pred_delta(F):
    return add(times((0,1),add(shift(F,1),scale(shift(F,-1),-1))),times((0,0,-1),add(shift(F,1),shift(F,-1),scale(F,-2))))

def quotient_rows(A):
    """Dynamic quotient-word count, with no paths inside port0."""
    n=len(A);out=[]
    for c in range(1,(n-1)//2+2):
        D={(tuple([0]*n),-1):1}
        for t in range(n+c-1):
            N={}
            for (counts,last),v in D.items():
                for a in range(n):
                    if counts[a]>=(c if a==0 else 1) or last==a or(last>=0 and not(A[last]>>a&1)):continue
                    cs=list(counts);cs[a]+=1;key=(tuple(cs),a);N[key]=N.get(key,0)+v
            D=N
        out.append(math.factorial(c)*sum(D.values()))
    return tuple(out)

def literal_words(A,c):
    n=len(A);out=[]
    # Choose port positions and permute each distinct exterior vertex.
    for spots in itertools.combinations(range(n+c-1),c):
        spots=set(spots)
        for perm in itertools.permutations(range(1,n)):
            it=iter(perm);w=tuple(0 if i in spots else next(it) for i in range(n+c-1))
            if all(A[a]>>b&1 for a,b in zip(w,w[1:])):out.append(w)
    return tuple(out)

def canon_pair(A,B):
    return min((relabel(A,(0,)+p),relabel(B,(0,)+p)) for p in itertools.permutations(range(1,len(A))))

def is_core(A):
    n=len(A);reach=[]
    for start in range(n):
        seen=1<<start;todo=seen
        while todo:
            b=todo&-todo;todo-=b;v=b.bit_length()-1
            new=A[v]&~seen;seen|=new;todo|=new
        reach.append(seen)
    return all(any(u!=v and reach[u]>>v&1 and reach[v]>>u&1 for u in range(n)) for v in range(n))

def core_budget(B):
    best=0
    for b in range(3):
        p=5**b
        if p>B:continue
        a=0
        while p*3<=B:p*=3;a+=1
        best=max(best,3*a+4*b)
    return best

def inverse_image(N):
    # All h>1 need h<=N/23; h1 uses the empty core and positive padding.
    B=N//23;K=core_budget(max(1,B));bank={(1,1):()}
    for n in range(1,K+1):
        for code in range(1<<(n*(n-1)//2)):
            C=tournament(n,code)
            if not is_core(C):continue
            h=H(C)
            if h<=B:bank.setdefault((h,E(C)),C)
    outputs={}
    for (h,e),C in bank.items():
        r=0 if C else 1
        while 6*e*2**r-h<=N:
            U=join(C,chain(r))
            for sg,Q in ((-1,QMIN),(1,QPLUS)):
                v=6*e*2**r+sg*h
                if v<=N:
                    T=substitute(Q,U)
                    need(H(T)==v and strong(T),'inverse certificate failure')
                    outputs.setdefault(v,{'core_rows':C,'h':h,'E_core':e,'padding':r,'sign':sg,'graph_rows':T})
            r+=1
    return K,bank,outputs

def main():
    need(flip(QMIN,1,3)==QPLUS,'single edge mismatch')
    need(strong(QMIN) and strong(QPLUS),'quotients not strong')
    words={str(i):{str(c):[''.join(map(str,w)) for w in literal_words(Q,c)] for c in (1,2,3)} for i,Q in enumerate((QMIN,QPLUS))}
    need(quotient_rows(QMIN)==(11,12,0) and quotient_rows(QPLUS)==(13,12,0),'isolator response rows')
    rowbanks={};ctx=0;literal=0;pairs=[]
    for n in range(1,6):
        bank={};rows=[];graphs=[]
        for code in range(1<<(n*(n-1)//2)):
            A=tournament(n,code);r=quotient_rows(A);graphs.append(A);rows.append(r);bank[r]=bank.get(r,0)+1;ctx+=1
            rr=tuple(math.factorial(c)*len(literal_words(A,c)) for c in range(1,len(r)+1))
            need(r==rr,'independent literal quotient word mismatch');literal+=len(r)
        rowbanks[str(n)]=[{'row':r,'labelled_count':c} for r,c in sorted(bank.items())]
        if n<=4:
            for r in bank:
                for s in bank:
                    need(not(s[0]-r[0]==2 and s[1:]==r[1:] and any(r[1:])),'smaller nontrivial isolator')
        if n==5:
            for A,r in zip(graphs,rows):
                for i in range(1,n):
                    for j in range(i+1,n):
                        B=flip(A,i,j);s=quotient_rows(B)
                        if s[0]-r[0]==2 and s[1:]==r[1:] and any(r[1:]):pairs.append((A,B))
    need(len(pairs)==48,'minimal pair multiplicity')
    need(all(quotient_rows(A)==(11,12,0) and quotient_rows(B)==(13,12,0) for A,B in pairs),'minimal pair rows')
    orbit=Counter(canon_pair(A,B) for A,B in pairs)
    need(sorted(orbit.values())==[24,24],'minimal pair orbit sizes')
    keys=list(orbit);need(canon_pair(reverse(keys[0][0]),reverse(keys[0][1]))==keys[1],'reversal orbit relation')
    modules=[tournament(n,c) for n in range(1,5) for c in range(1<<(n*(n-1)//2))]
    modules += [chain(5),chain(6),QMIN,QPLUS,join(tournament(3,1),chain(2))]
    # Include an explicit directed triangle, independent of code convention.
    modules.append((2,4,1))
    checked=0;literalgraphs=0;padchecks=0
    for U in modules:
        h=H(U);e=E(U);P=profile(U);F=cover_poly(P)
        need(e==2*(P[1]+(P[2] if len(P)>2 else 0)),'ordered cover evaluation')
        polys=[]
        for Q,sgn,pred in ((QMIN,-1,pred_minus),(QPLUS,1,pred_plus)):
            T=substitute(Q,U);need(strong(T),'non-strong actual output')
            need(H(T)==6*e+sgn*h,'universal scalar formula')
            poly=cover_poly(profile(T));need(poly==pred(F),'full polynomial operator');polys.append(poly);checked+=1
            if len(T)<=7:need(literal_H(T)==H(T),'literal graph H mismatch');literalgraphs+=1
        need(add(polys[1],scale(polys[0],-1))==pred_delta(F),'difference polynomial')
        for r in (0,1,2):
            V=join(U,chain(r));need(H(V)==h and E(V)==e*2**r,'padding sidecar factorization');padchecks+=1
    certificates=[]
    for k in range(1,9):
        U=chain(k);Ts=[substitute(Q,U) for Q in (QMIN,QPLUS)]
        hs=[H(T) for T in Ts]
        need(hs==[6*2**k-1,6*2**k+1],'exact dyadic pair')
        need(all(strong(T) for T in Ts),'dyadic strongness')
        need(sum((a^b).bit_count() for a,b in zip(*Ts))==2,'more than one arc changed')
        certificates.append({'k':k,'rows_minus':Ts[0],'rows_plus':Ts[1],'H':hs,'flipped_vertices':[k,k+2]})
    U=chain(4);p=[profile(substitute(Q,U)) for Q in (QMIN,QPLUS)]
    need(p[1][1]-p[0][1]==2 and p[1][2]-p[0][2]==-22,'profile monotonicity hostile')
    need(quotient_rows(flip(QMIN,1,2))!=(13,12,0),'wrong edge accepted')
    need(pred_plus((0,1))!=add(pred_minus((0,1)),(0,0,2)),'false 2xF multiplier accepted')
    inverseK,inversebank,inverseout=inverse_image(200)
    expected=[11,13,23,25,47,49,69,75,95,97,141,147,187,191,193,197]
    need(sorted(inverseout)==expected,'exact reflector image through200')
    need(sorted(inversebank)==[(1,1),(3,12),(5,32)] and inverseK==4,'exact small input core bank')
    print(json.dumps({'status':'PASS','context_universe_orders':[1,2,3,4,5],
        'labelled_contexts':ctx,'independent_literal_word_rows':literal,'row_banks':rowbanks,
        'minimal_exterior_budget':4,'minimal_one_edge_isolators':48,'port_preserving_orbit_sizes':sorted(orbit.values()),
        'minimal_orbits_interchanged_by_global_reversal':True,'exact_context_words':words,
        'actual_substitutions_with_full_profile_checks':checked,'independent_literal_graph_counts':literalgraphs,
        'actual_join_padding_sidecar_checks':padchecks,'dyadic_graph_certificates':certificates,
        'reflector_image_through200':sorted(inverseout),'inverse_core_order_cutoff':inverseK,
        'inverse_actual_core_bank':[{'h':h,'E':e,'rows':C} for (h,e),C in sorted(inversebank.items())],
        'inverse_graph_certificates':inverseout,'hostile_L4_H_increment':2,'hostile_L4_p2_increment':-22,
        'scope':'Universal identities are proved in PROOF.md. Finite minima refer only to the declared complete context universe. No global exact-spectrum completeness claim.'},sort_keys=True,indent=2))

if __name__=='__main__':main()
