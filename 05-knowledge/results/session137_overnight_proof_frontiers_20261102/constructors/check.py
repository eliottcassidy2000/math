#!/usr/bin/env python3
"""Portable exact consumers for actual-source geometry and right-local sealing."""
from collections import deque
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib,json
from migration import compiler,tournament,multiply

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def induced(A,keep):
    return tuple(sum(1<<j for j,v in enumerate(keep) if A[u]>>v&1) for u in keep)

def full(A,cut):
    n=len(A); return tuple(x|((not(cut>>i&1))<<n) for i,x in enumerate(A))+(cut,)

def as_source(T,x):
    keep=[i for i in range(len(T)) if i!=x]
    return induced(T,keep),sum(1<<j for j,v in enumerate(keep) if T[x]>>v&1)

@lru_cache(None)
def endpoint_table(A):
    n=len(A); rows=[[0]*n for _ in range(1<<n)]; H=[0]*(1<<n);H[0]=1
    incoming=[sum(1<<u for u in range(n) if A[u]>>v&1) for v in range(n)]
    for mask in range(1,1<<n):
        bits=mask
        while bits:
            bit=bits&-bits;bits-=bit;v=bit.bit_length()-1;rest=mask^bit
            if not rest: rows[mask][v]=1
            else:
                p=rest&incoming[v];s=0
                while p:
                    b=p&-p;p-=b;s+=rows[rest][b.bit_length()-1]
                rows[mask][v]=s
        H[mask]=sum(rows[mask])
    return tuple(H),tuple(tuple(r) for r in rows)

def HP(A):return endpoint_table(A)[0][-1]

@lru_cache(None)
def state(A,cut):
    n=len(A);H,end=endpoint_table(A);mask=(1<<n)-1
    rev=tuple(mask^A[i]^(1<<i) for i in range(n));start=endpoint_table(rev)[1]
    q=0
    for Q in range(1,mask+1):
        q+=H[mask^Q]*sum(start[Q][v] for v in range(n) if cut>>v&1)
    g=sum(end[mask][v] for v in range(n) if cut>>v&1)
    return H[-1],HP(full(A,cut)),q,g

def pos(v):
    h,F,q,g=v
    need((F-h)%2==0 and (q-F+h-g)%2==0,'positive-coordinate parity')
    return (F-h)//2,(q-F+h-g)//2,g,h-g

def native(v):
    C,L,g,t=v;return g+t,2*C+g+t,2*C+2*L+g,g

def positive_matrix(M):
    cols=[]
    for j in range(4):
        e=tuple(int(i==j) for i in range(4))
        cols.append(pos(multiply(M,native(e))))
    return tuple(tuple(cols[j][i] for j in range(4)) for i in range(4))

S0=((1,0,0,0),(0,1,1,1),(0,0,2,0),(0,0,0,0))
S1=((1,0,0,0),(0,1,0,0),(1,0,2,0),(1,0,0,0))

@lru_cache(None)
def matrix(recipe):
    if recipe=='S0':return S0
    if recipe=='S1':return S1
    n,code,z=recipe;return compiler(tournament(n,code),z)

def step(T,x,recipe):
    """Literal independent adjacency construction, keeping all original labels."""
    if isinstance(recipe,str):
        b=int(recipe[-1]); n=len(T);out=list(T)+[0]
        for i in range(n):
            if i!=x or b:out[i]|=1<<n
            else:out[n]|=1<<i
        return tuple(out),x
    n,code,z=recipe;G=tournament(n,code);N=len(T);out=list(T)+[0]*(n-1)
    labels=[x]+list(range(N,N+n-1))
    for i in range(N):
        if i!=x:
            for j in labels[1:]:out[i]|=1<<j
    for i in range(n):
        for j in range(i+1,n):
            u,v=(labels[i],labels[j]) if G[i]>>j&1 else (labels[j],labels[i])
            out[u]|=1<<v
    return tuple(out),labels[z]

def append_source(A,cut,recipe):
    T,x=step(full(A,cut),len(A),recipe)
    return as_source(T,x)

def chain(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def E(A):
    H=endpoint_table(A)[0];mask=(1<<len(A))-1
    return sum(H[s]*H[mask^s] for s in range(mask+1))

def local_coeff(G):
    """Independent native output coefficients, old port 0."""
    U=induced(G,list(range(1,len(G))));m=len(U);H,end=endpoint_table(U)
    allmask=(1<<m)-1;outcut=G[0]>>1
    rev=tuple(allmask^U[i]^(1<<i) for i in range(m));start=endpoint_table(rev)[1]
    s=sum(start[-1][v] for v in range(m) if outcut>>v&1)
    p=0
    for P in range(1,allmask+1):
        p+=sum(end[P][v] for v in range(m) if not(outcut>>v&1))*H[allmask^P]
    return H[-1],p,H[-1]-s,HP(G)-H[-1]

def build_bank():
    B={positive_matrix(S0):'S0',positive_matrix(S1):'S1'}
    for n in range(2,5):
        for code in range(1<<(n*(n-1)//2)):
            G=tournament(n,code)
            if sorted(sum(G[i]>>j&1 for j in range(1,n)) for i in range(1,n))!=list(range(n-1)):continue
            r=(n,code,n-1);B.setdefault(positive_matrix(matrix(r)),r)
    return list(B.values())

def saturation(bank,N):
    seed=(0,0,0,1);seen={seed};queue=deque([seed]);outputs={1}
    Ms=[positive_matrix(matrix(r)) for r in bank]
    while queue:
        v=queue.popleft()
        for M in Ms:
            w=tuple(min(N+1,a) for a in multiply(M,v));F=2*w[0]+w[2]+w[3]
            if F>N:continue
            outputs.add(F)
            if w not in seen:seen.add(w);queue.append(w)
    return outputs,len(seen)

def overflow_closure(bank,N):
    def pack(v):
        h,F,q,g=v
        if F>N:return None
        need(1<=h<=F and 0<=g<=h and q>=0,'nonactual fine output guard')
        return ('O',F) if q>=N else v
    seed=pack((1,1,0,0));seen={seed:None};queue=deque([seed]);outputs={1}
    while queue:
        v=queue.popleft()
        for r in bank:
            M=matrix(r)
            if v[0]=='O':
                F=v[1]
                if M[1][2]>0:continue
                if isinstance(r,str):
                    need(r=='S1','unexpected zero-cost keep-mark generator')
                    w=('O',F)
                else:
                    need(all(M[i][j]==0 for i in range(4) for j in (0,2,3)),'sealing not rank one')
                    w=pack(tuple(M[i][1]*F for i in range(4)))
            else:w=pack(multiply(M,v))
            if w is None:continue
            F=w[1] if w[0]=='O' else w[1];outputs.add(F)
            if w not in seen:seen[w]=(v,r);queue.append(w)
    return outputs,seen

def main():
    report={};axes=[]
    for n in range(1,11):
        A=chain(n);v=pos(state(A,(1<<n)-1))
        need(v==(0,2**(n-1)-1,1,0),'full-cut ray')
        if n>=2:
            first=pos(state(A,1));need(first==(2**(n-2),0,0,1),'first-cut ray')
        axes.append({'n':n,'full_cut':v,'first_cut':pos(state(A,1))})
    report['actual_axis_rays']=axes
    bases=[state(chain(1),0),state(chain(1),1),state(chain(2),1),state(chain(2),3)]
    need([pos(v) for v in bases]==[(0,0,0,1),(0,0,1,0),(1,0,0,1),(0,1,1,0)],'lattice witnesses')
    # Complete unit-height transitive cut universe through order 10.
    face=0
    for n in range(1,11):
        A=chain(n)
        for cut in range(1<<n):
            C,L,g,t=pos(state(A,cut))
            if L==g==0:
                need(t==1 and (C==0 or C&(C-1)==0),'unit-height face')
                face+=1
    need(pos(state(chain(4),3))==(5,1,0,1),'global H11 witness')
    report['face']={'all_cuts_through_chain_order':10,'face_instances':face,'fake_positive':[5,0,0,1],
                    'actual_H11':{'parent':chain(4),'cut':3,'state':state(chain(4),3)}}
    # A negative-column mutation passes the four calibration states as actual outputs,
    # but fails the next transitive ray. Positivity alone has the opposite failure.
    bad=((1,0,0,0),(0,1,0,0),(0,0,1,0),(-2,0,0,3))
    need(pos(state((2,4,1),0))==(0,0,0,3),'actual calibration image 3e_t')
    need([multiply(bad,pos(v)) for v in bases]==[(0,0,0,3),(0,0,1,0),(1,0,0,1),(0,1,1,0)],'calibration mutation')
    need(multiply(bad,(2,0,0,1))[-1]==-1,'negative-column mutation not detected')
    add4=((1,0,0,4),(0,1,0,0),(0,0,1,0),(0,0,0,1))
    need(multiply(add4,(1,0,0,1))==(5,0,0,1),'fake positive operation')
    report['hostiles']={'negative_column_passes_four_small_actual_tests':True,'next_actual_ray_rejects':True,
                        'positive_integer_parity_preserving_add4_rejected_by_face':True}
    # General local matrices: all labelled local tournaments through order 5.
    local_count=0;zero_b=0
    for n in range(2,6):
        for code in range(1<<(n*(n-1)//2)):
            r=(n,code,n-1);M=matrix(r);P=positive_matrix(M);local_count+=1
            need(min(x for row in P for x in row)>=0,'general compiler negative')
            need(tuple((P[2][j]+P[3][j])%2 for j in range(4))==(0,0,1,1),'compiler height parity')
            a,b,c,d=local_coeff(tournament(n,code))
            need(M[1]==(d,a,b,c),'native output law')
            if b==0:
                zero_b+=1;need(c==d==0,'sealing output law')
                need(all(M[i][j]==0 for i in range(4) for j in (0,2,3)),'sealing rank-one state')
    report['general_compiler']={'all_labelled_local_orders':[2,3,4,5],'matrices':local_count,'zero_cost_sealing_recipes':zero_b}
    # Every labelled source parent through order 3, every local through order 4.
    literal=0
    for n in range(1,4):
        for code in range(1<<(n*(n-1)//2)):
            A=tournament(n,code)
            for cut in range(1<<n):
                v=state(A,cut)
                for m in range(2,5):
                    for lc in range(1<<(m*(m-1)//2)):
                        r=(m,lc,m-1);B,newcut=append_source(A,cut,r)
                        need(state(B,newcut)==multiply(matrix(r),v),'literal compiler mismatch');literal+=1
    report['literal_compiler_checks']=literal
    # Collapse all words through length 3 in a genuinely mixed keep/migrate palette.
    palette=['S0','S1',(2,0,1),(3,2,2),(3,5,2)]
    collapse=0
    for ell in range(4):
        for word in product(palette,repeat=ell):
            for A,cut in [(chain(1),0),(chain(1),1),(chain(2),1),(chain(2),3)]:
                T=full(A,cut);x=len(A);root=x;v=state(A,cut)
                for r in word:T,x=step(T,x,r);v=multiply(matrix(r),v)
                need(state(*as_source(T,x))==v,'word actual state')
                need(all(T[i]>>j&1 for i in range(len(A)) for j in range(len(A)+1,len(T))),'original-parent dominance')
                if word:
                    G=induced(T,list(range(root,len(T))));a,b,c,d=local_coeff(G)
                    h,F,q,g=state(A,cut)
                    need(HP(T)==a*F+b*q+c*g+d*h,'global word collapse')
                    need((b==0)==all(G[0]>>j&1 for j in range(1,len(G))),'zero cost iff sealing')
                    if b:need(HP(T)>q,'productive cost gate')
                    else:need(HP(T)==a*F,'sealing product')
                collapse+=1
    report['collapse_word_source_checks']=collapse
    # The H-neutral padding hostile, including independent actual successful graph.
    A=chain(2);cut=1;v=state(A,cut)
    B,bc=append_source(A,cut,'S0');P,pc=append_source(A,cut,'S1')
    need(v==(1,3,2,0) and state(B,bc)[1]==5 and state(P,pc)==(1,3,5,1),'padding hostile')
    report['padding_hostile']={'source':{'parent':A,'cut':cut,'state':v},
       'successful_S0':{'parent':B,'cut':bc,'state':state(B,bc)},
       'neutral_S1':{'parent':P,'cut':pc,'state':state(P,pc)},'target_5_rejected_all_right_local_words':True}
    # Independent finite quotients, both all word lengths; exact scalar languages agree.
    bank=build_bank();need(len(bank)==32,'bank changed')
    N=173;out1,n1=saturation(bank,N);out2,seen=overflow_closure(bank,N)
    need(out1==out2,'overflow/saturation language mismatch')
    missing=[h for h in range(1,N+1,2) if h not in out2]
    need(missing==[7,21],'inherited positive language changed')
    report['all_length_quotient']={'cap':N,'bank_size':len(bank),'positive_saturation_states':n1,
       'source_cost_states':len(seen),'overflow_tags':sum(v[0]=='O' for v in seen),'missing_odds':missing,
       'support_sha256':hashlib.sha256(bytes(int(i in out2) for i in range(N+1))).hexdigest()}
    report['status']='PASS'
    report['scope']='Structural universal proofs use actual witnesses and inherited compiler; finite controls are not universal spectrum completeness.'
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__':main()
