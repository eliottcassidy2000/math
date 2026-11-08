#!/usr/bin/env python3
"""Portable exact controls for the Session227 graft boundary-escape theorem.
Only the standard library is used. Every guard survives python -O.
"""
import hashlib
import itertools
import json

L = ((1,-1,3,-2),(0,0,1,-1),(2,-1,4,-1),(1,-1,1,0))
M = ((0,1,-1),(3,4,1),(1,1,1))
MI = ((-3,2,-5),(2,-1,3),(1,-1,3))
SEED = (1,0,2,1)

def need(ok, why):
    if not ok:
        raise RuntimeError(why)

def mv(a, x):
    return tuple(sum(u*v for u,v in zip(row,x)) for row in a)

def mm(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))

def det(a):
    n=len(a)
    if n==0:
        return 1
    return sum((-1)**j*a[0][j]*det(tuple(tuple(row[k] for k in range(n) if k!=j)
                for row in a[1:])) for j in range(n))

def tournaments(n):
    pairs=list(itertools.combinations(range(n),2))
    for code in range(1<<len(pairs)):
        rows=[0]*n
        for k,(i,j) in enumerate(pairs):
            if code>>k&1:
                rows[i]|=1<<j
            else:
                rows[j]|=1<<i
        yield tuple(rows)

def valid(rows, root):
    n=len(rows)
    need(n>0 and 0<=root<n, 'root domain')
    for i,row in enumerate(rows):
        need(0<=row<1<<n and not row>>i&1,'row domain')
        for j in range(i):
            need(((row>>j)&1)+((rows[j]>>i)&1)==1,'tournament orientation')

def counts(rows, root):
    """Subset endpoint recurrences, independent of all graft formulas."""
    valid(rows,root)
    n=len(rows); size=1<<n
    dp=[[0]*n for _ in range(size)]
    starts=[[0]*n for _ in range(size)]
    for v in range(n):
        dp[1<<v][v]=1
    starts[1<<root][root]=1
    for mask in range(1,size):
        for last in range(n):
            p=dp[mask][last];q=starts[mask][last]
            if not (p or q):
                continue
            avail=rows[last]&~mask
            while avail:
                bit=avail&-avail;avail-=bit;v=bit.bit_length()-1
                dp[mask|bit][v]+=p
                starts[mask|bit][v]+=q
    hp=[sum(x) for x in dp];hp[0]=1
    sp=[sum(x) for x in starts]
    ep=[x[root] for x in dp]
    full=size-1
    state=(hp[full],ep[full],sum(sp[m]*hp[full^m] for m in range(size)),
           hp[full^(1<<root)])
    return state,sp[full],hp,sp,ep

def profile(rows,root):
    state,S,*_=counts(rows,root)
    return state,S

def defect(s):
    D,E,U,d=s
    return U-D-d+E

def direct(rows,root):
    """Permutation and split counting, including the missing-join defect."""
    n=len(rows)
    def path(word):
        return all(rows[u]>>v&1 for u,v in zip(word,word[1:]))
    D=E=U=S=d=f=0
    for w in itertools.permutations(range(n)):
        if path(w):
            D+=1;E+=int(w[-1]==root);S+=int(w[0]==root)
        for cut in range(n):
            if w[cut]==root and path(w[:cut]) and path(w[cut:]):
                U+=1
    rest=tuple(v for v in range(n) if v!=root)
    for w in itertools.permutations(rest):
        if path(w):
            d+=1
        for cut in range(1,len(rest)):
            if path(w[:cut]) and path(w[cut:]) and rows[root]>>w[cut-1]&1 and rows[root]>>w[cut]&1:
                f+=1
    return (D,E,U,d),S,f

def subsets(mask):
    s=mask
    while True:
        yield s
        if not s:
            break
        s=(s-1)&mask

def stats(rows,root):
    state,S,hp,sp,ep=counts(rows,root)
    D,E,U,d=state;full=(1<<len(rows))-1;rest=full^(1<<root)
    Q=sum(hp[m]*hp[full^m] for m in range(full+1))
    delta=sum(ep[m]*hp[full^m] for m in range(full+1))
    T=sum(sp[m]*sum(hp[z]*hp[(full^m)^z] for z in subsets(full^m)) for m in range(full+1))
    b=sum(hp[z]*hp[rest^z] for z in subsets(rest))
    return D,Q,S,E,U,delta,T,d,b

def operator(z,n):
    H,Q,s,e,a,de,T,q,b=z
    return ((H,-H,Q-H,2*H-Q),(e,-e,de-e,e-de),
            (a,s-a,2*T-a,3*a-s-2*T),
            (q,-q,b-q,2*q-b) if n>=2 else (1,0,0,0))

def graft(a,root,b,r):
    n=len(a);newbits=((1<<len(b))-1)<<n
    rows=tuple(row|(newbits if v!=root else 0) for v,row in enumerate(a))
    rows+=tuple((row<<n)|(1<<root) for row in b)
    return rows,n+r

def literal_op(rows,root,c):
    n=len(rows)
    if c=='O':
        return tuple(row|(1<<n) for row in rows)+(0,),root
    new=tuple(row|((1<<n) if v!=root else 0) for v,row in enumerate(rows))+(1<<root,)
    return new,(n if c=='J' else root)

def state_op(s,c):
    D,E,U,d=s
    if c=='O':return D,0,2*U,d
    if c=='Z':return D-E+U,d,2*U-d,d
    if c=='J':return D-E+U,D-E,D+U,D
    raise RuntimeError('bad operation')

def strong(rows):
    full=(1<<len(rows))-1
    for r in range(len(rows)):
        seen=todo=1<<r
        while todo:
            bit=todo&-todo;todo-=bit;v=bit.bit_length()-1
            nxt=rows[v]&~seen;seen|=nxt;todo|=nxt
        if seen!=full:return False
    return True

def main():
    digest=hashlib.sha256()
    def record(tag,*items):
        digest.update(json.dumps([tag,*items],separators=(',',':')).encode()+b'\n')
    bank=[];universe={};transition_checks=0;defect_checks=0
    for n in range(1,6):
        total=0
        for rows in tournaments(n):
            state,S,hp,_,_=counts(rows,0)
            pstate,pS,pf=direct(rows,0)
            need((state,S)==(pstate,pS),'permutation profile mismatch')
            need(defect(state)==pf,'missing-join cover interpretation')
            outs=[u for u in range(n) if rows[0]>>u&1]
            full=(1<<n)-1
            lower=2*sum(hp[full^(1<<0)^(1<<u)^(1<<v)] for u,v in itertools.combinations(outs,2))
            need(pf>=lower>=len(outs)*(len(outs)-1),'inherited defect lower bound')
            need((pf==0)==(len(outs)<=1),'zero-defect degree equivalence')
            D,E,U,d=state
            need(0<=E<=d and 0<=D-E<=U-d,'physical cone')
            need((D-E>0)==bool(rows[0]),'nonsink path-end criterion')
            for c in 'OZJ':
                out,r=literal_op(rows,0,c)
                snew,Snew=profile(out,r)
                need(snew==state_op(state,c),'literal operation formula')
                need(snew[1]==0 if c=='O' else snew[1]==snew[3] if c=='Z' else defect(snew)==0,
                     'terminal face')
                if n>=2:need(Snew==S,'root-start preservation')
                transition_checks+=1
            record('rooted',n,rows,state,S,pf,lower)
            if n<=4:bank.append((rows,state,S,stats(rows,0)))
            total+=1;defect_checks+=1
        universe[str(n)]=total
    need(universe=={'1':1,'2':2,'3':8,'4':64,'5':1024},'source universe')
    pair_checks=escape_checks=singleton_checks=0
    for a,sa,S_a,za in bank:
        for b,sb,S_b,zb in bank:
            rows,r=graft(a,0,b,0)
            sg,Sg=profile(rows,r)
            need(sg==mv(operator(zb,len(b)),sa),'general graft formula')
            H,Q,s,e,alpha,delta,T,q,bcover=zb
            need(Sg==(alpha*S_a if len(a)>=2 else s),'root-start multiplication')
            if len(b)>=2:
                D,E,U,d=sg
                need(d-E>=sa[3]*q>0,'strict deletion-end gap')
                need((E>0)==bool(a[0]),'old nonsink criterion')
                need((defect(sg)>0)==bool(b[0]),'gadget nonsink criterion')
                escape=E>0 and d>E and defect(sg)>0
                need(escape==(bool(a[0]) and bool(b[0])),'exact escape criterion')
                escape_checks+=int(escape)
            else:
                need(sg==state_op(sa,'J') and defect(sg)==0,'singleton gadget boundary')
                singleton_checks+=1
            record('graft',a,b,sg,Sg)
            pair_checks+=1
    need(pair_checks==5625 and singleton_checks==75,'pair universe')
    source_edge=(2,0);sink_edge=(0,1)
    need(stats(source_edge,0)==(1,4,1,0,2,1,3,1,2),'edge statistics')
    need(operator(stats(source_edge,0),2)==L,'edge matrix')
    for a,state,S,z in bank:
        sinkrows,sinkroot=graft(a,0,sink_edge,0)
        zrows,zroot=literal_op(a,0,'Z');zjrows,zjroot=literal_op(zrows,zroot,'J')
        need(profile(sinkrows,sinkroot)==profile(zjrows,zjroot),'sink-rooted edge word')
        n=len(a); perm=list(range(n))+[n+1,n]
        relabelled=[0]*(n+2)
        for u,row in enumerate(sinkrows):
            relabelled[perm[u]]=sum(1<<perm[v] for v in range(n+2) if row>>v&1)
        need(tuple(relabelled)==zjrows and perm[sinkroot]==zjroot,'sink edge graph relabelling')
    need(det(L)==0,'rank determinant')
    need(det(tuple(tuple(row[j] for j in (0,1,2)) for row in L[:3]))==-1,'rank-three minor')
    need(mv(L,(-2,-1,1,1))==(0,0,0,0),'edge kernel')
    need(all(L[0][j]-2*L[1][j]-L[3][j]==0 for j in range(4)),'edge image plane')
    I=((1,0,0),(0,1,0),(0,0,1))
    need(mm(M,MI)==I and mm(MI,M)==I and det(M)==-1,'restricted inverse')
    M2=mm(M,M);M3=mm(M2,M)
    need(all(M3[i][j]-5*M2[i][j]+M[i][j]+I[i][j]==0 for i in range(3) for j in range(3)),
         'cubic identity')
    x=SEED;H=[];N=[];states=[]
    for k in range(31):
        need(x[0]==2*x[1]+x[3],'iterate image plane')
        need(mv(M,(x[1],x[2],x[3]))==(mv(L,x)[1],mv(L,x)[2],mv(L,x)[3]),'restricted action')
        if k:
            need(x[1]>0 and x[3]>x[1] and defect(x)>0,'iterated escape')
        H.append(x[0]);N.append(state_op(x,'J')[0]);states.append(x)
        x=mv(L,x)
    for k in range(28):
        need(H[k+3]==5*H[k+2]-H[k+1]-H[k],'H recurrence')
        need(N[k+3]==5*N[k+2]-N[k+1]-N[k],'closed recurrence')
    need(H[:7]==[1,5,25,119,565,2681,12721],'H initial terms')
    need(N[:7]==[3,13,61,289,1371,6505,30865],'closed initial terms')
    rows=(2,0);root=0;literal=[]
    for k in range(5):
        state,S=profile(rows,root)
        need(state==states[k] and S==2**k,'literal edge iterates')
        if k:need(strong(rows),'literal iterate strongness')
        closed,cr=literal_op(rows,root,'J');closed_state,closed_S=profile(closed,cr)
        need(closed_state[0]==N[k] and closed_S==2**k and strong(closed),'literal closure')
        literal.append({'m':k,'order':len(rows),'state':state,'root_start':S,'closed_H':closed_state[0]})
        rows,root=graft(rows,root,source_edge,0)
    collision=state_op(state_op(SEED,'O'),'J')
    need(collision==(5,1,5,1) and states[1]==(5,1,9,3),'scalar collision')
    need(collision[0]==states[1][0] and defect(collision)==0 and defect(states[1])>0,'collision scope')
    witness=(18,28,9,17,4)
    need(direct(witness,4)[0][0]==13 and strong(witness),'literal H13 witness')
    print(json.dumps({'status':'PASS','labelled_root_zero_universe':universe,
        'defect_permutation_controls':defect_checks,'literal_one_vertex_transitions':transition_checks,
        'literal_graft_pairs':pair_checks,'strict_escaping_pairs':escape_checks,
        'singleton_gadget_controls':singleton_checks,'edge_matrix':L,
        'edge_rank':3,'image_relation':'D=2E+d','restricted_inverse':MI,
        'literal_iterates':literal,'H_first_seven':H[:7],'closed_H_first_seven':N[:7],
        'algebraic_recurrence_terms':len(H),'H13_rows':witness,
        'all_case_digest_sha256':digest.hexdigest()},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
