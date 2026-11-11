#!/usr/bin/env python3
"""Exact NEW controller tests and executable positive-word constructors.
No inherited census is rechecked. Graphs are literal, modular tests labelled.
"""
from collections import deque
from itertools import product
from math import gcd
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
S0=(-1,3,0,1) # (r,M,q,E), marked TT2


def step(s,c):
    r,M,q,E=s
    if c=='D':return r,2*M,q,0
    if c=='A':return -r,3*q+20*M,8*M,6*M+2*r-q
    a=int(c)
    return -r,(a*a+1)*M+a*q-E,(2*a-2)*M,2*r+a*M-q


def reduce_state(s,n):return (s[0],)+tuple(v%n for v in s[1:])


def full_bfs(n,initial=S0):
    root=reduce_state(initial,n);parent={root:None};queue=deque([root]);depth={root:0}
    while queue:
        s=queue.popleft()
        for c in '35':
            t=reduce_state(step(s,c),n)
            if t not in parent:
                parent[t]=(s,c);depth[t]=depth[s]+1;queue.append(t)
    return parent,depth


def quotient_bfs(n,sign):
    root=(sign,0,3%n) if sign==-1 else (sign,2%n,1%n)
    parent={root:None};queue=deque([root]);depth={root:0}
    while queue:
        r,q,M=queue.popleft();s=(r,q,M)
        for c,t in [('A',(-r,8*M%n,(3*q+20*M)%n)),('D',(r,q,2*M%n))]:
            if t not in parent:
                parent[t]=(s,c);depth[t]=depth[s]+1;queue.append(t)
    return parent,depth


def word_to(parent,target):
    if target not in parent:raise ValueError('Target not in this finite orbit')
    word=[]
    while parent[target] is not None:
        target,c=parent[target];word.append(c)
    return ''.join(reversed(word))


def construct_slope(n,a,sign=-1):
    """Actual source controller; n>1, gcd(n,30)=1. Returns word, graph, O, state."""
    if n<=1 or gcd(n,30)!=1:raise ValueError('Proved constructor requires n>1 coprime to30')
    parent,_=quotient_bfs(n,sign)
    word=word_to(parent,(sign,1,a%n))
    g,O,s=physical(word,sign)
    assert s[0]==sign and s[1]%n==a%n and s[2]%n==1
    assert len(g)<=8*n*n-2
    return word,g,O,s


def construct_full(p,target,initial=S0):
    """Residue-state positive word. Prime theorem except761; finite search may extend."""
    parent,_=full_bfs(p,initial)
    word=word_to(parent,reduce_state(target,p));s=initial
    for c in word:s=step(s,c)
    assert reduce_state(s,p)==reduce_state(target,p)
    return word,s


def close(g,O,t):
    n=len(g);out=[set(v) for v in g]+[set() for _ in range(t+1)];v=n+t
    for u in range(n):
        out[u].update(range(n,v))
        if u in O:out[v].add(u)
        else:out[u].add(v)
    for u in range(n,v):out[u].update(range(u+1,v+1))
    return out


def physical(word,sign=-1):
    g,O,s=([{1},set()],{0,1},S0) if sign==-1 else ([set()],{0},(1,1,2,1))
    for c in word:
        n=len(g)
        if c=='D':
            g=[set(v)|{n} for v in g]+[set()]
        else:
            t=2 if c=='3' else 3
            g=close(g,O,t)
            O={n+1,n+2} if c=='A' else {n,n+t-1}
        s=step(s,c)
    return g,O,s


def endpoint_tables(g):
    n=len(g);N=1<<n;incoming=[0]*n
    for u in range(n):
        for v in g[u]:incoming[v]|=1<<u
    end=[[0]*n for _ in range(N)]
    for v in range(n):end[1<<v][v]=1
    for mask in range(1,N):
        for v in range(n):
            if not mask>>v&1:continue
            rest=mask^(1<<v)
            if not rest:continue
            prev=rest&incoming[v]
            while prev:
                b=prev&-prev;prev-=b;end[mask][v]+=end[rest][b.bit_length()-1]
    return end


def recount(g,O):
    n=len(g);full=(1<<n)-1
    assert all((v in g[u]) != (u in g[v]) for u in range(n) for v in range(u+1,n))
    end=endpoint_tables(g);rev=[set() for _ in g]
    for u in range(n):
        for v in g[u]:rev[v].add(u)
    start=endpoint_tables(rev);H=[sum(row) for row in end];H[0]=1
    S=[sum(row[v] for v in O) for row in start]
    E=[sum(row[v] for v in O) for row in end]
    M=sum(H[full^mask]*S[mask] for mask in range(1,full+1))
    W=sum(E[mask]*S[full^mask] for mask in range(1,full))
    r=H[full]-W
    return r,M,H[full]+r,E[full]


def mm(A,B,n):
    return [[sum(a*b for a,b in zip(row,col))%n for col in zip(*B)] for row in A]


def eye(k):return [[int(i==j) for j in range(k)] for i in range(k)]


def mpow(A,k,n):
    B=eye(len(A))
    while k:
        if k&1:B=mm(B,A,n)
        A=mm(A,A,n);k//=2
    return B


def inv2(A,n):
    a,b=A[0];c,d=A[1];u=pow((a*d-b*c)%n,-1,n)
    return [[d*u%n,-b*u%n],[-c*u%n,a*u%n]]


def F(a):return [[-1,0,0,0],[0,a*a+1,a,-1],[0,2*a-2,0,0],[2,a,-1,0]]


def det3(cols):
    a,b,c=cols
    return a[0]*(b[1]*c[2]-b[2]*c[1])-b[0]*(a[1]*c[2]-a[2]*c[1])+c[0]*(a[1]*b[2]-a[2]*b[1])


def main():
    out={}
    # Closed quotient, physical graphs and closing counts, not arbitrary tuples.
    samples=[]
    words=['','A','D','AA','AD','DA','DD','ADD','DAD','DDA','ADA','DAA']
    words+=['3','5','33','35','53','55']
    for sign in [-1,1]:
        for word in words:
            g,O,s=physical(word,sign)
            assert recount(g,O)==s,(sign,word,s)
            G=close(g,O,1);H=sum(endpoint_tables(G)[-1])
            assert H==s[0]+2*s[1]
            samples.append({'sign':sign,'word':word,'order':len(g),'state_r_M_q_E':s,'one_tail_count':H})
    out['literal_graphs']=samples
    # Positive finite-word constructor gives small modular targets independently recounted.
    small=[]
    for a in [1,3]:
        word,g,O,s=construct_slope(7,a,-1)
        assert recount(g,O)==s
        small.append({'modulus':7,'target_slope':a,'word':word,'order':len(g),'state':s})
    out['literal_modular_constructor']=small
    # All sign/q/M residue vertices, plus explicit target words for each a.
    quotient=[]
    for n in [7,11,13,17,19,23,29,31,49,77,91,121]:
        D=[[1,0],[0,2]];A=[[0,8],[3,20]]
        L=mm(mm(A,D,n),inv2(A,n),n)
        comm=mm(mm(mm(D,L,n),inv2(D,n),n),inv2(L,n),n)
        assert comm==[[1,0],[5*pow(4,-1,n)%n,1]]
        for sign in [-1,1]:
            parent,depth=quotient_bfs(n,sign)
            assert all((sign,1,a) in parent for a in range(n))
            quotient.append({'modulus':n,'sign':sign,'orbit_size':len(parent),'largest_target_distance':max(depth[(sign,1,a)] for a in range(n))})
    out['quotient_orbits']=quotient
    # Hidden affine integer clock on every node in a finite mixed-word tree.
    clock_count=0
    for k in range(8):
        for w in product('35',repeat=k):
            s=S0
            for c in w:s=step(s,c)
            z=s[0]*(s[1]-4*s[2]+s[3]);t=step(s,'5');z2=t[0]*(t[1]-4*t[2]+t[3])
            assert z2==z-2
            u=step(t,'5');assert u[0]==s[0] and u[0]*(u[1]-4*u[2]+u[3])==z-4
            clock_count+=1
    out['clock_states']=clock_count
    # Pure full-state translation at semisimple primes; verify exceptional boundaries.
    translation=[];v=[1,-8,-13]
    for p in [3,7,11,13,17,19,23,29,31,43,79,101,127,257,769]:
        power=mpow(F(5),p*p-1,p);target=eye(4)
        for i,vi in enumerate(v,1):target[i][0]=vi*pow(10,-1,p)%p
        assert power==target,(p,power,target)
        translation.append(p)
    out['pure_translation_primes']=translation
    exceptional=[]
    for p in [5,761]:
        power=mpow(F(5),p*p-1,p)
        linear=[[power[i][j] for j in range(1,4)] for i in range(1,4)]
        assert linear!=eye(3)
        exceptional.append({'prime':p,'full_power_matrix':power,'pure_translation':False})
    out['semisimple_formula_exceptions']=exceptional
    d1=det3([v,[-1,4,11],[-9,-4,-7]])
    d2=det3([v,[-1,4,11],[-17,-8,-9]])
    assert (d1,d2,gcd(d1,d2))==(344,632,8)
    out['translation_span_minors']=[d1,d2]
    # These are EXACT residue-state universes, NOT literal tournament censuses.
    full=[]
    for n in range(3,32,2):
        parent,depth=full_bfs(n)
        assert len(parent)==2*n**3
        # Every retained BFS edge is checked, as is complete vertex count.
        for target,data in parent.items():
            if data is not None:
                prev,c=data;assert reduce_state(step(prev,c),n)==target
        full.append({'modulus':n,'states':len(parent),'largest_distance':max(depth.values()),'scope':'finite exact; composite cases not extrapolated'})
        if n==5:
            cert=[{'state':s,'parent':None if data is None else data[0],'letter':None if data is None else data[1]} for s,data in parent.items()]
            (HERE/'MOD5_TREE.json').write_text(json.dumps({'root':reduce_state(S0,5),'states':cert},indent=2)+'\n')
    out['full_state_orbits']=full
    sample=[]
    for p,t in [(3,(1,0,0,0)),(5,(-1,0,0,0)),(7,(-1,6,2,3)),(11,(1,0,1,10))]:
        w,s=construct_full(p,t)
        sample.append({'prime':p,'target':t,'word':w,'actual_state':s,'order':2+sum(3 if c=='3' else 4 for c in w)})
    out['full_constructor_samples']=sample
    (HERE/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'literal_graphs':len(samples)+len(small),'clock_states':clock_count,'quotient_orbits':len(quotient),'full_state_moduli':len(full),'full_state_vertices':sum(v['states'] for v in full),'translation_primes':translation,'span_minors':[d1,d2],'status':'PASS'},indent=2))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--construct-slope',nargs=3,type=int,metavar=('N','RESIDUE','SIGN'))
    group.add_argument('--construct-full',nargs=5,type=int,metavar=('P','SIGN','M','Q','E'))
    args=parser.parse_args()
    if args.construct_slope:
        n,a,sign=args.construct_slope
        if sign not in (-1,1):parser.error('SIGN must be -1 or1')
        word,g,O,state=construct_slope(n,a,sign)
        result={'word':word,'state_r_M_q_E':state,'selected_vertices':sorted(O),'outgoing':[sorted(row) for row in g]}
        print(json.dumps(result,indent=2))
    elif args.construct_full:
        p,sign,M,q,E=args.construct_full
        if sign not in (-1,1):parser.error('SIGN must be -1 or1')
        if p<3 or p==761 or any(p%d==0 for d in range(2,int(p**0.5)+1)):
            parser.error('The full constructor CLI requires an odd prime other than761')
        word,state=construct_full(p,(sign,M,q,E))
        g,O,check=physical(word)
        assert state==check
        print(json.dumps({'word':word,'state_r_M_q_E':state,'selected_vertices':sorted(O),'outgoing':[sorted(row) for row in g]},indent=2))
    else:main()
