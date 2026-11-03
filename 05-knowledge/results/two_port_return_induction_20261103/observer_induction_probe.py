#!/usr/bin/env python3
"""Two-port observer image guards and exact terminal-port liveness compiler.
No target-coverage scan. Matrix source is pinned; graph checks are independent.
"""
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
from itertools import permutations,product
from math import gcd,lcm
import hashlib,json

SOURCE=Path(__file__).resolve().parents[1]/'adaptive_state_induction_20261103'/'live_port.json'


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def dot(a,b):return sum(x*y for x,y in zip(a,b))


def rref(A):
    A=[list(map(Fraction,row)) for row in A];pivots=[];r=0
    for j in range(len(A[0])):
        i=next((i for i in range(r,len(A)) if A[i][j]),None)
        if i is None:continue
        A[r],A[i]=A[i],A[r];z=A[r][j];A[r]=[x/z for x in A[r]]
        for i in range(len(A)):
            if i!=r and A[i][j]:
                z=A[i][j];A[i]=[x-z*y for x,y in zip(A[i],A[r])]
        pivots.append(j);r+=1
        if r==len(A):break
    return A,pivots


def primitive(v):
    d=lcm(*(x.denominator for x in v));w=[int(d*x) for x in v];g=gcd(*w)
    w=[x//g for x in w]
    if next(x for x in w if x)<0:w=[-x for x in w]
    return w


def kernel(A):
    R,pivots=rref(A);free=[j for j in range(len(A[0])) if j not in pivots];out=[]
    for j in free:
        v=[Fraction(0)]*len(A[0]);v[j]=1
        for i,p in enumerate(pivots):v[p]=-R[i][j]
        out.append(primitive(v))
    need(all(dot(row,v)==0 for row in A for v in out),'exact primitive image equations')
    return len(pivots),out


def birth(T,v):
    n=len(T)
    return tuple(row|((1<<n) if i==v else 0) for i,row in enumerate(T))+(((1<<n)-1)^(1<<v),)


def run(T,ports,word):
    ports=list(ports)
    for a in word:
        r,d=divmod(a,3);n=len(T);T=birth(T,ports[r])
        if d<2:ports[d]=n
    return T,tuple(ports)


@lru_cache(None)
def H(T):
    n=len(T);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,value in enumerate(dp[mask]):
            choices=T[v]&~mask
            while value and choices:
                bit=choices&-choices;choices-=bit
                dp[mask|bit][bit.bit_length()-1]+=value
    return sum(dp[-1])


def strong(T):
    reach=[row|(1<<v) for v,row in enumerate(T)];full=(1<<len(T))-1
    for k in range(len(T)):
        for v in range(len(T)):
            if reach[v]>>k&1:reach[v]|=reach[k]
    return all(row==full for row in reach)


def profile(T,v):
    n=len(T);h=s=F=0
    for p in permutations(range(n)):
        bad=[i for i in range(n-1) if not(T[p[i]]>>p[i+1]&1)]
        if not bad:h+=1;s+=p[0]==v;F+=1
        elif len(bad)==1 and p[bad[0]]==v:F+=1
    q=sum(all(T[p[i]]>>p[i+1]&1 for i in range(n-2))
          for p in permutations([w for w in range(n) if w!=v]))
    return h,s,q,F


def extracted_profile(T,ports,r):
    h=H(T);keep=3*r+2;reset=3*r+(1-r);new_slot=1-r
    h1=H(run(T,ports,(keep,))[0]);hRR=H(run(T,ports,(keep,keep))[0])
    hRN=H(run(T,ports,(reset,3*new_slot+2))[0])
    F=Fraction(hRN-h,2);q=F-Fraction(hRR-h1,2);s=h+F-h1
    return h,s,q,F


def image_members(v,kernels):return [i for i,K in enumerate(kernels) if all(dot(v,z)==0 for z in K)]


INITIAL_PAIRS=((0,1),(1,2),(2,0))


def mandatory_sets(roots,terminal):
    n=3+len(roots)
    return [set(v for v in terminal if v<t)|set(v for v in roots[t-3:] if v<t) for t in range(3,n+1)]


def compile_marked(roots,terminal):
    """Two-color mandatory lifetimes, with initial arc and final slot guards."""
    live=mandatory_sets(roots,terminal)
    if any(len(S)>2 for S in live):return None
    vertices=set().union(*live);adj={v:set() for v in vertices};forced={terminal[0]:0,terminal[1]:1}
    for S in live:
        if len(S)==2:
            a,b=sorted(S);adj[a].add(b);adj[b].add(a)
    if len(live[0])==2:
        a,b=next(pair for pair in INITIAL_PAIRS if set(pair)==live[0])
        for v,color in ((a,0),(b,1)):
            if v in forced and forced[v]!=color:return None
            forced[v]=color
    colors={}
    for v in sorted(vertices):
        if v in colors:continue
        relative={v:0};todo=[v]
        while todo:
            a=todo.pop()
            for b in adj[a]:
                c=1-relative[a]
                if b in relative and relative[b]!=c:return None
                if b not in relative:relative[b]=c;todo.append(b)
        flips={relative[a]^forced[a] for a in relative if a in forced}
        if len(flips)>1:return None
        flip=next(iter(flips),0)
        colors.update({a:c^flip for a,c in relative.items()})
    initial=next((pair for pair in INITIAL_PAIRS if all(pair[colors[v]]==v for v in live[0])),None)
    need(initial is not None,'compatible oriented seed fill')
    ports=list(initial);word=[]
    for t,v in enumerate(roots,3):
        need(v in ports,'mandatory current root')
        root=ports.index(v);drop=colors[t] if t in colors else 2
        word.append(3*root+drop)
        if drop<2:ports[drop]=t
        need(all(ports[colors[a]]==a for a in live[t-2]),'retained mandatory names')
    need(tuple(ports)==tuple(terminal),'exact final slots')
    return {'initial_pair':initial,'commands':tuple(word),'terminal_ports':tuple(ports)}


def brute_terminals(roots):
    states=set(INITIAL_PAIRS)
    for t,v in enumerate(roots,3):
        out=set()
        for ports in states:
            if v not in ports:continue
            for drop in range(3):
                child=list(ports)
                if drop<2:child[drop]=t
                out.add(tuple(child))
        states=out
    return states


def main():
    source=json.loads(SOURCE.read_text());rank=source['rank'];M=[[[Fraction(z) for z in row] for row in A] for A in rank['query_matrices']]
    ranks=[];kernels=[]
    for A in M:
        r,K=kernel(A);ranks.append(r);kernels.append(K)
    need(ranks==[13,4,13,4,13,13],'command ranks')
    need(all(dot(row,z)==0 for row in M[1] for z in kernels[4]),'image1 subset image4')
    need(all(dot(row,z)==0 for row in M[3] for z in kernels[0]),'image3 subset image0')
    queries=rank['suffixes'];seed=tuple(rank['query_seed']);C3=(2,4,1)
    def query(T,ports):return tuple(H(run(T,ports,w)[0]) for w in queries)
    weak=(58,60,57,16,32,8);alternative=birth(birth(C3,0),0)
    T4=birth(C3,0);R5=tuple(sum(1<<j for j in range(5) if (j-i)%5 in(1,2)) for i in range(5))
    objects=[('weak_seed_scale3',weak,(0,1)),('strong_H9',alternative,(0,1)),('marked_H5_outsider',T4,(1,2)),('regular_H15',R5,(0,1))]
    receipt=[]
    for name,T,ports in objects:
        v=query(T,ports);ps=[]
        for r,p in enumerate(ports):
            z=profile(T,p);need(tuple(extracted_profile(T,ports,r))==z,'literal rooted profile extraction');ps.append(z)
        factors=[z[3]-z[0]-z[2]+z[1] for z in ps]+[z[2]-z[1] for z in ps]
        receipt.append({'name':name,'rows':T,'ports':ports,'H':H(T),'strong':strong(T),'query_vector':v,
                        'rooted_profiles':ps,'birth_face_factors':factors,'linear_last_images':image_members(v,kernels)})
    need(receipt[0]['query_vector']==tuple(3*x for x in seed) and not receipt[0]['strong'],'actual weak scalar ray')
    need(receipt[1]['H']==9 and receipt[1]['strong'] and receipt[1]['query_vector']!=receipt[0]['query_vector'],'strong same-H control')
    first=[]
    for a in range(6):
        T,ports=run(C3,(0,1),(a,));v=query(T,ports);first.append(v[:3])
    need(receipt[2]['query_vector'][:3]==(5,13,11) and (5,13,11) not in first,'first-birth marked outsider')
    need(receipt[2]['linear_last_images']==[5] and receipt[3]['birth_face_factors']==[4,4,2,2],'nested guard controls')
    weak_parent=(2,0,3);child,ports=run(weak_parent,(0,1),(5,));transport=(1,2,3,0)
    transported=[0]*4
    for i,row in enumerate(child):
        for j in range(4):
            if row>>j&1:transported[transport[i]]|=1<<transport[j]
    need(tuple(transported)==T4 and tuple(transport[v] for v in ports)==(1,2),'actual weak-parent inverse')
    need(H(weak_parent)==1 and not strong(weak_parent) and query(child,ports)==receipt[2]['query_vector'],'strong-child weak-parent distinction')
    word_layers=[()];tested=accepted=0
    for depth in range(4):
        next_words=[]
        for roots in word_layers:
            n=3+len(roots);brute=brute_terminals(roots)
            for terminal in permutations(range(n),2):
                certificate=compile_marked(roots,terminal)
                need((certificate is not None)==(terminal in brute),'independent marked-slot reachability')
                tested+=1
                if certificate is not None:
                    T=C3
                    for root in roots:T=birth(T,root)
                    U,ports=run(C3,certificate['initial_pair'],certificate['commands'])
                    need(U==T and ports==terminal,'literal compiled named source');accepted+=1
            next_words.extend(roots+(root,) for root in range(n))
        word_layers=next_words
    need(tested==2082,'small liveness universe')
    need(compile_marked((),(1,0)) is None and compile_marked((0,),(1,2)) is None,'orientation and width hostiles')
    return {'status':'exact rank/image certificates and independent named-source controls; no target-support scan',
            'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'command_ranks':ranks,
            'primitive_right_kernel_image_equations':kernels,'seed_ray_last_images':image_members(seed,kernels),
            'first_birth_H_and_two_next_counts':first,'actual_controls':receipt,
            'marked_H5_weak_inverse':{'parent_rows':weak_parent,'parent_ports':[0,1],'command':5,
                'parent_H':1,'child_H':5,'transport_to_physical_T4':transport},
            'terminal_liveness_checks':tested,'accepted_and_replayed':accepted,
            'liveness_scope':'all chronological root histories through3 births, all ordered terminal pairs; three oriented initial C3 gauges retained'}


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
