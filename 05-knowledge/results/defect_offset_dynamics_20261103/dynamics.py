"""Exact complementary two-ear responses on the finite defect atlases."""
from pathlib import Path
import sys,json
from collections import Counter
from functools import lru_cache
from itertools import product

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'ready_return_pentagon_20261103'))
from return_program import (graded,require,incoming_pair_source,induced,opposite,
                           rooted_data,tables,strong,field,observer,converse,project)
from verify import exceptional_signature
HERE=Path(__file__).resolve().parent


def named_kernel(skeleton,count,k,polarity,mask):
    state=list(skeleton)+[(-1,-1)]*count
    incoming=[None]+[i for i,(_,e) in enumerate(state) if (polarity if e<0 else (mask>>e)&1)]
    outgoing=[None]+[i for i,(s,_) in enumerate(state) if not (polarity if s<0 else (mask>>s)&1)]
    out=Counter()
    for i in incoming:
        for j in outgoing:
            if i is not None and i==j:continue
            parts=[p for t,p in enumerate(state) if t!=i and t!=j]
            parts.append((k if i is None else state[i][0],k if j is None else state[j][1]))
            out[graded.split(tuple(sorted(parts)))]+=1
    return out


@lru_cache(None)
def lift(k,polarity,mask):
    result={}
    for state in graded.skeletons(k):
        A=named_kernel(state,0,k,polarity,mask);B=named_kernel(state,1,k,polarity,mask)
        for (target,degree),value in A.items():B[target,degree+1]-=value
        require(all(value>=0 for value in B.values()),'nonnegative named-birth derivative')
        B=+B;result[state]=(A,B)
    return result


def apply_lift(vector,k,polarity,mask,cap=None):
    result=Counter();op=lift(k,polarity,mask)
    for (old,degree),value in vector.items():
        A,B=op[old]
        for (target,power),weight in A.items():
            d=degree+power
            if cap is None or len(target)+d<=cap:result[target,d]+=value*weight
        if degree:
            for (target,power),weight in B.items():
                d=degree-1+power
                if cap is None or len(target)+d<=cap:result[target,d]+=value*degree*weight
    return +result


def source(word):
    T=(2,0);ports=(0,1);v=graded.literal(T,ports)
    for a,b in word:
        command=(1,(1-a)+2*(1-b),2)
        T,ports=graded.birth(T,ports,*command)
        v=graded.apply(v,2,graded.operator(2,command))
    T,_=graded.birth(T,ports,0,3,2)
    v=apply_lift(v,2,0,3);ports+=(len(T)-1,)
    return T,ports,v


def last_h(vector,polarity,mask):
    h=0
    incoming=lambda e:polarity if e<0 else (mask>>e)&1
    for (skel,anon),weight in vector.items():
        parts=list(skel)+[(-1,-1)]*anon
        if len(parts)==1:
            s,e=parts[0];h+=weight*(incoming(e)+1-incoming(s))
        elif len(parts)==2:
            for i,j in ((0,1),(1,0)):h+=weight*incoming(parts[i][1])*(1-incoming(parts[j][0]))
    return h


@lru_cache(None)
def rows():
    basis=[(s,d) for s in graded.skeletons(3) for d in range(4-len(s)) if 1<=len(s)+d<=3]
    rows=[];commands=list(product(range(2),range(8),range(2),range(16)))
    for p,m in product(range(2),range(8)):
        first=[apply_lift(Counter({b:1}),3,p,m,2) for b in basis]
        for q,n in product(range(2),range(16)):
            rows.append([last_h(v,q,n) for v in first])
    return basis,commands,rows


def actual(T,ports,command):
    p,m,q,n=command
    T,_=graded.birth(T,ports,p,m,len(ports));ports+=(len(T)-1,)
    T,_=graded.birth(T,ports,q,n,len(ports));ports+=(len(T)-1,)
    return T,ports


def exact_vector(T,ports,v,command):
    p,m,q,n=command
    w=apply_lift(v,3,p,m);w=apply_lift(w,4,q,n)
    U,newports=actual(T,ports,command)
    return U,newports,w


def main():
    basis,commands,M=rows();print('basis',len(basis),'programs',len(M),flush=True)
    atlas=json.loads((HERE.parent/'ready_return_pentagon_20261103/languages.json').read_text())
    families=[(row['target'],tuple(map(tuple,f['tail'])),f['A'],f['B']) for row in atlas if row['target']<=6 for f in row['families']]
    all_hits=[]
    for C,tail,A,B in families:
        # Both actual source orientations are allowed, with transported port names.
        for dual in (False,True):
            candidates={(i,b) for i in range(len(commands)) for b in (2,4)}
            for k in range(1,len(tail)+7):
                T,ports,v=source(((1,1),)*k+tail)
                if dual:T,v=opposite(T),converse(v)
                vec=[v[b] for b in basis];h=A*2**k+B
                require(graded.scalar(v)==h,'source Hamilton formula')
                candidates={(i,b) for i,b in candidates if sum(x*y for x,y in zip(M[i],vec))==3*h+b}
                if not candidates:break
            for i,b in sorted(candidates):
                returns=[]
                for k in (1,2,3):
                    T,ports,v=source(((1,1),)*k+tail)
                    if dual:T,v=opposite(T),converse(v)
                    U,p,w=exact_vector(T,ports,v,commands[i]);hu=graded.scalar(w)
                    require(hu==3*(A*2**k+B)+b,'full vector versus finite-horizon row')
                    if k==1:require(tables(U)[0][-1]==hu,'independent literal child H');tables.cache_clear()
                    for con in (False,True):
                        V=opposite(U) if con else U
                        for slot,r in enumerate(p):
                            sig=exceptional_signature(V,r)
                            if sig is not None:returns.append(dict(k=k,converse=con,slot=slot,word=sig))
                all_hits.append(dict(C=C,tail=tail,A=A,B=B,dual=dual,offset=b,command=commands[i],strong=strong(U),returns=returns))
            if candidates:print('hit',C,tail,dual,[(commands[i],b) for i,b in sorted(candidates)],flush=True)
    (HERE/'search.json').write_text(json.dumps(dict(basis_size=len(basis),programs=len(commands),hits=all_hits),indent=2)+'\n')
    print('uniform hits',len(all_hits),'with all-k-small returns',sum(bool(r['returns']) for r in all_hits),flush=True)


if __name__=='__main__':main()
