#!/usr/bin/env python3
"""Exact controls for the S192 interval-cycle and directed-port constructions.
Only Python's standard library is used; checks survive python -O.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
import json

def require(ok, detail):
    if not ok: raise ValueError(detail)

def trim(p):
    while len(p)>1 and p[-1]==0: p.pop()
    return p

def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return trim(out)

def ev(p,z):return sum(x*z**i for i,x in enumerate(p))

def compositions(n, minimum=2):
    if n==0:yield ()
    for p in range(minimum,n+1):
        for tail in compositions(n-p,minimum):yield (p,)+tail

def chain(parts):
    n=sum(parts)+1; a=[sum(1<<j for j in range(i+1,n)) for i in range(n)];b=[0]
    for p in parts:b.append(b[-1]+p)
    for x,y in zip(b,b[1:]):a[x]^=1<<y;a[y]|=1<<x
    return a,b

def odd_cycles(a):
    n=len(a);out=[]
    for s in range(n):
        def rec(path,used):
            v=path[-1]
            if len(path)>=3 and len(path)%2 and a[v]>>s&1:out.append(tuple(path))
            for w in range(s+1,n):
                if not used>>w&1 and a[v]>>w&1:rec(path+[w],used|1<<w)
        rec([s],1<<s)
    return out

def interval_cycles(parts):
    a,b=chain(parts);out=[];bs=set(b)
    for i in range(len(parts)):
        for j in range(i+1,len(parts)+1):
            ints=[v for v in range(b[i],b[j]+1) if v not in bs]
            for mask in range(1<<len(ints)):
                sel=[v for t,v in enumerate(ints) if mask>>t&1]
                if (j-i+1+len(sel))%2:
                    cyc=tuple([b[i]]+sel+[b[t] for t in range(j,i,-1)])
                    if len(cyc)>=3:out.append((cyc,(i,j)))
    return out

def independence(cycles):
    masks=[sum(1<<v for v in c) for c in cycles];out=[1]
    def rec(start,used,k):
        if k:
            while len(out)<=k:out.append(0)
            out[k]+=1
        for i in range(start,len(masks)):
            if not used&masks[i]:rec(i+1,used|masks[i],k+1)
    rec(0,0,0)
    return out

def chain_polynomial(parts):
    previous=[1];current=[1]
    for index,p in enumerate(parts):
        q=2**(p-1)
        if index==0:nxt=[1,q//2]
        else:
            nxt=[0]*max(len(current),len(previous)+1)
            for i,x in enumerate(current):nxt[i]+=(1+q)*x
            for i,x in enumerate(previous):nxt[i]-=q*x;nxt[i+1]+=q*x//2
        previous,current=current,trim(nxt)
    return current

def matching_polynomial(parts):
    """Direct subsets of edges of a path, independently of interval recurrence."""
    if not parts:return [Fraction(1)]
    q=[2**(p-1) for p in parts];P=1
    for x in q:P*=1+x
    weights=[Fraction(q[0],1+q[0])]+[Fraction(q[i],(1+q[i-1])*(1+q[i])) for i in range(1,len(q))]
    coeff=[Fraction(0)]*(1+(len(q)+1)//2)
    for mask in range(1<<len(q)):
        if mask&(mask<<1):continue
        value=Fraction(P);d=mask.bit_count()
        for i,w in enumerate(weights):
            if mask>>i&1:value*=w
        coeff[d]+=value
    # Substitute t=z/2-1 exactly.
    from math import comb
    out=[Fraction(0)]*len(coeff)
    for d,x in enumerate(coeff):
        for j in range(d+1):out[j]+=x*comb(d,j)*Fraction((-1)**(d-j),2**j)
    require(all(v.denominator==1 for v in out),('matching nonintegral',parts,out))
    return trim([int(v) for v in out])

def general_poly(interiors):
    """Weighted interval scheduling, allowing empty blocks."""
    Z=[[1]]
    for k in range(1,len(interiors)+1):
        out=Z[-1].copy();m=0
        for i in range(k-1,-1,-1):
            m+=interiors[i]
            w=2**(m-1) if m else int((k-i)%2==0)
            if not w:continue
            prev=Z[max(0,i-1)]
            while len(out)<=len(prev):out.append(0)
            for j,x in enumerate(prev):out[j+1]+=w*x
        Z.append(trim(out))
    return Z[-1]

def induced_poly(parts,mask):
    _,b=chain(parts);result=[1];start=None
    for j in range(len(b)+1):
        live=j<len(b) and bool(mask>>b[j]&1)
        if live and start is None:start=j
        if not live and start is not None:
            ms=[sum(bool(mask>>v&1) for v in range(b[t-1]+1,b[t])) for t in range(start+1,j)]
            result=mul(result,general_poly(ms));start=None
    return result

def scalar_state(interiors):
    z,c,d=1,0,0
    for m in interiors:
        q=2**m;delta=int(m==0)
        nc=q*(z-d);nd=-delta*(z-c)
        z,c,d=z+nc+nd,nc,nd
    return z,c,d

def rooted_paths(a):
    n=len(a);dp=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        h[mask]=sum(dp[mask])
        for v,x in enumerate(dp[mask]):
            if not x:continue
            allowed=a[v]&~mask
            while allowed:
                bit=allowed&-allowed;allowed-=bit;w=bit.bit_length()-1
                dp[mask|bit][w]+=x
    return h,dp

def endpoint_matrix(a):
    n=len(a);M=[[0]*n for _ in range(n)]
    for word in permutations(range(n)):
        if all(a[x]>>y&1 for x,y in zip(word,word[1:])):M[word[0]][word[-1]]+=1
    return M

def tournament(n,code):
    a=[0]*n;t=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>t&1:a[i]|=1<<j
            else:a[j]|=1<<i
            t+=1
    return a

def direct_forest(a):
    n=len(a);edges=[(i,j) for i in range(n) for j in range(n) if a[i]>>j&1];out=defaultdict(int)
    for mask in range(1<<len(edges)):
        pred=[-1]*n;succ=[-1]*n;valid=True
        for t,(u,v) in enumerate(edges):
            if mask>>t&1:
                if succ[u]>=0 or pred[v]>=0:valid=False;break
                succ[u]=v;pred[v]=u
        if not valid:continue
        seen=set();pairs=[]
        for s in range(n):
            if pred[s]<0:
                t=s
                while t>=0:seen.add(t);end=t;t=succ[t]
                pairs.append((s,end))
        if len(seen)==n:out[tuple(sorted(pairs))]+=1
    return dict(out)

def arrival(profile,z,outgoing):
    out=defaultdict(int)
    for state,value in profile.items():
        def put(new):out[tuple(sorted(new))]+=value
        put(list(state)+[(z,z)])
        for i,(s,t) in enumerate(state):
            rest=list(state[:i]+state[i+1:])
            if not outgoing>>t&1:put(rest+[(s,z)])
            if outgoing>>s&1:put(rest+[(z,t)])
            if not outgoing>>t&1:
                for j,(u,v) in enumerate(state):
                    if j!=i and outgoing>>u&1:
                        rem=[pair for l,pair in enumerate(state) if l not in (i,j)]
                        put(rem+[(s,v)])
    return dict(out)

def recursive_forest(a):
    P={():1}
    for z in range(len(a)):P=arrival(P,z,a[z]&((1<<z)-1))
    return P

def main():
    report={'status':'PASS','scope':'exact bounded checks; all-order statements are proved in PROOF.md'}
    cases=cycles_checked=conflict_pairs=0
    for total in range(2,12):
        for parts in compositions(total):
            a,b=chain(parts);actual=odd_cycles(a);model=interval_cycles(parts)
            require(sorted(actual)==sorted(c for c,_ in model),('cycle normal form',parts))
            for i,(c,(u,v)) in enumerate(model):
                for d,(s,t) in model[i+1:]:
                    require(bool(set(c)&set(d))==(max(u,s)<=min(v,t)),('conflicts',parts,c,d))
                    conflict_pairs+=1
            expected=independence(actual);rec=chain_polynomial(parts);mat=matching_polynomial(parts)
            require(expected==rec==mat,('three polynomial calculations',parts,expected,rec,mat))
            require(len(rec)-1==(len(parts)+1)//2,('degree',parts))
            require(ev(rec,2)==__import__('math').prod(2**(p-1)+1 for p in parts),('Moon product',parts))
            cases+=1;cycles_checked+=len(actual)
    report['ordinary_chains']={'compositions':cases,'total_part_sums':[2,11],'cycles':cycles_checked,'conflict_pairs':conflict_pairs}
    decks=subsets=0
    for total in range(0,10):
        for parts in compositions(total):
            a,b=chain(parts);cs=odd_cycles(a);cm=[sum(1<<v for v in c) for c in cs];h,_=rooted_paths(a)
            for mask in range(1<<len(a)):
                literal=independence([c for c,m in zip(cs,cm) if m&mask==m]);abstract=induced_poly(parts,mask)
                require(literal==abstract,('induced polynomial',parts,mask,literal,abstract))
                require(ev(abstract,2)==h[mask],('induced OCF / path DP',parts,mask))
                subsets+=1
            decks+=1
    report['entire_induced_decks']={'chains':decks,'total_part_sums':[0,9],'subsets':subsets}
    general_count=0
    for total in range(0,9):
        for parts in compositions(total,1):
            a,b=chain(parts);p=general_poly([x-1 for x in parts]);literal=independence(odd_cycles(a))
            require(p==literal,('generalized cycles',parts,p,literal))
            h,dp=rooted_paths(a);z,c,d=scalar_state([x-1 for x in parts]);allmask=(1<<len(a))-1
            require(z==h[-1] and z-c==dp[-1][-1],('three-state path/end',parts,z,c,d))
            # B: first path arbitrary, second nonempty and starts at last vertex.
            root=len(a)-1;B=0
            for mask in range(1<<len(a)):
                if mask>>root&1:continue
                rem=allmask^mask
                # start-root count by reversing all arcs and ending at root.
                rev=[sum(1<<j for j in range(len(a)) if a[j]>>i&1) for i in range(len(a))]
                # DP precomputed once below would be faster, but the universe is small.
                if mask==0:rh,rdp=rooted_paths(rev)
                B+=h[mask]*rdp[rem][root]
            require(z-d==B,('three-state B response',parts,z,c,d,B))
            general_count+=1
    report['generalized_chains']={'compositions':general_count,'total_part_sums':[0,8]}
    forest_count=0
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            a=tournament(n,code);literal=direct_forest(a);transfer=recursive_forest(a)
            require(literal==transfer,('full directed forest',n,code))
            h,_=rooted_paths(a)
            require(sum(x for state,x in transfer.items() if len(state)==1)==h[-1],('forest/path',n,code))
            forest_count+=1
    report['directed_forest_transfer']={'all_labelled_tournaments_orders':[1,5],'tournaments':forest_count}
    A=[14,12,8,0];B=[14,8,10,0]
    AA=[14,12,24,16,3];BB=[14,8,26,16,3]
    require(endpoint_matrix(A)==endpoint_matrix(B), 'endpoint collision input')
    ha=rooted_paths(AA)[0][-1];hb=rooted_paths(BB)[0][-1]
    require((ha,hb)==(11,9),'endpoint collision outputs')
    report['endpoint_collision']={'A':A,'B':B,'endpoint_matrix':endpoint_matrix(A),'common_new_outneighbors':[0,1],'children':[AA,BB],'child_counts':[ha,hb]}
    flipbase,_=chain((2,3,3));flipped=flipbase.copy()
    flipped[3]^=1<<6;flipped[6]^=1<<3
    hostile=independence(odd_cycles(flipped))
    require(hostile==[1,48,11,1],('one-flip nonreal hostile',hostile))
    d,c,b,a=hostile
    discriminant=b*b*c*c-4*a*c*c*c-4*b*b*b*d-27*a*a*d*d+18*a*b*c*d
    require(discriminant==-159431,'negative cubic discriminant')
    require(rooted_paths(flipped)[0][-1]==149,'one-flip path count')
    report['one_flip_nonreal_hostile']={'base_parts':[2,3,3],'reversed_arc':[3,6],
        'rows':flipped,'polynomial':hostile,'discriminant':discriminant,'H':149}
    report['scalar_collision']={'parts':[[2,3,3],[3,2,3]],'polynomials':[chain_polynomial((2,3,3)),chain_polynomial((3,2,3))],'H':75}
    report['small_generalized_spectrum']={str(n+1):sorted({ev(general_poly([x-1 for x in p]),2) for p in compositions(n,1)}) for n in range(5)}
    report['clawed_real_rooted_example']={'parts':[2]*5,'polynomial':chain_polynomial((2,)*5),'leaves':[[0,1,2],[4,5,6],[8,9,10]],'center':[0,1,10,8,6,4,2]}
    for r in range(1,7):
        order=list(range(2*r+1));a=[sum(1<<j for j in range(i+1,len(order))) for i in order]
        require(rooted_paths(a)[0][-1]==1,('horizon transitive witness',r))
        require(sum(1 for i in order if i%2==0)==r+1,('horizon components',r))
    report['sharp_horizon_controls']={'r':[1,6],'construction':'transitive order old,new,old,...,new,old'}
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__':main()
