#!/usr/bin/env python3
"""Adjacency-only independent check of the proposed two-apex binary decoders."""
import itertools,json
from functools import lru_cache
from fractions import Fraction
S=('','0','1','2','00','10','20','30','01','11','21','31','000','100','200','010','210')
CX=(0,-5,-3,6,15,-9,-18,6,9,0,0,0,-4,6,6,-9,0)
CY=(0,-1,3,0,3,-3,0,0,-3,0,0,0,-2,0,0,3,0)
def require(c,msg):
    if not c:raise RuntimeError(msg)
@lru_cache(None)
def H(w):
    m=len(w);n=m+2;out=[0]*n
    for i,c in enumerate(w):
        for j in range(i+1,m):out[i]|=1<<j
        for k in range(2):
            if int(c)>>k&1:out[m+k]|=1<<i
            else:out[i]|=1<<(m+k)
    out[m]|=1<<(m+1);size=1<<n;d=[[0]*n for _ in range(size)]
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,size):
        for i,v in enumerate(d[mask]):
            if not v:continue
            todo=out[i]&~mask
            while todo:
                b=todo&-todo;todo-=b;j=b.bit_length()-1;d[mask|b][j]+=v
    return sum(d[-1])
def profile(w):return tuple(H(w+s) for s in S)
def decode(P):
    x=sum(a*b for a,b in zip(CX,P));y=sum(a*b for a,b in zip(CY,P))
    require(x%36==0 and y%12==0,'nonintegral code');x//=36;y//=12
    require(x>=0 and y>=0,'negative code')
    n=max(x.bit_length(),y.bit_length())
    return ''.join(str((x>>i&1)+2*(y>>i&1)) for i in range(n-1,-1,-1)),x,y

def closure_certificate():
    # Independently construct the three-run assignment state, with -1 = empty,
    # 0 = admissible last edge, 1 = inadmissible last edge. This gauge differs
    # from the producer and makes the no-middle x->y guard explicit.
    states=tuple(itertools.product((0,1),(-1,0,1),(-1,0,1),(False,True)))
    index={v:i for i,v in enumerate(states)};d=len(states)
    transitions=[]
    for letter in range(4):
        edges=[]
        for i,(order,P,Q,R) in enumerate(states):
            A,B=(0,1) if order==0 else (1,0)
            a=(letter>>A)&1;b=(letter>>B)&1
            edges.append((i,index[order,a,Q,R]))
            if Q!=-1 or a:edges.append((i,index[order,P,b,R]))
            if R or b:edges.append((i,index[order,P,Q,True]))
        transitions.append(edges)
    init=tuple(int(P==Q==-1 and not R) for order,P,Q,R in states)
    out=tuple(int(P!=1 and Q!=1 and (order==0 or Q==0)) for order,P,Q,R in states)
    def advance(v,a):
        z=[0]*d
        for i,j in transitions[a]:z[j]+=v[i]
        return tuple(z)
    def dual(v,a):
        z=[0]*d
        for i,j in transitions[a]:z[i]+=v[j]
        return tuple(z)
    def suffix(s):
        v=out
        for a in s[::-1]:v=dual(v,int(a))
        return v
    probes=tuple(suffix(s) for s in S)
    # Grow the invariant reachable span by exact rational elimination.
    pivots={};basis=[];queue=[]
    def insert(v):
        original=v;v=list(map(Fraction,v))
        for p,row in sorted(pivots.items()):
            if v[p]:
                c=v[p];v=[a-c*b for a,b in zip(v,row)]
        p=next((i for i,a in enumerate(v) if a),None)
        if p is None:return False
        scale=v[p];pivots[p]=tuple(a/scale for a in v)
        basis.append(original);queue.append(original);return True
    insert(init)
    while queue:
        v=queue.pop(0)
        for a in range(4):insert(advance(v,a))
    require(len(basis)==28,'unexpected reachable span')
    qx=tuple(int(o==0 and Q!=-1 and not R) for o,P,Q,R in states)
    qy=tuple(int(o==0 and Q==-1 and R) for o,P,Q,R in states)
    ell=tuple(int(o==0 and Q==-1 and not R) for o,P,Q,R in states)
    def dot(v,w):return sum(a*b for a,b in zip(v,w))
    require(dot(ell,init)==1 and all(dual(ell,a)==ell for a in range(4)),'constant clock')
    checks=0
    for bit,(coef,den,q) in enumerate(((CX,36,qx),(CY,12,qy))):
        residual=tuple(sum(coef[k]*probes[k][i] for k in range(len(S)))-den*q[i] for i in range(d))
        require(all(dot(residual,v)==0 for v in basis),'all-order extraction residual')
        require(dot(q,init)==0,'code initial value')
        for a in range(4):
            require(dual(q,a)==tuple(2*q[i]+((a>>bit)&1)*ell[i] for i in range(d)),'exact code recurrence')
            checks+=d
    # Reject the smallest coefficient mutation on the entire closed span.
    mutated=list(CX);mutated[1]+=1
    bad=tuple(sum(mutated[k]*probes[k][i] for k in range(len(S)))-36*qx[i] for i in range(d))
    require(any(dot(bad,v)!=0 for v in basis),'coefficient mutation escaped')
    return {'reachable_span_dimension':len(basis),'exact_digit_recurrence_coordinate_checks':checks,'all_order_extraction_residuals':2,'method':'independently rebuilt run transfer and invariant rational span; no producer import','mutated_coefficient_rejected':True}

def main():
    cases=0;profiles={};hostile=[]
    for n in range(5):
        for w0 in itertools.product('0123',repeat=n):
            w=''.join(w0);P=profile(w);v,x,y=decode(P)
            wantx=int(''.join(str(int(c)&1) for c in w) or '0',2)
            wanty=int(''.join(str(int(c)>>1) for c in w) or '0',2)
            require((x,y)==(wantx,wanty),'binary address mismatch '+w)
            require(v==w.lstrip('0'),'decode not original normalized word '+w)
            if P in profiles:require(profiles[P]==v,'profile collision beyond leading zeros')
            profiles[P]=v;cases+=1
            if w in ('','0','1','2','12','21'):hostile.append({'word':w,'H':P[0],'Qx':x,'Qy':y})
    P=profile('1');require(decode(P)[1:]!=(0,1),'swapped bits hostile did not reject')
    # Positive lattice multiples are not automatic words: decoder/replay must still reject them.
    negatives=[]
    for k in (7,21):
        P=tuple(k*v for v in profile(''));w,x,y=decode(P)
        require(profile(w)!=P,'spurious formal state accepted');negatives.append({'formal_multiple':k,'decoded_word':w,'Qx':x,'Qy':y,'replay_H':H(w)})
    print(json.dumps({'status':'PASS','all_order_certificate':closure_certificate(),'all_word_lengths':list(range(5)),'actual_profiles':cases,'normalized_profile_count':len(profiles),'independent_adjacency_H_evaluations':H.cache_info().currsize,'method':'generic Hamiltonian-path subset DP; no transfer or producer import','hostiles':hostile,'formal_state_rejections':negatives,'scope':'All-order code identities verified on a closed reachable span; actual-word membership requires exact replay. Scalar spectrum completeness remains open.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
