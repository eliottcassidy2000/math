"""Independent fixed-endpoint Fourier isolation; preserves prior frozen files."""
from common import core,response
from itertools import combinations

def fixed_endpoint_sum(Q,u,v,base_to_u=True,fixed=None):
    k=len(Q);m=k-1
    if k<2 or u==v:raise ValueError('distinct endpoints required')
    fixed=[0]*m if fixed is None else list(fixed)
    if len(fixed)!=m or any(not 0<=x<(1<<k) for x in fixed):raise ValueError('bad masks')
    R=[v,k]
    for i,w in enumerate(w for w in range(k) if w not in (u,v)):R.extend([w,k+i+1])
    cross=[]
    for a,b in zip(R,R[1:]):
        cross.append((b-k,a,1) if a<k else (a-k,b,0))
    second_bit_xor=0;total=0;residues={M:0 for M in [2,3,4,5,6,8,11]};queries=0
    for bits in range(1<<len(cross)):
        word=fixed[:];chi=1
        for j,(i,w,plusbit) in enumerate(cross):
            positive=bits>>j&1;bit=plusbit if positive else 1-plusbit
            word[i]=(word[i]&~(1<<w))|(bit<<w)
            if not positive:chi=-chi
        word[-1]=(word[-1]&~(1<<u))|((not base_to_u)<<u)
        z=response(Q,word)
        if z%2!=1:raise RuntimeError('Redei parity control failed')
        total+=chi*z;queries+=1;second_bit_xor^=((z-1)//2)&1
        for M in residues:residues[M]=(residues[M]+chi*(z%M))%M
    if total%2 or second_bit_xor!=(total//2)%2:raise RuntimeError("second-bit XOR failed")
    return total,residues,queries

if __name__=='__main__':
    tests=queries=0
    for k in [2,3,4]:
        for mask in range(1<<(k*(k-1)//2)):
            Q=core(k,mask)
            for u,v in combinations(range(k),2):
                fixed=[(11*i+mask+u+3*v)%(1<<k) for i in range(k-1)]
                for base_to_u in [True,False]:
                    total,res,n=fixed_endpoint_sum(Q,u,v,base_to_u,fixed)
                    expected=2*Q[u][v]-(0 if base_to_u else 2)
                    if total!=expected or any(res[M]!=expected%M for M in res):raise RuntimeError(('fixed endpoint failed',k,mask,u,v,base_to_u,total,res))
                    if base_to_u and any((res[M]==2)!=bool(Q[u][v]) for M in res if M>=3):raise RuntimeError('modular reconstruction failed')
                    tests+=1;queries+=n
        print('all labelled cores through',k,'PASS;',tests,'tests,',queries,'DP queries',flush=True)
    for mask in [0,2,41,341,1023]:
        Q=core(5,mask)
        for u,v in [(0,4),(3,1)]:
            for base_to_u in [True,False]:
                total,res,n=fixed_endpoint_sum(Q,u,v,base_to_u,[3,17,9,30]);expected=2*Q[u][v]-(0 if base_to_u else 2)
                if total!=expected or any(res[M]!=expected%M for M in res):raise RuntimeError('k5 sample failed')
                tests+=1;queries+=n
    print('TOTAL PASS:',tests,'tests,',queries,'DP queries; moduli2,3,4,5,6,8,11; both fixed-endpoint directions.')
