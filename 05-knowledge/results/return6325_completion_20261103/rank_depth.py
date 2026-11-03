"""Literal controls for the all-N, all-depth binomial response-rank theorem."""
from halving_response import wrapped_source,ear_values,pair_flip,birth,hp,need,HERE
from itertools import permutations
import sys,json,math
def fwht(a):
    a=list(a);h=1
    while h<len(a):
        for i in range(0,len(a),2*h):
            for j in range(i,i+h):a[j],a[j+h]=a[j]+a[j+h],a[j]-a[j+h]
        h*=2
    return a
def main():
    controls=[]
    for N in range(1,6):
        T=wrapped_source(N);n=len(T);f=ear_values(T,sys.argv[1]);seen=set()
        for cut in range(1,(1<<n)-1):
            vals=[]
            for x in range(1<<N):
                p=[2*b+(j^((x>>b)&1)) for b in range(N) for j in range(2)]+[2*N,2*N+1]
                vals.append(f[sum(1<<p[i] for i in range(n) if cut>>i&1)])
            seen.update(i for i,v in enumerate(fwht(vals)) if v)
        expected={x for x in range(1<<N) if x.bit_count()<=2}
        need(seen==expected,'complete one-birth Fourier support equals degree<=2')
        controls.append(dict(N=N,source_order=n,contexts=(1<<n)-2,rank=len(seen),expected=sum(math.comb(N,j) for j in range(min(2,N)+1))))
        print('rank-depth',N,len(seen),flush=True)
    gaps=[]
    for ell in range(7):
        n=4+ell;TT=tuple(sum(1<<j for j in range(i+1,n)) for i in range(n));counts=[]
        for x in range(4):
            p=list(range(n))
            if x&1:p[0],p[1]=p[1],p[0]
            if x&2:p[-2],p[-1]=p[-1],p[-2]
            U=tuple(sum(1<<j for j in range(n) if TT[p[i]]>>p[j]&1) for i in range(n))
            counts.append(hp(birth(U,1|(1<<(n-2)))))
        mixed=counts[0]-counts[1]-counts[2]+counts[3]
        need(mixed==2**ell+(ell==0),'spaced reader exact mixed Hamilton term')
        gaps.append(dict(middle_vertices=ell,counts=counts,mixed=mixed))
    (HERE/'rank_depth.json').write_text(json.dumps(dict(status='FINITE_EXACT_CONTROLS_FOR_ALL_SIZE_RANK_PROOF',one_birth=controls,spaced_readers=gaps),indent=2)+'\n')
if __name__=='__main__':main()
