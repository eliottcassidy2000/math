"""Exact packed cycle responses, with the terminal singleton tail batched by FLINT."""
from math import factorial,gcd
from time import perf_counter
from flint import fmpz,fmpz_poly
def need(b,s):
    if not b:raise ValueError(s)
def prefix(N,tail='convolution',backend='flint',audit=False):
    need(type(N) is int and N>=0,'nonnegative order bound')
    Z=fmpz if backend=='flint' else int
    fact=[factorial(k) for k in range(N+1)];v2=[k-k.bit_count() for k in range(N+1)];odd=[f>>v for f,v in zip(fact,v2)];D=Z(odd[N])
    width=max(1,N.bit_length());digit=(1<<width)-1;masks=[(1<<(width*k))-1 for k in range((N+1)//2+1)]
    states={(0,0):D};transitions=0;peak=1;rows=[];start=perf_counter()
    for ell in reversed(range(3,N+1,2)):
        row=sum(gcd(ell,j)<<(width*(j//2)) for j in range(1,ell,2));idx=(ell-1)//2;limit=N//ell
        den=[Z(ell**k*odd[k]) for k in range(limit+1)];base=[k*(ell-1)//2+ell*k*(k-1)//2-v2[k] for k in range(limit+1)]
        nxt={};cnt=0
        for (degree,response),weight in states.items():
            charge=(response>>(width*idx))&digit
            for k in range((N-degree)//ell+1):
                value=weight<<(base[k]+k*charge)
                if audit:need(value%den[k]==0,'exact raw-history division')
                value//=den[k];nd=degree+ell*k;cap=min(ell-2,N-nd);r=(response+k*row)&masks[(cap+1)//2]
                # The singleton response is still needed even when no space remains.
                if N==nd:r=(response+k*row)&digit
                key=(nd,r);nxt[key]=nxt.get(key,Z(0))+value;cnt+=1
        states=nxt;transitions+=cnt;peak=max(peak,len(states));rows.append([ell,len(states),cnt])
    core_seconds=perf_counter()-start;core_states=len(states);start=perf_counter();ans=[Z(0)]*(N+1)
    if tail=='loop':
        for (degree,r),w in states.items():
            for k in range(N-degree+1):
                z=w<<(k*(k-1)//2+k*r-v2[k]);need(not audit or z%odd[k]==0,'singleton exact division');ans[degree+k]+=z//odd[k]
        scale=D
    elif tail=='convolution':
        grouped={}
        for (degree,r),w in states.items():
            if r not in grouped:grouped[r]=[fmpz(0)]*(N+1)
            grouped[r][degree]+=fmpz(w)
        base=[(fmpz(D)<<(k*(k-1)//2-v2[k]))//odd[k] for k in range(N+1)]
        for r,U in grouped.items():
            low=next(i for i,x in enumerate(U) if x);L=N-low
            V=fmpz_poly([base[k]<<(k*r) for k in range(L+1)])
            prod=fmpz_poly(U[low:]).mul_low(V,L+1)
            for k in range(L+1):ans[low+k]+=prod[k]
        scale=D*D
    else:raise ValueError('tail mode')
    out=[]
    for z in ans:
        q,r=divmod(z,scale);need(r==0,'integral Burnside total');out.append(int(q))
    return out,dict(N=N,backend=backend,tail=tail,core_seconds=core_seconds,tail_seconds=perf_counter()-start,core_states=core_states,peak_states=peak,core_transitions=transitions,rows=rows)
if __name__=='__main__':
    import sys,json
    sys.set_int_max_str_digits(0)
    N=int(sys.argv[1]);mode=sys.argv[2] if len(sys.argv)>2 else 'convolution';a,s=prefix(N,mode);print(json.dumps(s));print(a[-1])
