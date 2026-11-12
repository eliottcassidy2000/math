"""Certified real/dyadic reconstruction of A000568. No floating rounding."""
from math import factorial,gcd
from collections import Counter
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'absolute'))
from certified_interval import integer_interval

def residue(n,b,*,single_inverse=True):
    """a(n) mod 2**b; exact odd-cycle sum, precision selected first.

    Recursive branches whose accumulated valuation is >= b vanish. The
    accumulated valuation is monotone at a multiplicity-group boundary.
    """
    if n<0 or b<0: raise ValueError('nonnegative n,b required')
    if b==0: return 0,dict(leaves=0,pruned=1,nodes=0,inverse_bits=0)
    facts=[factorial(k) for k in range(n+1)]
    vals=[(f&-f).bit_length()-1 for f in facts]
    odds=[f>>v for f,v in zip(facts,vals)]
    stats=Counter(); modulus=1<<b; common=odds[n]
    def visit(j,rem,v,den,previous):
        stats['nodes']+=1
        if v>=b:
            stats['pruned']+=1
            return 0
        if rem==0:
            stats['leaves']+=1
            precision=b-v
            stats['inverse_bits']+=precision
            if single_inverse:
                factor,remainder=divmod(common,den)
                assert remainder==0
                return (factor<<v)%modulus
            return (pow(den,-1,1<<precision)<<v)
        if j==1:
            k=rem
            nv=v+k*(k-1)//2-vals[k]+k*sum(m for ell,m in previous)
            return visit(0,0,nv,den*odds[k],())
        charge=sum(m*gcd(ell,j) for ell,m in previous)
        total=0
        for k in range(rem//j+1):
            nv=v+k*(j-1)//2+j*k*(k-1)//2-vals[k]+k*charge
            total+=visit(j-2,rem-j*k,nv,den*j**k*odds[k],
                         previous+((j,k),) if k else previous)
        return total%modulus
    answer=visit(n if n%2 else max(1,n-1),n,0,1,())%modulus
    if single_inverse: answer=answer*pow(common,-1,modulus)%modulus
    return answer,dict(stats)

def recover_from_interval(lo,hi,res,b):
    """Unique integer in [lo,hi] with specified exact residue, or fail."""
    if lo>hi or b<0: raise ValueError('invalid interval/precision')
    M=1<<b
    if hi-lo>=M: raise ValueError('interval not narrower than modulus')
    x=lo+(res-lo)%M
    if x>hi: raise ValueError('inconsistent residue and interval')
    return x

def complete(n,support_cut=3,*,backend="response"):
    lo,hi=integer_interval(n,support_cut)
    b=(hi-lo).bit_length() # 2**b > width including width=0
    if b==0: return lo,dict(bits=0,lo=str(lo),hi=str(hi),leaves=0)
    r,stats=(residue_response(n,b) if backend=="response" else residue(n,b))
    x=recover_from_interval(lo,hi,r,b)
    return x,dict(stats,bits=b,lo=str(lo),hi=str(hi))


def residue_response(n,b):
    """Session270 raw response recurrence over Z/(2**b), one final inverse.

    Deliberately does not use signed projective-gauge shifts: right shifts do
    not descend to fixed-precision residues without extra precision sidecars.
    """
    from collections import defaultdict
    if n<0 or b<0: raise ValueError('nonnegative n,b required')
    if b==0: return 0,dict(peak_states=0,transitions=0)
    mask=(1<<b)-1
    facts=[factorial(i) for i in range(n+1)]
    twos=[(f&-f).bit_length()-1 for f in facts]
    odds=[f>>v for f,v in zip(facts,twos)]
    width=max(1,n.bit_length()); digit=(1<<width)-1
    states={(0,0):1}; transitions=0; peak=1
    for ell in range(n if n%2 else n-1,2,-2):
        row=sum(gcd(ell,j)<<(width*i) for i,j in enumerate(range(1,ell,2)))
        idx=(ell-1)//2
        new=defaultdict(int)
        for (degree,response),weight in states.items():
            charge=(response>>(width*idx))&digit
            for k in range((n-degree)//ell+1):
                nd=degree+ell*k; rem=n-nd
                factor,rest=divmod(odds[nd],odds[degree]*ell**k*odds[k]); assert rest==0
                shift=k*(ell-1)//2+ell*k*(k-1)//2-twos[k]+k*charge
                assert shift>=0
                transitions+=1
                if shift>=b: continue
                val=((weight*factor)<<shift)&mask
                if not val: continue
                size=(min(ell-2,rem)+1)//2
                next_response=(response+k*row)&((1<<(width*size))-1)
                key=nd,next_response
                new[key]=(new[key]+val)&mask
        states={k:v for k,v in new.items() if v}; peak=max(peak,len(states))
    total=0
    for (degree,response),weight in states.items():
        rem=n-degree; charge=response&digit
        factor,rest=divmod(odds[n],odds[degree]*odds[rem]); assert rest==0
        shift=rem*(rem-1)//2+rem*charge-twos[rem]; assert shift>=0
        if shift<b: total=(total+((weight*factor)<<shift))&mask
    return total*pow(odds[n],-1,1<<b)&mask,dict(peak_states=peak,transitions=transitions)

if __name__=='__main__':
    sys.set_int_max_str_digits(0)
    import argparse,json,time
    p=argparse.ArgumentParser(); p.add_argument('n',type=int)
    p.add_argument('--support',type=int,default=3); args=p.parse_args()
    start=time.perf_counter(); x,stats=complete(args.n,args.support)
    stats['seconds']=time.perf_counter()-start
    print(json.dumps(dict(n=args.n,answer=str(x),stats=stats)))
