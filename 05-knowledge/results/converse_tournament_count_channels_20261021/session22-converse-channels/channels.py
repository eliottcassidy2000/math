"""Degree-local transfer for F_n(a,b), with signed core/restricted options."""
from collections import defaultdict
from fractions import Fraction
from math import factorial,gcd,comb

def channel_dp(bound, a=1, b=0, *, core=False,no_fixed=False):
    if type(bound) is not int or bound<0: raise ValueError('nonnegative literal integer bound required')
    if type(core) is not bool or type(no_fixed) is not bool: raise ValueError('literal Boolean modes required')
    if type(a) is not int or a<1 or type(b) is not int or b<0:
        raise ValueError("literal integers a>=1, b>=0 required")
    start=3 if no_fixed else 1
    width=max(1,bound.bit_length()); digit=(1<<width)-1
    masks=[(1<<(width*c))-1 for c in range(bound//2+2)]
    facts=[factorial(i) for i in range(bound+1)]
    twos=[(f&-f).bit_length()-1 for f in facts]
    odds=[f>>v for f,v in zip(facts,twos)]
    states={(0,0):1}; rows=[]; precomputed=0
    for ell in range(bound if bound%2 else bound-1,start-1,-2):
        idx=(ell-start)//2
        row=sum(gcd(ell,j)<<(width*i) for i,j in enumerate(range(start,ell,2)))
        tables={}
        # Transfer only depends on degree and multiplicities, not response.
        for degree,_ in states:
            if degree in tables: continue
            entries=[]
            for k in range((bound-degree)//ell+1):
                h=ell**k*odds[k]
                base=a*(k*(ell-1)//2+ell*k*(k-1)//2)-twos[k]+b*k
                for excess in (range((bound-degree)//ell-k+1) if core and k else (0,)):
                    nd=degree+ell*(k+excess)
                    factor,remainder=divmod(odds[nd],odds[degree]*h)
                    if remainder: raise ArithmeticError('transfer nonintegral')
                    if core and excess: factor*=(-1)**excess*comb(k+excess-1,excess)
                    count=max(0,(min(ell-2,bound-nd)-start)//2+1)
                    entries.append((k,nd,base,factor,masks[count],k*row))
                    precomputed+=1
            tables[degree]=entries
        new=defaultdict(int); transitions=0
        for (degree,response),weight in states.items():
            charge=(response>>(width*idx))&digit
            for k,nd,base,factor,mask,offset in tables[degree]:
                new[nd,(response+offset)&mask]+=(weight*factor)<<(base+a*k*charge)
                transitions+=1
        states={key:value for key,value in new.items() if value}
        rows.append({'cycle_length':ell,'states':len(states),'transitions':transitions})
    ans=[Fraction(0)]*(bound+1)
    for (d,r),w in states.items():
        if r: raise ArithmeticError('terminal response')
        ans[d]+=Fraction(w,odds[d])
    return ans,{'rows':rows,'transitions':sum(r['transitions'] for r in rows),'precomputed_transfers':precomputed,'peak_states':max([1]+[r['states'] for r in rows])}


def counts(bound):
    """All/SC/converse-merged and strong variants, including empty order zero.

    Singletons are strong; the empty tournament is not strong.
    """
    ordinary, stats = channel_dp(bound)
    even, _ = channel_dp(bound//2, 2, 0)
    odd, _ = channel_dp(max(0,(bound-1)//2), 2, 1)
    sc = [even[n//2] if n%2==0 else odd[n//2] for n in range(bound+1)]
    strong = [Fraction(0)]*(bound+1)
    strong_sc = [Fraction(0)]*(bound+1)
    for n in range(1,bound+1):
        strong[n] = ordinary[n]-sum(strong[k]*ordinary[n-k] for k in range(1,n))
        strong_sc[n] = sc[n]-(ordinary[n//2] if n%2==0 else 0)-sum(ordinary[j]*strong_sc[n-2*j] for j in range(1,(n-1)//2+1))
    result = {'ordinary':ordinary,'self_converse':sc,
              'merged':[(x+y)/2 for x,y in zip(ordinary,sc)],
              'strong':strong,'strong_self_converse':strong_sc,
              'strong_merged':[(x+y)/2 for x,y in zip(strong,strong_sc)]}
    for name, values in result.items():
        if any(v.denominator!=1 or v<0 for v in values):
            raise ArithmeticError('nonintegral or negative count: '+name)
        result[name] = [int(v) for v in values]
    return result
