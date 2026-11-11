"""Exact target-degree Burnside quotient, with mass-gauge response coalescence."""
from collections import defaultdict
from math import factorial, gcd

def count(n, *, gauge=True, switch=9, terminal3=True):
    if type(n) is not int or n < 0: raise ValueError('nonnegative literal integer required')
    facts=[factorial(i) for i in range(n+1)]
    twos=[(f&-f).bit_length()-1 for f in facts]
    odds=[f>>v for f,v in zip(facts,twos)]
    requested_gauge=gauge
    gauge=False
    width=max(1,n.bit_length()); digit=(1<<width)-1
    # Gauge coordinates b_j=(j*q_1-q_j)/2; b_1 vanishes identically.
    start=3 if gauge else 1
    masks=[(1<<(width*c))-1 for c in range(n//2+2)]
    states={(0,0):1}; rows=[]
    for ell in range(n if n%2 else n-1,2,-2):
        if requested_gauge and not gauge and ell<=switch:
            new=defaultdict(int)
            for (degree,response),weight in states.items():
                rem=n-degree; q1=response&digit; packed=0
                for j in range(3,min(ell,rem)+1,2):
                    qj=(response>>(width*((j-1)//2)))&digit
                    packed+=((j*q1-qj)//2)<<(width*((j-3)//2))
                new[degree,packed]+=weight << (rem*q1)
            states=dict(new); gauge=True; start=3
        if gauge and terminal3 and ell==3:
            total=0; tables={}; terms=0
            for (degree,response),weight in states.items():
                if degree not in tables:
                    R=n-degree; entries=[]
                    for k in range(R//3+1):
                        r=R-3*k
                        fac,rest=divmod(odds[n],odds[degree]*3**k*odds[k]*odds[r])
                        assert rest==0
                        base=R*(R-1)//2+3*k*k+(1-2*R)*k-twos[k]-twos[r]
                        entries.append((k,fac,base))
                    tables[degree]=entries
                charge=response&digit
                entries=tables[degree]
                low=min(base-2*k*charge for k,fac,base in entries)
                poly=sum(fac << (base-2*k*charge-low) for k,fac,base in entries)
                val=weight*poly
                if low>=0: val <<= low
                else:
                    assert val & ((1<<(-low))-1)==0
                    val >>= -low
                total+=val; terms+=len(entries)
            answer,rest=divmod(total,odds[n]); assert rest==0
            return answer,dict(rows=rows,transitions=sum(x['transitions'] for x in rows),
                              terminal_terms=terms,peak_states=max([1]+[x['states'] for x in rows]))
        idx=(ell-start)//2
        row=sum(((j-gcd(ell,j))//2 if gauge else gcd(ell,j))<<(width*i)
                for i,j in enumerate(range(start,ell,2)))
        tables={}
        for degree,_ in states:
            if degree in tables: continue
            entries=[]
            for k in range((n-degree)//ell+1):
                nd=degree+ell*k; rem=n-nd
                h=ell**k*odds[k]
                factor,remainder=divmod(odds[nd],odds[degree]*h)
                assert remainder==0
                base=k*(ell-1)//2+ell*k*(k-1)//2-twos[k]
                if gauge: base+=k*rem
                size=max(0,(min(ell-2,rem)-start)//2+1)
                entries.append((k,nd,base,factor,masks[size],k*row))
            tables[degree]=entries
        new=defaultdict(int); transitions=0
        for (degree,response),weight in states.items():
            charge=(response>>(width*idx))&digit
            for k,nd,base,factor,mask,offset in tables[degree]:
                shift=base-2*k*charge if gauge else base+k*charge
                val=weight*factor
                if shift>=0: val <<= shift
                else:
                    # The gauge can remove powers of two; exact by the proof.
                    assert val & ((1<<(-shift))-1)==0
                    val >>= -shift
                new[nd,(response+offset)&mask]+=val
                transitions+=1
        states=dict(new)
        rows.append(dict(cycle_length=ell,states=len(states),transitions=transitions))
    total=0
    for (degree,response),weight in states.items():
        rem=n-degree
        charge=0 if gauge else response&digit
        factor,rest=divmod(odds[n],odds[degree]*odds[rem]); assert rest==0
        total+=(weight*factor) << (rem*(rem-1)//2+rem*charge-twos[rem])
    answer,rest=divmod(total,odds[n]); assert rest==0
    return answer,dict(rows=rows,transitions=sum(x['transitions'] for x in rows),peak_states=max([1]+[x['states'] for x in rows]))
