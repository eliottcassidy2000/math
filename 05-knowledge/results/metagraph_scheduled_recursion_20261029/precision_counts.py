"""Exact normalized low bits of positive tournament Burnside channels.

Counts F_n(r,s), r>=1 independent tournament layers and s binary attributes.
The precision filter bounds the number of odd permutation cycles, not the
size of an individual cycle. Composite lengths and factorial powers survive.
"""
from collections import Counter
from math import factorial, gcd


def v2(n):
    if n <= 0:
        raise ValueError('positive input required')
    return (n & -n).bit_length()-1


def base_valuation(n,r,s):
    if n==0:
        return 0
    if n==2:
        return r+2*s-1
    return r*(n-1)//2+s if n%2 else r*n//2+2*s


def cycle_floor(k,r,s,odd):
    """Monotone valid lower bound for k>=3 odd or k>=4 even."""
    if odd:
        return r*(k-1)**2//2+s*(k-1)-(k-2)
    return r*k*(k-2)//2+s*(k-2)-(k-1)


def odd_parts(total,count,lower=1):
    """Each nondecreasing odd count-part partition occurs once."""
    if count==1:
        if total>=lower and total%2:
            yield (total,)
        return
    for a in range(lower,total//count+1,2):
        for tail in odd_parts(total-a,count-1,a):
            yield (a,)+tail


def term(parts,r,s,odd):
    k=len(parts)
    pair_gcd=sum(gcd(a,b) for i,a in enumerate(parts) for b in parts[i+1:])
    denominator=1
    for a,m in Counter(parts).items():
        denominator *= a**m*factorial(m)
    twos=v2(denominator)
    exponent=(r*(pair_gcd-(k-1)//2)+s*(k-1)-twos if odd else
              r*(pair_gcd-k//2)+s*(k-2)-twos)
    return exponent,denominator>>twos


def normalized_low_bits(n,bits,r=1,s=0,shortcut=True):
    if any(type(x) is not int for x in (n,bits,r,s)) or n<0 or bits<1 or r<1 or s<0:
        raise ValueError('integer n>=0, bits>=1, r>=1 and s>=0 required')
    modulus=1<<bits
    base=base_valuation(n,r,s)
    if n in (0,1,2):
        return 1,dict(base=base,cycle_cap=n,visited=1,retained=1,shortcut='boundary')
    delta=2*r+2*s-1
    if shortcut and n%2 and bits<=delta+1:
        return (pow(n,-1,modulus)+(1<<delta))%modulus,dict(
            base=base,cycle_cap=3,visited=0,retained=0,shortcut='sharp_first_correction')
    counts=[1] if n%2 else [2]
    k=counts[0]+2
    while k<=n and cycle_floor(k,r,s,bool(n%2))<bits:
        counts.append(k)
        k+=2
    total=0
    visited=retained=0
    by_cycles={}
    for k in counts:
        vk=rk=0
        for parts in odd_parts(n,k):
            vk+=1
            exponent,odd_denominator=term(parts,r,s,bool(n%2))
            if exponent<0:
                raise ArithmeticError('negative normalized valuation outside declared boundary')
            if exponent<bits:
                total += (1<<exponent)*pow(odd_denominator,-1,modulus)
                rk+=1
        by_cycles[k]=dict(visited=vk,retained=rk)
        visited+=vk;retained+=rk
    return total%modulus,dict(base=base,cycle_cap=counts[-1],visited=visited,
                              retained=retained,by_cycles=by_cycles,shortcut=None)


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('n',type=int)
    parser.add_argument('bits',type=int)
    parser.add_argument('--layers',type=int,default=1)
    parser.add_argument('--attributes',type=int,default=0)
    args=parser.parse_args()
    value,stats=normalized_low_bits(args.n,args.bits,args.layers,args.attributes)
    print(json.dumps(dict(n=args.n,bits=args.bits,layers=args.layers,
                         attributes=args.attributes,residue=value,stats=stats)))
