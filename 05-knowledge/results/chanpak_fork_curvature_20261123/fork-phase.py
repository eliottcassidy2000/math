"""Physical ordinal recipes for two-read fork phases; standard library only."""
from math import gcd

def ord2(m):
    if m < 3 or gcd(2, m) != 1:
        raise ValueError("odd modulus >=3 required")
    x = 2 % m
    a = 1
    while x != 1:
        a += 1
        x = 2*x % m
    return a

def totient(m):
    return sum(gcd(m, a) == 1 for a in range(1, m))

def family(k, m=None):
    if m is None:
        a,b = 2**k,3**k
        return a+1,4*a+2*b-2
    a,b = pow(2,k,m),pow(3,k,m)
    return (a+1)%m,(4*a+2*b-2)%m

def exit_reads(d,q,m=None):
    c=21*d+840*q
    v=c-2*d,c,c+2*d
    return v if m is None else tuple(x%m for x in v)

def inverses(m):
    """At primitive 2 and m coprime to 6, table of logs of all units."""
    e=ord2(m)
    if gcd(m,6)!=1 or e!=totient(m):
        raise ValueError("need 2 to generate the units and modulus coprime to 6")
    return e,{pow(2,a,m):a for a in range(e)}

def recipe(m,d,q,seed=(1,1)):
    """Ordinal sum: seed, a near-transitive blocks (order e+2), b singleton
    pads. Default (1,1) is an empty seed. Add e pads if recipe would be empty.
    """
    e,logs=inverses(m)
    d0,q0=seed
    a=logs[d*pow(d0,-1,m)%m]
    b=logs[q*pow(q0,-1,m)*pow(4,-a,m)%m]
    return e,a,b

BANK=((5,5),(19,8),(217,12),(223,12),(599,13))

def exact_bank(target, a=40, b=1):
    """All representations target = a*2^e + b*D by the bank, a,b>0."""
    for e in range(target.bit_length()):
        remainder=target-a*(1<<e)
        if remainder<1: break
        d,r=divmod(remainder,b)
        if r or not d: continue
        exps=[]
        rest=d
        for base, _ in BANK:
            power=0
            while rest%base==0:
                rest//=base
                power+=1
            exps.append(power)
        if rest!=1: continue
        pad=e-sum(a*b[1] for a,b in zip(exps,BANK))
        if pad>=0 and (sum(exps)+pad)>0:
            yield d,e,tuple(exps),pad

def audit():
    assert list(exact_bank(40*4096+217)) == [(217,12,(0,0,1,0,0),0)]
    assert list(exact_bank(40*4096+223)) == [(223,12,(0,0,0,1,0),0)]
    assert list(exact_bank(40*8192+599)) == [(599,13,(0,0,0,0,1),0)]
    assert list(exact_bank(40*16+1)) == [(1,4,(0,0,0,0,0),4)]
    for k in (-1,0,1):
        target=840*4096+(21+2*k)*223
        assert (223,12,(0,0,0,1,0),0) in exact_bank(target,840,21+2*k)
    for p in (11,13,19,29,37,59,61):
        e,logs=inverses(p)
        assert family(e,p)==(2,4)
        axes={}
        for k in range(1,e+1):
            d,q=family(k,p)
            if (d==0)!=(q==0):
                axes.setdefault((d==0,q==0),k)
        assert set(axes)=={(True,False),(False,True)}
        for d in range(p):
            for q in range(p):
                if d and q:
                    seeds=(1,1)
                elif not d and not q:
                    sd,sq=family(axes[True,False],p)
                    td,tq=family(axes[False,True],p)
                    seeds=(sd*td%p,sq*tq%p)
                else:
                    seeds=family(axes[d==0,q==0],p)
                d0,q0=seeds
                _,a,b=recipe(p,d or 1,q or 1,(d0 or 1,q0 or 1))
                actual=(d0*pow(2,a,p)%p,q0*pow(4,a,p)*pow(2,b,p)%p)
                assert actual==(d,q)
        print("full residue plane",p,axes)
    for m in (5,6,7,8):
        ways={
            5:lambda d:(-d,d,3*d),
            6:lambda d:(d,3*d,5*d),
            7:lambda d:(-2*d,0,2*d),
            8:lambda d:(3*d,5*d,7*d),
        }
        for d in range(1,m*2,2):
            for q in range(0,m*2,2):
                assert exit_reads(d,q,m)==tuple(x%m for x in ways[m](d))
    for k in range(1,16):
        d,q=family(k)
        assert q%3==((-1)**(k+2)*d)%3
if __name__=="__main__":
    audit()
