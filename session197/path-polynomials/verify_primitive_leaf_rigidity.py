#!/usr/bin/env python3
from verify_primitive_cancellation import *

def supported(n):
    v={}
    for p in (2,11,31,37,127):
        e=0
        while n%p==0:n//=p;e+=1
        v[p]=e
    return v,n

def run():
    residues=[pow(2,k,37) for k in range(1,37)]
    require(residues[-1]==1 and 1 not in residues[:-1],'exact order36')
    require(pow(2,18,37)==36 and pow(2,12,37)==26,'short exact order certificate')
    require(pow(2,6,9)==1,'nine-divisibility implication')
    accepted=[];cases=0;obstructions={'outside_support':0,'valuation_relation':0}
    for r in range(3,4097):
        n=2**r-1;d=2**(r-2)-1;g=gcd(n,d)
        require(g==(3 if r%2==0 else 1),'reduced Mersenne gcd')
        N,D=n//g,d//g
        require(N>4*D and N<=7*D and gcd(N,D)==1,'exact odd-ratio interval')
        if N%37==0:require(r%36==0 and N%3==0,'numerator37 forces unsupported3')
        if D%37==0:require((r-2)%36==0 and D%3==0,'denominator37 forces unsupported3')
        vn,rn=supported(N);vd,rd=supported(D)
        v={p:vn[p]-vd[p] for p in vn}
        ok=rn==rd==1 and v[127]==-v[31] and v[11]==-v[37]
        require(ok==(r==7),'variable-leaf exact support/valuation classification control')
        if ok:accepted.append(r)
        elif rn!=1 or rd!=1:obstructions['outside_support']+=1
        else:obstructions['valuation_relation']+=1
        cases+=1
    require(Fraction(22352,1147)/Fraction(352,37)==Fraction(127,62),'actual physical successful ratio')
    certificates=0
    for a,b,c in product(range(-5,6),repeat=3):
        # a,b,c are exponents of K,X,Y. Recover them from all valuations.
        v2=a+4*b+5*c;v11=b+c;v127=b;v31=-b;v37=-b-c
        B=v127;C=-v37-v127;A=v2-4*B-5*C
        require((A,B,C)==(a,b,c) and v127==-v31 and v11==-v37,'integer group membership certificate')
        certificates+=1
    require(127**2>7*31**2,'exponent at least2 exceeds interval')
    print(json.dumps({'status':'PASS','order37_residues':residues,'leaf_lengths_checked':cases,'accepted_leaf_lengths':accepted,'obstruction_controls':obstructions,'valuation_group_controls':certificates},sort_keys=True,indent=2))
if __name__=='__main__':run()
