from math import isqrt
from slot_profiles import profile,append,ordinal,external,strong

A=(158,164,56,242,98,1,167,52)
B0=(254,240,162,54,4,144,60,88)
B1=(142,136,18,20,131,31,191,44)
choices_a=((A,144),(A,6))
choices_b=((B0,61),(B1,79))
pa=[profile(*z) for z in choices_a]
pb=[profile(*z) for z in choices_b]
if pa!=[(389,5600,1168,1376,102,88,217),(389,5600,1168,1376,102,88,225)]:raise ValueError(pa)
if pb!=[(97,2560,1867,1715,69,97,-1093),(97,2560,1867,1715,69,97,-1095)]:raise ValueError(pb)
for i in range(2):
    for j in range(2):
        va,qa=choices_a[i];vb,qb=choices_b[j]
        r=ordinal(va,vb);q=qa|(qb<<len(va))
        if not strong(external(r,q,1)):raise ValueError('not strong')
        p=append(append((1,0,1,0),pa[i]),pb[j])
        if p!=(37733,3716343,-2454875+776*i-778*j,26841):raise ValueError(p)
        print(i,j,p)
def is_prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
for d in (0,1):
    n=3716343*(2**11)-2454875-2*d
    if not is_prime(n):raise ValueError(n)
    print('prime',n)
