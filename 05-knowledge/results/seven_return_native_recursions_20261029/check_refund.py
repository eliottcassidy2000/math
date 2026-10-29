"""Exhaustive small completion bounds, independently using factorials."""
from collections import Counter
from math import factorial
import json
def val(n):return (n&-n).bit_length()-1
def parts(n,low=1):
    if n==0:yield ();return
    for a in range(low,n+1,2):
        for tail in parts(n-a,a):yield (a,)+tail
def verify():
    checks=partitions=0;strict=0
    for n in range(1,31):
        for word in parts(n):
            partitions+=1;k=len(word);actual=sum(val(factorial(c)) for c in Counter(word).values())
            for length in range(k+1):
                alpha=word[:length];t=k-length
                if not alpha:basic=refined=val(factorial(k))
                else:
                    counts=Counter(alpha);a=alpha[-1];m=counts[a]
                    frozen=sum(val(factorial(c)) for x,c in counts.items() if x!=a)
                    basic=frozen+val(factorial(m+t));remaining=sum(word[length:]);excess=remaining-t*a
                    if excess:
                        cap=min(t,excess//2)
                        refined=frozen+max(val(factorial(m+t-v))+val(factorial(v)) for v in range(1,cap+1))
                    else:refined=basic
                if not actual<=refined<=basic<=val(factorial(k)):
                    raise ValueError((word,length,actual,refined,basic))
                checks+=1;strict+=refined<basic
    return dict(status='PASS',orders=[1,30],partitions=partitions,prefix_checks=checks,strict_sum_refinements=strict)
if __name__=='__main__':print(json.dumps(verify(),indent=2))
