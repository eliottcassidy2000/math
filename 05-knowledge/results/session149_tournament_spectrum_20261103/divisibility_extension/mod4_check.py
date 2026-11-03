import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'divisibility'))
from partition_count import count
from math import gcd

def phi(n):
    return sum(gcd(i,n)==1 for i in range(1,n))

for n in range(8,81,2):
    a,partitions,_=count(n)
    q,r=divmod(a,1<<(n//2))
    if r: raise ArithmeticError(f'n={n} not divisible')
    predicted=((1 if n%4==2 else -1)*(phi(n)//2))%4
    if q%4!=predicted: raise AssertionError(f'n={n}, got {q%4}, expected {predicted}')
    print(n,partitions,(a&-a).bit_length()-1,q%4,predicted,flush=True)
