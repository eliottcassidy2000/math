from itertools import permutations, combinations
from math import factorial
from fractions import Fraction
from partition_count import count

def direct(n):
    edges=list(combinations(range(n),2)); total=0
    for p in permutations(range(n)):
        seen=set(); types=[]
        for x in range(n):
            if x not in seen:
                i=x; size=0
                while i not in seen: seen.add(i);i=p[i];size+=1
                types.append(size)
        if any(k%2==0 for k in types):continue
        seen=set(); E=0
        for e in edges:
            if e not in seen:
                E+=1;f=e
                while f not in seen:
                    seen.add(f);f=tuple(sorted((p[f[0]],p[f[1]])))
        total+=1<<E
    if total % factorial(n) != 0:
        raise ArithmeticError("Nonintegral permutation Burnside sum")
    return total//factorial(n)
for n in range(1,9):
    a=direct(n); b=count(n)[0]
    if a != b:
        raise AssertionError(f"Mismatch at n={n}: {a} != {b}")
    print(n,a)
