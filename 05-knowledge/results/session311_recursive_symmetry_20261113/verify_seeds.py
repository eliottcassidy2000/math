from itertools import permutations
A=(2,28,25,81,97,79,7)
B=(2,80,3,7,173,143,189,15)
def iso_count(a,b):
    n=len(a)
    return sum(all((a[i]>>j&1)==(b[p[i]]>>p[j]&1)
                   for i in range(n) for j in range(i+1,n))
               for p in permutations(range(n)))
def opposite(a):
    return tuple(((1<<len(a))-1)^(1<<i)^r for i,r in enumerate(a))
if (iso_count(A,A),iso_count(A,opposite(A)),iso_count(B,B),iso_count(B,opposite(B)))!=(1,1,1,0):
    raise ValueError('seed symmetry failed')
print('PASS: both seeds rigid, A self-converse, B not self-converse')
