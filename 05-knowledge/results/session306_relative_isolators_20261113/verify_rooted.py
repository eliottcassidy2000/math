from itertools import permutations,product
from collections import Counter

def core(n, a):
    p=[(i,j) for i in range(n) for j in range(i+1,n)]
    def edge(i,j):
        if i==j:return 0
        if i>j:return 1-edge(j,i)
        return (a>>p.index((i,j)))&1
    return edge

def roots(n, edge, port):
    return [s for s in permutations(range(n)) if all(s[x]==x for x in port)
            and all(edge(i,j)==edge(s[i],s[j]) for i in range(n) for j in range(n))]
def nc(s):
    todo=set(range(len(s))); ans=0
    while todo:
        x=min(todo); ans+=1
        while x in todo:
            todo.remove(x); x=s[x]
    return ans

def orbit_count(n, g):
    return len({min(tuple(bits[s[i]] for i in range(n)) for s in g)
                for bits in product((0,1),repeat=n)})

checks=0
for n in range(1,5):
    for a in range(1<<(n*(n-1)//2)):
        edge=core(n,a)
        for port in [(),(0,)]:
            g=roots(n,edge,port)
            if orbit_count(n,g)*len(g)!=sum(2**nc(s) for s in g):
                raise ValueError('Burnside')
            checks+=1
def unmarked(edge,x):
    def arc(i,j):
        if i==j:return 0
        if i==2:return x[j]
        if j==2:return 1-x[i]
        return edge(i,j)
    return min(tuple(arc(s[i],s[j]) for i,j in [(0,1),(0,2),(1,2)])
               for s in permutations(range(3)))
e=core(2,1)
b=Counter(unmarked(e,x) for x in product((0,1),repeat=2))
if sorted(b.values())!=[1,3]:raise ValueError('unmarked hostile')
print(f'PASS {checks} Burnside checks and TT2 attachment collapse')
