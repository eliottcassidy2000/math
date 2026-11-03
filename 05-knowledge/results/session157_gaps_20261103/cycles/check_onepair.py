"""D=1 six-vertex classification with two directed triangles fixed."""
from itertools import combinations, permutations
from collections import defaultdict
E = list(combinations(range(6),2))
A = [(0,1),(1,2),(2,0),(3,4),(4,5),(5,3)]
fixed = {(min(i,j),max(i,j)): i<j for i,j in A}
free = [(i,j) for i,j in E if (i,j) not in fixed]
cycles = []
for size in (3,5):
    for subset in combinations(range(6),size):
        for perm in permutations(subset[1:]):
            cycles.append((subset[0],)+perm)
def signature(m):
    orient = dict(fixed)
    orient.update((e,bool(m & (1<<i))) for i,e in enumerate(free))
    def arc(a,b):
        return orient[min(a,b),max(a,b)] == (a<b)
    present = [x for x in cycles if all(arc(x[i],x[(i+1)%len(x)]) for i in range(len(x)))]
    triangles = [set(x) for x in present if len(x)==3]
    D = sum(not (a & b) for a,b in combinations(triangles,2))
    return len(triangles),len(present),D
def key(m):
    M = [[int(bool(m & (1<<free.index((i,j+3))))) for j in range(3)] for i in range(3)]
    choices=[]
    for a in range(3):
        for b in range(3):
            for swap in (False,True):
                for dual in (False,True):
                    Q = [[M[(i+a)%3][(j+b)%3] for j in range(3)] for i in range(3)]
                    if swap: Q = [[1-Q[j][i] for j in range(3)] for i in range(3)]
                    if dual: Q = [[1-Q[-i%3][-j%3] for j in range(3)] for i in range(3)]
                    choices.append(''.join(str(x) for row in Q for x in row))
    return min(choices)
groups = defaultdict(list)
for m in range(512): groups[key(m)].append(signature(m))
expected = {
'000000000':(2,2,1,2), '000000001':(4,6,1,18),
'000000011':(5,9,1,36), '000000111':(5,11,1,12),
'000001010':(6,12,1,18), '000001011':(6,13,1,36),
'000001110':(7,14,1,36), '000001111':(6,12,1,36),
'000011110':(7,16,1,36), '001010100':(8,20,1,6),
'001010101':(8,19,1,36),
}
actual={}
for k,values in groups.items():
    if len(set(values)) != 1: raise RuntimeError('invariant changed in orbit '+k)
    t,C,D=values[0]
    if D==1: actual[k]=(t,C,D,len(values))
if len(groups)!=22 or actual!=expected: raise RuntimeError('classification mismatch')
if any(C==8 for t,C,D,size in actual.values()): raise RuntimeError('counterexample')
print('PASS: 512 cross matrices, 22 orbits; 11 D=1 orbits match the table; C != 8.')
