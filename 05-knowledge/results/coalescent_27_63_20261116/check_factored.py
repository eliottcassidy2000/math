from itertools import combinations, permutations
from check_27_63 import profile, restriction
cores=((2,16,3,23,165,79,31,111),(16,5,17,7,74,31,175,63))
for a in cores:
    hist=[]
    for k in (3,5,7):
        cycles=[]
        for ss in combinations(range(8),k):
            cycles.extend(ss for tail in permutations(ss[1:])
                if all(a[x]>>y&1 for x,y in zip((ss[0],)+tail,tail+(ss[0],))))
        if k==3:
            assert all(set(s)&set(t) for s,t in combinations(cycles,2))
        hist.append(len(cycles))
    assert hist==[8,17,6] and profile(a)[-1]==63
b=restriction(cores[0],[6,3,0,1,4])
assert b==(30,28,8,16,4)
assert profile(b)[-1]==3
print('Both packing polynomials and the 3-new-vertex two-edit pattern verified.')
