from check_27_63 import profile, scc_sizes
a=(62,124,120,112,32,64,17)
assert profile(a)[-1]==27
for h,u,v,m in [(63,3,6,63),(189,0,2,24),(567,0,6,288),(1701,4,8,304),(5103,3,6,16),(15309,0,1,56),(45927,0,1,8056),(137781,3,6,2576)]:
    b=list(a);b[u]^=1<<v;b[v]^=1<<u
    a=tuple(z|(((m>>j)&1)^1)<<len(b) for j,z in enumerate(b))+(m,)
    assert scc_sizes(a)==(len(a),) and profile(a)[-1]==h
print('Connected strong rewrites verified from 27 through 137781.')
def ordinal(a,b):
    n=len(a)
    return tuple(z|(((1<<len(b))-1)<<n) for z in a)+tuple(z<<n for z in b)
def grow(a,x):
    n=len(a)
    return tuple(z|(((x>>j)&1)^1)<<n for j,z in enumerate(a))+(x,)
tri=(2,4,1)
for a in [(2,32,11,19,7,29),(16,5,9,3,46,15)]:
    assert scc_sizes(a)==(6,) and profile(a)[-1]==27
    for b in [ordinal(a,tri),ordinal(tri,a)]:
        assert profile(b)[-1]==81
        assert all(profile(grow(b,x))[-1]!=189 for x in range(512))
print('H189 residual excluded by 2048 literal orientation checks.')
