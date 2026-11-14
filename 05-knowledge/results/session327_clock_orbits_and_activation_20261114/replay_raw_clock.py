from profile_tools import profile, strong
def chain(n):
    return tuple(sum(1<<j for j in range(n) if j==i+1 or j<=i-2) for i in range(n))
def values(z):
    out=[1,z,z*z]
    for n in range(3,13):
        out.append(z*out[-1]+out[-2]+z*out[-3])
    return out
a,b=values(1),values(2)
for n in range(3,13):
    T=chain(n)
    if not strong(T) or profile(T)!=(a[n],b[n]):
        raise RuntimeError('raw clock recurrence failed')
if [tuple(zip(a,b))[n] for n in (10,11,12)]!=[(193,11200),(355,29780),(653,79184)]:
    raise RuntimeError('raw clock values failed')
print('literal diameter chain: independent H,Q2 through order 12')
