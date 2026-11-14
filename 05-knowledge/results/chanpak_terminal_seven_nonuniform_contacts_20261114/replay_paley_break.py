from orbital import rows7, prime, induced_table, triangle_modules

def attached(mask_v,mask_w):
    u=rows7((1,2,4))
    return tuple(a|((not(mask_v>>j&1))<<7)|((not(mask_w>>j&1))<<8)
                 for j,a in enumerate(u))+(mask_v|256,mask_w)

for v,w,H,F in ((96,115,2621,41140),(36,63,2471,38668)):
    rows=attached(v,w)
    T=induced_table(rows)
    assert (T[-1],sum(a*b for a,b in zip(T,T[::-1])))==(H,F)
assert prime(2621) and prime(353) and 2471==7*353
assert prime(326479) and 5171*(1<<6)-4465==326479
# All 128 apex contact subsets on the base orbit, with no other vertices:
for mask in range(128):
    u=rows7((1,2,4))
    w=tuple(v|((not(mask>>j&1))<<7) for j,v in enumerate(u))+(mask,)
    k=mask.bit_count()
    H=induced_table(w)[-1]
    assert H%3==0
    if mask==0: A=H
    if mask==1: B=H
    # For this bare example both extreme counts are 189.
    assert H==189+((B-189)//6)*k*(7-k) if mask else H==189
print('one-contact divisibility and two prime-lane witnesses checked')
from orbital import color_poly_value
T=induced_table(attached(36,63))
assert tuple(z%7 for z in color_poly_value(T,3))==(0,0,1)
for r in range(1,32):
    from collections import Counter
    assert max(Counter(z.bit_count() for z in triangle_modules(list(attached(36,63)),r)).values())<7
