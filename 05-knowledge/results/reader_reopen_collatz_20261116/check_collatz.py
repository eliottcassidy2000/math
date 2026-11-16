from verify import profile
def oddpart(x):
    return x//(x&-x)
def sel(h):
    return 1 if h%4==1 else -1
def step(h,sig):
    return oddpart(3*h+sig)
def hword(w):
    A=p=0;h=1
    for c in map(int,w):
        if not c:h+=A+p
        A=2*A+c;p=c
    return h
seed=(4,5,352,263,271,347,411,319,3)
h,g,R,S=profile(seed,8)
assert (h,g,R,S)==(181,97,161,188)
assert (R+S-g,h-2*S+2*g)==(252,-1)
for h in range(3,4000,2):
    y=step(h,sel(h))
    assert y<h and y%3
for e in range(6,181):
    sig=(-1)**e
    h=(7*2**e-sig)//3
    assert sel(h)==sig and step(h,sig)==7
    u=step(h,-sig);w=step(u,-sig);z=step(w,-sig)
    assert (u,w,z)==(7*2**(e-1)-sig,21*2**(e-3)-sig,63*2**(e-5)-sig)
    assert all(x not in (7,21) for x in (u,w,z)) and z<h
    assert 27*h==32*z+23*sig
    if sig==-1:assert hword('1001'+'0'*(e-5)+'100')==z
for path,signs in [
    ([9,13,5],[-1,1]),([19,29,11],[1,1]),
    ([37,55,41,31],[-1,-1,1]),([75,113,85,1],[1,1,1])]:
    assert all(step(v,s)==w for v,w,s in zip(path,path[1:],signs))
assert [hword(w) for w in ('100','1100','101010')]==[5,11,31]
print('Signed 7-bypass and exact endpoint readers verified.')
