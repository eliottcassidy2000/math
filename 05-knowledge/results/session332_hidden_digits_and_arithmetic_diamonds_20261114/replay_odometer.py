from itertools import product

F=((2,3**7*5**31),(5,3**7),(3,1))
FC=((1,3**7*5**31),)+F[1:]

def val(x,p):
    j=0
    while x%p==0:
        x//=p;j+=1
    return j

def unpack(x):
    h,a,b=(val(x,p) for p in (2,3,5))
    assert x==2**h*3**a*5**b
    return h,a,b

def step(x,program=F):
    for a,b in program:
        if x%b==0:return a*(x//b)

def pack(k):
    h,k=divmod(k,256)
    b,a=divmod(k,8)
    return 2**h*3**a*5**b

def index(x):
    q,a,b=unpack(x)
    return 256*q+8*b+a
for k in range(4096):
    assert step(pack(k))==pack(k+1)
    assert step(pack(k%256),FC)==pack((k+1)%256)
for a,b in product(range(40),range(40)):
    x=3**a*5**b
    for k in range(50):
        y=step(x)
        assert index(y)==index(x)+1
        x=y
        _,a,b=unpack(x)
        if a<8 and b<32:
            break
    else:raise AssertionError('normalisation exceeded test bound')
for k in range(256):
    _,a,b=unpack(pack(k)); h=3**a*5**b; q=12**a*32**b
    assert val(h,3)==a and val(h,5)==b and val(q,2)==2*a+5*b
assert [pack(k) for k in range(8,12)]==[5,15,45,135]
assert max(pow(a,128,512) for a in range(1,512,2))==1
assert pow(5,64,512)!=1
print('three-fraction 3+5-bit odometer, self-normalisation, and 512 contrast: PASS')
