from itertools import product
P='133220032220'; Q='311103311002'
U='1332202101012022331';V='022230022331';W='200113301113';Z='2001131232321311002'
def rc(w):return ''.join(str(3-int(z)) for z in reversed(w))
def row(w):return int(''.join(str(int(z)//2) for z in w),2)
def amp(w):
    assert row(w)%13==0 and row(rc(w))%13==0
    return row(w)//13%120,row(rc(w))//13%120
assert pow(2,12,195)==1
assert all(row(g)%195==0 for g in [V,W,rc(V),rc(W)])
assert [amp(w) for w in [P,Q,U+Z]]==[(30,45),(45,30),(69,69)]
def width(n,d):
    m=n-12*abs(d)
    if m<0 or m%2:return None
    b=(pow(19,-1,12)*(m//2))%12
    if 38*b>m:return None
    return (m-38*b)//24,b
def build(n,d):
    k,b=width(n,d)
    return (P if d>=0 else Q)*abs(d)+(P+Q)*k+(U+Z)*b
assert width(24,1) is None and width(120,1) is None
for d in range(-10,11):
    assert width(12*abs(d)+394,d) is None
    for n in range(12*abs(d)+396,12*abs(d)+516,2):
        assert len(build(n,d))==n
    for n in range(12*abs(d)+444,12*abs(d)+564,2):
        found=set()
        for p,dp in [(P,1),(Q,-1)]:
            for q,dq in [(P,1),(Q,-1)]:
                tape=p+build(n-24,d-dp-dq)+q
                found.add(amp(tape))
        a=[30,45] if n%4==0 else [54,69]
        assert found==set(product(a,a))
print('Exact all-even conductor, modulus-120 width wheel, and independent marked phases pass.')
