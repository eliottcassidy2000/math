from fractions import Fraction
def one(s):
    h,x,d,g=s
    return h+x+d,2*x-g,-g,g
def tri(s):
    h,x,d,g=s
    return 3*h+3*x+d,12*x-3*g,-g,3*g
for h,x,d,g in [(1,1,-1,1),(21,51,-9,19),(5365,13425,-153,459)]:
    s=(h,x,d,g)
    a=one(s)[0]-h;b=tri(s)[0]-3*h;c=one(one(s))[0]-h
    assert Fraction(b-a,2)==x and a-x==d and Fraction(b-c,2)==g
    v=s
    for j,coef in enumerate([(3,1,0),(45,3,-10),(567,9,-168)],1):
        v=tri(v)
        assert v[0]-3**j*h==sum(i*z for i,z in zip(coef,s[1:]))
assert 3*(3*(-168)+90)-(45*(-168)+5670)==648
print('Three normalized probes, exact recovery, and triangle-only rank certificate passed.')
