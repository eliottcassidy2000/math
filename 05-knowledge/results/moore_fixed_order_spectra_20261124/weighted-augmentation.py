from math import factorial

def min_sq(odd, even, increments):
    total, i = odd, 0
    while increments:
        b = min(increments, even)
        total += (4+8*i)*b
        increments -= b
        b = min(increments, odd)
        total += (8+8*i)*b
        increments -= b
        i += 1
    return total

def test():
    n,k,r,q=3250,1625,22,70
    near,exceptions=0,[]
    for s in range(r,k+1):
        for t in range(max(s,k+q-s+1),k+1):
            d=r*(s+t-k-q)
            z=115*s*t-2*n*d
            if z>=0 and z*z>=225*s*t*(n-s)*(n-t):
                continue
            near+=1
            c=d-1
            for x,y in ((s,t),(t,s)):
                rem=n-x
                if 57*y<c or (115*y-2*c)*x-2*c*rem<=2*x*rem:
                    continue
                if min_sq(0,x,c)+min_sq(y,rem-y,57*y-c)>225*y+4*y*y:
                    break
            else:
                exceptions.append((s,t,d))
    assert near==246 and exceptions==[(228,1625,3476)]
    x,y,c=1625,228,3474
    rem=n-x
    assert (115*y-2*c)*x-2*c*rem>2*x*rem
    lhs=min_sq(0,x,c)+min_sq(y,rem-y,57*y-c)
    assert lhs==225*y+4*y*y==259236
    assert 4*1625+2*224==2*c
    assert 11*228+12*1397==115*y-2*c
    assert 4**2*1401+6**2*224+11**2*228+12**2*1397==lhs
    S=224
    increments=(115*S-4*n-S)//8
    assert 115*S==4*n+S+8*increments and increments==1567<=n-S
    assert 4**2*(n-S)+5**2*S+(12**2-4**2)*increments==254592
    assert 225*S+4*S*S==251104<254592
    m,b,E=k-q,k-q-238,52929
    assert 2*E<b*b
    num=(21**(21*m)*69**(69*2*q)*factorial(m)*b**E)**22
    den=(22**(20*m)*70**(68*2*q)*factorial(q)**2*(b+1)**E)**22*(2*m)
    assert (den<<(18927*22))<num<(den<<(18928*22))
    print('246 near-cuts, a single unit-defect size, and uniform exponent 18927')

if __name__=='__main__':
    test()
