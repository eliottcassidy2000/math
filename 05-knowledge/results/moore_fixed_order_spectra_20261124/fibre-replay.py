from math import factorial

def min_sq(odd, even, increments):
    result, i = odd, 0
    while increments:
        b = min(increments, even)
        result += (4+8*i)*b
        increments -= b
        b = min(increments, odd)
        result += (8+8*i)*b
        increments -= b
        i += 1
    return result

def verify_factor(r=22, q=71, k=1625):
    n, spectral, integer = 2*k, 0, []
    for s in range(r, k+1):
        for t in range(max(s,k+q-s+1),k+1):
            demand = r*(s+t-k-q)
            z = 115*s*t-2*n*demand
            if z>=0 and z*z>=225*s*t*(n-s)*(n-t):
                spectral += 1
                continue
            cap = demand-1
            for X,Y in ((s,t),(t,s)):
                rest = n-X
                if 57*Y < cap or (115*Y-2*cap)*X - 2*cap*rest <= 2*X*rest:
                    continue
                margin = min_sq(0,X,cap)+min_sq(Y,rest-Y,57*Y-cap)-225*Y-4*Y*Y
                if margin>0:
                    integer.append(margin)
                    break
            else:
                raise AssertionError((s,t))
    assert (spectral,len(integer),min(integer))==(604424,82,44)
    m=k-q
    b,E=m-238,52929
    num=(r-1)**((r-1)*(k+q)*11)*factorial(m)**11*b**(E*11)
    den=r**((r-2)*(k+q)*11)*factorial(r)**q*(b+1)**(E*11)
    assert den<<(18903*11)<num<den<<(18904*11)
    print(spectral,len(integer),min(integer),18903)

def height(w):
    return 1+sum(2 if j==i+1 else 2**(j-i-1)
                 for i in range(len(w)) for j in range(i+1,len(w))
                 if w[i]=='1' and w[j]=='0')

def binary(a):
    return format(a,'b') if a else ''

def flipped_words(a,b):
    A,B=binary(a),binary(b)
    refl=''.join(str(1-int(x)) for x in B[::-1])
    return A+'0'+refl,A+'1'+refl

def h(a):
    return height(binary(a))

def lower(a,b):
    return 2*a*b+h(a)+h(b)-1+2*((min(a,b)+1)//2)

def solve(H,d):
    if d<1 or H<1:
        return []
    # Solutions up to swapping the two flanks / taking converse.
    out=[(0,a) for a in (2*d-1,2*d) if h(a)==H]
    for eps in (0,1):
        for eta in (0,1):
            lo,hi=1,H+1
            while lo<hi:
                u=(lo+hi)//2
                a,b=2*u-eps,2*(u+d)-eta
                if lower(a,b)<H:
                    lo=u+1
                else:
                    hi=u
            a,b=2*lo-eps,2*(lo+d)-eta
            if lower(a,b)==H:
                out.append((a,b))
    return sorted(out)

if __name__=='__main__':
    verify_factor()
    for a in range(256):
        for b in range(a,256):
            u,v=((a+1)//2),((b+1)//2)
            H,d=lower(a,b),v-u
            w0,w1=flipped_words(a,b)
            assert (height(w0),height(w1))==(H,H+2*d)
            if d:
                sol=solve(H,d)
                assert (a,b) in sol
                assert all(lower(c,e)==H and (e+1)//2-(c+1)//2==d for c,e in sol)
    print('all 32896 nondecreasing flank pairs with both flanks below 256')
