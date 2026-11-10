# THM-5109 (D) finite part: every two-step sheet-compatible chain has t' <= 27 (Ellison + Herschfeld).
# Enumerate all (t, a) with 3^a < 2^t, t <= 16, and all a' with a+t+1 <= a' and 3^a' < 2^27;
# report every (t, a, a', t') with h' + 3^a' = 2^t', where h = 2^t - 3^a, h' = h + 2^(2a+1).
ispow=lambda r: r>0 and r&(r-1)==0
sols=[]; checked=0
for t in range(1,40):
    a=0
    while 3**a < 2**t:
        h=2**t-3**a; hp=h+2**(2*a+1)
        ap=a+t+1
        while 3**ap < 2**60:
            checked+=1
            s=hp+3**ap
            if ispow(s): sols.append((t,a,h,hp,ap,s.bit_length()-1))
            ap+=1
        a+=1
print('pairs checked',checked,'solutions',sols)
# also: one-step solutions of the law h + 3^a = 2^t for small h (the only possible dyadic one-step values)
