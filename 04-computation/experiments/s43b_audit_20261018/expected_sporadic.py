from math import comb, isqrt
FAM = {'A1': lambda k: comb(2*k,k)**3, 'A2': lambda k: comb(2*k,k)**2*comb(3*k,k),
       'A3': lambda k: comb(2*k,k)**2*comb(4*k,2*k), 'A4': lambda k: comb(2*k,k)*comb(3*k,k)*comb(6*k,3*k)}
CMP = {'A1': [-512,-64,-8,1,16,64,256,4096], 'A2': [-8640,-1728,-192,-27,8,64,108,216,1458,3375],
       'A3': [-82944,-12288,-3969,-1024,-144,81,256,648,2304,20736], 'A4': [-32768,-3375,1728,8000,54000]}
primes = [p for p in range(7, 301) if all(p % d for d in range(2, isqrt(p)+1))]
E = 0.0; N = 0
for f, Af in FAM.items():
    for p in primes:
        A = [Af(k) % p for k in range(p)]
        for m in CMP[f]:
            if (6*m) % p == 0: continue
            w = pow(m % p, -1, p); s = 0; pw = 1
            for k in range(p):
                s = (s + A[k]*pw) % p; pw = pw*w % p
            if s == 0:
                E += 1.0/p; N += 1
print('supersingular instances', N, 'expected sporadic mod p^3 coincidences sum 1/p =', round(E, 2))
