"""Since b-a = 1/2 for 2F1(1/4,3/4;1;z), a quadratic transformation exists.
Test candidates numerically to find the exact one."""
from mpmath import mp, mpf, hyp2f1, sqrt, ellipk, pi, nstr
mp.dps = 40

def F(z): return hyp2f1(mpf(1)/4, mpf(3)/4, 1, z)

tests = {
 "(1+sqrt z)^{-1/2} 2F1(1/2,1/2;1; 2 sqrt z/(1+sqrt z))":
    lambda z: (1+sqrt(z))**mpf(-0.5)*hyp2f1(mpf(1)/2, mpf(1)/2, 1, 2*sqrt(z)/(1+sqrt(z))),
 "(1-sqrt z)^{-1/2} 2F1(1/2,1/2;1; -2 sqrt z/(1-sqrt z))":
    lambda z: (1-sqrt(z))**mpf(-0.5)*hyp2f1(mpf(1)/2, mpf(1)/2, 1, -2*sqrt(z)/(1-sqrt(z))),
 "2F1(1/8,3/8;1;4z(1-z))":
    lambda z: hyp2f1(mpf(1)/8, mpf(3)/8, 1, 4*z*(1-z)),
 "(1-z)^{-1/4} 2F1(1/2,1/2;1; (1-sqrt(1-z))/2 )":
    lambda z: (1-z)**mpf(-0.25)*hyp2f1(mpf(1)/2, mpf(1)/2, 1, (1-sqrt(1-z))/2),
 "((1+sqrt(1-z))/2)^{-1/2} 2F1(1/2,1/2;1;(1-sqrt(1-z))/(1+sqrt(1-z)))":
    lambda z: ((1+sqrt(1-z))/2)**mpf(-0.5)*hyp2f1(mpf(1)/2, mpf(1)/2, 1,
              (1-sqrt(1-z))/(1+sqrt(1-z))),
}
zs = [mpf("0.1"), mpf("0.37"), mpf("0.8")]
for name, f in tests.items():
    errs = []
    for z in zs:
        try: errs.append(abs(F(z)-f(z)))
        except Exception as e: errs.append(None)
    good = all(e is not None and e < mpf(10)**-30 for e in errs)
    print(("MATCH  " if good else "no     ") + name,
          "" if good else f"  err@0.37={nstr(errs[1],3) if errs[1] is not None else 'err'}")
