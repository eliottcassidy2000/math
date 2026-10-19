"""The golden Collatz circle (THM-4717): x(theta) = (1/2) sum_{r>=0} 3^r 2^(-floor(r phi + theta)), theta in [0, phi).
Rotation R(theta) = theta + 1 mod phi; letter s(theta) = 1 iff theta >= 1/phi.  Claim: x(R theta) = (3x(theta)+1)/2 if
s(theta) = 1, else x(theta)/2.  Also: c_min = x(phi^-) = sum 3^r 2^(-ceil((r+1)phi)), c_max = x(0); Cantor structure."""
from mpmath import mp, mpf, floor, ceil, sqrt
import random
mp.dps = 50
phi = (1 + sqrt(5))/2
R = 6000
def x(theta): return sum(mpf(3)**r / mpf(2)**int(floor(r*phi + theta)) for r in range(R)) / 2
random.seed(5); bad = 0
for _ in range(300):
    th = mpf(random.random()) * phi
    s = 1 if th >= 1/phi else 0
    th2 = (th + 1) % phi
    y = (3*x(th) + 1)/2 if s == 1 else x(th)/2
    if abs(y - x(th2)) > mpf(10)**(-40): bad += 1
print(f"conjugacy x(theta+1 mod phi) = T_s(x(theta)), s = [theta >= 1/phi]: failures {bad}/300 at 1e-40")
cmin = sum(mpf(3)**r / mpf(2)**int(ceil((r+1)*phi)) for r in range(R))
cmax = x(mpf(0))
print("c_min =", mp.nstr(cmin, 45)); print("c_max =", mp.nstr(cmax, 45)); print("c_max/c_min =", mp.nstr(cmax/cmin, 30))
print("x(phi^-) - c_min =", mp.nstr(x(phi - mpf(10)**(-30)) - cmin, 5))
# monotone decreasing step function with dense jumps: sample 20000 thetas, count distinct values, max gap in the value set
ths = sorted(mpf(random.random())*phi for _ in range(5000))
vals = [x(t) for t in ths]
assert all(vals[i] >= vals[i+1] for i in range(len(vals)-1))
print("x is non-increasing on 5000 samples; distinct values:", len(set(mp.nstr(v, 30) for v in vals)), "of 5000")
# jumps: the r-th term drops by half when theta crosses (j - r phi) mod ...; total variation check
print("total variation c_max - c_min =", mp.nstr(cmax - cmin, 20), " = sum over r of term drops:", mp.nstr(sum(mpf(3)**r*(mpf(2)**(-int(floor(r*phi))) - mpf(2)**(-int(ceil((r+1)*phi))+1)) for r in range(R))/2, 20))
