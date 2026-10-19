"""The golden Collatz limit cycle: Fibonacci-shape balanced rational 3x+1 cycles x_t = B_t/D converge.
Claims: x_min -> c_min = sum_{r>=0} 3^r 2^(-ceil((r+1) phi)),  x_max -> c_max = (1/2) sum_{r>=0} 3^r 2^(-floor(r phi)),
and the limit set is x(theta) = (1/2) sum_r 3^r 2^(-floor(r phi + theta)), theta in [0,1), on which the Collatz step
dictated by the Sturmian letter s_0(theta) = floor(phi + theta) - floor(theta) acts as the rotation theta -> theta + phi mod 1.
Convergence needs 3 < 2^phi, i.e. log2 3 < phi (1.58496 < 1.61803)."""
from mpmath import mp, mpf, floor, ceil, sqrt, log, nsum, inf
from math import gcd
import sys
mp.dps = 40
phi = (1 + sqrt(5))/2
def c_min(R=4000): return sum(mpf(3)**r / mpf(2)**int(ceil((r+1)*phi)) for r in range(R))
def c_max(R=4000): return sum(mpf(3)**r / mpf(2)**int(floor(r*phi)) for r in range(R)) / 2
def x_theta(theta, R=4000): return sum(mpf(3)**r / mpf(2)**int(floor(r*phi + theta)) for r in range(R)) / 2
cm, cM = c_min(), c_max()
print(f"ratio of consecutive terms 3/2^phi = {3/mpf(2)**phi}")
print(f"c_min = {cm}")
print(f"c_max = {cM}, c_max/c_min = {cM/cm}")
# compare with the exact Fibonacci-shape cycles
def mech(p, q, t): return [(((k+1)*p + t)//q) - ((k*p + t)//q) for k in range(q)]
def carry(w):
    P, Bc, Q = 1, 0, 1
    for b in w:
        if b: P, Bc, Q = 3*P, 3*Bc + Q, 2*Q
        else: Q = 2*Q
    return Bc
F = [1, 1]
for _ in range(20): F.append(F[-1] + F[-2])
for n in range(4, 15):
    p, q = F[n], F[n+1]
    D = 2**q - 3**p
    elems = [carry(mech(p, q, t)) for t in range(q)]
    xmin, xmax = mpf(min(elems))/D, mpf(max(elems))/D
    print(f"(q,p)=({q},{p}): x_min = {mp.nstr(xmin, 12)}, x_max = {mp.nstr(xmax, 12)}, x_min - c_min = {mp.nstr(xmin - cm, 4)}, x_max - c_max = {mp.nstr(xmax - cM, 4)}")
# circle-rotation conjugacy on the limit set: x(theta + phi mod 1) == T_(s_0(theta))(x(theta))
import random
random.seed(3); bad = 0
for _ in range(200):
    th = mpf(random.random())
    s0 = int(floor(phi + th)) - int(floor(th))          # first letter of the Sturmian word with intercept theta, slope phi?? see note
    x0 = x_theta(th); x1 = x_theta((th + phi) % 1)
    y = (3*x0 + 1)/2 if s0 == 1 else x0/2
    if abs(y - x1) > mpf(10)**(-25): bad += 1
print(f"conjugacy x(theta + phi) = T_s0(x(theta)) with s0 = floor(phi+theta) - floor(theta): failures {bad}/200 at 1e-25")
# alternative convention: letter = floor((theta + phi)) - floor(theta) uses slope phi > 1 -> letters in {1,2}; try the slope 1/phi convention
bad2 = 0
a = 1/phi
def x_alt(theta, R=4000): return sum(mpf(3)**r / mpf(2)**int(floor(r*phi + theta)) for r in range(R)) / 2
for _ in range(50):
    th = mpf(random.random()); x0 = x_alt(th)
    for s0, nxt in ((1, (3*x0+1)/2), (0, x0/2)):
        pass
print("done")
