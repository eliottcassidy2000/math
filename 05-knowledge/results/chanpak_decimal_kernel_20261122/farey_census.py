"""Exact finite address bound and numerical phase constants for the E61 charts."""
from math import lgamma, log
from fractions import Fraction
from functools import lru_cache
from fermat_receiver import TARGET
from hyperbola_controls import chart


def tower_bound(low, high):
    """Upper bound for address pairs on [low,high] in the canonical nested tower."""
    total = 0
    for n in range(4, 1000, 4):
        h, q, c, t, e, r, k = chart(n)
        if k > high:
            return total, n
        assert c+h <= q and t+h <= q and e+r <= q
        s = (high-k)//q
        u = max(0, (low-1-k)//q)
        # D(s) <= s(1+log s), while D(u) >= u(log u-1).
        upper = 1+(high-k)//t+(high-k)//c+s*(100+70*s.bit_length())//100
        lower = (u*(69*(u.bit_length()-1)-100)+99)//100 if u else 0
        total += upper-max(0, lower)
    raise ArithmeticError("No terminal bare height found")


def phase_parameters(data, m, i, j):
    """Two intercepts and the unfixed radical-inverse width."""
    h, q, c, t, e, r, k = data
    B = m & -m

    def unrevealed(z):
        reverse = 0
        for _ in range(B.bit_length()-1):
            reverse = 2*reverse + (z&1)
            z //= 2
        return Fraction(B-1-reverse, B)

    d = Fraction(h, q*m*B)
    a = (i+Fraction(c,q)+Fraction(h,q)*unrevealed(j))/m
    b = (j+Fraction(t,q)+Fraction(h,q)*unrevealed(i))/m
    return a, b, d


def phase_constant(data, m, i, j):
    """Coefficient of S/m^2 after the term (S/m^2) log(S/m^2)."""
    a,b,d=map(float,phase_parameters(data,m,i,j))
    return (lgamma(a)-lgamma(a+d)+lgamma(b)-lgamma(b+d))/d - 1


UNIT = 10**10


@lru_cache(None)
def digamma_interval(z):
    """Integer lower/upper bounds for UNIT times -psi(z)-gamma."""
    assert 0 < z < 1
    lo = hi = 0
    for l in range(256):
        term = UNIT*(1-z)/((l+1)*(l+z))
        lo += term.numerator//term.denominator
        hi += -(-term.numerator//term.denominator)
    tail = UNIT*(1-z)/(255+z)
    return lo, hi-(-tail.numerator//tail.denominator)


def gamma_interval(data, m, i, j):
    a,b,d = phase_parameters(data,m,i,j)
    lo = digamma_interval(a+d)[0] + digamma_interval(b+d)[0]
    hi = digamma_interval(a)[1] + digamma_interval(b)[1]
    return Fraction(lo, UNIT), Fraction(hi, UNIT)


def log_interval(m):
    """Bound log(m) by atanh after splitting off powers of two."""
    def core(x):
        t=(x-1)/(x+1)
        lower=2*sum(t**(2*k+1)/(2*k+1) for k in range(20))
        return lower,lower+2*t**41/(41*(1-t*t))
    power=m.bit_length()-1
    a,b=core(Fraction(2))
    c,d=core(Fraction(m,2**power))
    return power*a+c,power*b+d


if __name__ == "__main__":
    low, high = 2**120, 2**121-1
    upper, cutoff = tower_bound(low, high)
    odd_targets = (high-low+1)//2
    print("first excluded bare spine", cutoff, "upper bound on presentations", upper)
    print("number of odd targets", odd_targets, "inverse fraction >", odd_targets//upper)

    data = chart(12)
    h, q, c, t, e, r, k = data
    m = 34
    selected = []
    selected_bounds = []
    constants = []
    for i in range(m):
        for j in range(m):
            C = phase_constant(data, m, i, j) - 2*log(m)
            constants.append(C)
            if (k+t*i+c*j+q*i*j+r*(i&1)+e*(j&1)-TARGET) % 17 == 0:
                selected.append(C)
                selected_bounds.append(gamma_interval(data,m,i,j))
    full = phase_constant(data, 1, 0, 0)
    assert abs(sum(constants)/m**2-full) < 1e-6
    print("mod 17 cells", len(selected), "/", m**2,
          "relative leading weight", Fraction(17*len(selected), m**2))
    print("linear coefficients relative to own leading weight",
          "mod 17", sum(selected)/len(selected), "all phases", full)

    logs = log_interval(m)
    base = gamma_interval(data,1,0,0)
    lower = sum(a for a,b in selected_bounds)/len(selected)-2*logs[1]-base[1]
    upper = sum(b for a,b in selected_bounds)/len(selected)-2*logs[0]-base[0]
    assert Fraction(-19,10)<lower<=upper<Fraction(-17,10)
    print("certified intervals imply -1.9 < relative correction < -1.7")
