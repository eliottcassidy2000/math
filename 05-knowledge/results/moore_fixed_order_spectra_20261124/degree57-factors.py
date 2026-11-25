from math import comb, factorial

N = 3250
K = N // 2

# Exact parity-constrained spectral minimum for the only candidate set sizes.
m = 1577
cases = []
for x in range(K-m, (N-m)//2+1):
    y = N-m-x
    if 520*x*y > 9*N*m:
        continue
    costs = sorted([4+8*i for i in range(3)]*m
                   + [8+8*i for i in range(3)]*x)
    assert 57*x <= len(costs)
    margin = x + sum(costs[:57*x]) - (225*x+4*x*x)
    assert margin > 0
    cases.append((x,margin))
assert cases == list(zip(range(48,55), [56,124,184,236,280,316,344]))

# Random cut splitting an existing perfect matching: union plus Chebyshev.
den = 2**56
bad = N * sum(comb(56,j) for j in range(15))
assert bad*500**2 + 22750*den < den*500**2

# For every positive-demand cut, at least one of the following two analytic
# inequalities verifies the r-factor condition. Neither uses float rounding.
r = 13
checks = [0,0]
for s in range(1,K+1):
    for t in range(max(s,K-s+1),K+1):
        a, u = s+t-K, K-t
        w = 16*s-r*a
        if 2*w >= u and w*(w-u) >= u*s*(s-1):
            checks[0] += 1
            continue
        z = 115*s*t - 2*N*r*a
        assert z >= 0 and z*z >= 225*s*t*(N-s)*(N-t), (s,t)
        checks[1] += 1

def check_bound(m, degree, edges, exp, r=None):
    b = m+1-4*degree-3
    num, den = factorial(m)*b**edges, (b+1)**edges
    if r is not None:
        num *= (r-1)**((r-1)*m)
        den *= r**((r-2)*m)
    assert num > den << exp
    assert num < den << (exp+1)
    return exp

print('spectral parity excesses:', cases)
print('cut tests using girth / spectrum:', checks)
print('uniform floor exponent:',check_bound(1578,59,52929,14437))
print('selected spine floor exponent:',check_bound(1625,58,47626,18706,r=13))
