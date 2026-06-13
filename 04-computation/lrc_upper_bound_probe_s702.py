# Probe the upper-bound witness rule for the coprime lower-half law.
# For gcd(a,n)=1, a<n/2, N=n+a: claim min_{v in R} ||v θ|| <= 2/N for ALL θ.
# Case 1 (easy): ||aθ|| <= 1/N  => runner 2a has ||2aθ|| <= 2/N.  (2a in R since a<n/2)
# We probe the remaining regime ||aθ|| > 1/N to see which runner saves it.
from math import gcd
from fractions import Fraction

def nn(x):
    r = x % 1.0
    return min(r, 1-r)

worst_case = []
for n in range(5, 30):
    for a in range(1, (n+1)//2):
        if gcd(a,n)!=1 or 2*a>=n: continue
        N = n+a
        R = [v for v in range(1,n) if v!=a]
        # sample theta densely on rationals p/Q
        viol=0; case2_needs=set(); worst=0
        for Q in range(2, 4*N):
            for p in range(1,Q):
                th=p/Q
                if nn(a*th) <= 1.0/N + 1e-12:    # case 1 covered by 2a
                    continue
                # find closest runner in R
                ds=[(nn(v*th),v) for v in R]
                dmin,vmin=min(ds)
                if dmin > 2.0/N + 1e-9:
                    viol+=1
                # record which runner type saves it (relative to a)
                case2_needs.add(vmin - a if vmin>a else ('low',vmin))
                worst=max(worst,dmin)
        if viol:
            print(f"n={n} a={a}: VIOLATION count={viol}")
print("probe done (no VIOLATION lines = upper bound holds in case-2 regime too)")

# Now: in the case-2 regime, is the saving runner always '2a' fails but some structured one works?
# Targeted: when ||aθ||>1/N, check if EITHER (n+j closest is high) and runner (something) works.
