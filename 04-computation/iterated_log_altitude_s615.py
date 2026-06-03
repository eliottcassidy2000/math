import math
# Shortcut Syracuse map
def T(n): return n//2 if n%2==0 else (3*n+1)//2

# ---- Abstraction test 1: altitude theorem ----
# If log R(N) <= rho*log N + C (contraction at LEVEL 1, i.e. in log-space),
# then iterations to reach O(1) ~ (1/log(1/rho)) * loglog N.  Sanity: R=sqrt -> rho=1/2 -> loglog base 2.
def tau(R, N, thresh=4):
    c=0
    while N>thresh and c<10000: N=R(N); c+=1
    return c
import statistics
for name,R,exp in [("N/2 (level0->1 log)", lambda x:x//2, "log2 N"),
                   ("sqrt N (level1->2)", lambda x:int(math.isqrt(x))+1 if x>4 else 1, "log2log2 N")]:
    print(name)
    for N in [2**8,2**16,2**32,2**64,2**128]:
        ll=math.log2(N); lll=math.log2(ll)
        print(f"  N=2^{int(math.log2(N)):<4} tau={tau(R,N):<4} log2N={ll:<7.2f} log2log2N={lll:.2f}")

# ---- The gem: Collatz EPOCH count is a genuine double-log ----
# An "epoch" = run the shortcut map for (current bit-length) steps. Per-step drift = sqrt(3)/2 ~0.866,
# so log2(value) contracts by factor rho~0.79 per epoch  =>  epoch-count ~ C * loglog n.
def epochs(n, thresh=4):
    e=0
    while n>thresh and e<10000:
        b=n.bit_length()
        for _ in range(b):
            n=T(n)
            if n<=thresh: break
        e+=1
    return e

import random
print("\nCollatz epoch-count vs loglog n  (epoch = bit-length-many shortcut steps):")
print(f"{'bits':>5} {'samples':>8} {'mean_epochs':>12} {'log2log2 n':>11} {'ratio':>7}")
random.seed(1)
rows=[]
for bits in [16,32,64,128,256,512,1024]:
    samp=[random.randrange(2**(bits-1),2**bits)|1 for _ in range(400)]
    es=[epochs(n) for n in samp]
    m=statistics.mean(es)
    lll=math.log2(bits)  # log2 log2 n  (n~2^bits so log2 n ~ bits)
    rows.append((bits,m,lll))
    print(f"{bits:>5} {len(samp):>8} {m:>12.3f} {lll:>11.3f} {m/lll:>7.3f}")

# Fit mean_epochs = a*log2log2 n + b  (linear in loglog => confirms double-log)
xs=[r[2] for r in rows]; ys=[r[1] for r in rows]
n=len(xs); sx=sum(xs); sy=sum(ys); sxx=sum(x*x for x in xs); sxy=sum(x*y for x,y in zip(xs,ys))
a=(n*sxy-sx*sy)/(n*sxx-sx*sx); b=(sy-a*sx)/n
ss_tot=sum((y-sy/n)**2 for y in ys); ss_res=sum((y-(a*x+b))**2 for x,y in zip(xs,ys))
print(f"\nFit epochs = {a:.3f}*log2log2(n) + {b:.3f}   R^2={1-ss_res/ss_tot:.4f}")
print("predicted slope 1/log2(1/0.79) =", 1/math.log2(1/0.7925))

# ---- Contrast: raw STEP count is single-log (linear in log n / bits) ----
def steps(n,thresh=4):
    c=0
    while n>thresh and c<10**7: n=T(n); c+=1
    return c
print("\nContrast: raw step-count vs log2 n (single log):")
for bits in [16,32,64,128,256]:
    samp=[random.randrange(2**(bits-1),2**bits)|1 for _ in range(200)]
    m=statistics.mean(steps(n) for n in samp)
    print(f"  bits={bits:>4} mean_steps={m:>9.1f}  steps/log2n={m/bits:.3f}")
