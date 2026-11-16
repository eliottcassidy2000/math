import sys
from itertools import combinations
from check_27_63 import restriction, profile, key, scc_sizes
from rewires27 import entries

def pad(a,p,s):
    n=p+len(a)+s
    m=len(a)
    return tuple(
        (((((1<<s)-1)<<(p+m))|(a[i-p]<<p)) if p<=i<p+m else ((1<<n)-(1<<(i+1))))
        for i in range(n)
    )

if __name__=='__main__':
    sources, hosts=entries(sys.argv[1])
    strong={key(a):a for a in sources.values() if scc_sizes(a)==(len(a),)}
    seen={}
    for host in hosts:
        for k in range(8):
            for S in combinations(range(8),k):
                a=restriction(host,S)
                for n in (6,7):
                    if k>n:continue
                    for p in range(n-k+1):
                        b=pad(a,p,n-k-p)
                        for e,f in combinations(range(n),2):
                            c=list(b);c[e]^=1<<f;c[f]^=1<<e
                            if profile(c)[-1]!=27:continue
                            z=key(c)
                            if z not in strong:continue
                            if z not in seen or k>seen[z][0]:seen[z]=k,host,S,p,(e,f),tuple(c)
    for z,a in strong.items():print('Source',a,seen.get(z))
