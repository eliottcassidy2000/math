"""Full eta-ratio coefficient recurrences; independent finite product check."""
from math import comb
from pathlib import Path
import json

CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)

def recurrence(p,h,N):
    c=[1]
    for n in range(1,N+1):
        z=-h*sum(sum(d for d in range(1,j+1) if j%d==0 and d%p)*c[n-j] for j in range(1,n+1))
        need(z%n==0,'integral logarithmic recurrence'); c.append(z//n)
    return c

def literal_product(p,h,N):
    c=[1]+[0]*N
    for d in range(1,N+1):
        if d%p==0: continue
        f={0:1}
        for k in range(1,N//d+1):
            if h>0:
                if k<=h: f[k*d]=(-1)**k*comb(h,k)
            else: f[k*d]=comb(-h+k-1,k)
        c=[sum(c[n-j]*v for j,v in f.items() if j<=n) for n in range(N+1)]
    return c

def main():
    rows=[]; N=40
    for p in [2,3,5,7,13]:
        h=24//(p-1); c=recurrence(p,h,N+1); d=recurrence(p,-h,N)
        need(c==literal_product(p,h,N+1),'eta product')
        need(d==literal_product(p,-h,N),'inverse eta product')
        need(all(sum(c[j]*d[n-j] for j in range(n+1))==0 for n in range(1,N+1)),'reciprocal product')
        coeff=[c[n+1]+p**(h//2)*d[n-1] for n in range(1,N+1)]
        need(c[1]+h==0,'constant normalization')
        rows.append({'level':p,'eta_exponent':h,'Fricke_scale':p**(h//2),'positive_coefficients':coeff})
    need([r['positive_coefficients'][0] for r in rows]==[4372,783,134,51,12],'leading coefficients')
    result={'status':'FINITE-EXACT q^40 checks of proved recurrences and cited modular functions',
            'levels':rows,'checks':CHECKS}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
