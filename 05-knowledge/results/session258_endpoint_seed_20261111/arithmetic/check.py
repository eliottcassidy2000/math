#!/usr/bin/env python3
"""Exact, targeted controls; no assertion statements; Python -O safe."""
import json, math, pathlib

def require(ok, message):
    if not ok: raise RuntimeError(message)
def prime(n):
    return n >= 2 and all(n%d for d in range(2,math.isqrt(n)+1))
def factor(n):
    out=[]; p=2
    while p*p<=n:
        e=0
        while n%p==0: n//=p; e+=1
        if e: out.append([p,e])
        p+=1
    if n>1:out.append([n,1])
    return out

def move(s,l):
    D,E,U,d=s
    return {'O':(D,0,2*U,d),'Z':(D-E+U,d,2*U-d,d),'J':(D-E+U,D-E,D+U,D)}[l]
def run(w):
    s=(1,0,2,1)
    for l in w:s=move(s,l)
    return s
def F(s):return s[0]-s[1]+s[2]
def main():
    ps=[p for p in range(2,100) if prime(p)]
    require(len(ps)==25 and ps[-1]==97,'first25')
    rows=[]; period=1
    for p in ps[1:]:
        d=1
        while pow(2,d,p)!=1:d+=1
        residues=[r for r in range(d) if (99*pow(2,r,p)+23)%p==0]
        if residues:period=math.lcm(period,d)
        for k in range(1,4*d+1):
            require(((99*2**k+23)%p==0)==(k%d in residues),'modular phase')
        rows.append(dict(prime=p,order=d,forbidden=residues))
    require(period==428468040,'period')
    table=[]
    for k in range(1,17):
        n=99*2**k+23; fac=factor(n)
        require(math.prod(p**e for p,e in fac)==n,'factor reconstruction')
        require(all(prime(p) for p,e in fac),'factor primality')
        table.append(dict(k=k,value=n,factors=fac,wheel_survivor=all(n%p for p in ps)))
    false=[r for r in table if r['wheel_survivor'] and not prime(r['value'])]
    require([r['k'] for r in false]==[16],'first wheel false positive through16')
    # Proved obstruction gets a small complete exponent-box control, not a proof substitute.
    two_run_checks=0
    for k in range(1,65):
        n=99*2**k+23
        require(n%3==2 and n%9==5,'power gates')
        require((n+1)%3==0 and (n+1)%9!=0,'successor power gate')
        if k>=2:require((n-1)%4==2 and (n-1)//2==99*2**(k-1)+11,'unit-offset slope')
        for r in range(2,n.bit_length()+1):
            q,rem=divmod(n-2,2**r-1)
            is_mersenne=(q>=3 and (q+1)&q==0)
            require(rem!=0 or not is_mersenne,'two-run exclusion')
            two_run_checks+=1
    # Inherited physical O/J twin, plus explicit ordered continuation distinction.
    twin_controls=0
    for a in range(0,7):
      for b in range(0,7):
        s=run('O'*a+'Z'*b); D,E,U,d=s
        o,j=move(s,'O'),move(s,'J')
        require(o!=j and F(o)==F(j)==D+2*U,'physical twin')
        for k in range(0,7):
            x,y=o,j
            for _ in range(k):x,y=move(x,'Z'),move(y,'Z')
            require(F(x)-F(y)==(2**(k+1)-2)*(U-d),'twin splitting')
            twin_controls+=1
    require(F(run('OOZZZZZ'))==443 and F(run('OOOOOZZ'))==443,'443 sources')
    # Independent affine composition checks on Catalan same-clock carry hostile.
    def carry(word):
        B=0; K=0
        for e in word:B=3*B+2**K; K+=e
        return B,2**K-3**len(word)
    require(carry([1,3])==(5,7) and carry([2,2])==(7,7),'carry hostile')
    ready=[]
    for k in range(1,13):
        D=2**k; H=99*D+23; U=450*3**k+204*D-117
        ready.append(dict(k=k,H=H,E=3,S=27*D+2,U=U,V=H+1,d=27*D+3,K=450*3**k+240*D-126))
    return dict(status='FINITE-EXACT controls accompanying all-k proofs; no spectrum-completeness claim',
        first25_primes=ps,phase_rows=rows,valid_common_index_period=period,
        factor_table=table,two_run_controls=two_run_checks,twin_controls=twin_controls,
        inherited_ready_states=ready,hostile=dict(k=16,n=6488087,factors=[1301,4987]))
if __name__=='__main__':
    value=main(); out=json.dumps(value,indent=2,sort_keys=True)+'\n'
    target=pathlib.Path(__file__).with_name('receipt.json')
    import sys
    if '--write' in sys.argv:target.write_text(out)
    else:require(target.read_text()==out,'receipt mismatch')
    print('PASS: exact phases, factors, two-run exclusions, twin splitting and carry controls')
