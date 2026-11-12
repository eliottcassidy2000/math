#!/usr/bin/env python3
"""Independent literal selected-cover DP and exact prime-adic twin-clock controls."""
import json

def check(x, msg):
    if not x: raise RuntimeError(msg)

def inflate(rows, selected, v, m):
    old=list(range(len(rows)))
    types=old[:v]+[v]*m+old[v+1:]
    out=[]
    for i,a in enumerate(types):
        out.append(sum(1<<j for j,b in enumerate(types) if i!=j and ((a==b and i<j) or (a!=b and (rows[a]>>b)&1))))
    marks=sum(1<<i for i,a in enumerate(types) if (selected>>a)&1)
    return out,marks

def end_tables(rows):
    n=len(rows); N=1<<n; dp=[[0]*n for _ in range(N)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,N):
        for i,z in enumerate(dp[s]):
            if not z: continue
            rem=rows[i]&~s
            while rem:
                bit=rem&-rem;rem-=bit;j=bit.bit_length()-1
                dp[s|bit][j]+=z
    return dp

def literal(rows,O):
    n=len(rows);full=(1<<n)-1
    e=end_tables(rows)
    rev=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)]
    s=end_tables(rev)
    eo=[sum(r[i] for i in range(n) if O>>i&1) for r in e]
    so=[sum(r[i] for i in range(n) if O>>i&1) for r in s]
    H=sum(e[full]);W=sum(eo[a]*so[full^a] for a in range(1,full))
    U=so[full]+sum(sum(e[a])*so[full^a] for a in range(1,full))
    return {'n':n,'H':H,'W':W,'delta':H-W,'U':U}

def defect(m): return 1-4*(3**m-3*2**m+3)
def val(n,p):
    if not n:return None
    a=0
    while n%p==0:n//=p;a+=1
    return a

def main():
    cores=[((90,28,9,16,36,79,30),72,2),((54,28,24,17,32,78,31),68,0)]
    controls=[]
    for rows,O,v in cores:
        seq=[]
        for m in range(1,6):
            r=literal(*inflate(rows,O,v,m))
            check(r['delta']==defect(m),'literal defect mismatch')
            seq.append({'m':m,**r})
        controls.append({'rows':rows,'selected':O,'port':v,'values':seq})
    for m in range(3,501):
        d=defect(m)-1
        check(d<0,'nonunit sign')
        check(val(d,2)==(3 if m%2 else 4),'dyadic valuation')
        check(val(d,3)==(1 if m%2 else 2+val(m,3)),'ternary valuation')
    near=[]
    for a in range(1,8):
        m=2*3**a;d=defect(m)-1
        check(val(d,3)==a+2,'near-seed valuation')
        near.append({'a':a,'m':m,'v3':val(d,3),'decimal_digits_abs_defect':len(str(abs(d)))})
    # Good-prime physical lengths reproduce a unit residue to requested precision.
    congruences=[]
    for p in (5,7,11,13):
        for a in range(1,7):
            m=2+(p-1)*p**(a-1)
            dmod=(-4*(pow(3,m,p**a)-3*pow(2,m,p**a)+3))%(p**a)
            check(dmod==0,'good-prime residue')
            congruences.append({'p':p,'precision':a,'length':m})
    print(json.dumps({'literal_controls':controls,'valuation_range':[3,500], 'ternary_near_seeds':near,'good_prime_congruences':congruences},sort_keys=True,indent=2))
if __name__=='__main__':main()
