"""Independent arithmetic sieve on EVERY saved H bin; recursive H on 16 samples.
Full presentation enumeration and the full partner-H table remain C++ only.
"""
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
import json

def need(b,s):
    if not b:raise RuntimeError(s)

def count_paths(rows):
    @lru_cache(None)
    def end(last,left):
        if not left:return 1
        opts=rows[last]&left;total=0
        while opts:
            bit=opts&-opts;opts-=bit
            total+=end(bit.bit_length()-1,left^bit)
        return total
    return sum(end(v,((1<<len(rows))-1)^(1<<v)) for v in range(len(rows)))

def source(minus):
    rows=[0]*12
    for i in range(12):
        for j in range(i+1,12):
            if i//4!=j//4:forward=(j//4-i//4)%3==1
            elif j%4==3:forward=not(minus and i//4==2)
            else:forward=(j-i)%3==1
            a,b=(i,j) if forward else (j,i);rows[a]|=1<<b
    return rows

def partner(code):
    rows=[0]*12;k=0
    for i in range(12):
        for j in range(i+1,12):
            forward=j==i+1 or not(code>>k&1)
            if j>i+1:k+=1
            a,b=(i,j) if forward else (j,i);rows[a]|=1<<b
    return rows

def linear_sieve(N):
    primes=[];lp=[0]*(N+1);mu=[0]*(N+1);lam=mu.copy();mu[1]=lam[1]=1
    for n in range(2,N+1):
        if not lp[n]:lp[n]=n;primes.append(n);mu[n]=lam[n]=-1
        for p in primes:
            if p*n>N:break
            lp[p*n]=p;lam[p*n]=-lam[n]
            if n%p==0:mu[p*n]=0;break
            mu[p*n]=-mu[n]
    return mu,lam

def main():
    data=json.loads(Path(__file__).with_name('hostile.json').read_text(encoding='utf-8-sig'))
    mu,lam=linear_sieve(max(h for r in data for h,m in r['partner_H_histogram']))
    result=[]
    for i,r in enumerate(data):
        need(count_paths(source(i))==r['H']==462915,'independent source')
        need(r['presentations']*r['Aut']==r['H'],'free presentation fibres')
        need(len({h for h,m in r['partner_H_histogram']})==len(r['partner_H_histogram']),'distinct H bins')
        mh=[0]*3;lh=[0]*2
        for h,m in r['partner_H_histogram']:
            need(h%2==1 and m>0 and m%r['Aut']==0,'parity and multiplicity')
            mh[mu[h]+1]+=m;lh[(lam[h]+1)//2]+=m
        need(mh==r['mu_hist_minus_zero_plus'] and lh==r['lambda_hist_minus_plus'],'all-bin arithmetic sieve')
        need(sum(mh)==sum(lh)==r['H'],'total mass')
        need(mh[2]-mh[0]==r['L_muH'] and lh[1]-lh[0]==r['L_lambdaH'],'signed readouts')
        need(sum(h*m for h,m in r['partner_H_histogram'])==r['L2_one'],'inherited second moment')
        for code,h in r['sample']:need(count_paths(partner(code))==h,'independent recursive sample')
        result.append({'word':r['word'],'distinct_H_bins':len(r['partner_H_histogram']),
                       'mean_mu':str(Fraction(r['L_muH'],r['H'])),
                       'mean_lambda':str(Fraction(r['L_lambdaH'],r['H']))})
    need([r['L_muH'] for r in data]==[-2916,297],'hostile Mobius sign reversal')
    need([r['L_lambdaH'] for r in data]==[-3321,-1323],'hostile Liouville separation')
    need([r['L2_one'] for r in data]==[202898835819,202869738135],'inherited full-sum controls')
    print(json.dumps({'independent_source_counts':2,'independent_partner_counts':16,
                      'all_H_bins_arithmetic_audited':True,'results':result},sort_keys=True,indent=2))

if __name__=='__main__':main()
