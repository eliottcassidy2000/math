#!/usr/bin/env python3
from check_integer_termination import step, require
from check_exceptional_cylinders import WORDS,chart,fast_run
import json

def v2(n): return (n&-n).bit_length()-1

def compressed_run(n):
    # Maximal valuation1 run, followed by its first non1 valuation.
    t=v2(n+1);u=(n+1)>>t;c=v2(3**t*u-1)
    endpoint=(3**t*u-1)>>c
    j=t;d=t+c;carry=3**t-2**t
    require((3**j*n+carry)==(1<<d)*endpoint,'compressed carry')
    # Exact source guard on u; c must not be replaced by lower bound.
    mod=1<<(c+1)
    require(u%mod==pow(3**t,-1,mod)*(1+(1<<c))%mod,'compressed guard')
    return endpoint,j,d,carry,t,c

def main():
    charts=[chart(w) for w in WORDS[1:]]
    repaired=[]
    for row in charts:
        r,M=row['residue'],row['modulus'];A=3**row['J'];D=row['D'];C=row['C']
        gap=(1<<(D+2))-3*A;debt=3*C+(1<<D)
        require(gap>0 and gap*r>debt,'uniform originalsource certificate')
        for k in range(0,2000,2):
            n=r+M*k;y,w,states=fast_run(n)
            require(w==tuple(row['word']) and y>n and y%4==1,'chart')
            z,a=step(y)
            require(a>=2 and z<n,'originalsource return')
            require((3**(row['J']+1)*n+debt)==(1<<(D+a))*z,'retry carry')
        repaired.append({'residue':r,'modulus':2*M,'gap':gap,'debt':debt,'max_odd_steps':row['J']+1})
    for n in range(3,200000,16):
        x,a=step(n);y,b=step(x)
        require(a==1 and b>=3 and y<n,'3 mod16 repair')
    for n in range(1,100002,2):
        endpoint,j,d,c,t,last=compressed_run(n)
        x=n;word=[]
        for _ in range(j):x,a=step(x);word.append(a)
        require(x==endpoint and word==[1]*(t-1)+[last+1],'native compression')
        require(sum(word)==d,'clock')
    power_returns=[]
    for t in range(2,301):
        n=(1<<t)-1;y,j,d,c,tt,last=compressed_run(n)
        require((y<n)==(t in (2,4,8)),'powerfamily classification')
        if y<n:power_returns.append({'t':t,'source':n,'endpoint':y})
    require(step(27)==(41,1) and step(41)==(31,2) and 31>27,'local not original descent')
    print(json.dumps({'repaired_rays':repaired,'compressed_sources':50001,
      'power_family_test_t_max':300,'power_family_descents':power_returns,
      'local_descent_hostile':[27,41,31]},sort_keys=True))
if __name__=='__main__':main()
