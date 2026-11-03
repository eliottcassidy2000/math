"""All-slot edge response and its transport through triangle substitution."""
from itertools import product
from pathlib import Path
import json
import sympy as s
from mechanisms import *
HERE=Path(__file__).resolve().parent

def charge(N,r):
    digits=bits(N);r=s.sympify(r);k=len(digits)
    return s.expand(sum(d*r**i*((r+1)**(k-1-i)+(r-1)**(k-1-i))/2 for i,d in enumerate(digits)))

def response(w,u,r):
    m=len(w)
    if not m or r==0:return 0
    # Rational symbolic arithmetic, never floating point.
    b=[s.sympify(r)**j*((s.sympify(r)+1)**(m-1-j)+(s.sympify(r)-1)**(m-1-j))/2 for j in range(m)]
    area=sum(s.sympify(r)**(m-1-(j-i))*((s.sympify(r)+1)**(j-i-1)-(s.sympify(r)-1)**(j-i-1))/2
             *(w[i]*u[j]-w[j]*u[i]) for i in range(m) for j in range(i+1,m))
    return s.expand(2*r*(sum(v*(y-x) for v,x,y in zip(b,w,u))+area))

def all_subset_h(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        h[mask]=sum(dp[mask])
        for last,value in enumerate(dp[mask]):
            available=rows[last]&~mask
            while value and available:
                bit=available&-available;available-=bit
                dp[mask|bit][bit.bit_length()-1]+=value
    return h

def slot_counts(rows,maximum):
    h=all_subset_h(rows);N=len(h);old=[0]*N;old[0]=1;values=[0]
    for r in range(1,maximum+1):
        new=[0]*N
        for mask in range(N):
            sub=mask
            while True:
                new[mask]+=h[sub]*old[mask^sub]
                if sub==0:break
                sub=(sub-1)&mask
        old=new;values.append(old[-1])
    return values

def triangle(rows):
    n=len(rows)
    # X -> c -> d -> X, with the actual old labels retained.
    return tuple(r|(1<<n) for r in rows)+(1<<(n+1),(1<<n)-1)

def main():
    checks=0
    for m in range(5):
        for w in product((0,1),repeat=m):
            for u in product((0,1),repeat=m):
                plus=slot_counts(source(w,u),4);minus=slot_counts(source(w,u,False),4)
                for r in range(5):
                    need(plus[r]-minus[r]==response(w,u,r),'all-slot contraction versus ordered subset convolution')
                    checks+=1
    r=s.Symbol('r');examples=[];charge_checks=0
    for P in range(1,17):
        for R in range(1,17):
            w,u=chart(P,R)
            expected=2*r*(r**P.bit_length()*charge(R,r)-r**R.bit_length()*charge(P,r))
            need(s.expand(response(w,u,r)-expected)==0,'weighted binary cancellation');charge_checks+=1
    for j in range(2,13):
        w,u=chart(2**j,2**j+2)
        need(s.expand(response(w,u,r)-2*r**(2*j+2))==0,'whole-profile offset-zero family')
        if j>=3:
            w,u=chart(2**j+3,2**j+5)
            need(s.expand(response(w,u,r)-2*r**(2*j))==0,'whole-profile offset-three family')
    for P in (1,2,3,4,5,6):
        w,u=chart(P,P+2);polynomial=response(w,u,r)
        lower=source(w,u,False);upper=source(w,u)
        after=hp(triangle(upper))-hp(triangle(lower))
        expected=int(polynomial.subs(r,1)+polynomial.subs(r,2))
        need(after==expected,'literal triangle substitution response')
        if P in (1,2,3):need(s.expand(polynomial-2*r**3)==0,'strong realization of padded triangle difference')
        if P in (4,5):need(s.expand(polynomial-2*r**6)==0,'higher response degree with same scalar jump')
        examples.append({'P':P,'R':P+2,'order':len(lower),'response_polynomial':str(s.factor(polynomial)),
                         'ordinary_jump':int(polynomial.subs(r,1)),'after_triangle_jump':after,
                         'upper_H':hp(upper),'lower_H':hp(lower)})
    need([x['after_triangle_jump'] for x in examples]==[18,18,18,130,130,-126],'context-sensitive sign control')
    out={'status':'PROVED_SCOPED','subset_convolution_checks':checks,
         'symbolic_charge_checks':charge_checks,'profile_family_parameters_through':12,
         'universe':'all pairs of binary words with lengths0 through4; r=0 through4',
         'examples':examples,'inherited_operator':'session159_constructions_20261103/profile_operator/NOTE.md'}
    (HERE/'profile.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print('PASS: all-slot response;',checks,'literal profile checks; triangle jumps 18/130/-126')
if __name__=='__main__':main()
