"""Independent of the decoder: check44 words and every candidate block modulus."""
from pathlib import Path
from math import gcd,isqrt
import json

def check(condition,message):
    if not condition:raise RuntimeError(message)

def phi_read(word,n,r,modulus,reflection=False):
    check(isinstance(word,int) and 0<=word<2**n,'word width')
    digits=[(word>>i)&1 for i in range(n)]
    check(all(not(digits[i] and digits[(i+1)%n]) for i in range(n)),'cyclic guard')
    return sum(digit*pow(r,(-i if reflection else i)%n,modulus) for i,digit in enumerate(digits))%modulus

src=json.loads(Path(__file__).with_name('period55_certificate.json').read_text())
n,N,r=src['period'],src['N'],src['root']
check(n==55,'period')
f=[0,1]
for i in range(2,n+2):f.append(f[-1]+f[-2])
a,b=f[n-1]-1,f[n]
check(gcd(a,b)==1 and f[n-1]+f[n+1]==N,'cyclic ideal quotient')
check(a*a+a*b-b*b==-N,'norm')
check((a+b*r)%N==0 and (r*r-r-1)%N==0,'chosen root')
check(pow(r,n,N)==1 and all(pow(r,d,N)!=1 for d in (1,5,11)),'clock')
divs={1};product=1
for p,e in src['factorization']:
    check(p>=2 and all(p%q for q in range(2,isqrt(p)+1)),'prime factor')
    check(e>0,'factor exponent')
    product*=p**e
    divs={d*p**j for d in divs for j in range(e+1)}
check(product==N,'complete factorization')
proper=divs-{1,N}
rows=src['certificate']
check(len(rows)==len(proper) and {row['d'] for row in rows}==proper,'candidate block coverage')
for row in rows:
    d=row['d']
    for label in ('x','y'):
        w=row['word_'+label]
        check(phi_read(w,n,r,N)==row[label],'source reader value')
        check(phi_read(w,n,r,N,True)==row['R'+label],'reflection reader value')
    check((row['x']-row['y'])%d==0,'same block')
    check((row['Rx']-row['Ry'])%d!=0,'separated by reflection')
print(json.dumps({'period':n,'states':N,'proper_nontrivial_divisors':len(proper),'word_occurrences_checked':2*len(rows),'enumeration':False,'decoder_imported':False,'conclusion':'primitive translation/reflection action, conditional only on stated all-odd reader theorem'},sort_keys=True))
