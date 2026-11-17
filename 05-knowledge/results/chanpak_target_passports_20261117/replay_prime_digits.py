from math import prod,gcd
import json
from pathlib import Path
from replay_ordinal_digits import BANK,profile,state,crt
j=json.loads(Path(__file__).with_name('prime_cert.json').read_text())
registry={int(k):v for k,v in j['certificate'].items()}

done={2}
def prove(n):
    if n in done:return
    if n<3:raise ValueError(n)
    v=registry[n]
    if prod(p**e for p,e,_ in v)!=n-1:raise ValueError(n)
    for p,e,a in v:
        prove(p)
        if pow(a,n-1,n)!=1 or gcd(pow(a,(n-1)//p,n)-1,n)!=1:raise ValueError((n,p,a))
    done.add(n)
A=[(profile(a,b),profile(a,c)) for a,b,c,*_ in BANK]
o=(0,1,3,2,5,6,4)
for k,p,b in j['outputs']:
    J,W,K0,E0=state(A,o,(0,)*len(o))
    _,_,K,_=state(A,o,b)
    if W*2**k+K!=7*p or crt(A,o,J,W,K0,7*p,k)!=b:raise ValueError((k,p,b))
    prove(p)
    print(k,p,b)
