import json
from math import gcd,prod
from pathlib import Path

certificate=json.loads(Path(__file__).with_name("prime_repair_cert.json").read_text())
checked={2}
def prime(p):
    if p in checked:return
    row=certificate[str(p)]
    if p<3 or prod(q**e for q,e,_ in row)!=p-1:raise ValueError(p)
    for q,_,a in row:
        prime(q)
        if pow(a,p-1,p)!=1 or gcd(pow(a,(p-1)//q,p)-1,p)!=1:raise ValueError(p)
    checked.add(p)

F,Z=20+64*470,32+64*470
for b,p in ((0,4239511297982033),(1,4239511297981997)):
    w=bin(F)[2:]+str(b)+''.join(str(1-int(c)) for c in reversed(bin(Z)[2:]))
    L=len(w); h=1; P=1<<(L+1)
    for i,bit in enumerate(w):
        if bit!='1':continue
        for j in range(i+1,L):
            if w[j]!='0':continue
            d=j-i
            h+=(1<<(d-1))+(d==1)
            P+=(1<<(L-d))*(3**(d-1)+1)
    if p!=5*P-3*h:raise ValueError((p,h,P))
    prime(p);print(b,h,P,p,7*p)
