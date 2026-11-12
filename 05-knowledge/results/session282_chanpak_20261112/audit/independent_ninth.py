#!/usr/bin/env python3
from independent_carrier import packet,check
import json
D=(242,244,249,243,224,64,128,32)
def tt(n):return tuple(sum(1<<j for j in range(i+1,n))for i in range(n))
def join(*parts):
    out=[];shift=0;n=sum(map(len,parts))
    for A in parts:
        future=((1<<n)-1)^((1<<(shift+len(A)))-1)
        out +=[(x<<shift)|future for x in A];shift+=len(A)
    return tuple(out)
a=5;b=3;deletion=join(tt(a),D,tt(b));inc=25|(122<<a)|(5<<(a+len(D)));n=len(deletion)
T=tuple([(((1<<n)-1)^inc)<<1]+[(x<<1)|((inc>>i)&1)for i,x in enumerate(deletion)])
p=packet(T)
check(p==dict(H=54381,h=15,Q=196608,U=76389,V=139878,d=15,e=0,k=0),'ninth independent packet')
V=16*139878-15;H=54365+8*139878;z=(V-15-18)//15;P=2*z
check(z==149200 and P==298400,'exact carry initialization')
def g(x):
    if x==0:return 1
    if x&1:return g(x//2)
    v=x//2;return g(v)+v+(v%2)
N=H+15*(g(P)-1)+(V-15)*P
check(g(P)==291769 and N==667830121109,'eighteen baseline')
check((V-15)-15*z==18,'nine to eighteen shift')
print(json.dumps(dict(status='PASS',rows=T,packet=p,prepared_H=H,prepared_V=V,P=P,gP=g(P),eighteen_start=N,eighteen_end=N+34),indent=2,sort_keys=True))
