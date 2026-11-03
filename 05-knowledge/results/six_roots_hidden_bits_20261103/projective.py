#!/usr/bin/env python3
"""Exact 63/65 actions on the 65 path addresses of E8, not vertex actions."""
import json
from pathlib import Path
P=Path(__file__).resolve().parent
def need(v,msg):
    if not v:raise RuntimeError(msg)
def mul(a,b):
    z=0
    while b:
        if b&1:z^=a
        b>>=1;a<<=1
        if a&64:a^=67 # X^6+X+1
    return z
def power(a,n):
    z=1
    while n:
        if n&1:z=mul(z,a)
        a=mul(a,a);n>>=1
    return z
def trace(a):
    z=0
    for _ in range(6):z^=a;a=mul(a,a)
    return z
def main():
    need(len({power(2,k) for k in range(63)})==63 and power(2,63)==1,'F64 primitive')
    c=next(a for a in range(64) if trace(a)==1)
    need(all(mul(y,y)^y^c for y in range(64)),'quadratic has no root')
    def times(x,y):
        u,v=x;s,t=y
        return mul(u,s)^mul(c,mul(v,t)),mul(u,t)^mul(v,s)^mul(v,t)
    def xp(x,n):
        z=(1,0)
        for _ in range(n):z=times(z,x)
        return z
    units=[(u,v) for u in range(64) for v in range(64) if mul(u,u)^mul(u,v)^mul(c,mul(v,v))==1]
    need(len(units)==65,'norm-one size')
    gen=next(x for x in units if xp(x,5)!=(1,0) and xp(x,13)!=(1,0))
    need(xp(gen,65)==(1,0),'order65')
    def act(x,z):
        u,v=x
        if z is None:return None if v==0 else mul(u,power(v,62))
        den=mul(v,z)^u^v;num=mul(u,z)^mul(c,v)
        return None if den==0 else mul(num,power(den,62))
    orbit=[];z=None
    for _ in range(65):orbit.append(z);z=act(gen,z)
    need(z is None and len(set(orbit))==65,'regular projective action')
    rows=[sum(1<<j for j in range(i+1,8)) for i in range(8)];rows[0]^=128;rows[7]|=1
    def decode(t):
        if t is None:return list(range(8))
        return [i+1 for i in range(6) if t>>i&1]+[7,0]+[i+1 for i in range(6) if not t>>i&1]
    paths=[decode(t) for t in orbit]
    need(all(all(rows[a]>>b&1 for a,b in zip(w,w[1:])) for w in paths),'literal paths')
    need(len({tuple(w) for w in paths})==65,'all path addresses')
    # Complement of the six subset bits is translation by 63, fixing infinity.
    need(all(set(decode(t)[:decode(t).index(7)])=={i+1 for i in range(6) if t>>i&1} for t in range(64)),'inverse decoder')
    # Ranking excess: a different observable from Hamilton-path count.
    excess={p:sum((p-d)*(1 if pow(d,(p-1)//2,p)==1 else -1) for d in range(1,p)) for p in (7,23,67)}
    need(excess=={7:7,23:69,67:67},'Paley ranking')
    out={'status':'FINITE-EXACT explicit address action; see general proof in PORT_MEMORY','F64_polynomial':'X^6+X+1','quadratic_parameter':c,'norm_one_generator':gen,'cycle65':orbit,'Paley_natural_order_excess':excess,'lost':'This action permutes Hamilton paths, not tournament vertices or rooted ear contexts.'}
    (P/'projective.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS: 63-cycle, norm-one 65-cycle, all actual E8 paths; Paley excess 7/69/67')
if __name__=='__main__':main()
