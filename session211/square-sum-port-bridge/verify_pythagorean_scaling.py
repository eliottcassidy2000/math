#!/usr/bin/env python3
"""Fresh finite replay of inherited THM-4510; not a new all-order proof.
For every primitive 1<=s<t<=60 with s²+t²=u², follow the unique path
in the three-target graph on n=t²-1 and t². These are also Q_n paths,
not counts of all Q_n paths.
"""
from hashlib import sha256
from math import gcd,isqrt
from pathlib import Path
import json

def need(c,msg):
    if not c:raise RuntimeError(msg)

def main():
    out=[]
    for t in range(2,61):
        for s in range(1,t):
            u=isqrt(s*s+t*t)
            if gcd(s,t)!=1 or u*u!=s*s+t*t:continue
            target={s*s,t*t,u*u};checks=[]
            for n in (t*t-1,t*t):
                adj={x:sorted(z-x for z in target if 1<=z-x<=n and z!=2*x) for x in range(1,n+1)}
                leaves=[x for x in adj if len(adj[x])==1]
                need(len(leaves)==2 and all(1<=len(adj[x])<=2 for x in adj),'degree structure')
                w=[leaves[0]]
                while True:
                    nxt=[z for z in adj[w[-1]] if len(w)<2 or z!=w[-2]]
                    if not nxt:break
                    need(len(nxt)==1,'unique successor')
                    w.append(nxt[0]);need(len(w)<=n,'no cycle')
                need(sorted(w)==list(range(1,n+1)),'spanning permutation')
                need(all(x+y in target for x,y in zip(w,w[1:])),'target sums')
                checks.append({'n':n,'endpoints':leaves,'path_sha256':sha256(','.join(map(str,w)).encode()).hexdigest()})
            out.append({'triple':[s,t,u],'checks':checks})
    need(len(out)==11,'universe size')
    r={'status':'PASS','scope':'all 11 primitive integer triples 1<=s<t<=60; two prescribed orders each',
       'note':'Checks inherited construction, not all square-sum paths or a new theorem.','rows':out}
    Path(__file__).with_name('MACRO_RESULT.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2))
if __name__=='__main__':main()
