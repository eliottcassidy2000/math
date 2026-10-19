#!/usr/bin/env python3
"""Complete bounded-integer-clock inverse tree, with independent forward replay.

The cutoff in N is proved in the companion note, not fitted to scan output.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
from root_clock_channels_20261018o import edge, check


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--clock',type=int,default=12)
    ap.add_argument('--routes',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();T=args.clock
    check(0<=T<=16,'declared certificate clock range')
    path=Path(__file__).parent/'first_arrival_20261018/V_closed_Y1000.txt'
    V={int(s.split()[0]) for s in path.read_text().splitlines() if s and not s.startswith('#')}
    vs=sorted(V);vi={s:i for i,s in enumerate(vs)}
    seen={s:(s,0,0,0) for s in V};stack=sorted(V)
    while stack:
        y=stack.pop();u,t,d,p=seen[y];b=y;j=0
        while 2*j+1<=T-t:
            if b%3:
                a=2 if b%3==1 else 1;dt=2*j+a;x=(2**a*b-1)//3
                if dt<=T-t and x not in V:
                    check(x%2==1 and x%8!=5,'inverse source primitive')
                    target=(u,t+dt,d+1,x if d==0 else p)
                    check(edge(x)==(y,dt),'inverse edge forwards')
                    if x in seen:check(seen[x]==target,'unique stopped path')
                    else:seen[x]=target;stack.append(x)
            b=4*b+1;j+=1
    for s,record in seen.items():
        x=s;t=d=0;p=0
        while x not in V:
            p=x;x,dt=edge(x);t+=dt;d+=1
            check(t<=T,'forward replay within clock budget')
        check((x,t,d,p)==record,'independent literal first-arrival replay including entry port')
    counts=Counter(u for u,t,d,p in seen.values());K2=sum(k*k for k in counts.values())
    bound=(8*K2+32+4)//5
    check(5*bound-32>=8*K2,'uniform contractive tail cutoff integer inequality')
    rows=[[s,vi[u],t,d,p] for s,(u,t,d,p) in sorted(seen.items())]
    args.routes.write_text(f'{len(V)} {T} {bound} {K2}\n'+''.join(' '.join(map(str,r))+'\n' for r in rows))
    args.out.write_text(json.dumps(dict(status='complete finite inverse tree with forward replay',
        clock=T,sources=len(seen),largest_source=max(seen),terminal_fibre_square_sum=K2,
        N_tail_cutoff=bound,ports=vs,routes=rows,
        proof_tail='W <= 3 sum K_u^2; S_N >= 15N/8-12; 3 sum K_u^2 <= S_N above cutoff'),indent=2)+'\n')
    print(T,len(seen),K2,bound)


if __name__=='__main__':main()
