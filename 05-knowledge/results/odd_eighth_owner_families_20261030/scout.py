"""Exact exploratory maximum search for residue-complete, two-band speed cores."""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import gcd
from pathlib import Path
import json
import time

HERE=Path(__file__).resolve().parent


def maximum(speeds):
    speeds=sorted(set(speeds))
    denominators=sorted({a+b for a,b in combinations_with_replacement(speeds,2)})
    best=F(0); witnesses=set()
    for den in denominators:
        for num in range(1,den//2+1):
            val=den
            for v in speeds:
                r=v*num%den; val=min(val,r,den-r)
                if val*best.denominator<best.numerator*den: break
            if val*best.denominator>best.numerator*den:
                best=F(val,den); witnesses={F(num,den)}
            elif val*best.denominator==best.numerator*den:
                witnesses.add(F(num,den))
    return best,sorted(witnesses)


def owner(q,extra):
    return min(max([u]+[v for v in extra if v%q==u]) for u in range(1,q) if gcd(u,q)==1)


def main():
    start=time.monotonic(); rows=[]
    for mask in range(128):
        extra=[9+i for i in range(7) if mask>>i&1]
        s=owner(8,extra)
        for m in range(1,11):
            V=list(range(1,8))+extra+[8*m]
            val,times=maximum(V)
            if m>s and val!=F(m,8*m+s): raise RuntimeError((mask,m,s,val))
            rows.append({'mask':mask,'extra':extra,'m':m,'s':s,'value':str(val),
                         'half_witnesses':[str(t) for t in times]})
    (HERE/'scout.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'rows':len(rows),'seconds':time.monotonic()-start,
                      'distinct_values':sorted({r['value'] for r in rows},key=F)}))


if __name__=='__main__':main()
