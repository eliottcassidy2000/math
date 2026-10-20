#!/usr/bin/env python3
"""Finite probe of a deliberately weaker positive-charge certificate."""
import json
from fractions import Fraction

def main():
    best=Fraction(-1);records=[];cases=0;max_depth=0
    for n in range(3,1000002,2):
        x=n;v=n-1;K=0;m=0;word=[]
        while x>=n:
            raw=3*x+1;a=(raw&-raw).bit_length()-1;x=raw>>a
            v=3*v+((1<<(K+a)) if a==1 else 0)
            K+=a;m+=1;word.append(a)
            if m>2000:raise RuntimeError('declared first-exit cap')
        if v >= (n-1)*(1<<K):raise RuntimeError(('unpaid positive-charge bound',n))
        cases+=1;max_depth=max(max_depth,m)
        lost=v+(1-x)*(1<<K);margin=(n-x)*(1<<K)
        if lost*best.denominator>best.numerator*margin:
            best=Fraction(lost,margin)
            records.append({'source':n,'steps':m,'clock':K,'endpoint':x,
                            'lost_reset_dividend_over_actual_descent':str(best),'word':word})
    print(json.dumps({'status':'FINITE-EXACT; all-order first-exit assertion OPEN',
                      'universe':[3,1000001,2],'sources':cases,'route_cap':2000,
                      'maximum_first_exit_steps':max_depth,'records':records},indent=2,sort_keys=True))

if __name__=='__main__':main()
