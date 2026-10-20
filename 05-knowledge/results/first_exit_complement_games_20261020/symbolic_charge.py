#!/usr/bin/env python3
"""Exhaust exact first-exit cylinders, not a bounded source range.

Every positive odd source n has x=(A*n+B)/d and relaxed z upper
(A*(n-1)+C)/d. The residue has modulus 2*d. None means no upper bound.
All decisions use integer arithmetic; exceptions remain active under -O.
"""
import argparse
import json
from collections import Counter

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def least(r, mod, lo):
    return lo + (r-lo) % mod

def run(depth):
    counts=Counter(); frontier=Counter(); failures=[]; closest=None
    # A, B, C, d, residue, minimum source, maximum source, word
    stack=[(1,0,0,1,1,3,None,())]
    while stack:
        A,B,C,d,r,lo,hi,w=stack.pop()
        m=len(w); counts['prefixes']+=1
        # At a>=2 next upper is 3*U/2^a. Find an entire terminal tail.
        amax=2
        while (d<<amax)*(lo-1)<=3*(A*(lo-1)+C):
            amax+=1
        counts['automatic_exponent_tails']+=1
        need((d<<amax)*(lo-1)>3*(A*(lo-1)+C), 'tail bound')
        for a in range(1,amax):
            AA=3*A; BB=3*B+d; dd=d<<a
            CC=3*C+(dd if a==1 else 0)
            rr=((dd-BB)*pow(AA,-1,2*dd))%(2*dd)
            ll=least(rr,2*dd,lo)
            if hi is not None and ll>hi:
                continue
            counts['nonempty_children']+=1
            need((rr-r)%(2*d)==0,'child cylinder not subset')
            delta=dd-AA
            hh=hi
            if delta>0:
                exitlo=max(ll,BB//delta+1)
                ex=least(rr,2*dd,exitlo)
                if hi is None or ex<=hi:
                    counts['paid_exit_cylinders']+=1
                    if delta*(ex-1)<=CC:
                        failures.append({'source':ex,'word':w+(a,),
                                         'endpoint':(AA*ex+BB)//dd,
                                         'A':AA,'B':BB,'C':CC,'denominator':dd})
                        return dict(counts), dict(frontier), failures
                    # Payment increases with source in the same cylinder.
                hh=min(hi,BB//delta) if hi is not None else BB//delta
            if hh is not None and ll>hh:
                continue
            if m+1<depth:
                stack.append((AA,BB,CC,dd,rr,ll,hh,w+(a,)))
            else:
                frontier['cylinders']+=1
                frontier['bounded_cylinders']+=hh is not None
                frontier['least_source']=min(frontier.get('least_source',ll),ll)
    return dict(counts),dict(frontier),failures

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--depth',type=int,default=16)
    args=ap.parse_args();need(args.depth>=1,'positive depth')
    counts,frontier,failures=run(args.depth)
    print(json.dumps({'status':'FINITE-EXACT symbolic cylinder certificate',
        'scope':'all positive odd sources whose first strict descent has at most depth odd steps',
        'depth':args.depth,'counts':counts,'unresolved_frontier':frontier,
        'failures':failures},indent=2,sort_keys=True))
    if failures:raise SystemExit(2)

if __name__=='__main__':main()
