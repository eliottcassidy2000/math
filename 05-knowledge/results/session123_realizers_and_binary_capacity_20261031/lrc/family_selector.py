"""Construct witnesses for the certified bounded-eight-core/five-tail family.

Requires the Session122 universal family theorem. No new family is asserted.
Huge speeds are never expanded into their full periodic combs.
"""
from fractions import Fraction as F
from math import ceil,floor
from witness import conductor,norm,need
R=F(1,14)
C={1:F(1,7),2:F(2,7),3:F(3,7),4:F(4,7),5:F(10,7)}

def attach(intervals,v):
    out=[]
    for lo,hi in intervals:
        for n in range(ceil(v*lo+R-1),floor(v*hi-R)+1):
            a=max(lo,(n+R)/v);b=min(hi,(n+1-R)/v)
            need(a<=b,'intersection endpoint order')
            out.append((a,b))
    return out

def choose_family(core,tails):
    original_core=tuple(core);original_tails=tuple(tails)
    need(all(type(x) is int and 1<=x<=13 for x in original_core),'bounded core values')
    need(len(set(original_core))==len(original_core) and 8<=len(original_core)<=13,'distinct bounded core')
    need(len(original_core)+len(original_tails)==13,'thirteen constraints')
    need(all(type(x) is int and x>0 for x in original_tails),'positive integer tails')
    core=set(original_core);tails=sorted(original_tails)
    stats={'promoted_or_padded':0,'finite_tail_intersections':0,'max_expanded_speed':max(core),'cutoffs':[]}
    while tails and tails[0]<=13:
        v=tails.pop(0)
        d=v if v not in core else next(d for d in range(1,14) if d not in core)
        core.add(d);stats['promoted_or_padded']+=1
    def finish(t,method):
        need(all(norm(v*t)>=R for v in (*original_core,*original_tails)),'all original constraints')
        stats['method']=method
        return t,stats
    if not tails:return finish(R,'bounded_pool')
    J=[(F(0),F(1))]
    for v in sorted(core):J=attach(J,v)
    while tails:
        k=len(tails);lo,hi=max(J,key=lambda z:z[1]-z[0],default=(F(0),F(0)));W=hi-lo
        need(W>0,'positive width from universal completion')
        cutoff=ceil(C[k]/W);stats['cutoffs'].append(cutoff)
        if tails[0]>=cutoff:
            if k==1:
                v=tails[0];p=(v*lo)%1
                delta=R-p if p<R else F(0) if p<=1-R else 1+R-p
                t=lo+delta/v;need(t<=hi,'one-tail inclusive cutoff')
                return finish(t,'one_tail')
            t,st=conductor((lo,hi),tails,R);stats['terminal_stats']=st
            return finish(t,'conductor')
        v=tails.pop(0);J=attach(J,v)
        stats['finite_tail_intersections']+=1;stats['max_expanded_speed']=max(stats['max_expanded_speed'],v)
    need(bool(J),'completed finite row nonempty')
    return finish(J[0][0],'finite_leaf')
