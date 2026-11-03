"""Common-frame complement/birth programs on the observer-neutral switch."""
from escape_lifts import *


def hpaths(T,prefix=(),left=None):
    if left is None:left=set(range(len(T)))
    if not left:yield prefix
    else:
        for v in sorted(left):
            if not prefix or T[prefix[-1]]>>v&1:
                yield from hpaths(T,prefix+(v,),left-{v})


def run(m):
    A,_=source4(m,1,2,2);B,_=source4(m,1,3,1)
    h=hwing(m,1,2,2);u=m+5;v=m+2
    shared=0;hits=[];single=[set(),set()];totals=[0,0]
    for p in hpaths(A):
        if any(not B[x]>>y&1 for x,y in zip(p,p[1:])):continue
        shared+=1
        X=native.native(A,p);Y=native.native(B,p)
        xx=ear_values(X);yy=ear_values(Y)
        for mask,(x,y) in enumerate(zip(xx,yy)):
            for i,z in enumerate((x,y)):
                if z in (3*h+2,3*h+4):single[i].add(z);totals[i]+=1
            if set((x,y))=={3*h+2,3*h+4}:
                U=birth(X,mask);V=birth(Y,mask)
                hits.append(dict(path=p,mask=mask,hA=x,hB=y,strongA=strong(U),strongB=strong(V),targetA=U,targetB=V))
        tables.cache_clear()
    row=dict(m=m,h=h,flip=(u,v),shared_paths=shared,single_target_values=[sorted(x) for x in single],literal_target_counts=totals,paired=hits)
    print('m',m,'h',h,'shared',shared,'totals',totals,'paired',len(hits),flush=True)
    return row


if __name__=='__main__':
    rows=[run(m) for m in (1,2,3)]
    (HERE/'switch_search.json').write_text(json.dumps(rows,indent=2)+'\n')
