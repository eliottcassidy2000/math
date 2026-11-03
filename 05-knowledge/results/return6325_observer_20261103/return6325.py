"""Literal strong-endpoint repair of the actual retained-root 6325 source."""
from functools import lru_cache
import json
from support import HERE,source,need,save,junction,flip,paths,phases,observer,values
from observer import records

def guard(T,a,b):
    if not T[a]>>b&1:a,b=b,a
    for P in paths(T):
        for i in range(len(P)-2):
            if P[i:i+2]==(a,b):return P,i
    return None
def main():
    U,r,b,names=source();d=json.loads((HERE/'target6325.json').read_text());edits=d['flips'];full=(1<<len(edits))-1
    @lru_cache(None)
    def seek(mask):
        T=U
        for j,(a,c) in enumerate(edits):
            if mask>>j&1:T=flip(T,a,c)
        if not junction.strong(T):return None
        if mask==full:return ()
        for j,(a,c) in enumerate(edits):
            if mask>>j&1 or guard(T,a,c) is None:continue
            tail=seek(mask|1<<j)
            if tail is not None:return (j,)+tail
        return None
    order=seek(0);need(order is not None,'strong-endpoint guarded schedule')
    T=U;trace=[];fresh=max(names)+1;n=len(T)
    for j in order:
        a,c=edits[j];P,i=guard(T,a,c);left,right=P[i:i+2]
        X=junction.ear(T,((1<<n)-1)^(1<<left));p1=P+(n,);p2=P[:i+1]+(n,)+P[i+1:]
        Y=junction.native(X,p1);Z=junction.native(Y,p2);W=tuple(x&((1<<n)-1) for x in Z[:n])
        need(W==flip(T,a,c) and junction.strong(W),'isolated arc repair with strong endpoint')
        trace.append({'helper':fresh,'arc':[names[left],names[right]],'paths':[[names[v] if v<n else fresh for v in p] for p in (p1,p2)],
                      'H_before':junction.H(T),'H_after':junction.H(W),'rows':W})
        T=W;fresh+=1
    z=junction.profile(T,r);need(z==(6325,891,893,6327),'literal ready return')
    need(values(T,r)[1]==z,'independent endpoint profile')
    V=junction.ear(junction.ear(T,1<<r),1<<n);need(junction.H(V)==18979,'next ternary offset')
    save('return6325.json',{'status':'FINITE_EXACT_ALL_OLD_NAMES_RETAINED','source_rows':U,'names':names,'root':r,'spare':b,
      'trace':trace,'final_rows':T,'profile':z,'native_moves':2*len(order),'helper_births':len(order),'old_deletions':0,
      'peak_order':n+1,'initial_phases':phases(U,r),'final_phases':phases(T,r),'spare_phases':phases(T,b),
      'two_port_observer':records(observer(T,(r,b))),'next_rows':V,'next_names':names+[fresh,fresh+1],
      'next_root':n,'next_H':18979,'next_phases':phases(V,n)})
    print('PASS return6325',2*len(order),'native moves; no old deletions; next18979',flush=True)
if __name__=='__main__':main()
