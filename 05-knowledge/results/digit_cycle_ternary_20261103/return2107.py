"""Guarded native-pair repair retaining every old vertex at 2107."""
import json
from functools import lru_cache
from support import HERE,need,save,source,junction,paths,flip,phases,observer,values
from observer import records

def guard(T,a,b):
    if not T[a]>>b&1:a,b=b,a
    for P in paths(T):
        for i in range(len(P)-2):
            if P[i:i+2]==(a,b):return P,i
    return None
def main():
    U,r,b,names=source();d=json.loads((HERE/'target2107.json').read_text());chosen=None
    for witness in d['witnesses']:
        edits=witness['flips'];full=(1<<len(edits))-1
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
        schedule=seek(0)
        if schedule is not None:chosen=witness;break
    need(chosen is not None,'guarded strong repair exists')
    T=U;trace=[];fresh=max(names)+1;n=len(T)
    for j in schedule:
        a,c=chosen['flips'][j];P,i=guard(T,a,c);left,right=P[i:i+2]
        X=junction.ear(T,((1<<n)-1)^(1<<left));p1=P+(n,);p2=P[:i+1]+(n,)+P[i+1:]
        need(all(X[v]>>w&1 for v,w in zip(p1,p1[1:])),'first native guard')
        Y=junction.native(X,p1)
        need(all(Y[v]>>w&1 for v,w in zip(p2,p2[1:])),'second native guard')
        Z=junction.native(Y,p2);W=tuple(x&((1<<n)-1) for x in Z[:n])
        need(W==flip(T,a,c) and junction.strong(W),'isolated repair and strong endpoint')
        trace.append({'helper':fresh,'arc':[names[left],names[right]],'H_before':junction.H(T),'H_after':junction.H(W),
            'paths':[[names[v] if v<n else fresh for v in path] for path in (p1,p2)],'rows':W})
        T=W;fresh+=1
    profile=junction.profile(T,r);need(profile[0]==2107 and profile[3]==2109,'actual retained2-ready root')
    need(values(T,r)[1]==profile,'independent endpoint count')
    V=junction.ear(junction.ear(T,1<<r),1<<n);need(junction.H(V)==6325,'literal next offset-four child')
    save('return2107.json',{'status':'FINITE_EXACT_ALL_OLD_NAMES_RETAINED','source_rows':U,'names':names,
        'root':r,'spare':b,'trace':trace,'final_rows':T,'profile':profile,'native_moves':2*len(schedule),
        'helper_births':len(schedule),'old_deletions':0,'peak_order':n+1,'initial_phases':phases(U,r),
        'final_phases':phases(T,r),'spare_phases':phases(T,b),'two_port_observer':records(observer(T,(r,b))),
        'next_rows':V,'next_names':names+[fresh,fresh+1],'next_root':n,'next_H':6325,'next_phases':phases(V,n)})
    print('PASS native return2107',2*len(schedule),'native moves; all names; profile',profile,'next6325',flush=True)
if __name__=='__main__':main()
