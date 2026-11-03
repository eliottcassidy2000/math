"""A literal named return to the2-ready701 lane, with native guards."""
import json
from functools import lru_cache
from common import HERE,need,save,junction,source701,observer
from return_native import paths,flip,phases
from observer import records

def guard(T,a,b):
    if not T[a]>>b&1:a,b=b,a
    for P in paths(T):
        for i in range(len(P)-2):
            if P[i:i+2]==(a,b):return P,i
    return None

def main():
    U,r,b,names=source701();d=json.loads((HERE/'target_match.json').read_text());selected=None
    for witness in d['witnesses_first80']:
        start=tuple(witness['intermediate_rows']);flips=witness['flips'];full=(1<<len(flips))-1
        if not junction.strong(start):continue
        @lru_cache(None)
        def seek(mask):
            T=start
            for j,(a,bb) in enumerate(flips):
                if mask>>j&1:T=flip(T,a,bb)
            if not junction.strong(T):return None
            if mask==full:return ()
            for j,(a,bb) in enumerate(flips):
                if mask>>j&1 or guard(T,a,bb) is None:continue
                tail=seek(mask|1<<j)
                if tail is not None:return (j,)+tail
            return None
        schedule=seek(0)
        if schedule is not None:selected=witness;break
    need(selected is not None,'strong guarded six-edge repair schedule exists')
    P=tuple(selected['native_path']);T=junction.native(U,P);trace=[{'kind':'native','path':[names[i] for i in P],
          'H_before':701,'H_after':junction.H(T),'rows':T}]
    deleted=selected['deleted_index'];keep=tuple(i for i in range(10) if i!=deleted);T=junction.sub(T,keep)
    live=selected['remaining_names'];trace.append({'kind':'delete','name':names[deleted],'H_after':junction.H(T),'rows':T})
    fresh=15;states=[{'rows':T,'root':selected['root'],'spare':selected['spare'],'stage':'after native and deletion'}]
    for j in schedule:
        a,bb=selected['flips'][j];P,i=guard(T,a,bb);n=len(T);left,right=P[i:i+2]
        X=junction.ear(T,((1<<n)-1)^(1<<left));p1=P+(n,);p2=P[:i+1]+(n,)+P[i+1:]
        need(all(X[v]>>w&1 for v,w in zip(p1,p1[1:])),'first actual native guard')
        Y=junction.native(X,p1);need(all(Y[v]>>w&1 for v,w in zip(p2,p2[1:])),'second actual native guard')
        Z=junction.native(Y,p2);W=tuple(row&((1<<n)-1) for row in Z[:n])
        need(W==flip(T,a,bb),'isolated old-edge reversal after deleting helper')
        trace.append({'kind':'native_pair','helper':fresh,'arc':[live[left],live[right]],
           'paths':[[live[v] if v<n else fresh for v in path] for path in (p1,p2)],
           'H_before':junction.H(T),'H_after':junction.H(W),'rows':W})
        T=W;fresh+=1
    rr=selected['root'];bb=selected['spare'];prof=junction.profile(T,rr)
    need(tuple(T)==tuple(selected['target_rows']) and prof[0]==701 and prof[3]==703,'actual2-ready701 endpoint')
    V=observer(T,(rr,bb));child=junction.ear(junction.ear(T,1<<rr),1<<len(T))
    need(junction.H(child)==2107,'literal next offset-four child')
    save('native_return.json',{'status':'FINITE_EXACT_NATIVE_RETURN_AT701_IN_TWO_READY_LANE','source_rows':U,'source_names':names,
      'input_root':r,'input_spare':b,'trace':trace,'final_rows':T,'final_names':live,'root':rr,'spare':bb,'profile':prof,
      'native_moves':1+2*len(schedule),'helper_births':len(schedule),'old_deleted_names':[names[deleted]],'peak_order':10,
      'protected_names':[names[i] for i in d['protected_indices']],'repair_schedule':schedule,
      'source_root_phases':phases(U,r),'final_root_phases':phases(T,rr),'final_spare_phases':phases(T,bb),
      'full_two_port_observer':records(V),'next_rows':child,'next_names':live+[fresh,fresh+1],
      'next_root':len(T),'next_H':2107,'next_phases':phases(child,len(T))})
    print('PASS named701 return:13 native moves, one old deletion, peak10, profile',prof,'next2107',flush=True)
if __name__=='__main__':main()
