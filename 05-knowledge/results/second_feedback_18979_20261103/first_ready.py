"""A size-dependent return at the first stable clone size, with literal guards."""
import json,sys,subprocess
from itertools import combinations
from functools import lru_cache
from support import HERE,PREV,INCOMING,need,save,junction,clone,phases
sys.path.insert(0,str(PREV))
import moment_search as mom

def hpath(T,V):
    P=[]
    for v in V:
        i=next((i for i,u in enumerate(P) if T[v]>>u&1),len(P));P.insert(i,v)
    need(all(T[a]>>b&1 for a,b in zip(P,P[1:])),'actual insertion Hamilton guard');return P

def compile_flip(T,a,b,r,serial,names):
    n=len(T)
    if not T[a]>>b&1:a,b=b,a
    P0=hpath(T,[i for i in range(n) if i not in (a,b)]);z,s,w=n,n+1,n+2
    cuts=[1<<b,(1<<(n+1))-1,((1<<(n+2))-1)^(1<<a)];X=T
    for cut in cuts:X=junction.ear(X,cut)
    P=(a,b,z,*P0,s,w);Q=(a,w,b,z,*P0,s)
    frame=lambda path:(path[path.index(r)-1] if path.index(r) else None,path[path.index(r)+1] if path.index(r)+1<len(path) else None)
    need(frame(P)==frame(Q),'common root frame')
    Z=junction.native(junction.native(X,P),Q);U=junction.sub(Z,tuple(range(n)))
    diff={(i,j) for i in range(n) for j in range(i+1,n) if ((U[i]>>j)^(T[i]>>j))&1}
    need(diff=={tuple(sorted((a,b)))} and r not in (a,b),'exact isolated exterior flip')
    labels=names+[f'helper{serial}:z',f'helper{serial}:s',f'helper{serial}:w']
    return U,{'arc':[names[a],names[b]],'cuts':cuts,'helpers':labels[n:],'paths':[[labels[i] for i in path] for path in (P,Q)],'rows':U,'completed_strong':junction.strong(U)}

def main():
    code=(INCOMING/'funnel_bank.cpp').read_text();needle='vector<I>targets={18977,18979};';need(code.count(needle)==1,'inherited target site')
    file=HERE/'first_ready.generated.cpp';file.write_text(code.replace(needle,'vector<I>targets={11735,11737};'),encoding='utf8',newline='\n')
    exe=HERE/'first_ready.exe';subprocess.run(['g++','-std=c++20','-O3',str(file),'-o',str(exe)],check=True)
    result=subprocess.run([str(exe),'17'],text=True,capture_output=True,check=True);bank=[json.loads(x) for x in result.stdout.splitlines()]
    targets=[]
    for d in bank:
        if d['kind']!='hit':continue
        P=mom.unrank(d['n'],d['rank']);T=junction.ear(P,1);r=len(P);p=junction.profile(T,r)
        need(p[0]==d['target'] and p[3]-p[0]==1 and junction.strong(T),'literal complete first ready fibre')
        targets.append({'parent_order':d['n'],'rank':d['rank'],'rows':T,'root':r,'profile':p})
    seed=json.loads((PREV/'grandparent_repair.json').read_text());T0=tuple(seed['repaired_rows']);q=seed['root']
    T,ids=clone(T0,q,2);names=[str(seed['source_names'][a])+(f'.{j}' if j else '') for a,j in ids];r=ids.index((q,0))
    clone_cut=sum(1<<i for i in range(len(T0)) if T0[i]>>q&1)|(1<<q)
    appended=junction.ear(T0,clone_cut);address_to_append=[a if j==0 else len(T0) for a,j in ids]
    need(junction.sub(appended,tuple(address_to_append))==T,'literal one-birth clone with old identities preserved')
    need(junction.profile(T,r)[0::3]==(3911,3912),'actual first stable repaired clone')
    source=junction.ear(junction.ear(T,1<<r),1<<len(T));names+=['reset-internal','reset-tip'];root=len(T);incoming=r
    need(junction.H(source)==11735 and junction.strong(source),'actual 3H+2 child of repaired clone')
    target=next(x for x in targets if x['profile'][0]==11735 and len(x['rows'])==14);A=tuple(target['rows']);n=len(A)
    protected={'6','21','31','31.1','reset-internal','reset-tip'};best=None
    for deleted in combinations([i for i,v in enumerate(names) if v not in protected],2):
        keep=[i for i in range(len(source)) if i not in deleted];U=junction.sub(source,tuple(keep))
        if not junction.strong(U):continue
        rr=keep.index(root);ww=keep.index(incoming);mapping=[ww]+[i for i in range(n) if i not in (rr,ww)]+[rr]
        def distance(M):return sum(((A[i]>>j)^(U[M[i]]>>M[j]))&1 for i in range(n) for j in range(i+1,n))
        d=distance(mapping)
        while True:
            changed=False
            for i,j in combinations(range(1,n-1),2):
                mapping[i],mapping[j]=mapping[j],mapping[i];dd=distance(mapping)
                if dd<d:d=dd;changed=True;break
                mapping[i],mapping[j]=mapping[j],mapping[i]
            if not changed:break
        item=(d,tuple(deleted),tuple(mapping),tuple(keep))
        if best is None or item<best:best=item
    need(best is not None,'source-preserving deletion embedding exists')
    d,deleted,mapping,keep=best;ns=[names[i] for i in keep];U=junction.sub(source,keep);rr=keep.index(root);sp=ns.index('6')
    inverse=[mapping.index(i) for i in range(n)];target_rows=junction.sub(A,tuple(inverse));trace=[]
    need(U[rr]==target_rows[rr],'same literal root star')
    edges=[(i,j) for i in range(n) for j in range(i+1,n) if ((U[i]>>j)^(target_rows[i]>>j))&1]
    initial=U
    @lru_cache(None)
    def order(mask):
        if mask==(1<<len(edges))-1:return ()
        for k,(a,b) in enumerate(edges):
            if mask>>k&1:continue
            V=list(initial);nextmask=mask|1<<k
            for j,(u,v) in enumerate(edges):
                if nextmask>>j&1:V[u]^=1<<v;V[v]^=1<<u
            if junction.strong(tuple(V)):
                tail=order(nextmask)
                if tail is not None:return (k,)+tail
        return None
    chosen=order(0);need(chosen is not None,'strong completed-repair ordering')
    edges=[edges[k] for k in chosen]
    for serial,(a,b) in enumerate(edges):
        U,step=compile_flip(U,a,b,rr,serial,ns);need(step['completed_strong'],'every completed repair strong');trace.append(step)
    need(U==target_rows and junction.H(U)==11735 and junction.profile(U,rr)[3]-11735==1,'actual returned source')
    next_source=junction.ear(junction.ear(U,1<<rr),1<<n);need(junction.H(next_source)==35207 and junction.strong(next_source),'literal next ternary child')
    save('first_ready.json',{'status':'FINITE_EXACT_COMPLETE_READY_FIBRES_AND_LITERAL_NATIVE_RETURN','bank':bank,'ready_targets':targets,
        'repaired_seed_rows':T0,'repaired_seed_names':seed['source_names'],'clone_addresses':ids,'clone_birth_cut':clone_cut,'address_to_append':address_to_append,'cloned_source':T,'cloned_root':r,
        'source':source,'source_names':names,'source_root':root,'source_inport':incoming,'source_spare':names.index('6'),'source_phases':phases(source,root),
        'deleted_names':[names[i] for i in deleted],'retained_names':ns,'target_root':rr,'target_spare':sp,'target_to_retained_map':mapping,
        'after_deletion':junction.sub(source,keep),'trace':trace,'target':U,'target_profile':junction.profile(U,rr),'target_root_phases':phases(U,rr),'target_spare_phases':phases(U,sp),
        'native_events':2*d,'helper_births':3*d,'helper_deletions':3*d,'peak_order_during_repairs':n+3,'next_rows':next_source,'next_H':35207,
        'scope':'one actual m=2 return on the named repaired-source branch, retaining both clone addresses and old pins6/21/31; two old deletions are necessary by the complete one-ready11735 fibre; no uniform return or native-cost optimality'})
    print('PASS actual 11735 return',2*d,'native events; deleted',[names[i] for i in deleted],'; completed strong',sum(x['completed_strong'] for x in trace),'of',len(trace),flush=True)
if __name__=='__main__':main()
