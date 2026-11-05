import sys,time,json
from collections import deque
from transfer import *

def closure(w,M,max_states=1000000):
    states,mat=automaton(w); final=[i for i,s in enumerate(states) if len(s)==1]
    v=(1,)+(0,)*(len(states)-1); seen={v:''}; q=deque([v]); counts={}; steps=0;t=time.time()
    while q:
        v=q.popleft(); word=seen[v]
        h=sum(v[i] for i in final)
        if h and h not in counts:counts[h]=word
        for mask in range(1<<w):
            nv=advance(v,mat,mask,M+1)
            if sum(nv[i] for i in final)>M:continue
            if nv not in seen:
                seen[nv]=word+str(mask); q.append(nv)
        steps+=1
        if len(seen)>max_states:return {'complete':False,'states':len(seen),'processed':steps,'time':time.time()-t,'counts':counts}
        if steps%100000==0:print('progress',steps,len(seen),len(q),time.time()-t,flush=True)
    return {'complete':True,'states':len(seen),'processed':steps,'time':time.time()-t,'counts':counts}

if __name__=='__main__':
    w=int(sys.argv[1]);M=int(sys.argv[2]);r=closure(w,M,int(sys.argv[3]) if len(sys.argv)>3 else 1000000)
    print({k:v for k,v in r.items() if k!='counts'})
    print('gaps',[i for i in range(1,M+1,2) if i not in r['counts']])
    with open(f'session197/ribbon/closure-{w}-{M}.json','w') as f:json.dump(r,f,indent=2)
