"""Literal native implementation of the marked 77 -> 233 source return."""
from itertools import permutations
import json
from return_probe import HERE,old,junction,need,observer,values
from observer import records
from golden_transport import phase,golden

def phases(T,r):
    v,profile=values(T,r);h,ss,ee,u,p,d=v
    need((u-ee)%2==0,'integral golden coordinates')
    g=(h,(u-ee)//2,(u+ee)//2)
    return {'profile':profile,'golden':g,'golden_phase':phase(g),
            'all_five_Q4_phases':{name:phase(golden.step(g,golden.BANK[name])) for name in ('q0','q1','q2','q4','q5')}}

def paths(T):
    full=(1<<len(T))-1
    def walk(path,mask):
        if mask==full:yield tuple(path);return
        avail=(T[path[-1]]&~mask) if path else full
        while avail:
            bit=avail&-avail;avail-=bit;v=bit.bit_length()-1
            yield from walk(path+[v],mask|bit)
    yield from walk([],0)
def guard(T,a,b):
    if not T[a]>>b&1:a,b=b,a
    return next((p for p in paths(T) if p[:2]==(a,b)),None)
def flip(T,a,b):
    out=list(T);out[a]^=1<<b;out[b]^=1<<a;return tuple(out)
def main():
    d=json.loads((HERE/'return_probe.json').read_text());U=tuple(d['child233']);candidate=d['minimizers'][0]
    flips=[tuple(x) for x in candidate['flips']];deleted=candidate['deleted'][0];schedules=[]
    for sequence in permutations(range(5)):
        T=U;names=list(range(9));hs=[233];valid=True;peak=9
        for op in sequence:
            if op==4:
                keep=[i for i,n in enumerate(names) if n!=deleted];T=junction.sub(T,tuple(keep));names=[names[i] for i in keep]
            else:
                a,b=(names.index(n) for n in flips[op])
                if guard(T,a,b) is None:valid=False;break
                T=flip(T,a,b);peak=max(peak,len(T)+1)
            if not junction.strong(T):valid=False;break
            hs.append(junction.H(T))
        if valid:schedules.append({'order':sequence,'H':hs,'peak':peak})
    need(schedules,'at least one strong legal native schedule')
    chosen=min(schedules,key=lambda x:(x['peak'],max(x['H']),x['order']))
    T=U;names=list(range(9));trace=[];fresh=9
    for op in chosen['order']:
        before=junction.H(T)
        if op==4:
            keep=[i for i,n in enumerate(names) if n!=deleted];T=junction.sub(T,tuple(keep));names=[names[i] for i in keep]
            trace.append({'operation':'delete','name':deleted,'H_before':before,'H_after':junction.H(T)});continue
        a,b=(names.index(n) for n in flips[op]);P=guard(T,a,b);need(P is not None,'actual Hamilton guard')
        n=len(T);first,second=P[:2]
        incoming=(1<<n)-1;incoming^=1<<first
        X=junction.ear(T,incoming)
        p1=P+(n,);p2=P[:1]+(n,)+P[1:]
        need(all(X[x]>>y&1 for x,y in zip(p1,p1[1:])),'first native path')
        Y=junction.native(X,p1)
        need(all(Y[x]>>y&1 for x,y in zip(p2,p2[1:])),'second native path')
        Z=junction.native(Y,p2);new=tuple(row&((1<<n)-1) for row in Z[:n])
        need(new==flip(T,a,b),'exact isolated old-edge reversal')
        trace.append({'operation':'native_arc','arc':[names[first],names[second]],'helper':fresh,
                      'paths':[[names[v] if v<n else fresh for v in p] for p in (p1,p2)],'H_before':before,'H_after':junction.H(new)})
        fresh+=1;T=new
    r=names.index(7);b=names.index(6);V=observer(T,(r,b));z=junction.profile(T,r)
    need(z==(233,36,37,234),'actual ready endpoint at retained child root')
    next_child=junction.ear(junction.ear(T,1<<r),1<<len(T))
    need(junction.H(next_child)==701,'the returned source actually feeds its next ternary output')
    out={'status':'FINITE_EXACT_NATIVE_RETURN','legal_schedules':len(schedules),'chosen_schedule':chosen,
         'trace':trace,'final_names':names,'final_rows':T,'retained_original_names':[5,6,4],
         'root_name':7,'other_port_name':6,'root_profile':z,'full_observer':records(V),
         'native_moves':8,'temporary_helpers':4,'old_vertex_deletions':[deleted],'peak_order':chosen['peak'],'next_ternary_H':701}
    out['sourced_phases']={'source77':phases(tuple(d['source77']),5),'child233':phases(U,7),
                          'ready233':phases(T,r),'child701':phases(next_child,len(T))}
    (HERE/'return_native.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: exact native return',chosen,'profile',z,'next H701',flush=True)
if __name__=='__main__':main()
