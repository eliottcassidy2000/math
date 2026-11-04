"""Independent permutation controls for the incoming S182/C2 comparison."""
from pathlib import Path
import itertools,json
HERE=Path(__file__).resolve().parent

def need(ok,message):
    if not ok:raise RuntimeError(message)

def profile(rows):
    n=len(rows);h=0;f=[0]*n;e=[0]*n
    for path in itertools.permutations(range(n)):
        bad=[i for i,(u,v) in enumerate(zip(path,path[1:])) if not(rows[u]>>v&1)]
        if not bad:
            h+=1;e[path[-1]]+=1
            for v in path:f[v]+=1
        elif len(bad)==1:f[path[bad[0]]]+=1
    return h,f,e

def converse(rows):
    full=(1<<len(rows))-1
    return tuple(full^row^(1<<v) for v,row in enumerate(rows))

def family(t,m,ports):
    n=t+m+1;r=n-1;rows=[0]*n
    for i in range(m):
        rows[i]=sum(1<<j for j in range(i+1,m))
        if i in ports:rows[i]|=1<<r
    for i in range(t):rows[m+i]=((1<<m)-1)|sum(1<<(m+j) for j in range(i+1,t))
    rows[r]=sum(1<<v for v in range(r) if v not in ports)
    return tuple(rows),[*(f'b{i}' for i in range(m)),*(f'a{i}' for i in range(t)),'r']

def unit_start_family(m):
    a=m;w=m+1
    rows=[sum(1<<j for j in range(i+1,m))|(1<<a) for i in range(m)]
    return tuple(rows+[(1<<w),(1<<m)-1]),w

def birth(rows,cut):
    n=len(rows)
    return tuple(row|(((cut>>v)&1)<<n) for v,row in enumerate(rows))+(((1<<n)-1)^cut,)

def c2(rows,a):
    full=(1<<len(rows))-1
    return birth(birth(rows,full^(1<<a)),full)

def main():
    ours=(20,1,2,7,14);our_target=c2(ours,2)
    incoming,names=family(1,3,{0,1,2});incoming_target,target_names=family(1,5,{0,2,4})
    comparisons=[]
    for role,a,b in (('source',ours,incoming),('target',our_target,incoming_target)):
        ha,fa,unused=profile(a);hb,fb,unused=profile(b)
        need(ha==hb,'equal H control')
        need(sorted(fa)!=sorted(fb),'unrooted nonisomorphism certificate')
        comparisons.append(dict(role=role,H=ha,our_rows=a,incoming_rows=b,
            our_F=fa,incoming_F=fb,our_sorted_F=sorted(fa),incoming_sorted_F=sorted(fb)))
    clocks=[]
    for t in range(1,5):
        rows,names=family(t,3,{0,1,2});h,f,e=profile(rows);reverse_h,s,unused=profile(converse(rows))
        need(reverse_h==h==7*(1<<t)-5,'source clock count')
        defects=[x-h for x in s]
        need(defects[-1]==(1<<t)-2 and defects[2]==1,'root and terminal auxiliary defects')
        need(defects[1]>=3 and all(defects[v]>=6 for v in [0,*range(3,3+t)]),'other auxiliary types')
        plus2=c2(rows,2);h2,f2,unused=profile(plus2)
        need(h2==3*h+2 and f2[-3]==3*f[-1],'actual C2 plus-two root response')
        clocks.append(dict(t=t,names=names,rows=rows,H=h,F=f,start_defects=defects,
            plus2_rows=plus2,plus2_H=h2,plus2_retained_root_E2=f2[-3]))
    rows,names=family(2,3,{0,1,2});h,f,unused=profile(rows)
    r=len(rows)-1;child=c2(rows,r);hh,ff,unused=profile(child)
    tau=[(ff[w]-3*f[w])//2 for w in range(len(rows))]
    need(h==23 and hh==73,'isolated root has sigma two')
    need(tau==[4,8,19,1,1,3] and all(v>0 for v in tau),'all physical marks have positive endpoint debt')
    need(not set(ff)&{3*x for x in f},'even arbitrary output re-marking misses the endpoint target')
    alternate,anames=family(2,5,{0,2,4});ah,af,unused=profile(alternate)
    need(ah==hh and af[-1]==3*f[r]-2==85 and ff[r]==93,'same-source same-count different root response')
    unit_controls=[]
    for m in range(1,7):
        u,w=unit_start_family(m);uh,uf,unused=profile(u);unused,unused,start=profile(converse(u))
        need(uh==uf[w]==(1<<m)+1 and start[w]==1,'zero-ready unit-start family')
        unit_controls.append(dict(m=m,rows=u,w=w,H=uh,F=uf[w],start=start[w]))
    triangle=(2,4,1);triangle_child=c2(triangle,2)
    th,tf,unused=profile(triangle_child);unused,unused,ts=profile(converse(triangle_child))
    need((th,tf[0],ts[0])==(9,9,3),'literal tripler loses unit start')
    output=dict(status='FINITE_EXACT_CONTROLS_FOR_PROVED_BRIDGE_BOUNDARY',
        incoming_commit='b486271d24',comparisons=comparisons,clock_controls=clocks,
        isolated_t2=dict(source_rows=rows,source_names=names,H=h,F=f,C2_rows=child,
            C2_H=hh,C2_F=ff,tau_old_marks=tau,incoming_H=ah,
            incoming_root_E2=af[-1],C2_old_root_E2=ff[r]),
        zero_ready_unit_start_controls=unit_controls,
        C3_tripler_boundary=dict(rows=triangle_child,w=0,H=th,F=tf[0],start=ts[0]),
        scope='All-order no-sigma2-tau0 statement uses the proof in INCOMING_BRIDGE.md, not extrapolation from four clocks.')
    (HERE/'incoming_bridge.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'nonisomorphic_equal_H_pairs':[x['H'] for x in comparisons],
        'isolated_tau':tau,'root_E2_same_H73':[af[-1],ff[r]],'clock_controls':len(clocks)},sort_keys=True))

if __name__=='__main__':main()
