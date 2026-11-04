"""Independent matrix replay; full relations plus a tagged joint response prefix."""
from pathlib import Path
import sys,json,subprocess
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'return6325_completion_20261103'))
from audit import need,matrix,birth,native,strong,hp
from transport import chart
def induced(T,V):return tuple(sum(((T[a]>>b)&1)<<j for j,b in enumerate(V)) for a in V)
def main():
    data=json.loads((HERE/'returns.json').read_text());states=[];checks=[]
    for z in data:
        T=tuple(z['source']);matrix(T);old=z['source_names'];names=z['names'];r=z['root'];sp=z['spare']
        filename='return6325_plus2.json' if z['source_H']==18977 else 'return6325.json'
        inherited=json.loads((HERE.parent/'return6325_completion_20261103'/filename).read_text())
        need(list(T)==inherited['next_rows'] and old==inherited['next_names'] and z['source_root']==inherited['next_root'],'exact inherited named source')
        V=[old.index(x) for x in names];T=induced(T,V)
        need(T==tuple(z['after_deletion']) and strong(T),'explicit deletion endpoint')
        need(set(old)-set(names)==set(z['deleted_old_names']),'exact discarded names')
        need(names[r]==old[z['source_root']] and names[sp]==6 and 31 in names,'retained root, inport, spare')
        states.append(('source',z['source_H'],z['budget'],tuple(z['source']),old,z['source_root'],z['source_spare']))
        for step in z['trace']:
            before=T;N=len(T);local=names+[step['helper']];X=birth(T,step['cut'])
            P,Q=[[local.index(a) for a in path] for path in step['paths']]
            def frame(path):
                j=path.index(r);return (path[j-1] if j else None,path[j+1] if j+1<len(path) else None)
            need(frame(P)==frame(Q),'common live-root entry/exit')
            Z=native(native(X,P),Q);T=induced(Z,list(range(N)))
            diff={(i,j) for i in range(N) for j in range(i+1,N) if ((T[i]>>j)^(before[i]>>j))&1}
            arc=tuple(sorted(names.index(a) for a in step['arc']))
            need(diff=={arc} and r not in arc,'exact isolated arc, unchanged root star')
            need(T==tuple(step['rows']) and strong(T),'literal strong completed repair')
        need(T==tuple(z['target']),'literal target')
        U=birth(birth(T,1<<r),1<<len(T));need(U==tuple(z['next_rows']) and strong(U),'literal strong next child')
        states.extend([('target',z['source_H'],z['budget'],T,names,r,sp),('next',z['source_H'],z['budget'],U,z['next_names'],z['next_root'],z['next_spare'])])
        checks.append(dict(H=z['source_H'],offset=2*z['budget'],native_guards=2*len(z['trace']),deleted_old=z['deleted_old_names']))
    text=''.join(' '.join(map(str,[len(T),r,s,*T]))+'\n' for _,_,_,T,_,r,s in states)
    p=subprocess.run([sys.argv[1]],input=text,text=True,capture_output=True,check=True);result=[]
    for state,line in zip(states,p.stdout.splitlines(),strict=True):
        stage,H,b,T,names,r,s=state;v=list(map(int,line.split()));h,P=v[:2];ports=[v[2:7],v[7:12]];E1,S1,q,F,S2=ports[0]
        need((S2-E1)%2==0,'integer golden observer');g=(h,(S2-E1)//2,(S2+E1)//2)
        if stage=='target':need(h==H and F-H==b,'independent returned readiness')
        if stage=='next':need(h==3*H+2*b,'independent literal ternary count')
        result.append(dict(stage=stage,source_H=H,budget=b,rows=T,names=names,ports=[names[r],names[s]],H=h,P2=P,
            port_fields=['E1','S1','deletion_H','E2','S2'],port_data=ports,W2_end_start=v[12:16],golden=g,charts=chart(g)))
    # Existing literal hostile: independent source count and one distinguishing cut.
    control=json.loads((HERE.parent/'return6325_completion_20261103/halving_response.json').read_text())['hostile']
    for T,h in zip(control['rows'],control['responses']):
        need(hp(tuple(T))==462915 and hp(birth(tuple(T),256))==h,'literal 81/27 named-cut hostile')
    from phases_bridge import five_charts
    for z in result:
        z['five_Q4']=five_charts(z['golden']);T=z['rows'];r=z['names'].index(z['ports'][0]);I=[i for i in range(len(T)) if T[i]>>r&1]
        if len(I)==1:
            keep=[i for i in range(len(T)) if i!=r];p=keep.index(I[0]);U=induced(T,keep)
            request=' '.join(map(str,[len(U),p,p,*U]))+'\n'
            v=list(map(int,subprocess.run([sys.argv[1]],input=request,text=True,capture_output=True,check=True).stdout.split()))
            a,s,d,F=v[0],v[3],v[4],v[5];need(a+F-s==z['H'],'parent clone interface')
            z['clone_law']=dict(parent=[a,s,d,F],coefficient=2*(F-d),constant=z['H']-2*(F-d),defect_first=s,defect_after_first=d)
    report=dict(status='FINITE_EXACT_INDEPENDENT_MATRIX_AND_FORWARD_ENDPOINT_REPLAY',branches=checks,states=result,
        joint_scope='Exact ordered two-port two-bin prefix plus full literal relation; not claimed closed or minimal',
        F343_scope='Current observer vector codec, not inherited interval ancestry; old clocks remain frozen',hostile=control)
    (HERE/'audit_returns.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS four native returns, joint ports and field charts',flush=True)
if __name__=='__main__':main()
