"""Independent literal replay, endpoint DP, three named ports and field codecs."""
from pathlib import Path
import sys,json,subprocess,importlib.util
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'return6325_completion_20261103'))
from audit import need,matrix,birth,native,strong,hp
from transport import chart
spec=importlib.util.spec_from_file_location('phase_bridge',HERE.parent/'return18979_palindrome_20261103/phases_bridge.py');bridge=importlib.util.module_from_spec(spec);spec.loader.exec_module(bridge)
def induced(T,V):return tuple(sum(((T[a]>>b)&1)<<j for j,b in enumerate(V)) for a in V)
def query(T,r,s):
    return list(map(int,subprocess.run([sys.argv[1]],input=' '.join(map(str,[len(T),r,s,*T]))+'\n',text=True,capture_output=True,check=True).stdout.split()))
def main():
    old={z['next_H']:z for z in json.loads((HERE.parent/'return18979_palindrome_20261103/returns.json').read_text())}
    data=json.loads((HERE/'moving_returns.json').read_text());states=[];branches=[]
    for z in data:
        inherited=old[z['source_H']];T=tuple(z['source']);names=z['source_names'][:];matrix(T)
        need(list(T)==inherited['next_rows'] and names==inherited['next_names'] and z['source_root']==inherited['next_root'],'exact inherited source and mark')
        oldmark=names[z['source_root']]
        if z['permanent_birth']:
            b=z['permanent_birth'];T=birth(T,b['cut']);names.append(b['name']);need(T==tuple(b['rows']) and strong(T),'literal permanent birth')
        r=names.index(z['moved_mark_name']);need(names[r]!=oldmark,'actual moved mark')
        need([names[i] for i in range(len(T)) if T[i]>>r&1]==[oldmark],'old mark retained as unique entry port')
        V=[names.index(a) for a in z['names']];T=induced(T,V);names=z['names'];r=z['root'];sp=z['spare']
        need(T==tuple(z['after_deletion']) and strong(T),'literal deletion endpoint')
        need(set(z['source_names'])-set(names)==set(z['deleted_old_names']),'all old deletions accounted')
        need(names[r]==z['moved_mark_name'] and names[sp]==6 and oldmark in names,'three physical named ports retained')
        states.append(('source',z['source_H'],z['budget'],tuple(z['source']),z['source_names'],[oldmark,6]))
        for step in z['trace']:
            before=T;N=len(T);local=names+[step['helper']];X=birth(T,step['cut']);P,Q=[[local.index(a) for a in p] for p in step['paths']]
            def frame(p):
                j=p.index(r);return (p[j-1] if j else None,p[j+1] if j+1<len(p) else None)
            need(frame(P)==frame(Q),'common new-root entry/exit')
            Z=native(native(X,P),Q);T=induced(Z,list(range(N)))
            diff={(i,j) for i in range(N) for j in range(i+1,N) if ((T[i]>>j)^(before[i]>>j))&1};arc=tuple(sorted(names.index(a) for a in step['arc']))
            need(diff=={arc} and r not in arc,'one isolated edge, preserved marked star')
            need(T==tuple(step['rows']) and strong(T),'literal strong completed repair')
        need(T==tuple(z['target']),'actual returned relation');U=birth(birth(T,1<<r),1<<len(T));need(U==tuple(z['next_rows']) and strong(U),'literal strong child')
        states.extend([('target',z['source_H'],z['budget'],T,names,[names[r],6,oldmark]),('next',z['source_H'],z['budget'],U,z['next_names'],[z['next_names'][z['next_root']],6,oldmark])])
        branches.append(dict(H=z['source_H'],offset=2*z['budget'],next_H=z['next_H'],native_guards=2*len(z['trace']),old_mark=oldmark,new_mark=z['moved_mark_name'],permanent_birth=z['permanent_birth'] is not None,deleted_old=z['deleted_old_names']))
    result=[]
    for stage,H,b,T,names,ports in states:
        pairs=[];g=None;h=None
        for i in range(len(ports)):
          for j in range(i+1,len(ports)):
            r,s=names.index(ports[i]),names.index(ports[j]);v=query(T,r,s);need(h is None or h==v[0],'consistent full count');h=v[0]
            if i==0 and j==1:
                need((v[6]-v[2])%2==0,'integer golden vector');g=[h,(v[6]-v[2])//2,(v[6]+v[2])//2]
                if stage=='target':need([h,v[3],v[4],v[5]]==next(z['profile'] for z in data if z['source_H']==H and z['budget']==b) and v[5]-h==b,'independent ready profile')
                if stage=='next':need(h==3*H+2*b,'independent ternary count')
                if stage=='source':need(h==H,'source count')
            pairs.append(dict(ports=[ports[i],ports[j]],P2=v[1],fields=['E1','S1','deletion_H','E2','S2'],port_data=[v[2:7],v[7:12]],ordered_W2=v[12:16]))
        result.append(dict(stage=stage,source_H=H,budget=b,H=h,rows=T,names=names,ports=ports,joint_pairs=pairs,golden=g,charts=chart(g),five_Q4=bridge.five_charts(g)))
    control=json.loads((HERE.parent/'return6325_completion_20261103/halving_response.json').read_text())['hostile']
    for T,h in zip(control['rows'],control['responses']):need(hp(tuple(T))==462915 and hp(birth(tuple(T),256))==h,'81/27 named-cut hostile')
    result=dict(status='FINITE_EXACT_INDEPENDENT_NATIVE_AND_FORWARD_PORT_REPLAY',branches=branches,states=result,hostile=control,
        scope='All pairwise ordered two-bin responses on three physical ports, plus full relation; not a closed/minimal tensor quotient.',
        clock_scope='Inherited ancestral clocks frozen; current F49/F343 observer codecs and five Q4 presentations recomputed. No free phase reset.')
    (HERE/'audit_moving.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS seven actual moving-mark returns and hostile control')
if __name__=='__main__':main()
