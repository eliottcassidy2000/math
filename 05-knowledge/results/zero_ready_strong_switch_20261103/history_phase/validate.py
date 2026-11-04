"""Read-only sources; independent endpoints; explicit named gauge; literal 81/27.

python validate.py [--input records.json] [--output receipt.json]
records is a list (or {states:[...]}) of {tag,rows,names,root_name,ports}.
The default replays the four literal old stages and both zero-ready gauges.
Pairwise W2 is a bounded observer, NOT the complete graded port forest.
Literal rows/names remain in the receipt as the lossless sidecar.
"""
from pathlib import Path
from itertools import combinations,product
import argparse,json,subprocess,tempfile,hashlib

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
OLD=RESULTS/'triangle_return_170807_20261103'
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def norm(x):return json.loads(json.dumps(x))
def digest(p):return hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def induced(T,V):return [sum(((T[a]>>b)&1)<<j for j,b in enumerate(V)) for a in V]
def birth(T,cut):
    n=len(T)
    return [a|(((cut>>i)&1)<<n) for i,a in enumerate(T)]+[((1<<n)-1)^cut]
def strong(T):
    n=len(T);f=(1<<n)-1
    for r in range(n):
        seen=1<<r
        while True:
            more=seen
            for v in range(n):
                if seen>>v&1:more|=T[v]
            if more==seen:break
            seen=more
        if seen!=f:return False
    return True
def mul(a,b):
    c=[0]*5
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    for i in (4,3):c[i-1]-=c[i];c[i-3]-=c[i]
    return [x%7 for x in c[:3]]
def power(a,k):
    z=[1,0,0]
    while k:
        if k&1:z=mul(z,a)
        a=mul(a,a);k//=2
    return z
def chart(g):
    h,x,y=g;v=[a%7 for a in g];a,b,c=(x+4*(h-y))%7,y%7,(h-y)%7
    need(power(v,343)==v,'cubic codec closure')
    need([(b+c)%7,(a-4*c)%7,b]==v,'invertible phase coordinates')
    return dict(observer_mod7=v,F343_vector_codec=[v,power(v,7),power(v,49)],
                F343_codec_degree=1 if power(v,7)==v else 3,F49_charge=[a,b,c],
                F49_frobenius=[(a+4*b)%7,(-b)%7,c])
ATLAS=read(RESULTS/'golden_seam_recurrence_20261103/atom_atlas.json')
def five(g):
    out={}
    for name in ('q0','q1','q2','q4','q5'):
        M=ATLAS[name]['matrix'];v=[sum(g[i]*M[i][j] for i in range(3)) for j in range(3)]
        h,x,y=v;out[name]=dict(golden=v,phase=[4*(h+2*x-y)%7,y%7,(h-y)%7])
    return out
class Engine:
    def __init__(self,exe):self.exe=exe;self.cache={}
    def query(self,T,ports):
        key=(tuple(T),tuple(ports))
        if key not in self.cache:
            n=len(T);k=len(ports);txt=' '.join(map(str,[n,k,*T,*ports]))+'\n'
            v=list(map(int,subprocess.run([str(self.exe)],input=txt,text=True,capture_output=True,check=True).stdout.split()))
            need(len(v)==2+5*k+k*k,'endpoint output arity')
            self.cache[key]=dict(H=v[0],P2=v[1],port_data=[v[2+5*i:7+5*i] for i in range(k)],ordered_W2=[v[2+5*k+i*k:2+5*k+(i+1)*k] for i in range(k)])
        return self.cache[key]
    def hp(self,T):return self.query(T,[0])['H']
def observe(E,z):
    T=z['rows'];ns=z['names'];ports=z['ports'];root=z['root_name']
    need(len(ns)==len(T) and len(set(ns))==len(ns),'distinct physical names')
    need(len(set(ports))==len(ports) and root in ports,'distinct retained ports including root')
    q=E.query(T,[ns.index(p) for p in ports]);j=ports.index(root);e=q['port_data'][j]
    need((e[4]-e[0])%2==0,'integral golden extraction')
    g=[q['H'],(e[4]-e[0])//2,(e[4]+e[0])//2]
    pairs=[]
    for i,j in combinations(range(len(ports)),2):
        pairs.append(dict(ports=[ports[i],ports[j]],P2=q['P2'],port_data=[q['port_data'][i],q['port_data'][j]],ordered_W2=[q['ordered_W2'][a][b] for a in (i,j) for b in (i,j)]))
    out=dict(**z,H=q['H'],strong=strong(T),joint_pairs=pairs,all_ports=q,golden=g,charts=chart(g),five_Q4=five(g))
    if 'expected_H' in z:need(out['H']==z['expected_H'],'new record expected H')
    return out
def defaults():
    d=read(OLD/'interface.json');saved=read(OLD/'audit.json')['states'];states=[]
    for s in saved:
        states.append(dict(tag=s['stage'],rows=s['rows'],names=s['names'],root_name=50,ports=[50,6,51],expected_H=s['H']))
    zeros=read(OLD/'zero_ready.json')
    for j,z in enumerate(zeros):
        side='source' if j==0 else 'target'
        states.append(dict(tag='zero_'+side,rows=d[side],names=d[side+'_names'],root_name=z['root'],ports=list(dict.fromkeys([z['root'],6,50,51,z['inport']])),expected_H=z['H']))
    return states,saved,zeros
def hostile(E):
    h=read(RESULTS/'return6325_completion_20261103/halving_response.json')['hostile'];out=[]
    need(h['cut']==256,'literal named cut')
    for side,T in enumerate(h['rows']):
        need(E.hp(T)==462915,'hostile common H')
        response=E.hp(birth(T,256));need(response==h['responses'][side],'hostile cut response')
        deletions=[E.hp(induced(T,[v for v in range(12) if v not in (a,b)])) for a,b in combinations(range(12),2)]
        need(min(deletions)==7821,'literal pair-deletion floor')
        group=[]
        for t in range(3 if side==0 else 1):
            for shifts in product(range(3),repeat=3):
                perm=tuple(4*((b+t)%3)+(0 if q==0 else 1+(q-1+shifts[b])%3) for b in range(3) for q in range(4))
                need(all(((T[i]>>j)&1)==((T[perm[i]]>>perm[j])&1) for i in range(12) for j in range(12)),'literal symmetry action');group.append(perm)
        pool=set(group);need(len(pool)==(81 if side==0 else 27),'distinct subgroup elements')
        need(all(tuple(g[p[i]] for i in range(12)) in pool for g in pool for p in pool),'subgroup closure')
        out.append(dict(rows=T,H=462915,cut=256,response=response,pair_deletion_counts=deletions,subgroup_order=len(pool),guaranteed_positive_increment=2*min(len(pool),min(deletions))))
    need(out[0]['response']!=out[1]['response'],'equal-H quotient is hostile to named cut')
    return out
def complement(T,P):
    n=len(T);need(sorted(P)==list(range(n)),'native guard visits each actual vertex')
    need(all(T[a]>>b&1 for a,b in zip(P,P[1:])),'native guard is directed at time of use')
    fixed=[0]*n
    for a,b in zip(P,P[1:]):fixed[a]|=1<<b;fixed[b]|=1<<a
    f=(1<<n)-1
    return [T[a]^(f^(1<<a)^fixed[a]) for a in range(n)]
def native_replay(E,path,indices,approved):
    data=read(path);out=[]
    for index in indices:
        z=data[index];T=z['source_rows'][:];ns=z['source_names'][:]
        key=(tuple(T),tuple(ns));need(key in approved,'source is an exact previously validated named target')
        source_h=approved[key];need(E.hp(T)==source_h==z['source_H'],'source count and lineage')
        frames=[]
        for step in z['trace']:
            n=len(T);labels=ns+[step['helper']];need(step['helper'] not in ns,'fresh temporary helper')
            P,Q=[[labels.index(a) for a in p] for p in step['paths']]
            def frame(p,name):
                j=p.index(labels.index(name))
                return [labels[p[j-1]] if j else None,labels[p[j+1]] if j+1<len(p) else None]
            for name in z['common_frame_ports']:need(frame(P,name)==frame(Q,name),'common physical guard frame')
            X=birth(T,step['cut']);Y=complement(X,P);Z=complement(Y,Q);U=induced(Z,list(range(n)))
            diff=[{ns[a],ns[b]} for a,b in combinations(range(n),2) if ((U[a]^T[a])>>b)&1]
            need(diff==[set(step['arc'])],'exact single old edge changed')
            need(U==step['rows'] and strong(U),'recorded strong old-vertex endpoint')
            T=U;frames.append({str(name):frame(P,name) for name in z['common_frame_ports']})
        for b in z['prefix_births']:
            need(b['name'] not in ns,'fresh permanent prefix name');T=birth(T,b['cut']);ns.append(b['name'])
            need(T==b['rows'] and strong(T),'literal strong permanent prefix birth')
        need(T==z['target_rows'] and ns==z['target_names'],'exact final rows and names')
        need(set(z['source_names']).issubset(ns) and not z['deleted_old_names'],'all17old names retained')
        rr=ns.index(z['target_root']);ww=ns.index(z['target_inport']);aa=[ns.index(a) for a in z['prefix_names']];bb=[ns.index(a) for a in z['suffix_names']]
        need(sorted([rr,*aa,*bb])==list(range(len(T))),'literal zero-ready partition')
        need([j for j in range(len(T)) if T[j]>>rr&1]==[ww],'new unique inport')
        need(all(T[a]>>b&1 for a in aa for b in bb),'prefix precedes suffix')
        need(all(T[a]>>b&1 for j,a in enumerate(aa) for b in aa[j+1:]),'ordered transitive prefix')
        B=induced(T,bb);need(B==z['suffix_rows'] and strong(B),'actual strong marked suffix')
        q=E.query(B,[bb.index(ww)]);need([q['H'],q['port_data'][0][3]]==[z['suffix_h'],z['suffix_F']],'new suffix independent h/F')
        qr=E.query(T,[rr]);h=qr['H'];need(qr['port_data'][0][3]==h,'target root zero-ready')
        need(h==z['suffix_h']+(1<<len(aa))*z['suffix_F']==3*source_h+z['offset']==z['target_H'],'actual 3H+c and component law')
        need(z['native_complements']==2*len(z['trace']),'native command accounting')
        approved[(tuple(T),tuple(ns))]=h
        out.append(dict(source_file=path.name,index=index,source_H=source_h,H=h,offset=z['offset'],order=len(T),root=z['target_root'],inport=z['target_inport'],prefix_length=len(aa),h_B=q['H'],F_w=q['port_data'][0][3],native_complements=2*len(z['trace']),temporary_births=len(z['trace']),permanent_births=len(z['prefix_births']),common_frame_ports=z['common_frame_ports'],common_frames=frames,old_names_retained=len(z['source_names'])))
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path);ap.add_argument('--native',type=Path);ap.add_argument('--native-plan',type=Path);ap.add_argument('--branches',type=int,nargs='+',default=[0,1,4]);ap.add_argument('--output',type=Path,default=HERE/'validation.json');args=ap.parse_args()
    states,saved,zeros=defaults()
    custom=[]
    if args.input:
        custom=read(args.input)
        if isinstance(custom,dict):custom=custom['states']
    plans=read(args.native_plan) if args.native_plan else ([dict(file=str(args.native),branches=args.branches)] if args.native else [])
    with tempfile.TemporaryDirectory(prefix='marked-port-audit-') as td:
        exe=Path(td)/'endpoints.exe';subprocess.run(['g++','-O2','-std=c++20',str(HERE/'endpoint_backward.cpp'),'-o',str(exe)],check=True)
        E=Engine(exe);observations=[observe(E,z) for z in states]
        for got,old in zip(observations[:4],saved):
            for key in ('H','joint_pairs','golden','charts','five_Q4'):need(norm(got[key])==norm(old[key]),'old literal field '+key)
        for got,old in zip(observations[4:],zeros):
            for key in ('H','joint_pairs','golden','charts','five_Q4'):need(norm(got[key])==norm(old[key]),'zero-ready gauge field '+key)
            B=old['suffix_rows'];bn=old['suffix_names'];q=E.query(B,[bn.index(old['inport'])])
            need([q['H'],q['port_data'][0][3]]==[old['h_B'],old['F_w']],'suffix h/F')
            need(got['H']==old['h_A']*old['h_B']+old['P_A2']*old['F_w'],'zero-ready formula')
            rootindex=got['ports'].index(got['root_name']);need(got['all_ports']['port_data'][rootindex][3]==got['H'],'zero-ready defect zero')
        controls=hostile(E);extra=[observe(E,z) for z in custom]
        old=read(OLD/'interface.json');approved={(tuple(old['target']),tuple(old['target_names'])):170807};native=[]
        for plan in plans:native.extend(native_replay(E,Path(plan['file']),plan['branches'],approved))
    pins=[OLD/'interface.json',OLD/'audit.json',OLD/'zero_ready.json',RESULTS/'return6325_completion_20261103/halving_response.json',RESULTS/'golden_seam_recurrence_20261103/atom_atlas.json',RESULTS/'return6325_completion_20261103/transport.py']
    pins.extend(p for p in (args.input,args.native_plan) if p is not None);pins.extend(Path(p['file']) for p in plans)
    result=dict(status='FINITE_EXACT_INDEPENDENT_BACKWARD_ENDPOINT_REPLAY',states=observations,new_states=extra,native_interfaces=native,hostile_81_27=controls,source_pins={str(p.resolve().relative_to(RESULTS.parent.parent)).replace('\\','/'):digest(p) for p in pins},scope='Named rows plus all pairwise endpoint/two-slot data and current finite-field codecs. Pairwise data do not assert closure of the full graded joint-port forest. No tournament automorphism is inferred from Frobenius.')
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: six named baseline gauges; all joint pairs; F49/F343 and five Q4 charts; literal 81/27, cut256 and 132 deletion controls; '+str(len(extra))+' new named states; '+str(len(native))+' native interfaces')
if __name__=='__main__':main()
