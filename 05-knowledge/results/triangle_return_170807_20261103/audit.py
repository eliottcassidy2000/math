"""Independent backward Hamilton DP, native replay and retained port charts."""
from pathlib import Path
from itertools import combinations,permutations,product
import importlib.util,json,sys,subprocess,os
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'return56941_two_feedback_20261103'))
from audit_moving import need,matrix,birth,native,strong,chart,bridge
BIN=Path(sys.argv[1])
def run(exe,text):return subprocess.run([str(BIN/exe)],input=text,text=True,capture_output=True,check=True).stdout
def hps(Ts):return list(map(int,run('hp.exe',''.join(' '.join(map(str,[len(T),*T]))+'\n' for T in Ts)).split()))
def inspect(Ts):return [json.loads(x) for x in run('inspect.exe',''.join(' '.join(map(str,[len(T),*T]))+'\n' for T in Ts)).splitlines()]
def endpoint(T,r,s):return list(map(int,run('endpoint.exe',' '.join(map(str,[len(T),r,s,*T]))+'\n').split()))
def substitute(T,blocks,As):
    groups=[(vs,A) for vs,A in zip(blocks,As)]+[([v],(0,)) for v in range(len(T)) if all(v not in b for b in blocks)]
    loc=[(i,j) for i,(vs,A) in enumerate(groups) for j in range(len(A))]
    return tuple(sum(1<<b for b,(k,l) in enumerate(loc) if (i==k and groups[i][1][j]>>l&1) or (i!=k and T[groups[i][0][0]]>>groups[k][0][0]&1)) for i,j in loc)
def nonempty(T):
    n=len(T);h=c=0
    def path(p):return all(T[a]>>b&1 for a,b in zip(p,p[1:]))
    for p in permutations(range(n)):
        h+=path(p)
        for k in range(1,n):c+=path(p[:k]) and path(p[k:])
    return h,c
def main():
    z=json.loads((HERE/'interface.json').read_text());T=tuple(z['source']);names=z['source_names'];r=z['source_root'];sp=z['source_spare'];states=[('source',T,names)]
    for step in z['trace']:
        N=len(T);local=names+[step['helper']];X=birth(T,step['cut']);P,Q=[[local.index(a) for a in p] for p in step['paths']]
        Y=native(X,P);V=native(Y,Q);U=tuple(x&((1<<N)-1) for x in V[:N]);matrix(U)
        changed=[{names[i],names[j]} for i in range(N) for j in range(i+1,N) if ((T[i]>>j)^(U[i]>>j))&1]
        need(changed==[set(step['arc'])] and U==tuple(step['rows']) and strong(U),'native literal edge and strong repair')
        def frame(p,v):
            j=p.index(v);return (p[j-1] if j else None,p[j+1] if j+1<len(p) else None)
        need(all(frame(P,v)==frame(Q,v) for v in (r,sp)),'joint root/spare path frame')
        T=U;states.append(('edge'+str(len(states)),T,names))
    T=birth(T,z['birth_cut']);need(T==tuple(z['target']) and strong(T),'strong literal final birth');states.append(('target',T,z['target_names']))
    counts=hps([T for _,T,_ in states]);need(counts[0]==56935 and counts[-2:]==[67599,170807],'independent backward counts')
    joint=[]
    for (stage,T,labels),H in zip(states,counts):
        ports=[50,6,51];pairs=[];g=None
        for a,b in combinations(ports,2):
            v=endpoint(T,labels.index(a),labels.index(b));need(v[0]==H,'independent endpoint H agrees')
            pairs.append(dict(ports=[a,b],P2=v[1],port_data=[v[2:7],v[7:12]],ordered_W2=v[12:]))
            if (a,b)==(50,6):
                need((v[6]-v[2])%2==0,'integral golden coordinates');g=[H,(v[6]-v[2])//2,(v[6]+v[2])//2]
        joint.append(dict(stage=stage,rows=T,names=labels,H=H,joint_pairs=pairs,golden=g,charts=chart(g),five_Q4=bridge.five_charts(g)))
    # Preserve the original named hostile, not only its equal unrooted H.
    hostile=json.loads((HERE.parent/'return6325_completion_20261103/halving_response.json').read_text())['hostile']
    need(hps(hostile['rows'])==[462915,462915],'81/27 common source H')
    need(hps([birth(tuple(T),256) for T in hostile['rows']])==hostile['responses'],'81/27 distinguished actual cut')
    deletions=[]
    for side,U in enumerate(hostile['rows']):
        cases=[]
        for a,b in combinations(range(len(U)),2):
            V=[i for i in range(len(U)) if i not in (a,b)]
            cases.append(tuple(sum(((U[v]>>w)&1)<<j for j,w in enumerate(V)) for v in V))
        hh=hps(cases);need(min(hh)==7821,'81/27 pair-deletion floor')
        group=[]
        for t in range(3 if side==0 else 1):
            for shifts in product(range(3),repeat=3):
                perm=tuple(4*((b+t)%3)+(0 if q==0 else 1+(q-1+shifts[b])%3) for b in range(3) for q in range(4))
                need(all(((U[i]>>j)&1)==((U[perm[i]]>>perm[j])&1) for i in range(12) for j in range(12)),'actual hostile group action');group.append(perm)
        pool=set(group);need(len(pool)==(81 if side==0 else 27),'distinct inherited symmetry elements')
        need(all(tuple(g[h[i]] for i in range(12)) in pool for g in group for h in group),'literal group closure')
        deletions.append(dict(minimum=min(hh),counts=hh,subgroup_order=len(pool),guaranteed_positive_increment=2*min(len(pool),min(hh))))
    (HERE/'hostile_pair_deletions.json').write_text(json.dumps(deletions,indent=2)+'\n')
    zero=[]
    for U,labels,root,inport,prefix,expected in [(tuple(z['source']),names,50,39,[51],[18977,18979]),(tuple(z['target']),z['target_names'],39,31,[50,51],[2195,42153])]:
        rr=labels.index(root);ww=labels.index(inport);aa=[labels.index(v) for v in prefix];bb=[i for i in range(len(U)) if i!=rr and i not in aa]
        def induced(V):return tuple(sum(((U[i]>>j)&1)<<b for b,j in enumerate(V)) for i in V)
        A=induced(aa);B=induced(bb)
        need([i for i in range(len(U)) if U[i]>>rr&1]==[ww],'unique actual in-port')
        need(strong(B) and all(U[a]>>b&1 for a in aa for b in bb),'canonical SCC suffix and prefix order')
        qa=endpoint(A,0,0);qb=endpoint(B,bb.index(ww),bb.index(labels.index(6)));qr=endpoint(U,rr,labels.index(6))
        need([qb[0],qb[5]]==expected and qr[5]==qr[0],'zero-ready component fields')
        need(qr[0]==qa[0]*qb[0]+qa[1]*qb[5],'literal Z1 formula')
        ports=list(dict.fromkeys([root,6,50,51,inport]));pairs=[]
        for a,b in combinations(ports,2):
            vv=endpoint(U,labels.index(a),labels.index(b));pairs.append(dict(ports=[a,b],P2=vv[1],port_data=[vv[2:7],vv[7:12]],ordered_W2=vv[12:]))
        g=[qr[0],(qr[6]-qr[2])//2,(qr[6]+qr[2])//2]
        zero.append(dict(root=root,inport=inport,prefix_names=prefix,suffix_names=[labels[i] for i in bb],suffix_rows=B,H=qr[0],h_A=qa[0],P_A2=qa[1],h_B=qb[0],F_w=qb[5],joint_pairs=pairs,golden=g,charts=chart(g),five_Q4=bridge.five_charts(g)))
    (HERE/'zero_ready.json').write_text(json.dumps(zero,indent=2)+'\n')
    # Every pruned-in second-birth prefix gets independent literal H replay.
    pref=json.loads((HERE/'three_births.json').read_text());source=tuple(z['source'])
    prefTs=[birth(birth(source,a),b) for a,b,_,_ in pref['prefixes']]
    need(hps(prefTs)==[v[3] for v in pref['prefixes']],'all192 prefix counts independently replayed')
    need(not pref['found'] and pref['weak_hits']==0 and len(prefTs)==192,'declared complete birth receipt')
    # Exact run coefficients from literal TT1/TT2 inputs, plus C3 controls.
    atoms=[(0,),(2,0),(6,4,0),(2,4,1)];runrows=[];cases=[]
    for A,B in product(atoms,repeat=2):
        S=substitute(source,[[3,5],[6,4]],[A,B]);V=substitute(tuple(z['target']),[[3,5],[6,4]],[A,B]);cases.extend([S,V]);runrows.append(dict(A=A,B=B,hcA=nonempty(A),hcB=nonempty(B)))
    hs=hps(cases)
    for i,row in enumerate(runrows):
        a,c=row['hcA'];b,d=row['hcB'];expected=[12631*a*b+4449*c*b+8805*a*d+4449*c*d,37943*a*b+13297*c*b+26541*a*d+13297*c*d]
        need(hs[2*i:2*i+2]==expected,'independent all-module coefficients');row.update(counts=expected,offset=expected[1]-3*expected[0])
    # Forward endpoint replay is independent for the four coefficient bases.
    low=[j for j,row in enumerate(runrows) if len(row['A'])<=2 and len(row['B'])<=2]
    forward=inspect([cases[2*j+k] for j in low for k in (0,1)])
    need([r['H'] for r in forward]==[hs[2*j+k] for j in low for k in (0,1)],'second coefficient path')
    # General triangle obstruction: every labelled strong quotient of order3,4.
    modules=[]
    for n in (3,4):
        edges=list(combinations(range(n),2))
        for mask in range(1<<len(edges)):
            U=[0]*n
            for k,(a,b) in enumerate(edges):
                if mask>>k&1:U[a]|=1<<b
                else:U[b]|=1<<a
            if strong(tuple(U)):
                for v in range(n):modules.append(substitute(U,[[v]],[(2,4,1)]))
    observed=inspect(modules);need(all(x['plus2']==x['plus4']==0 for x in observed),'proper strong triangle universal finite controls')
    boundary=inspect([(2,4,1)])[0];need(boundary['plus2']>0,'proper-module hypothesis is necessary')
    weak_boundary=inspect([(10,12,9,0)])[0];need(weak_boundary['plus2']>0,'strong-source hypothesis is necessary')
    saved=json.loads((HERE.parent/'return56941_two_feedback_20261103/feedback_ports.json').read_text())['states'];Urows=[];samples=[]
    for record in saved:
        if record['stage']!='second_port':continue
        U=record['rows'];ns=record['names'];tri=[ns.index(a) for a in (32,22,'feedback_q')];outside=[i for i in range(len(U)) if i not in tri]
        need(all(U[tri[i]]>>tri[(i+1)%3]&1 for i in range(3)),'literal triangle cycle')
        need(all(len({(U[v]>>j)&1 for v in tri})==1 for j in outside),'literal module incidence')
        if 2<=record['m']<=4:Urows.append(U);samples.append(record['m'])
    ur=inspect(Urows);need(all(min(v-x['H'] for v in x['E2'])==3 and not x['plus2'] and not x['plus4'] for x in ur),'all-root floor and actual birth obstruction')
    report=dict(status='FINITE_EXACT_SAME_AGENT_INDEPENDENT_ALGORITHMS',counts=counts,states=joint,hostile=hostile,
        prefix_counts_replayed=len(prefTs),run_law=runrows,proper_triangle_controls=len(modules),whole_triangle_hostile=boundary,uniform_samples=[dict(m=m,**v) for m,v in zip(samples,ur)],
        phase_scope='Current F49/charge and F343 vector codec, all five Q4 charts. Frozen ancestral clocks remain linked; no free source reset.',
        scope='Independent source and prefix counts, literal guards, exact coefficient bases and stated finite controls. Full discovery menus are in separate receipts.')
    (HERE/'audit.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS native170807, joint phases,81/27,192 prefixes,run law and triangle controls')
if __name__=='__main__':main()
