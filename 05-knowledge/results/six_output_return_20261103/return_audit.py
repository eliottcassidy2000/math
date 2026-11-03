"""Independent replay and joint all-size coefficients for the next actual returns."""
from substitution import *
from audit import clone,lift_path,raw,coefficients,evaluate,observer as clone_observer
from port_rank import strong,one_separator
def run():
    previous=json.loads((HERE/'sources.json').read_text())
    byhead={3*d['source_H']+d['offset']:d for d in previous}
    programs=json.loads((HERE/'continuations.json').read_text());reports=[]
    for d in programs:
        h=d['source_H'];p0=byhead[h];T=tuple(p0['target']);q=d['port'];W=tuple(d['target']);p=d['first_path'];s=d['second_path'];c=d['offset']
        need(q==p0['port'] and triangle_matching(T,q)[0]==triangle_matching(W,q)[0]==1,'same named capacity-one port')
        need(strong(T) and strong(W),'strong endpoints');legal(T,p);U=native.native(T,p);legal(U,s)
        need(birth(native.native(U,s),d['mask'])==W,'actual second-level literal program')
        def frame(path):
            i=path.index(q);return (path[i-1] if i else None,path[i+1] if i+1<len(path) else None)
        need(frame(p)==frame(s),'common module entry/exit')
        sampled=[[],[]]
        for m in range(1,5):
            C,ad=clone(T,q,m);D,wd=clone(W,q,m);pp=lift_path(p,ad);ss=lift_path(s,ad)
            legal(C,pp);V=native.native(C,pp);legal(V,ss);V=native.native(V,ss)
            mask=sum(1<<i for i,a in enumerate(ad) if d['mask']>>a[0]&1)
            need(birth(V,mask)==D,'full lifted target relation')
            for ix,(X,xd) in enumerate(((C,ad),(D,wd))):
                sampled[ix].append(raw(X,xd.index((2,0))));tables.cache_clear()
        cs=[]
        for side in sampled:
            cc=[coefficients([side[i][j] for i in range(3)]) for j in range(3)]
            need(all(evaluate(cc[j],4)==side[3][j] for j in range(3)),'fourth-clone independent control');cs.append(cc)
        need(cs[1][0][0]==3*cs[0][0][0]+c and cs[1][0][1]==3*cs[0][0][1] and cs[1][0][2]==0,'uniform all-size offset')
        joint=[]
        for m in range(1,7):
            a,b=clone_observer(cs[0],m),clone_observer(cs[1],m)
            need(b==(3*a[0]+c,3*a[1],3*a[2]+c),'full golden return on the two witnessed mixed words')
            joint.append(dict(m=m,input=field_frames(a),output=field_frames(b)))
        need(tables(W)[0][-1]==3*h+c==sum(1 for _ in hpaths(W)),'independent unpruned target paths');tables.cache_clear()
        report=dict(source=h,offset=c,target=3*h+c,word=[p0['offset'],c],port=q,guardian=one_separator(W,q),run_coefficients=cs,
                    frames=joint,ancestry=ancestry(p0['source_H']))
        reports.append(report);print('next uniform return',h,c,report['target'],flush=True)
    out=dict(reports=reports)
    (HERE/'return_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print('return audit',hashlib.sha256((HERE/'return_audit.json').read_bytes()).hexdigest(),flush=True);return out
if __name__=='__main__':run()
