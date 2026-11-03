"""Universal child-profile transport, with actual nontransitive children."""
import sys,json,hashlib,importlib.util
from pathlib import Path
HERE=Path(__file__).resolve().parent;OLD=HERE.parent/'adaptive_four_returns_20261103'
sys.path.insert(0,str(OLD))
from audit import profile,tables,native,birth,need,legal,field_frames,ancestry,hpaths
from port_rank import triangle_matching

def substitute(T,q,S):
    ad=[(v,j) for v in range(len(T)) for j in range(len(S) if v==q else 1)]
    rows=tuple(sum(1<<b for b,(w,k) in enumerate(ad)
        if (v==w==q and S[j]>>k&1) or (v!=w and T[v]>>w&1)) for v,j in ad)
    return rows,ad
def expand(p,ad,q,R):
    return tuple(ad.index((v,j)) for v in p for j in (R if v==q else (0,)))
def child_profile(S):
    H=tables(S)[0];full=(1<<len(S))-1;previous=[0]*(full+1);previous[0]=1;answer=[]
    for k in range(1,4):
        current=[0]*(full+1)
        for mask in range(1,full+1):
            sub=mask
            while sub:
                current[mask]+=H[sub]*previous[mask^sub];sub=(sub-1)&mask
        previous=current;answer.append(current[full])
    tables.cache_clear();return answer
def observer(coeff,Q):
    h,e,s=[sum(a*b for a,b in zip(row,Q)) for row in coeff]
    need((s-e)%2==0,'integer golden observer');return h,(s-e)//2,(s+e)//2
def run():
    programs=json.loads((OLD/'lift_certificates.json').read_text());reports=json.loads((OLD/'audit.json').read_text())['lifts']
    sources={row['H']:tuple(row['source']) for row in json.loads((OLD/'branches.json').read_text())}
    children={'singleton':(0,),'TT2':(2,0),'C3':(2,4,1),'strong4':(2,4,9,3),'source_C3':(14,4,8,2)}
    results=[]
    for name,S in children.items():
        Q=child_profile(S);R=next(hpaths(S));need(Q[1]%2==0 and Q[2]%6==0,'ordered-component factorial guards')
        for d,report in zip(programs,reports):
            h=d['source_H'];q=d['port'];c=d['offset'];T=sources[h];W=tuple(d['target'])
            C,ad=substitute(T,q,S);D,wd=substitute(W,q,S)
            p=expand(d['first_path'],ad,q,R);s=expand(d['second_path'],ad,q,R)
            legal(C,p);U=native.native(C,p);legal(U,s);U=native.native(U,s)
            mask=sum(1<<i for i,a in enumerate(ad) if d['mask']>>a[0]&1)
            need(birth(U,mask)==D,'literal arbitrary-child complement square')
            actual=[]
            for X,xd in ((C,ad),(D,wd)):
                actual.append(profile(X,xd.index((2,0)))['observer']);tables.cache_clear()
            predicted=[observer(coeff,Q) for coeff in report['run_coefficients']]
            need(actual==predicted,'independent full nontransitive-child joint responses')
            need(actual[1][0]==3*actual[0][0]+c*Q[0],'universal homogeneous correction')
            results.append(dict(child=name,Q=Q,base=h,offset=c,observers=actual,constant_offset_error=c*(Q[0]-1)))
        print('arbitrary child',name,Q,'six literal squares passed',flush=True)
    # Constructor imported from the established packet; group orders/ranks are inherited.
    sys.path.insert(0,str(HERE.parent/'overnight_rank_return_20261029'))
    from atom_precision import hostile
    controls=[]
    for flip in (False,True):
        S=tuple(hostile(flip));Q=child_profile(S);ranks=[triangle_matching(S,q)[0] for q in range(12)]
        controls.append(dict(flip=flip,Q=Q,triangle_matching_ranks=ranks))
    need(controls[0]['Q']==controls[1]['Q'],'81/27 have equal first three total child responses')
    need(controls[0]['triangle_matching_ranks']==controls[1]['triangle_matching_ranks'],'labelled-triangle rank hostile')
    need(controls[0]['Q'][1:]==[9011898,80614764],'inherited two/three-run masses')
    embeddings=[]
    for report in reports:
        pair=[observer(coeff,controls[0]['Q']) for coeff in report['run_coefficients']]
        embeddings.append(dict(base=report['source'],offset=report['offset'],observers=pair,
                               frames=[field_frames(v) for v in pair],ancestry=ancestry(report['source'])))
    out=dict(literal=results,hostile=controls,hostile_symbolic_embeddings=embeddings,
             scope='24-vertex target values use the proved profile substitution, not direct enumeration')
    (HERE/'substitution.json').write_text(json.dumps(out,indent=2)+'\n')
    print('substitution audit',hashlib.sha256((HERE/'substitution.json').read_bytes()).hexdigest(),controls,flush=True)
    return out
if __name__=='__main__':run()
