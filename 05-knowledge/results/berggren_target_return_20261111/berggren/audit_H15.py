"""Independent finite-core inverse, module slopes, and named H5993 maps."""
from pathlib import Path
import hashlib,json
from check import need,H,strong,induced,birth,root_profile,state
from audit_reset import coefficient_polynomial,inflate
from audit_defect_escape import BASE,replacement_trace,g

HERE=Path(__file__).resolve().parent;RESET=HERE.parent/'reset'


def literal(core,cut,A,B):
    left=list(map(int,bin(A)[2:])) if A else []
    right=[1-int(x) for x in reversed(bin(B)[2:])] if B else []
    l=len(left);c=len(core);n=l+c+len(right);T=[0]*(n+1)
    for i in range(n):
        for j in range(i+1,n):
            forward=(core[i-l]>>(j-l)&1) if l<=i<j<l+c else 1
            if forward:T[i]|=1<<j
            else:T[j]|=1<<i
    for i,bit in enumerate(left+[(cut>>i)&1 for i in range(c)]+right):
        if bit:T[n]|=1<<i
        else:T[i]|=1<<n
    return tuple(T)


def recover(core,cut):
    p={(A,B):H(literal(core,cut,A,B)) for A,B in
       ((0,0),(1,0),(2,0),(3,0),(0,1),(0,2),(0,3),(1,1),(2,2))}
    h=p[0,0];i=(p[3,0]-p[1,0])//2;o=(p[0,3]-p[0,1])//2
    s=p[1,0]-h-i;e=p[0,1]-h-o;d=(p[2,0]-h-2*i)//2
    L=(p[2,2]-p[2,0]-p[0,2]+h)//4
    need(p[0,2]==h+2*d+2*o,'transverse calibration')
    need(p[1,1]==h+L+i+o+s+e,'zero cross-parity coefficient')
    return(h,d,L,i,o,s,e)


def inverse(N,co):
    h,d,L,i,o,s,e=co
    need(d>0 and L>0 and min(i,o,s,e)>=0,'positive complete chart')
    def value(A,B):return h+d*(g(A)+g(B)-2)+L*A*B+i*A+o*B+s*(A&1)+e*(B&1)
    found=set();cap=max(0,(N-h)//L)
    for A in range(1,cap+1):
        for B in range(1,cap//A+1):
            if value(A,B)==N:found.add((A,B))
    for left,slope,parity in((True,i,s),(False,o,e)):
        if slope:
            for u in range(max(0,(N-h)//slope)+1):
                pair=(u,0) if left else(0,u)
                if value(*pair)==N:found.add(pair)
        else:
            need(h==d and parity==0 and N%d!=0,'entire unbounded uniform axis excluded')
    return sorted(found)


def expanded_map(kept,mapping,k):
    # Both marked placeholders are vertex 1. Expand the checked base map.
    out=[];target=[]
    for u,v in zip(kept,mapping):
        if u==1:
            need(v==1,'same marked placeholder')
            out.extend(range(1,k+1));target.extend(range(1,k+1))
        else:
            out.append(u+(k-1 if u>1 else 0))
            target.append(v+(k-1 if v>1 else 0))
    return out,target


def main():
    data=json.loads((RESET/'H15_inverse.json').read_text());calibrated=[];small=[]
    for rec in data['cores']:
        core=tuple(rec['core']);need(H(core)==(5 if rec['name']=='S4' else 15),'declared core count')
        need(len(rec['inverse'])==1<<len(core),'every core cut')
        rows=[];hits=[]
        for row in rec['inverse']:
            cut=row['cut'];co=recover(core,cut);ans=inverse(5993,co)
            need(list(co)==row['coefficients'],'independently recovered chart')
            need([list(x) for x in ans]==row['solutions'],'complete independent rectangle and axis inverse')
            rows.append(dict(cut=cut,coefficients=co,solutions=ans))
            for A,B in ans:
                T=literal(core,cut,A,B)
                need(H(T)==5993 and strong(T),'literal positive hit')
                hits.append((cut,A,B))
                if rec['name']=='S4':small.append((cut,A,B,T))
        need(len(hits)==(6 if rec['name']=='S4' else 0),'whole-menu answers')
        calibrated.append(dict(name=rec['name'],rows=rows,hits=hits))
    modules=[]
    for cut,A,B,T in small:
        producer=next(r for r in data['module_checks'] if (r['cut'],r['A'],r['B'])==(cut,A,B))
        polys=[]
        for v in range(len(T)):
            co=coefficient_polynomial(T,v)
            need(co==producer['vertices'][v]['runs'][1:],'independent complete backward run polynomial')
            need(not(co[1]==3402 and not any(co[2:])),'required ternary slope absent in entire placeholder universe')
            polys.append(co)
        modules.append(dict(cut=cut,A=A,B=B,polynomials=polys))
    need(sum(len(x['polynomials']) for x in modules)==84,'all84 marked vertices')
    target=next(T for cut,A,B,T in small if (cut,A,B)==(1,3,57))
    rec=data['actual_direct_plus2'];trace=replacement_trace(BASE,target,rec['kept'],rec['map'],rec['insertion_target_labels'])
    need(len(BASE)==12 and len(target)==13 and len(rec['kept'])==10,'exact retention scope')
    need(trace==[1997,409,1599,3067,5993],'direct target trace')
    need(rec['map'][rec['kept'].index(0)]==0 and rec['map'][rec['kept'].index(11)]==12,'same physical root and port')
    p=root_profile(target,0);need(p['H']==5993 and p['F']==5998 and p['indegree']==1 and target[12]&1,'defect5 endpoint')
    co=coefficient_polynomial(target,1);need(co[:2]==[5993,3684] and not any(co[2:]),'complete dyadic target family')
    for m in range(1,5):
        S=inflate(BASE,1,m);T=inflate(target,1,m);kept,mapping=expanded_map(rec['kept'],rec['map'],m)
        need(induced(S,kept)==induced(T,mapping),'uniform retained actual arcs')
        need(len(kept)==m+9 and len(S)==m+11 and len(T)==m+12,'uniform retention counts')
        ps,pt=root_profile(S,0),root_profile(T,0)
        need(ps['F']-ps['H']==3 and pt['F']-pt['H']==5,'uniform interface types')
        need((H(S),H(T))==(1134*(1<<m)-271,3684*(1<<m)-1375),'uniform count controls')
        need(H(T)-3*H(S)-2==282*((1<<m)-2),'exact target-preservation failure coordinate')
    ready=json.loads((RESET/'H15_ready.json').read_text());ready_controls=[]
    need(tuple(ready['source_rows'])==target,'same actual direct endpoint')
    for rec in ready['targets']:
        T=tuple(rec['rows']);rewrite=rec['witness'];r,p=rec['root'],rec['port']
        tr=replacement_trace(target,T,rewrite['kept'],rewrite['map'],rewrite['insertion_target_labels'])
        need(tr==rewrite['trace'],'literal ready rebase trace')
        need(rewrite['map'][rewrite['kept'].index(0)]==r and rewrite['map'][rewrite['kept'].index(12)]==p,'persistent root and sole old port')
        q=root_profile(T,r);need(q['H']==5993 and q['F']==5995 and q['indegree']==1 and T[p]>>r&1,'literal defect2 ready target')
        need(H(birth(birth(T,r),len(T)))==17983,'actual subsequent+4')
        ready_controls.append(dict(retained=len(rewrite['kept']),trace=tr,next_H=17983))
    paths=[RESET/name for name in('H15_inverse.py','H15_inverse.json','H15_ready.py','H15_ready.json','H15_MENU.md')]
    out=dict(status='PASS',scope='Independent2448literal signature probes; complete product-region andaxis inverses; all84 full placeholder polynomials; literal retained maps and ready positives. Ready census completeness and retention optimality not recomputed.',
             charts=calibrated,module_checks=modules,direct_trace=trace,ready_controls=ready_controls,
             source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.exists()})
    (HERE/'H15_audit.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print('PASS: complete H15/S4 inverse;84 independent module polynomials;1997->5993->ready5993->17983 literal persistent-root maps')


if __name__=='__main__':main()
