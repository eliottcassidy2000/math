"""Exact bounded controls. Run normally and with python3 -O; no assert gates."""
from pathlib import Path
from itertools import permutations
import hashlib,json,subprocess,tempfile
from gap import count,nine_fields,pair_fields,induced,compositions,need
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent

def literal(rows,bins=1,starts=(),ends=()):
    total=0
    for order in permutations(range(len(rows))):
        for sizes in compositions(len(rows),bins):
            words=[];pos=0
            for size in sizes:words.append(order[pos:pos+size]);pos+=size
            if any(not words[b] or words[b][0]!=v for b,v in starts):continue
            if any(not words[b] or words[b][-1]!=v for b,v in ends):continue
            if all(all(rows[a]>>b&1 for a,b in zip(w,w[1:])) for w in words):total+=1
    return total

def rows_code(n,code):
    rows=[0]*n;bit=0
    for a in range(n):
        for b in range(a+1,n):
            if code>>bit&1:rows[a]|=1<<b
            else:rows[b]|=1<<a
            bit+=1
    return tuple(rows)

def birth(rows,incoming):
    n=len(rows);out=list(rows);r=0
    for v in range(n):
        if incoming>>v&1:out[v]|=1<<n
        else:r|=1<<v
    return tuple(out+[r])

def wrapper(rows,w,k):
    n=len(rows);out=list(rows)
    for j in range(k):out.append((1<<n)-1|sum(1<<(n+l) for l in range(j+1,k)))
    return birth(tuple(out),1<<w)

def nine_literal(rows,w):
    vs=[v for v in range(len(rows)) if v!=w];q=induced(rows,vs)
    return [literal(rows),literal(rows,starts=((0,w),)),literal(q),literal(rows,2,ends=((0,w),)),
            literal(rows,ends=((0,w),)),literal(rows,2),literal(rows,2,starts=((0,w),)),
            literal(rows,3,ends=((0,w),)),literal(q,2)]

def eight(rows,w):
    z=nine_literal(rows,w)
    h,s,q,F,e,P2,S2,E3,Q2=z
    return [h,e,S2,F,Q2,q,P2,literal(rows,3,starts=((0,w),))]

def phase(g):
    def mul(a,b):
        c=[0]*5
        for i,x in enumerate(a):
            for j,y in enumerate(b):c[i+j]+=x*y
        for i in (4,3):c[i-1]-=c[i];c[i-3]-=c[i]
        return tuple(x%7 for x in c[:3])
    def power(a,k):
        out=(1,0,0)
        while k:
            if k&1:out=mul(out,a)
            a=mul(a,a);k//=2
        return out
    h,x,y=g;z=tuple(v%7 for v in g);a,b,c=(x+4*(h-y))%7,y%7,(h-y)%7
    need(power(z,343)==z,'cubic Frobenius')
    return dict(observer_mod7=z,F49_charge=[a,b,c],F49_frobenius=[(a+4*b)%7,(-b)%7,c],
                F343_vector_codec=[z,power(z,7),power(z,49)],
                F343_codec_degree=next(i for i in (1,3) if power(z,7**i)==z))

def lift(v,k):
    h,s,q,F,e,P2,S2,E3,Q2=v[:9]
    H=h+(1<<k)*F;St=(1<<k)*P2-S2;Pt=2*(3**k*E3+(1<<k)*P2-S2)
    need((St-e)%2==0,'integral golden coordinates')
    return dict(H=H,P2=Pt,root_fields=[e,h,h,H,St],golden=[H,(St-e)//2,(St+e)//2])

def spare_lift(v,k):
    a2=1<<k;a3=3**k
    eu,su,qu,Fu,S2u,Awu,Auw,Bwu,G,D,Auw3,SE,Fdel=v[9:]
    return [eu+a2*D,Auw,qu+a2*Fdel,2*a2*Fu-Awu-Bwu+2*a3*G,S2u+a2*(Auw3+SE)]

def compile_cpp(source,path):
    subprocess.run(['c++','-std=c++20','-O2',str(source),'-o',str(path)],check=True,capture_output=True)

def cpp(exe,jobs):
    data=''.join(' '.join(map(str,[len(rows),w,u,*rows]))+'\n' for rows,w,u in jobs)
    return [list(map(int,line.split())) for line in subprocess.run([str(exe)],input=data,text=True,
             capture_output=True,check=True).stdout.splitlines()]

def main():
    small=[];queries=0
    for code in range(64):
        rows=rows_code(4,code);fb=(0,1)
        for m in (1,2,3):
            qs=[((),()),((),((0,0),)),(((0,0),),()),(((0,0),),((0,1),))]
            if m>1:qs.extend([(((1,1),),((0,0),)),((),((0,0),(1,1)))])
            for starts,ends in qs:
                need(count(rows,fb,m,starts,ends)==literal(rows,m,starts,ends),'gap versus literal path-bin enumeration')
                queries+=1
        z=nine_literal(rows,0);need(nine_fields(rows,fb,0)==z,'nine baseline')
        h,s,q,F,e,P2,S2,E3,Q2=z
        child=birth(rows,1)
        expect_R=[h+F-s,q,q,2*F-q,e,2*(E3+P2-S2),S2+Q2,3*E3-Q2,2*Q2]
        expect_N=[h+F-s,h-s,h,h+F,e,2*(E3+P2-S2),P2-S2,E3+P2,P2]
        need(nine_fields(child,(0,1,4),0)==expect_R,'native R nine-scalar closure')
        need(nine_fields(child,(0,1,4),4)==expect_N,'native N nine-scalar closure')
        z8=eight(rows,0);h,e,S2,F,Q2,q,P2,S3=z8
        childO=birth(rows,((1<<4)-1)^1);childSink=birth(rows,(1<<4)-1)
        expectO=[h+S2-e,q,2*S2-q,F+Q2,2*Q2,q,2*(P2-F+S3),3*S3-Q2]
        expectSink=[h,0,2*S2,F,2*Q2,q,2*P2,3*S3]
        for child,ex in ((childO,expectO),(childSink,expectSink)):
            v=nine_fields(child,(0,1,4),0)
            got=[v[0],v[4],v[6],v[3],v[8],v[2],v[5],count(child,(0,1,4),3,starts=((0,0),))]
            need(got==ex,'opposite singleton/sink eight-scalar closure')
        small.append((rows,0,1))
    with tempfile.TemporaryDirectory(prefix='zero-ready-interface-') as td:
        exe=Path(td)/'bin_fields';compile_cpp(HERE/'bin_fields.cpp',exe)
        vals=cpp(exe,small)
        for (rows,w,u),v in zip(small,vals):
            need(v[:9]==nine_literal(rows,w),'C++ endpoint DP nine fields')
            mixed=[literal(rows,ends=((0,u),)),literal(rows,starts=((0,u),)),literal(induced(rows,[0,2,3])),
                   literal(rows,2,ends=((0,u),)),literal(rows,2,starts=((0,u),)),
                   literal(rows,2,starts=((0,w),),ends=((0,u),)),literal(rows,2,starts=((0,u),),ends=((0,w),)),
                   literal(rows,2,starts=((1,w),),ends=((0,u),)),literal(rows,3,ends=((0,w),(1,u))),
                   literal(rows,2,ends=((0,w),(1,u))),literal(rows,3,starts=((0,u),),ends=((0,w),)),
                   literal(rows,3,starts=((0,u),),ends=((1,w),)),
                   literal(induced(rows,[0,2,3]),2,ends=((0,0),))]
            need(v[9:]==mixed,'C++ mixed fields versus literal counts')
        actual=json.loads((RESULTS/'triangle_return_170807_20261103/zero_ready.json').read_text());receipts=[]
        for obj in actual:
            rows=tuple(obj['suffix_rows']);names=obj['suffix_names'];w=names.index(obj['inport']);u=names.index(6);k=len(obj['prefix_names'])
            v=cpp(exe,[(rows,w,u)])[0];got=lift(v,k);saved=obj['joint_pairs'][0]
            need(got['H']==obj['H'] and got['golden']==obj['golden'],'actual marked wrapper H and golden')
            need(got['root_fields']==saved['port_data'][0] and got['P2']==saved['P2'],'actual root profile and P2')
            need(spare_lift(v,k)==saved['port_data'][1],'actual mixed spare lift')
            need(json.loads(json.dumps(phase(got['golden'])))==obj['charts'],'actual field charts')
            receipts.append(dict(root=obj['root'],inport=obj['inport'],spare=6,nine=v[:9],mixed=v[9:],
                                 **got,spare_fields=spare_lift(v,k),charts=phase(got['golden'])))
        collision=[];base=(20,1,2,7,14)
        for root in range(5):
            reached=1<<root
            while True:
                expanded=reached
                for v in range(5):
                    if reached>>v&1:expanded|=base[v]
                if expanded==reached:break
                reached=expanded
            need(reached==31,'phase collision component is strong')
        for w in (0,3):
            order=[0,1,2,3,4] if w==0 else [3,1,2,0,4]
            rows=induced(base,order);v=cpp(exe,[(rows,0,1)])[0];T=wrapper(rows,0,2)
            got=lift(v,2);profile=pair_fields(T,tuple(range(5)),7,1)
            need(profile[0]==got['H'] and profile[2:7]==got['root_fields'],'literal same-named collision wrapper')
            collision.append(dict(component_rows=rows,mark=0,spare=1,nine=v[:9],**got,charts=phase(got['golden'])))
        need([(x['nine'][0],x['nine'][3]) for x in collision]==[(9,11),(9,11)],'collision same h,F')
        need([x['golden'] for x in collision]==[[53,162,163],[53,156,157]],'collision unequal golden')
        hostile=json.loads((RESULTS/'return6325_completion_20261103/halving_response.json').read_text())['hostile'];hr=[]
        ep=Path(td)/'endpoint';compile_cpp(RESULTS/'return18979_palindrome_20261103/endpoint.cpp',ep)
        query_words=('', '1','2','3','01','11','21','31','02','12','22','32','201','301','021','221','321')
        query_jobs=[]
        for code in (0,13,42,63):
            for word in query_words:
                T=rows_code(4,code)
                for letter in word:T=birth(T,int(letter))
                query_jobs.append((T,0,1))
        for (T,r,s),profile in zip(query_jobs,cpp(ep,query_jobs)):
            need(count(T,(0,1))==profile[0],'universal17 query graph independent subset DP')
        for rows,response in zip(hostile['rows'],hostile['responses']):
            v=cpp(exe,[(rows,8,6)])[0];need(v[0]==hostile['baseline'],'81/27 common baseline')
            need(v[0]+v[3]-v[1]==response,'81/27 same distinguished native cut')
            T=wrapper(rows,8,1);profile=cpp(ep,[(T,13,6)])[0];got=lift(v,1)
            need(profile[0]==got['H'] and profile[1]==got['P2'] and profile[2:7]==got['root_fields'],'81/27 literal new wrapper endpoint DP')
            need(profile[7:12]==spare_lift(v,1),'81/27 new wrapper spare mixed lift')
            hr.append(dict(mark=8,spare=6,nine=v[:9],old_cut_response=response,**got,spare_fields=profile[7:12],charts=phase(got['golden'])))
    pins=['gap.py','bin_fields.cpp','check.py','NOTE.md']
    source_paths=[HERE/p for p in pins]+[RESULTS/'triangle_return_170807_20261103'/p for p in ('ZERO_READY.md','zero_ready.json')]+[
        RESULTS/'return6325_completion_20261103/halving_response.json',RESULTS/'return18979_palindrome_20261103/endpoint.cpp',
        RESULTS/'polarized_path_calculus_20261103/CALCULUS.md',RESULTS/'adaptive_state_induction_20261103/COMMON_OBSERVER.md']
    out=dict(scope='FINITE-EXACT bounded controls; polynomial identities proved in NOTE.md',
        small=dict(tournaments=64,orders=[4,5],gap_literal_queries=queries,native_R_N_transforms=128,opposite_singleton_sink_transforms=128,
                   cpp_mixed_literal_controls=64,static17_query_controls=68),actual=receipts,phase_collision=collision,hostile81_27=hr,
        source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths})
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
