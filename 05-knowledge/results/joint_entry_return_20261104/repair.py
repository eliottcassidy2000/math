"""Exact controls for the paid triangle repair and its joint-entry boundary."""
from pathlib import Path
from itertools import permutations,combinations
import json,sys,subprocess,os,hashlib,shutil,tempfile,importlib.util
from entry import tables,need
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parent
spec=importlib.util.spec_from_file_location('native_model',RESULTS/'adaptive_metagraph_reformulation_20261103/metagraph_probe.py')
native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native)

def strong(T):
    for v in range(len(T)):
        reached=1<<v
        while True:
            nxt=reached
            for i in range(len(T)):
                if reached>>i&1:nxt|=T[i]
            if nxt==reached:break
            reached=nxt
        if reached!=(1<<len(T))-1:return False
    return True

def profile(T,w):
    H,J,S=tables(T);full=(1<<len(T))-1
    return [H[full],S[w],H[full^(1<<w)],sum(row[w] for row in J)]

def dial(T,w):
    X=tuple(T);last=w
    for _ in range(3):
        X=native.extend(X,set(range(len(X)))-{last});last=len(X)-1
    return X

def guard(T,w):
    n=len(T)
    cycle=next((w,)+p for p in permutations([v for v in range(n) if v!=w])
               if native.legal(T,(w,)+p) and T[p[-1]]>>w&1)
    X=dial(T,w);a,b,c=n,n+1,n+2
    path=(b,c,a,cycle[-1],w)+cycle[1:-1]
    U,V,W,p,q,case=native.traced_gadget(X,c,a,path)
    Y=native.restrict(W,n+3)
    need(Y==native.flip(X,c,a) and strong(Y),'Literal strong paid repair')
    def frame(word):
        i=word.index(w);return (word[i-1] if i else None,word[i+1] if i+1<len(word) else None)
    need(frame(p)==frame(q),'Old component mark common path frame')
    return X,Y,dict(source=list(T),mark=w,source_cycle=cycle,helper_out=U[-1],paths=[p,q],target=list(Y))

def endpoint_exe():
    source=RESULTS/'return56941_two_feedback_20261103/endpoint_checked.cpp'
    suffix='.exe' if os.name=='nt' else ''
    exe=Path(tempfile.gettempdir())/('joint-entry-endpoint-'+hashlib.sha256(source.read_bytes()).hexdigest()[:16]+suffix)
    if not exe.exists():
        compiler=shutil.which('c++') or shutil.which('g++') or str(Path.home()/'scoop/apps/gcc/current/bin/g++.exe')
        subprocess.run([compiler,'-O2','-std=c++20',str(source),'-o',str(exe)],check=True)
    return exe

def endpoint(T,w,u):
    data=' '.join(map(str,[len(T),w,u,*T]))+'\n'
    return list(map(int,subprocess.run([str(endpoint_exe())],input=data,text=True,capture_output=True,check=True).stdout.split()))

def main():
    tested=0;guards=0;unique_dials=0
    for n in range(3,5):
        pairs=list(combinations(range(n),2))
        for code in range(1<<len(pairs)):
            T=native.decode(n,code)
            if not strong(T):continue
            for w in range(n):
                h,s,q,F=profile(T,w);X,Y,program=guard(T,w)
                need(profile(Y,w)==[h+4*F,q,q,8*F-q],'Full paid-repair profile')
                need(profile(X,0)[0]==3*h+3*F-s,'Inherited N3 count')
                _,J,_=tables(T);HX,JX,SX=tables(X)
                need(SX[:n]==[row[w] for row in J],'Joint entry equals old-start counts after N3')
                if s==1:
                    need((1 in SX)==(q==1),'Reentry boundary')
                    unique_dials+=1
                tested+=1;guards+=2
    # Repeated q=1 repair at the same physical old mark.
    T=(2,4,1);chain=[];w=0
    for step in range(3):
        old=profile(T,w);X,Y,program=guard(T,w);new=profile(Y,w)
        need(new[1:3]==[1,1],'Reusable unique-start entry')
        need(new==[old[0]+4*old[3],1,1,8*old[3]-1],'Exact recurrence')
        for s in range(1,12):
            for t in range(1,12):
                need(new[0]+(1<<t)*new[3] not in (3*(old[0]+(1<<s)*old[3])+2,3*(old[0]+(1<<s)*old[3])+4),'Ternary clock hostile')
        chain.append(dict(step=step,source_profile=old,dial_H=3*old[0]+3*old[3]-1,target_profile=new,
                          count_cost=2*old[0]-old[3]-1,program=program));T=Y
    # Sharp two-applications boundary, retaining the old mark after the first.
    X=dial((2,4,1),0);XX=dial(X,0)
    need(profile(X,0)[1:3]==[1,3] and min(tables(XX)[2])>1,'Two dial applications are sharp')
    # The current source is exactly absorption of its marked component.
    zero=json.loads((RESULTS/'triangle_return_170807_20261103/zero_ready.json').read_text())[1]
    D=tuple(zero['suffix_rows']);w=zero['suffix_names'].index(31);X=dial(D,w);Y=native.flip(X,len(D)+2,len(D))
    old=endpoint(D,w,zero['suffix_names'].index(6));new=endpoint(Y,w,zero['suffix_names'].index(6))
    need([old[0],old[3],old[4],old[5]]==[2195,1,53,42153],'Actual inherited component')
    need([new[0],new[3],new[4],new[5]]==[170807,53,53,337171],'Actual absorbed old-mark profile')
    absorption=dict(source_rows=list(D),target_rows=list(Y),mark_index=w,old_fields=old,new_fields=new)
    # Literal 81/27 sources plus the inherited common cut, independently counted.
    hostile=json.loads((RESULTS/'return6325_completion_20261103/halving_response.json').read_text())['hostile']
    counts=[];birth_counts=[];joint=[]
    for U in hostile['rows']:
        U=tuple(U);v=endpoint(U,0,8);counts.append(v[0]);joint.append(v)
        # Inherited birth helper uses incoming mask256: only vertex8 beats it.
        V=native.extend(U,set(range(len(U)))-{8});birth_counts.append(endpoint(V,0,8)[0])
    need(counts==[462915,462915] and birth_counts==hostile['responses'] and birth_counts[0]!=birth_counts[1],'81/27 transporter hostile')
    # The only two-entry escape, and the sharp all-size three-entry family.
    sys.path.insert(0,str(RESULTS/'zero_ready_component_returns_20261103/interface'))
    from check import phase
    D=(2,4,9,3);bridges=[]
    for m in range(10):
        X=dial(D,3);h,s,q,F=profile(D,3);col=[row[3] for row in tables(D)[1]]
        need(col[:3]==[1,1,1],'Three core ports have unit joint entry')
        expected=[5,2,3,6] if m==0 else [3+3*(1<<m),3,3,3+3*(1<<m)]
        need([h,s,q,F]==expected,'Exact bridge input profile')
        pairs=[]
        for v in range(3):
            vals=endpoint(X,v,3);need(vals[0]==3*h+3*F-s and vals[3]==1 and vals[4]>1,'Actual restored entry and retained debt')
            golden=[vals[0],(vals[6]-vals[2])//2,(vals[6]+vals[2])//2]
            pairs.append(dict(ports=[v,3],fields=vals,golden=golden,charts=phase(golden)))
        bridges.append(dict(m=m,source_rows=list(D),source_profile=[h,s,q,F],joint_entry_column=col,target_rows=list(X),pairs=pairs))
        D=native.extend(D,{0,1,2})
    out=dict(status='FINITE_EXACT',small_marked_sources=tested,native_guards=guards,unique_start_inputs=unique_dials,
             repair_chain=chain,two_dial_final_min_start=min(tables(XX)[2]),
             actual_absorption=absorption,
             sharp_entry_bridges=bridges,
             hostile=dict(source_rows=hostile['rows'],source_H=counts,joint_pairs_0_8=joint,incoming_cut=256,responses=birth_counts))
    (HERE/'repair.json').write_text(json.dumps(out,indent=2)+'\n')
    print('marked small controls',tested,'native guards',guards,'chain',[(z['source_profile'][0],z['target_profile'][0]) for z in chain],
          'hostile',birth_counts)

if __name__=='__main__':main()
