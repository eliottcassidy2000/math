"""Reproduce exact source-reset certificates. No removable assertions.

All exhaustive universes and theorem scopes are specified in AUDIT.md.
C++ counts literal adjacency subsets; Python uses independent recurrences
and explicit named path operations. The inherited arbitrary arc compiler
is used only to implement the new cheaper target permutation.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations,permutations
from math import gcd,lcm
from pathlib import Path
import argparse,hashlib,importlib.util,json,subprocess,sys,tempfile
from clone_probe import coefficients

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent
TMP=Path(tempfile.gettempdir())/'math-source-reset-branches-20261103'

def need(b,m):
    if not b:raise RuntimeError(m)

def load(name,path):
    sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m

old=load('reset_source_inverse',RESULTS/'golden_inverse_defects_20261103/inverse.py')
native=load('reset_arc_compiler',RESULTS/'ternary_repair_clock_20261103/native.py')

def run(cmd):return subprocess.run(list(map(str,cmd)),check=True,capture_output=True,text=True).stdout.replace('\r\n','\n').strip()
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def converse(T):return tuple(((1<<len(T))-1)^r^(1<<i) for i,r in enumerate(T))
def ear(T,cut):
    n=len(T);return tuple(r|((1<<n) if cut>>i&1 else 0) for i,r in enumerate(T))+(((1<<n)-1)^cut,)
def birth(T,r):return ear(T,1<<r)
def legal(T,p):return sorted(p)==list(range(len(T))) and all(T[a]>>b&1 for a,b in zip(p,p[1:]))
def complement(T,p):
    need(legal(T,p),'literal Hamilton guard');U=list(converse(T))
    for a,b in zip(p,p[1:]):U[a]|=1<<b;U[b]&=~(1<<a)
    return tuple(U)
def strong(T):
    for R in (T,converse(T)):
        mask=1
        while True:
            oldmask=mask
            for v in range(len(R)):
                if mask>>v&1:mask|=R[v]
            if mask==oldmask:break
        if mask!=(1<<len(R))-1:return False
    return True

def table(T):
    """Independent incoming-edge DP, indexed by subset and final vertex."""
    n=len(T);D=[[0]*n for _ in range(1<<n)]
    incoming=converse(T)
    for mask in range(1,1<<n):
        todo=mask
        while todo:
            bit=todo&-todo;todo-=bit;v=bit.bit_length()-1;rest=mask^bit
            if rest==0:D[mask][v]=1;continue
            prev=incoming[v]&rest
            while prev:
                b=prev&-prev;prev-=b;D[mask][v]+=D[rest][b.bit_length()-1]
    return D

def H(T):return 1 if not T else sum(table(T)[-1])
def profiles(T):
    end=table(T);start=table(converse(T));hs=[sum(r) for r in end];hs[0]=1;full=(1<<len(T))-1
    return [(hs[full],start[full][v],hs[full^(1<<v)],sum(end[m][v]*hs[full^m] for m in range(1,full+1))) for v in range(len(T))]
def relabel(T,order):return tuple(sum(((T[a]>>b)&1)<<j for j,b in enumerate(order)) for a in order)
def observer(T,r=0):
    # These are the actual golden coordinates of T, not of its converse.
    C=profiles(converse(T))[r];h,s,_,F=C
    need((F-s)%2==0,'integral golden observer')
    return h,(F-s)//2,(F+s)//2
def fmul(a,b):
    x,y=a;z,w=b;return ((x*z+2*y*w)%7,(x*w+y*z+4*y*w)%7)
POWER=[(1,0)]
for _ in range(47):POWER.append(fmul(POWER[-1],(0,1)))
def field(v):
    h,x,y=v;alpha=(4*(h+2*x-y)%7,y%7)
    return dict(alpha=alpha,charge=(h-y)%7,theta_exponent=POWER.index(alpha) if alpha!=(0,0) else None)
def append_chart(T,r,c):
    """Append a physical block; old -> new except new-last -> old port."""
    B=old.literal((c,));n=len(T);d=len(B);allnew=((1<<d)-1)<<n
    U=[row|allnew for row in T]+[row<<n for row in B]
    U[r]&=~(1<<(n+d-1));U[n+d-1]|=1<<r
    return tuple(U),n

def q4_record(T,r=0):
    v=observer(T,r);pred={}
    for c in ('q0','q1','q2','q4','q5'):
        U,s=append_chart(T,r,c);target=old.step(v,old.BANK[c]);actual=observer(U,s)
        need(actual==target,('all five literal Q4 charts',c,v,target,actual));pred[c]=dict(state=target,field=field(target))
    return dict(state=v,field=field(v),Q4=pred)

def read_ancestry(path):
    rows=[]
    for line in Path(path).read_text().splitlines():
        a=list(map(int,line.split()));h,n=a[:2];T=tuple(a[2:2+n]);k=a[2+n];w=tuple(a[3+n:3+n+k]);kids=a[4+n+k:]
        rows.append(dict(H=h,rows=T,roots=w,kids=kids))
    return rows

def census(exe):
    raw=TMP/'ancestry.txt';stdout=run([exe,295,raw]);rows=read_ancestry(raw)
    need(len(rows)==2820,'all-word unrestricted adaptive class count')
    hist={};controls=[]
    for h in (31,95,97):
        mask=[0]*len(rows)
        for i in reversed(range(len(rows))):
            for j in range(3):
                if rows[i]['H']==3*h+2*j:mask[i]|=1<<j
            for k in rows[i]['kids']:mask[i]|=mask[k]
        hist[h]=dict(sorted(Counter(mask[i] for i,r in enumerate(rows) if r['H']==h).items()))
    need(hist[31]=={2:1,5:2},'three H31 classes cannot branch to both')
    selected=[]
    for i,s in enumerate(rows):
        if s['H']!=31:continue
        need(H(s['rows'])==31,'independent H31 literal count');queue=[i];hits={}
        for j in queue:
            t=rows[j]
            if t['H'] in (93,95,97):hits.setdefault(t['H'],j)
            if t['H']<97:queue.extend(t['kids'])
        chosen=[rows[j] for j in hits.values()]
        for t in chosen:need(H(t['rows'])==t['H'],'independent descendant literal count')
        selected.append(dict(source=s,children=chosen))
    (TMP/'selected.json').write_text(json.dumps(selected),encoding='utf-8')
    return rows,dict(stdout=stdout,sha256=digest(raw),histograms=hist,H31_sources=selected)

def parse_native(text):
    lines=text.splitlines();need(lines[0].startswith('SUCCESS'),'actual reset found')
    paths=[tuple(map(int,l.split()[1:])) for l in lines if l.startswith('path ')]
    target=next(tuple(map(int,l.split()[1:])) for l in lines if l.startswith('target '))
    cuts={d:next(list(map(int,l.split()[1:])) for l in lines if l.startswith(name)) for d,name in ((2,'plus2'),(4,'plus4'))}
    return paths,target,cuts

def finite_resets(ear_exe,native_exe,family_exe,phase_exe,native_audit,phase_audit):
    T=(2,4,41,115,103,67,7);need(T==old.literal(('a','a','r5')),'inherited source phase')
    need(observer(T)==(31,14,25),'inherited full seed observer');first=TMP/'ear_tree.txt'
    out=run([ear_exe,4,first,*T]);tree=[list(map(int,l.split())) for l in first.read_text().splitlines()]
    need(len(tree)==13 and tree[0][5]==3,'complete four-generation one-ear obstruction')
    records=[]
    for child in tree[1:5]:
        h=child[4];U=tuple(child[7:]);need(U==ear(T,child[2]),'named initial ear')
        transcript=run([native_exe,2,*U]);need(run([native_audit,2,*U])==transcript,'independent permutation universe: reset')
        paths,V,cuts=parse_native(transcript)
        need(len(paths)==2,'one-step minimality in searched native language')
        chain=[U]
        for p in paths:chain.append(complement(chain[-1],p))
        need(chain[-1]==V and H(V)==h,'neutral literal reset endpoint')
        names=list(range(len(V)));core=V;deleted=[]
        while not strong(core):
            n=len(core);candidates=[i for i,row in enumerate(core) if row in (0,((1<<n)-1)^(1<<i))]
            need(candidates,'only certified universal-vertex deletion');j=candidates[0]
            deleted.append(dict(name=names.pop(j),kind='sink' if core[j]==0 else 'source'))
            core=relabel(core,[i for i in range(n) if i!=j])
        need(len(deleted)==(1 if h==97 else 0) and H(core)==h and 0 in names,'strong-core reset endpoint and retained root')
        strong_search=run([native_exe,'--strong',2,*U])
        need(run([native_audit,'--strong',2,*U])==strong_search,'independent permutation universe: strong target')
        need(strong_search.startswith('NO_SUCCESS') if h==97 else strong_search.startswith('SUCCESS'),'strong-target radius-two boundary')
        branches={}
        for d in (2,4):
            cut=cuts[d][0];W=ear(V,cut);corecut=sum(((cut>>v)&1)<<i for i,v in enumerate(names));WC=ear(core,corecut)
            need(H(W)==H(WC)==3*h+d and strong(WC),'strong affine child after certified padding deletion')
            branches[d]=dict(cut=cut,rows=W,core_cut=corecut,strong_rows=WC,H=3*h+d,observer=q4_record(WC,names.index(0)))
        records.append(dict(initial_cut=child[2],H=h,source=U,paths=paths,target=V,
                            intermediate_counts=[H(R) for R in chain],deleted_universal_vertices=deleted,retained_names=names,strong_core=core,
                            transport=[q4_record(R) for R in chain]+[q4_record(core,names.index(0))],branches=branches,search=transcript,strong_search=strong_search))
    family={};phases={}
    for cut in (32,40):
        U=ear(T,cut);family[cut]=run([family_exe,*U])
    need('clone_lifts 0' in family[40] and 'clone_lifts 2' in family[32],'complete untouched-clone scope')
    for cut in (32,40,61,63):
        U=ear(T,cut);h=H(U);x,y=(41,75) if h==95 else (40,75)
        phases[cut]=run([phase_exe,x,y,*U])
        need(run([phase_audit,x,y,*U])==phases[cut],'independent permutation universe: full phase')
        if h==95:need('phase_hits 0' in phases[cut],'strict original95 field phase obstructed')
        lines=phases[cut].splitlines()
        for i,line in enumerate(lines):
            if not line.startswith('HIT'):continue
            a=line.split();r=int(a[2]);v=tuple(map(int,a[4:]));P=tuple(map(int,lines[i+1].split()[1:]));Q=tuple(map(int,lines[i+2].split()[1:]));V=tuple(map(int,lines[i+3].split()[1:]))
            need(complement(complement(U,P),Q)==V and observer(V,r)==v,'independent full phase witness replay')
            need(h==97 and not strong(V) and 0 in V and v[1]==v[2],'weak zero-endpoint phase boundary')
    return dict(source=T,source_transport=q4_record(T),one_ear_census=out,ear_tree_sha256=digest(first),resets=records,clone_template_scans=family,strict_phase_scans=phases)

def block_profile(z,k,A):
    h,s,q,F=z;a=h-s;d=F-q;B=k-A-1
    return (a+q+(2**k-1)*d,
            2**B*d if A else a+(2**B-1)*d,
            a+q+(2**(k-1)-1)*d,
            2**A*a+(2**(B+1)*3**A-2**A)*d+2**(A+1)*q+(s-q if B==0 else 0))

def reset_divisor(z,k,A):
    h,s,q,F=z;a=h-s;d=F-q
    return gcd(d,a+q,s+(2**(k-1)-1)*q) if A==k-1 else gcd(a,q,d)

def rotate_to_oldest(T,leaves,marked):
    U=list(T)
    # The tested mark is the newest source of the leaf block.
    for v in leaves:
        if v==marked:continue
        need(T[marked]>>v&1,'source mark is newest');U[marked]^=1<<v;U[v]^=1<<marked
    return tuple(U)

def block_controls():
    checked=0;strong_positions=0
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            T=[0]*n;bit=0
            for i in range(n):
                for j in range(i+1,n):
                    a,b=(i,j) if code>>bit&1 else (j,i);T[a]|=1<<b;bit+=1
            T=tuple(T);zs=profiles(T)
            for r,z in enumerate(zs):
                U=T
                for k in range(1,5):
                    U=birth(U,r)
                    if k<2:continue
                    ps=profiles(U);states=[]
                    for A in range(k):
                        v=n+k-1-A;expected=block_profile(z,k,A)
                        need(ps[v]==expected,'independent all-position block formula')
                        need(gcd(*expected)==reset_divisor(z,k,A),'oldest/non-oldest divisor law')
                        states.append((expected[0],expected[1],expected[3]));checked+=1
                    if strong(T):need(len(set(states))==k,'all strong ages distinguished by golden observer');strong_positions+=k
    z=(5,2,3,6);T=birth((2,4,1),0);need(profiles(T)[3]==z,'inherited factor family seed')
    rotations=[]
    for k in range(2,8):
        U=T
        for _ in range(k):U=birth(U,3)
        marked=len(U)-1;V=rotate_to_oldest(U,range(4,len(U)),marked)
        need(H(U)==H(V)==3*(2**k+1),'actual neutral stable rotation')
        need(profiles(V)[marked]==block_profile(z,k,k-1),'oldest selected root')
        work=U;program=[]
        for v in range(4,marked):
            before=work;work,c=native.arc_command(work,marked,v);program.append(c)
            enlarged=ear(before,((1<<len(before))-1)^c['birth_out'])
            for spine in c['paths']:enlarged=complement(enlarged,spine)
            replay=tuple(row&((1<<len(before))-1) for row in enlarged[:-1])
            need(replay==work,'independent stabilized native replay')
        need(work==V and len(program)==k-1,'linear stabilized native clock')
        v=(block_profile(z,k,k-1)[0],(block_profile(z,k,k-1)[3]-3)//2,(block_profile(z,k,k-1)[3]+3)//2)
        # v belongs to the converse rooted graph.
        need(observer(converse(V),marked)==v,'root profile to golden field transport')
        rotations.append(dict(k=k,source=U,target=V,marked=marked,profile=block_profile(z,k,k-1),
                              state=v,field=field(v),native_moves=2*(k-1),peak_order=len(U)+1,program=program))
    need(rotations[0]['field']['alpha']!=rotations[2]['field']['alpha'],'gcd clock does not determine F49 phase')
    need(rotations[4]['field']['alpha']==(0,0) and gcd(*rotations[4]['profile'])==1,'zero field is not nonprimitive')
    return dict(literal_profile_checks=checked,strong_position_checks=strong_positions,rotations=rotations)

def precision_profiles(profile_exe,rows):
    path=TMP/'profiles.txt';out=run([profile_exe,TMP/'ancestry.txt',path]);records=[list(map(int,l.split())) for l in path.read_text().splitlines()]
    maxima=(0,None);bad=[];primitive=0;twobad=[];transitions=0
    address={tuple(row['roots']):i for i,row in enumerate(rows)}
    for id,r,h,s,q,F in records:
        if id or r==0:
            transitions+=1;child_h=h+F-s;child=address.get(tuple(rows[id]['roots'])+(r,))
            need(child_h>h and (child is not None)==(child_h<=295),'independent full transition closure via rooted insertion')
            if child is not None:need(rows[child]['H']==child_h and rows[child]['rows']==birth(rows[id]['rows'],r),'all literal ancestry and parent response checks')
        if gcd(h,s,q,F)>1:continue
        primitive+=1;g=gcd(F-q,h-s+q);need(g>0 and g%2==1,'positive odd precision modulus')
        power=2%g;initial=power;values=[]
        while not values or power!=initial:
            values.append(gcd(g,s+(power-1)*q));power=2*power%g
        if min(values)>1:bad.append([id,r,h,s,q,F,g,values]);continue
        first=values.index(1)+2
        if first>maxima[0]:maxima=(first,[id,r,h,s,q,F,g,values])
        if len(values)>1 and values[0]>1 and values[1]>1:twobad.append([id,r,h,s,q,F,g,values])
    need(not bad and maxima[0]==7 and len(twobad)==7,'finite primitive reset phase classification')
    example=maxima[1];id,r,h,s,q,F,g,cycle=example;T=rows[id]['rows']
    need(profiles(T)[r]==(h,s,q,F),'independent phase hostile source')
    need((h,s,q,F,g)==(179,10,41,566,105),'clock hostile parameters')
    need(transitions==25196,'complete independent birth universe')
    return dict(stdout=out,sha256=digest(path),rooted_profiles=len(records),primitive_profiles=primitive,independent_birth_transition_checks=transitions,
                no_primitive_age=bad,first_two_ages_both_nonprimitive=twobad,
                latest_first_primitive_age=maxima[0],phase_hostile=dict(rows=T,roots=rows[id]['roots'],root=r,profile=(h,s,q,F),modulus=g,period=len(cycle),divisors=cycle))

def hostile_control():
    graphs=[]
    for flip in (False,True):
        T=[0]*12
        for u,v in combinations(range(12),2):
            b,i=divmod(u,4);c,j=divmod(v,4)
            arc=(c-b)%3==1 if b!=c else i==0 or (i>0 and (j-i)%3==1)
            if flip and b==c==2:arc=not arc
            T[u if arc else v]|=1<<(v if arc else u)
        graphs.append(tuple(T))
    cyclic=lambda T:{t for t in combinations(range(12),3) if all(sum(T[v]>>w&1 for w in t if w!=v)==1 for v in t)}
    need(cyclic(graphs[0])==cyclic(graphs[1]),'same labelled cyclic triples')
    iso=load('source_reset_hostile_iso',RESULTS/'regular_parent_symmetry_20261029/verify.py')
    groups=[len(iso.iso_maps(T,T)) for T in graphs]
    counts=[H(T) for T in graphs]
    # The inherited named probe is the dual, outdegree-one ear, after a sink.
    responses=[H(ear(tuple(row|(1<<12) for row in T)+(0,),((1<<13)-1)^(1<<8))) for T in graphs]
    need(groups==[81,27] and counts==[462915,462915] and responses==[2390013,1754919],'full inherited hostile retained')
    return dict(rows=graphs,groups=groups,H=counts,cyclic_triples=len(cyclic(graphs[0])),port8_outdegree_one_after_sink=responses)

def reflection_controls():
    # New specialization of incoming Session158's reflection-sandwich theorem.
    def ending_path(T,root):
        full=(1<<len(T))-1
        @lru_cache(None)
        def visit(mask,last):
            if mask==full:return (last,)
            for u in range(len(T)):
                if not mask>>u&1 and T[u]>>last&1:
                    tail=visit(mask|1<<u,u)
                    if tail is not None:return tail+(last,)
            return None
        p=visit(1<<root,root);need(p is not None and legal(T,p),'strong ending-path guard');return p
    result=[]
    for k in (2,3):
        X=birth((2,4,1),0)
        for _ in range(k):X=birth(X,3)
        r=len(X)-1;X=rotate_to_oldest(X,range(4,len(X)),r);p=ending_path(X,r);A=complement(X,p);n=len(A)
        right=converse(A);tailmask=((1<<(n+1))-1)<<n
        S=tuple(row|tailmask for row in A)+(((1<<n)-1)<<(n+1),)+tuple(row<<(n+1) for row in right)
        P=p+(n,)+tuple(n+1+j for j in reversed(p));Y=complement(S,P)
        sigma=tuple(n+1+i if i<n else n if i==n else i-n-1 for i in range(2*n+1))
        need(all(((Y[i]>>j)&1)==((Y[sigma[j]]>>sigma[i])&1) for i in range(len(Y)) for j in range(len(Y))),'literal self-converse reflection')
        h,s,q,F=profiles(X)[r];pred=36*4**k+4*9**k+24*6**k+6*2**k+8*3**k-5
        need(pred==F*F+2*h*F+(h-s)**2==H(Y) and strong(Y),'independent reflected-source count')
        need(H(S)==H(A)**2,'literal first reflection endpoint')
        result.append(dict(k=k,source=X,marked=r,source_path=p,first_endpoint=S,reflected_path=P,strong_endpoint=Y,H=pred))
    polynomial=[1]
    for eigen in (1,2,3,4,6,9):
        out=[0]*(len(polynomial)+1)
        for i,a in enumerate(polynomial):out[i]-=eigen*a;out[i+1]+=a
        polynomial=out
    return dict(controls=result,exponential_coefficients={1:-5,2:6,3:8,4:36,6:24,9:4},minimal_rational_recurrence_order=6,annihilator_ascending=polynomial)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();TMP.mkdir(exist_ok=True)
    executables={}
    for name in ('ancestry_probe','ear_probe','native_probe','family_probe','profile_probe','phase_probe'):
        exe=TMP/(name+'.exe');run(['g++','-O3','-std=c++17',HERE/(name+'.cpp'),'-o',exe]);executables[name]=exe
    for name in ('native_probe','phase_probe'):
        exe=TMP/(name+'_audit.exe');run(['g++','-O2','-DPERMUTATION_PATHS','-std=c++17',HERE/(name+'.cpp'),'-o',exe]);executables[name+'_audit']=exe
    rows,c=census(executables['ancestry_probe'])
    clones=json.loads(run([sys.executable,'-B',HERE/'clone_probe.py',TMP/'selected.json']))
    # The all-five-Q4 checks are deliberately literal, not matrix-only.
    result=dict(scope='scoped all-size oldest-port divisor and clock theorems; exact bounded native branch repairs; global spectrum OPEN',
                ancestry=c,cloning=clones,native=finite_resets(executables['ear_probe'],executables['native_probe'],executables['family_probe'],executables['phase_probe'],executables['native_probe_audit'],executables['phase_probe_audit']),
                module_age=block_controls(),precision=precision_profiles(executables['profile_probe'],rows),hostile81_27=hostile_control(),reflection=reflection_controls())
    canonical=json.dumps(result,indent=2,sort_keys=True)+'\n';target=HERE/'receipt.json'
    if args.write:target.write_text(canonical,encoding='utf-8',newline='\n')
    else:need(target.read_text(encoding='utf-8')==canonical,'frozen receipt equality')
    print('PASS 2820 complete adaptive classes; four two-native branch repairs; all five Q4 phases; module ages and exact divisor clocks')

if __name__=='__main__':main()
