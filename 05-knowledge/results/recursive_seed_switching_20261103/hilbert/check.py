"""Small exact Hilbert/Gray sealed-reader and named-frame controls."""
from pathlib import Path
from functools import lru_cache
from itertools import product,combinations,permutations
from fractions import Fraction
import importlib.util,sys,json,hashlib
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
OLD=RESULTS/'zero_ready_component_returns_20261103'
sys.path.insert(0,str(OLD/'inverse'))
import reader
sys.path.insert(0,str(OLD/'interface'))
from gap import count,pair_fields,induced,need
spec=importlib.util.spec_from_file_location('old_phase',OLD/'interface/check.py');oldphase=importlib.util.module_from_spec(spec);spec.loader.exec_module(oldphase)
spec=importlib.util.spec_from_file_location('native_compiler',OLD/'inverse/native.py');native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native)

@lru_cache(None)
def point(d,i):
    if d==0:return (0,0)
    s=1<<(d-1);q,j=divmod(i,4**(d-1));x,y=point(d-1,j)
    return ((y,x),(x,y+s),(x+s,y+s),(2*s-1-y,s-1-x))[q]

def gray_word(d,i):
    return format(i^(i>>1),'0'+str(2*d)+'b') if d else ''

def step(z,bit):
    a,b,c,d=z
    return (a,0,c,2*d) if bit else (a+d-b,c,c,2*d-c)

def state(seed,w,k,word):
    h,e,t,F,p,q,_,_=reader.all_fields(seed)[w]
    z=(h+(1<<k)*F,e,q,t+(1<<k)*p)
    for bit in word:z=step(z,int(bit))
    return z

def output(z):a,b,c,d=z;return a+d-b

def make(seed,seed_names,w,word,free,seal,prefix,root,order):
    need(len(word)==len(free),'free word/names length')
    at={v:i for i,v in enumerate(order)};edges=set()
    def arc(a,b):edges.add((a,b))
    for i,a in enumerate(seed_names):
        for j,b in enumerate(seed_names):
            if seed[i]>>j&1:arc(a,b)
    for j,(v,bit) in enumerate(zip(free,word)):
        for a in seed_names:
            arc(v,a) if a==w and bit=='0' else arc(a,v)
        for a in free[:j]:arc(a,v)
    U=seed_names+free
    for a in U:arc(seal,a) if a==w else arc(a,seal)
    for j,a in enumerate(prefix):
        for b in prefix[j+1:]+U+[seal]:arc(a,b)
    for a in order:
        if a!=root:arc(a,root) if a==w else arc(root,a)
    rows=tuple(sum(1<<at[b] for x,b in edges if x==a) for a in order)
    need(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(len(rows)) for j in range(i)),'literal tournament')
    return rows

def atlas(word,seed=(2,4,1),w=0,k=2):
    n=len(seed);seal=n;prefix=list(range(n+1,n+k+1));root=n+k+1;free=list(range(n+k+2,n+k+2+len(word)));order=list(range(n+k+2+len(word)))
    return make(seed,list(range(n)),w,word,free,seal,prefix,root,order),dict(mark=w,seal=seal,prefix=prefix,root=root,free=free,order=order)

def converse(T):return tuple(((1<<len(T))-1)^row^(1<<i) for i,row in enumerate(T))
def flip(T,a,b):
    U=list(T);U[a]^=1<<b;U[b]^=1<<a;return tuple(U)

def literal_native(T,path):
    need(sorted(path)==list(range(len(T))) and all(T[a]>>b&1 for a,b in zip(path,path[1:])),'independent directed Hamilton guard')
    keep={frozenset((a,b)) for a,b in zip(path,path[1:])};U=T
    for a,b in combinations(range(len(T)),2):
        if frozenset((a,b)) not in keep:U=flip(U,a,b)
    return U

def gadget(T,a,b):
    U,V,W,p,q,case,inverse=native.compile_arc(T,a,b)
    need(literal_native(U,p)==V and literal_native(V,q)==W,'independent two-native replay')
    restricted=tuple(row&((1<<len(T))-1) for row in W[:len(T)])
    need(restricted==flip(T,a,b),'exact one-arc endpoint')
    return restricted,dict(pair=[a,b],helper_out_mask=U[-1],paths=[p,q],inverse=inverse,case=case)

def chart(T,root,feedback):
    other=next(v for v in range(len(T)) if v!=root)
    f=pair_fields(T,feedback,root,other);e,S2=f[2],f[6]
    g=[f[0],(S2-e)//2,(S2+e)//2]
    need((S2-e)%2==0,'integral rooted golden')
    return dict(H=f[0],golden=g,charts=oldphase.phase(g))

def main():
    # Exact observer intertwining: after converse, bit1 is command8 and bit0 command11.
    common=json.loads((RESULTS/'adaptive_state_induction_20261103/two_port_envelope.json').read_text())
    F=json.loads((OLD/'interface/natural_pair.json').read_text())['factors_in_profile_order_h_s_q_F']['mark_root']
    mats={1:[[1,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,2]],
          0:[[1,0,0,0],[-1,0,0,0],[0,1,1,-1],[1,0,0,2]]}
    for bit,cmd in ((1,8),(0,11)):
        A=[[Fraction(x) for x in row] for row in common['rank']['query_matrices'][cmd]]
        left=[[sum(mats[bit][i][l]*F[l][j] for l in range(4)) for j in range(17)] for i in range(4)]
        right=[[sum(F[i][l]*A[l][j] for l in range(17)) for j in range(17)] for i in range(4)]
        need(left==right,'exact4-to17 intertwining under converse command '+str(cmd))
    levels=[];profiles=0;births=0;native_checks=[]
    for d in range(4):
        records=[]
        for i in range(4**d):
            word=gray_word(d,i);T,m=atlas(word);z=state((2,4,1),0,2,word)
            H=count(T,(0,m['root']));need(H==output(z),'literal full H versus four-register sealed output')
            need(reader.strong(T),'all literal atlas full graphs strong')
            B=reader.graph((2,4,1),0,word+'0');v=reader.run(reader.all_fields((2,4,1))[0],word+'0')
            need(H==v[0]+4*v[3] and reader.strong(B),'seed-reader and strong-component control')
            vs=[v for v in range(len(T)) if v!=m['seal']];P=induced(converse(T),vs);w=vs.index(0);r=vs.index(m['root'])
            f=pair_fields(P,(w,r),w,r)
            need([f[0],f[3],f[4],f[5]]==list(z),'four state is actual converse-parent rooted profile')
            profiles+=1
            if d:
                parent,parentmeta=atlas(gray_word(d-1,i//4));need(induced(T,list(range(len(parent))))==parent,'actual induced parent at address prefix')
                need(tuple(v//2 for v in point(d,i))==point(d-1,i//4),'Hilbert cell ancestry')
                # Replay the two interior births as actual fixed-port cuts.
                U=parent
                for j in range(len(parentmeta['free']),len(m['free'])):
                    bit=int(word[j]);out=(1<<m['seal'])|(0 if bit else 1<<m['mark'])
                    U=reader.birth(U,out)
                    rev=converse(U);previous=converse(tuple(row&((1<<(len(U)-1))-1) for row in U[:-1]))
                    incoming=(1<<m['seal'])|(0 if bit else 1<<m['mark'])
                    expected=tuple(row|((1<<(len(U)-1)) if incoming>>v&1 else 0) for v,row in enumerate(previous))+(sum(1<<v for v in range(len(previous)) if not(incoming>>v&1)),)
                    need(rev==expected,'literal converse sparse-in command8/11')
                    births+=1
                need(U==T,'birth replay equals atlas child')
            if i:
                old,oldmeta=atlas(gray_word(d,i-1));diff=[(a,b) for a,b in combinations(range(len(T)),2) if ((old[a]^T[a])>>b)&1]
                need(len(diff)==1 and 0 in diff[0],'Gray neighbors differ by one marked-tail arc')
                need(sum(abs(a-b) for a,b in zip(point(d,i),point(d,i-1)))==1,'unit Hilbert scan join')
                if d==2:
                    fixed,proof=gadget(old,*diff[0]);need(fixed==T,'native Hilbert scan endpoint')
                    native_checks.append(dict(from_index=i-1,to_index=i,**proof))
            records.append([i,word,*point(d,i),H])
        levels.append(dict(depth=d,cells=len(records),value_min=min(r[-1] for r in records),value_max=max(r[-1] for r in records),records=records if d<=1 else None,full_record_sha256=hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest()))
    # Every actual small seed/mark, including weak and singleton seeds.
    strong_controls=0;weak_component_example=None
    for n in (1,2,3):
        for code in range(1<<(n*(n-1)//2)):
            seed=[0]*n
            for j,(a,b) in enumerate(combinations(range(n),2)):
                u,v=(a,b) if code>>j&1 else (b,a);seed[u]|=1<<v
            for w in range(n):
                for length in range(3):
                    for bits in product('01',repeat=length):
                        word=''.join(bits)
                        for k in (1,2):
                            T,m=atlas(word,tuple(seed),w,k);need(reader.strong(T),'universal collar full strongness')
                            component=reader.graph(tuple(seed),w,word+'0')
                            if not reader.strong(component) and weak_component_example is None:weak_component_example=dict(seed=seed,mark=w,word=word,prefix=k,full_rows=T,component_rows=component)
                            strong_controls+=1
    # Local edge response: suffix includes free letters only, then implicit zero-seal.
    gradients=0
    seed=(2,4,1)
    for prelen in range(4):
        for pre in map(''.join,product('01',repeat=prelen)):
            a,b,c,d=state(seed,0,2,pre)
            for length in range(4):
                for suffix in map(''.join,product('01',repeat=length)):
                    delta=output(state(seed,0,2,pre+'0'+suffix))-output(state(seed,0,2,pre+'1'+suffix))
                    predicted=d-b-2*c if not suffix else d-b-c*((1<<len(suffix))+sum(1<<j for j,x in enumerate(suffix) if x=='0')+(suffix[0]=='0'))
                    need(delta==predicted,'exact local Gray edge-gradient formula');gradients+=1
    jumps=[]
    for d in range(1,5):
        m=2*d;x='1'*m;y='1'*(m-1)+'0'
        delta=output(state(seed,0,2,y))-output(state(seed,0,2,x))
        need(delta==19*(1<<(m-1))-2,'unbounded H-jump on a unit Hilbert scan edge')
        inv=lambda g: next(i for i in range(1<<m) if i^(i>>1)==g)
        ix,iy=inv(int(x,2)),inv(int(y,2));need(abs(ix-iy)==1,'hostile is consecutive in Gray scan')
        jumps.append(dict(depth=d,indices=[ix,iy],words=[x,y],delta_H=delta))
    # One real seed promotion, then mark switch to the retained old spare1.
    source,m=atlas('0');U=reader.graph(seed,0,'0');U_names=[0,1,2,7]
    promoted_seed=reader.birth(U,3);promoted_names=U_names+[8]
    promoted=reader.birth(source,(1<<3)|1|2)
    want=make(promoted_seed,promoted_names,0,'',[],3,[4,5],6,list(range(9)))
    need(promoted==want,'one-vertex promotion preserves all old names/arcs and collar')
    target=make(promoted_seed,promoted_names,1,'',[],3,[4,5],6,list(range(9)))
    changes=[(3,0),(3,1),(6,0),(6,1)];current=promoted;reset=[]
    for pair in changes:current,proof=gadget(current,*pair);reset.append(proof)
    need(current==target and reader.strong(target),'four-edge frame retarget exact endpoint')
    need(induced(target,[0,1,2,7,8])==promoted_seed,'all absorbed unsealed-seed arcs preserved')
    # Hostile to reinterpreting the old word at the new mark.
    reset_without_birth=make(U,U_names,1,'',[],3,[4,5],6,list(range(8)))
    naive=make(seed,[0,1,2],1,'0',[7],3,[4,5],6,list(range(8)))
    lost=[list(p) for p in combinations(range(8),2) if ((reset_without_birth[p[0]]^naive[p[0]])>>p[1])&1]
    need(lost==[[0,7],[1,7]],'old reader word must be absorbed, not reinterpreted')
    promotion=dict(source_rows=source,promoted_rows=promoted,target_rows=target,old_spare_new_mark=1,
        permanent_birth_name=8,permanent_birth_out_mask=11,absorbed_seed_names=promoted_names,absorbed_seed_rows=promoted_seed,
        retarget_pairs=changes,native_program=reset,source_chart=chart(source,6,(0,6)),target_chart=chart(target,6,(0,1,6)),
        naive_word_reinterpretation_differences=lost)
    # Two precise historical hostiles; retain proof-corrected statements in NOTE.
    triangle=(2,4,1);flipped=flip(triangle,0,1)
    need(sorted(x.bit_count() for x in triangle)!=sorted(x.bit_count() for x in flipped),'equal scores do not preserve sorted score under an arc flip')
    transitive=tuple(sum(1<<j for j in range(i+1,4)) for i in range(4));rev=flip(transitive,0,3)
    h0=reader.tables(transitive)[0][-1];h1=reader.tables(rev)[0][-1]
    need((h0,h1)==(1,5),'historical contraction-bound hostile')
    # Four Hilbert linear frames, with the Gray carry, really give eight controller states.
    # Frame encoded (sign,swap); q0 swaps, q3 negative-swaps, q1/q2 identity.
    frames={(1,0,0)};todo=list(frames)
    while todo:
        sign,swap,carry=todo.pop()
        for q in range(4):
            x,y=q//2,q%2;out=(x^carry,y^x)
            new=(sign*(-1 if q==3 else 1),swap^(q in (0,3)),y)
            if new not in frames:frames.add(new);todo.append(new)
    need(len(frames)==8,'finite orientation/carry controller has eight reachable states')
    paths=[HERE/'check.py',HERE/'NOTE.md',OLD/'inverse/reader.py',OLD/'inverse/native.py',OLD/'interface/gap.py',OLD/'interface/check.py',OLD/'interface/natural_pair.json',
       RESULTS/'adaptive_state_induction_20261103/two_port_envelope.json',RESULTS/'return701_algebra_20261103/SPACE_FILLING.md',
       RESULTS/'return701_algebra_20261103/hilbert_collatz.py',RESULTS/'gray_selector_recursion_20261018e.md',
       RESULTS/'adaptive_metagraph_reformulation_20261103/METAGRAPH.md',RESULTS/'adaptive_metagraph_reformulation_20261103/metagraph_probe.py',
       RESULTS.parent.parent/'01-canon/theorems/THM-452-doubling-tiling-coset-law.md',RESULTS.parent.parent/'01-canon/theorems/THM-830-b3-deletion-deck-mirror-current-calculus.md',
       RESULTS.parent.parent/'04-computation/gray_codes_space_filling_s136.py',RESULTS.parent.parent/'04-computation/gray_code_s160.py']
    out=dict(status='PROVED identities with FINITE-EXACT controls; no coverage or recursive bank closure',levels=levels,
        exact_observer_intertwining_commands={'bit1':8,'bit0':11},four_state_matrices=mats,controller_states=sorted(frames),
        literal_profile_controls=profiles,interior_births=births,all_seed_mark_full_strong_controls=strong_controls,
        weak_component_strong_wrapper_hostile=weak_component_example,edge_gradient_controls=gradients,
        gray_native_programs=native_checks,unbounded_jump_controls=jumps,promotion=promotion,
        first_level_missing_returns=dict(parent_H=33,targets=[101,103],actual_four_children=[x[-1] for x in levels[1]['records']]),
        historical_hostiles=dict(equal_score_triangle_before=triangle,after=flipped,transitive4_rows=transitive,reversed_long_arc_rows=rev,H_values=[h0,h1]),
        source_sha256={str(p.relative_to(RESULTS.parent.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
