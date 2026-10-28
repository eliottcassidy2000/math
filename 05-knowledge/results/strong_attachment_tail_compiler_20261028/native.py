"""Actual tournament ports, transporter frames, and bounded unnamed closure."""
from collections import Counter,defaultdict,deque
from itertools import combinations,permutations,product
from pathlib import Path
import argparse,importlib.util,json

ROOT=Path(__file__).resolve().parents[3]
def source(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
core=source('native_core','05-knowledge/results/strong_atom_continuation_vitali_20261026/core.py')
oracle=source('native_oracle','05-knowledge/results/relation_tower_frobenius_ports_20261028/check_root_switch.py')
def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def pair_model(n,p,converse=False):
    free=[(i,j) for i in range(n) for j in range(i+2,n)];index={e:k for k,e in enumerate(free)}
    base=0;channels=[]
    for k,(i,j) in enumerate(free):
        a,b=p[i],p[j];lo,hi=sorted((a,b));flip=(a>b)^converse
        if hi==lo+1:
            if not flip:base|=1<<k
        else:
            if flip:base|=1<<k
            channels.append((index[lo,hi],k))
    tables=[]
    for offset in range(0,len(free),8):
        tables.append(tuple(sum(((x>>(old-offset))&1)<<new for old,new in channels if offset<=old<offset+8) for x in range(256)))
    def transform(x):
        value=base
        for offset,table in enumerate(tables):value^=table[(x>>(8*offset))&255]
        return value
    return transform

def native_classes(n):
    free=[(i,j) for i in range(n) for j in range(i+2,n)];index={e:k for k,e in enumerate(free)}
    N=1<<len(free);parent=list(range(N));sizes=[1]*N
    def find(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    for a in range(n):
        for b in range(a+2,n):
            p=tuple(range(a))+tuple(range(a+1,b+1))+(a,)+tuple(range(b+1,n))
            move=pair_model(n,p);zeros=1<<index[a,b];ones=0
            if a:ones|=1<<index[a-1,a+1]
            if b+1<n:ones|=1<<index[a,b+1]
            for x in range(N):
                if x&zeros or x&ones!=ones:continue
                u,v=find(x),find(move(x))
                if u!=v:
                    if sizes[u]<sizes[v]:u,v=v,u
                    parent[v]=u;sizes[u]+=sizes[v]
    classes=[find(x) for x in range(N)];ids={a:i for i,a in enumerate(sorted(set(classes)))}
    classes=[ids[a] for a in classes]
    need(len(ids)==[1,1,1,2,4,12,56,456][n],'native ordinary class count')
    if n<=5:
        canonical=[core.canonical(core.presentation(n,x)) for x in range(N)]
        equivalence={}
        for k,c in zip(classes,canonical):equivalence.setdefault(k,set()).add(c)
        need(all(len(v)==1 for v in equivalence.values()) and len({next(iter(v)) for v in equivalence.values()})==len(ids),'permutation class oracle')
    reflect=pair_model(n,tuple(reversed(range(n))),True)
    blue=[n>=3 and reflect(x)==x for x in range(N)]
    return classes,blue

def join_mask(n,x,left):
    t=core.presentation(n,x)
    if left:t=( (1<<(n+1))-2,)+tuple(row<<1 for row in t)
    else:t=tuple(row|(1<<n) for row in t)+(0,)
    return core.free_mask(t)

def closure(maxn=7):
    banks={n:native_classes(n) for n in range(1,maxn+1)}
    summary=[];last=None;partitions={}
    def separate(n,x,y):
        q,blue=banks[n];N=len(q)
        if (q[x],blue[x])!=(q[y],blue[y]):
            return dict(word=[],terminal_n=n,observations=[[q[x],blue[x]],[q[y],blue[y]]])
        for flipped in (False,True):
            a,b=(x^(N-1),y^(N-1)) if flipped else (x,y)
            prefix=['complement'] if flipped else []
            if (q[a],blue[a])!=(q[b],blue[b]):
                return dict(word=prefix,terminal_n=n,observations=[[q[a],blue[a]],[q[b],blue[b]]])
            if n<maxn:
                for left in (True,False):
                    u,v=join_mask(n,a,left),join_mask(n,b,left)
                    if partitions[n+1][u]!=partitions[n+1][v]:
                        answer=separate(n+1,u,v)
                        answer['word']=prefix+['singleton-left' if left else 'singleton-right']+answer['word']
                        return answer
        raise RuntimeError('unequal cells have no separating context')
    for n in range(maxn,0,-1):
        q,blue=banks[n];N=len(q)
        signatures=[]
        for x in range(N):
            signature=(q[x],blue[x])
            if last is not None:signature+=(last[join_mask(n,x,True)],last[join_mask(n,x,False)])
            signatures.append(signature)
        keys=[(signatures[x],signatures[x^(N-1)]) for x in range(N)]
        ids={key:i for i,key in enumerate(sorted(set(keys)))}
        cells=[ids[k] for k in keys]
        partitions[n]=cells
        reset_only=len({(q[x],q[x^(N-1)],blue[x]) for x in range(N)})
        strong=[x for x in range(N) if core.is_strong(core.presentation(n,x))]
        witness=None;groups=defaultdict(list)
        for x in range(N):groups[q[x],q[x^(N-1)],blue[x]].append(x)
        for group in groups.values():
            if len({cells[x] for x in group})>1:
                a=group[0];b=next(x for x in group if cells[x]!=cells[a]);witness=[a,b];break
        summary.append(dict(n=n,presentations=N,ordinary_classes=len(set(q)),reset_complement_colour_cells=reset_only,
                            bounded_attachment_cells=len(ids),strong_presentations=len(strong),strong_cells=len({cells[x] for x in strong}),
                            split_witness=witness,separating_context=separate(n,*witness) if witness else None))
        if n==6:
            remaining=defaultdict(list)
            for x,c in enumerate(cells):remaining[c].append(x)
            residual=[]
            for members in remaining.values():
                for x,y in combinations(members,2):
                    found=None
                    # Extend only unresolved pairs to order eight; no order-eight bank assumed.
                    for f0,L0,f1,L1,f2 in product((False,True),repeat=5):
                        u,v=x,y;word=[];size=n
                        for flip,left in [(f0,L0),(f1,L1)]:
                            if flip:
                                mask=(1<<((size-1)*(size-2)//2))-1;u^=mask;v^=mask;word.append('complement')
                            u,v=join_mask(size,u,left),join_mask(size,v,left);size+=1
                            word.append('singleton-left' if left else 'singleton-right')
                        if f2:
                            mask=(1<<((size-1)*(size-2)//2))-1;u^=mask;v^=mask;word.append('complement')
                        if not oracle.isomorphisms(core.presentation(size,u),core.presentation(size,v),False):
                            found=dict(word=word,terminal_n=size,terminal_masks=[u,v]);break
                    residual.append(dict(pair=[x,y],separating_context=found))
            summary[-1]['order_eight_residual_probes']=residual
        last=cells
    return banks,list(reversed(summary))

def tuple_key(group,vertices):return min(tuple(g[v] for v in vertices) for g in group)
def compose(g,h):return tuple(g[h[i]] for i in range(len(g)))
def inverse(g):return tuple(g.index(i) for i in range(len(g)))

def frame_test(t):
    n=len(t);G=oracle.isomorphisms(t,t,True);identity=tuple(range(n))
    frame=None
    for size in range(n+1):
        frame=next((b for b in combinations(range(n),size) if sum(all(g[v]==v for v in b) for g in G)==1),None)
        if frame is not None:break
    need(frame is not None,'faithful permutation action')
    orbit={tuple(g[v] for v in frame):g for g in G};need(len(orbit)==len(G),'free frame orbit')
    for g,h in product(G,repeat=2):
        relative=compose(inverse(g),h)
        need(compose(g,relative)==h and relative in G,'LCA relative transporter')
    counts={}
    for r in (2,3):
        full=set();star=set();pairs=set();fibres=defaultdict(set)
        for ports in permutations(range(n),r):
            key=tuple_key(G,ports);full.add(key)
            singles=tuple(tuple_key(G,(v,)) for v in ports)
            ss=singles+tuple(tuple_key(G,(ports[0],ports[j])) for j in range(1,r))
            pp=singles+tuple(tuple_key(G,(ports[i],ports[j])) for i,j in combinations(range(r),2))
            star.add(ss);pairs.add(pp);fibres[pp].add(key)
        counts[str(r)]=dict(joint_orbits=len(full),star_records=len(star),all_pair_records=len(pairs),
                           largest_pair_fibre=max(map(len,fibres.values())))
    return dict(n=n,Aut=len(G),minimum_free_frame_size=len(frame),first_frame=frame,tuple_counts=counts,
                two_relative_frame_coordinates=len(G)**2)

def marker_attachment(t,ports,order):
    n=len(t);r=len(ports);rows=list(t)+[0]*r
    for j,v in enumerate(ports):
        rows[n+j]|=1<<v
        for u in range(n):
            if u!=v:rows[u]|=1<<(n+j)
    for i,j in combinations(range(r),2):
        a,b=(i,j) if order.index(i)<order.index(j) else (j,i)
        rows[n+a]|=1<<(n+b)
    return tuple(rows)

def atom_experiments(banks):
    c3=(2,4,1);t5=tuple(sum(1<<((i+d)%5) for d in (1,2)) for i in range(5))
    p7=tuple(sum(1<<((i+d)%7) for d in (1,2,4)) for i in range(7))
    s9=core.substitute(c3,[c3]*3)[0]
    A=tuple(row|8 for row in c3)+(0,);op=core.converse(A)
    t12=core.substitute(c3,[A]*3)[0];u12=core.substitute(c3,[op,A,A])[0]
    G12=set(oracle.isomorphisms(t12,t12,True));H12=set(oracle.isomorphisms(u12,u12,True))
    cyclic=lambda t: [c for c in combinations(range(12),3) if all(sum(t[u]>>v&1 for v in c if u!=v)==1 for u in c)]
    need(cyclic(t12)==cyclic(u12) and len(cyclic(t12))==67,'81/27 cyclic triples')
    need(H12<G12 and all(all(g[v]!=v for v in range(12)) for g in G12-H12),'all nonempty stabilizers coincide')
    transporter=[sum(g[3]==7 for g in group) for group in (G12,H12)]
    need(transporter==[27,0],'81/27 cross-block transporter')
    h12=[oracle.hamilton(t) for t in (t12,u12)]
    need(h12==[462915]*2,'81/27 Hamilton mass')
    frames={name:frame_test(t) for name,t in [('C3',c3),('T5',t5),('P7',p7),('C3[C3,C3,C3]',s9),('T81',t12),('T27',u12)]}
    G=oracle.isomorphisms(p7,p7,True)
    need(tuple(2*i%7 for i in range(7)) in G,'literal F8 Frobenius automorphism')
    a,b=(0,1,3),(0,1,5)
    need(all(tuple_key(G,(a[i],a[j]))==tuple_key(G,(b[i],b[j])) for i,j in combinations(range(3),2)),'P7 pair hostile')
    need(tuple_key(G,a)!=tuple_key(G,b),'P7 joint hostile')
    markers=[]
    for order in permutations(range(3)):
        h=[oracle.hamilton(marker_attachment(p7,ports,order)) for ports in (a,b)]
        markers.append(dict(marker_order=order,H=h))
    common=[sum(all(p7[u]>>v&1 for u in ports) for v in range(7)) for ports in (a,b)]
    need(common==[0,1],'native ternary guard')
    complements={}
    for name,t in [('C3',c3),('T5',t5),('P7',p7)]:
        n=len(t);paths=list(core.paths(t));presentations={core.free_mask(core.relabel(t,p)) for p in paths}
        q,blue=banks[n];N=len(q);source_class=q[next(iter(presentations))]
        neighbours={q[x^(N-1)] for x in presentations}
        edges=defaultdict(set)
        for x in range(N):edges[q[x]].add(q[x^(N-1)])
        reached={source_class};todo=deque(reached)
        while todo:
            for v in edges[todo.popleft()]:
                if v not in reached:reached.add(v);todo.append(v)
        reached_states=[x for x in range(N) if q[x] in reached]
        partners=[]
        for cls in sorted(neighbours):
            x=next(x for x in range(N) if q[x]==cls);child=core.presentation(n,x)
            partners.append(dict(class_id=cls,Aut=len(oracle.isomorphisms(child,child,True)),H=oracle.hamilton(child),strong=core.is_strong(child)))
        complements[name]=dict(H=len(paths),path_orbits=len(presentations),first_partner_classes=len(neighbours),
            partners=partners,closure_ordinary_classes=len(reached),closure_presentations=len(reached_states),
            minimum_reset_complement_colour_cells=len({(q[x],q[x^(N-1)],blue[x]) for x in reached_states}))
    fields=source('native_fields','05-knowledge/results/frobenius_lca_metagraph_20261028/fields.py')
    field_controls=[]
    for r,poly in [(6,67),(8,283)]:
        f=fields.Field(r,poly)
        normal=next(f.orbit(a) for a in range(1,1<<r) if len(f.orbit(a))==r and fields.rank(f.orbit(a))==r)
        admissible=[];conflicts=[]
        for k in range(r):
            # An unordered-pair orbit is orientable iff it never returns reversed.
            conflict=None
            for u,v in combinations(range(r),2):
                x,y=(u+k)%r,(v+k)%r
                while (x,y)!=(u,v):
                    if (x,y)==(v,u):conflict=[u,v];break
                    x,y=(x+k)%r,(y+k)%r
                if conflict is not None:break
            if conflict is None:admissible.append(k)
            else:conflicts.append(dict(phase=k,reversed_pair=conflict))
        need(admissible==([0,2,4] if r==6 else [0]),'normal-frame orientability')
        field_controls.append(dict(field=1<<r,normal_basis=normal,tournament_admissible_phases=admissible,conflicts=conflicts))
    # Realize F64's surviving C3 on six distinct vertices of an actual strong atom.
    ports=(0,3,1,4,2,5)
    shift=tuple((i//3)*3+(i+1)%3 if i<6 else i for i in range(9))
    need(shift in oracle.isomorphisms(s9,s9,True),'F64 odd phase on strong atom')
    need(all(shift[ports[j]]==ports[(j+2)%6] for j in range(6)),'normal-basis phase intertwiner')
    field_controls[0]['strong_atom_realization']=dict(atom='C3[C3,C3,C3]',ports=ports,phase_two_permutation=shift)
    return dict(frames=frames,P7_triple_hostile=dict(ports=[a,b],common_out=common,marker_attachments=markers),
                complement_closures=complements,field_frame_controls=field_controls,
                transporter_hostile=dict(Aut=[len(G12),len(H12)],H=h12,labelled_cyclic_triples=67,
                    transporter_3_to_7=transporter,identical_nonempty_pointwise_stabilizers=True))

def main():
    p=argparse.ArgumentParser();p.add_argument('--output');args=p.parse_args()
    banks,closed=closure(7)
    result=dict(bounded_unnamed=closed,actual_atoms=atom_experiments(banks),scope_max_order=7)
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:Path(args.output).write_text(text,encoding='utf-8',newline='\n')
    else:print(text,end='')

if __name__=='__main__':main()
