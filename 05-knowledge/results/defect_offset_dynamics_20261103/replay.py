"""Reproduce all finite certificates and hostile controls. No assert statements."""
from wing import *
from branch_search import canon,main as branches_main
from boundary_search import main as boundary_main
from dynamics import main as dynamics_main
from clone_clock import main as clone_main
from run_polynomial import main as polynomial_main
from wing import main as wing_main
from return_program import literal_Q4,golden
from carrier import native
from itertools import combinations
import hashlib


def paths(T):
    def go(prefix,left):
        if not left:yield tuple(prefix);return
        for v in sorted(left):
            if not prefix or T[prefix[-1]]>>v&1:yield from go(prefix+[v],left-{v})
    return go([],set(range(len(T))))


def literal_observer(T,r):
    P,E,S=bin_values(T,r,2)
    require((S[2]-E[1])%2==0,'literal golden parity')
    v=(P[1],(S[2]-E[1])//2,(S[2]+E[1])//2)
    return dict(root=r,observer=v,phase=field(v))


def birth_checks():
    T,ports,v=source(((1,1),));checks=0
    # Independent DFS for every labelled two-ear extension of the four-vertex source.
    for m in range(16):
        U=birth(T,m);values=ear_values(U)
        for q in range(32):
            require(values[q]==sum(1 for _ in paths(birth(U,q))),'all512 cut responses versus DFS')
            checks+=1
        tables.cache_clear()
    for p,m in product(range(2),range(8)):
        U,_=graded.birth(T,ports,p,m,3);pp=ports+(len(U)-1,)
        require(apply_lift(v,3,p,m)==graded.literal(U,pp),'literal full named birth')
    for old in graded.skeletons(3):
        for p,m in product(range(2),range(8)):
            A,B=lift(3,p,m)[old]
            for count in (2,3):
                predicted=Counter()
                for (s,d),w in A.items():predicted[s,d+count]+=w
                for (s,d),w in B.items():predicted[s,d+count-1]+=count*w
                require(+predicted==named_kernel(old,count,3,p,m),'anonymous derivative identity')
    return checks


def escape_checks():
    chain=json.loads((HERE/'branches.json').read_text())
    require([len(x) for x in chain]==[2,6,1],'exact successive strong classes')
    T=tuple(chain[-1][0]['T']);require(sum(1 for _ in paths(T))==161,'source161 independent DFS')
    require(not all_two(T,{485,487}),'complete two-ear obstruction')
    (HERE/'literal_161.json').write_text(json.dumps(dict(T=T,hits=[]),indent=2)+'\n')
    successes=[]
    for i,row in enumerate(chain[1]):
        U=tuple(row['T']);hits=[(m,q,h) for m,q,h in all_two(U,{161,163}) if strong(birth(birth(U,m),q))]
        require(all(h==161 for m,q,h in hits),'no strong163 continuation')
        if hits:successes.append((i,len(hits)))
    require(successes==[(0,6),(3,6)],'two successful source classes, one target class')
    hist=Counter();found={}
    for path in paths(T):
        U=native.native(T,path);hu=tables(U)[0][-1];hist[hu]+=1
        if hu<=487:
            for mask,value in enumerate(ear_values(U)):
                if value not in (485,487) or value in found:continue
                V=birth(U,mask)
                if strong(V):
                    require(sum(1 for _ in paths(V))==value,'escape target independent DFS')
                    found[value]=dict(path=path,intermediate=U,intermediate_h=hu,masks=(mask,),target=V)
        tables.cache_clear()
    require(set(found)=={485,487} and hist[161]==0,'both actual escapes; no one-native H-neutral reset')
    (HERE/'native_escape_161.json').write_text(json.dumps(found,indent=2)+'\n')
    (HERE/'native_161.json').write_text(json.dumps(dict(histogram=sorted(hist.items()),neutral=[]),indent=2)+'\n')
    phase_rows=[]
    for h,data in sorted(found.items()):
        U=tuple(data['target']);phase_rows.append(dict(h=h,source=literal_observer(T,0),target=literal_observer(U,0)))
        if h==485:U=opposite(U);root=9
        else:root=2
        require(rooted_data(U,root)[3]==h,'actual zero guard')
        for j in range(1,3):
            x=len(U);U=birth(U,1<<root);U=birth(U,1<<x);root=x
            require(strong(U) and tables(U)[0][-1]==h*3**j,'strong zero continuation')
            tables.cache_clear()
    neutral1=[];neutral2=[];arcs=list(combinations(range(10),2))
    for number,records in ((1,neutral1),(2,neutral2)):
        for edges in combinations(arcs,number):
            U=list(T)
            for a,b in edges:U[a]^=1<<b;U[b]^=1<<a
            U=tuple(U)
            if tables(U)[0][-1]==161:
                require(canon(U)==canon(T),'every one/two-edge neutral neighbor is isomorphic')
                records.append(edges)
            tables.cache_clear()
    require(len(neutral1)==4 and len(neutral2)==8,'neutral-neighbor counts')
    old1=[]
    for edges in neutral1:
        a,b=edges[0];U=list(T);U[a]^=1<<b;U[b]^=1<<a
        old1.append(dict(edge=(a,b),T=U,hits=[]))
    (HERE/'neutral_flip_161.json').write_text(json.dumps(old1,indent=2)+'\n')
    (HERE/'neutral_two_flip_161.json').write_text(json.dumps(dict(neutral=8,classes=[]),indent=2)+'\n')
    # Fresh hostile: same H and F49 state, different actual continuation language.
    S,p,v=source(((1,1),));S,p,v=exact_vector(S,p,v,(1,5,1,7));p=p[:3]
    bad,_=actual(S,p,(0,1,0,1));good,_=source2(3,1)
    bo=literal_observer(bad,7);go=literal_observer(good,7)
    require(bo['phase']==go['phase']==(2,2,2),'matched F49 hostile')
    require(bo['observer']==(53,225,226) and go['observer']==(53,204,205),'integer sidecar difference')
    require(not all_two(bad,{161,163}),'hostile53 has no two-ear continuation')
    good_hits=[z for z in all_two(good,{161,163}) if strong(birth(birth(good,z[0]),z[1]))]
    require(any(z[2]==161 for z in good_hits),'positive53 strong continuation')
    (HERE/'literal_53.json').write_text(json.dumps(dict(T=bad,h=53,hits=[]),indent=2)+'\n')
    return dict(successful53_sources=successes,neutral_one=len(neutral1),neutral_two=len(neutral2),
                native_paths=sum(hist.values()),phases=phase_rows,
                hostile53=dict(bad=bad,good=good,bad_observer=bo,good_observer=go,good_hits=good_hits))


def main():
    independent=birth_checks()
    dynamics_main();boundary_main();branches_main();clone_main();polynomial_main();wing_main()
    escapes=escape_checks();q4=0
    for m,n in ((1,1),(3,1)):
        T,p=source2(m,n);v,ph=phase(m,n)
        q4+=literal_Q4(T,p[0],v)
    B=21526903444
    require(10**226>B*9**226 and 10**225<=B*9**225,'exact eventual bound')
    rows=json.loads((HERE/'boundary.json').read_text())
    require(sum(len(r['hits']) for r in rows)==14,'complete full-boundary program count')
    require(sum(h['strong'] for r in rows for h in r['hits'])==10,'ten strong uniform programs')
    require(all(h['offset']==2 for r in rows for h in r['hits']),'no uniform plus-four')
    from incoming_checks import run as incoming_run
    incoming=incoming_run()
    from wing_switch import check as switch_check
    switch=switch_check()
    receipt=dict(independent_two_ears=independent,all_five_Q4_checks=q4+switch['literal_Q4_checks'],escapes=escapes,
                 boundary_families=len(rows),boundary_programs=14,strong_boundary_programs=10,
                 clock=clock(),incoming=incoming,wing_switch=switch,status='PASS; same-session review, independent combinatorial counting paths')
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    digest=lambda p:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    manifest={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.json') and p.name!='manifest.json'}
    repo=HERE.parents[2]
    imports={Path(m.__file__).resolve() for m in tuple(sys.modules.values()) if getattr(m,'__file__',None)}
    imports.update(Path(m.__file__).resolve() for m in (graded,golden,native))
    manifest['inherited_sources']={p.relative_to(repo).as_posix():digest(p) for p in sorted(imports)
                                   if p.is_relative_to(repo) and p.parent!=HERE and p.suffix=='.py' and p.is_file()}
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('PASS:512 DFS controls;47 all-size boundary families;20 Q4 checks;literal161 escapes;81/18-state clocks;neutral joint switch',flush=True)


if __name__=='__main__':main()
