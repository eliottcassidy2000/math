#!/usr/bin/env python3
"""Independent finite audit; no producer code, labels, or result table imported."""
from array import array
from collections import Counter
from hashlib import sha256
from itertools import combinations, product, permutations
from pathlib import Path
import argparse, csv, gzip, json, random, subprocess, sys, tempfile
P=Path(__file__).resolve().parent
I=P/'inputs'
CHECKS=Counter()
def need(value,why):
    CHECKS[why]+=1
    if not value: raise RuntimeError(why)
def read(rel): return json.loads((I/rel).read_text())
def add(a,b): return ((((a>>4)+(b>>4))%128)<<4)|((a^b)&15)
def sub(a,b): return ((((a>>4)-(b>>4))%128)<<4)|((a^b)&15)
def neg(a): return sub(0,a)
def group_sum(xs):
    q=0
    for x in xs:q=add(q,x)
    return q

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    provenance=json.loads((P/'input_provenance.json').read_text())
    for z in provenance:
        b=(P/z['path']).read_bytes();need(len(b)==z['bytes'] and sha256(b).hexdigest()==z['sha256'],'pinned input digest')
    old=read('native_packet_decomposition_20261023/native/domains.json')['sources']
    electrical=read('metagraph_order8_electrical_20261023/receipt.json')['order8']
    X=electrical['higher_layers'][-1]['lifts'][0];Y=electrical['radical_lifts_mod4'];potentials=[X]+Y[:5];mods=[256]+[4]*5
    label=lambda s:16*(X[s]%128)+sum((Y[j][s]%2)<<(j-1) for j in range(1,5))
    need(all((X[v]-Y[0][v]-Y[4][v])%2==0 for v in range(len(X))),'original redundant parity relation')
    U,V=1965,2234
    need(all(2*(f[U]-f[V])%m==0 for f,m in zip(potentials,mods)),'source-row dyadic U/V equality')
    supply=read('compatible_charged_helpers_20261023/supply/selected_supply.json')
    new=read('guarded_component_completion_20261024/sections/native_2125_receipt.json')
    data={}
    for z in old:
        by={r['target']:r for r in z['exact_all_shared_clocks_target_domain']}
        data[z['source']]={'records':[by[U],by[V]],'fixed':z['fixed_clocks'],'blocks':[]}
    core=sorted(data)
    for z in supply['helpers']:
        ss=z['selected'];data[z['source']]={'records':[ss['before'],ss['after']],'fixed':ss['common_clocks'],'blocks':[]}
    by={r['target']:r for r in new['records']};data[2125]={'records':[by[U],by[V]],'fixed':[],'blocks':new['retained_blocks']}
    costs={s:sum(a!=b for a,b in zip(*(r['path'] for r in z['records']))) for s,z in data.items()}
    cheap=sorted(s for s in core if costs[s]==2);outside=sorted(set(core)-set(cheap))
    free_helpers=sorted(set(data)-set(core)-{368,776,3139})
    need(len(data)==82 and len(cheap)==20 and len(outside)==8 and len(free_helpers)==51,'exact original native populations')
    need([costs[s] for s in outside]==[4,3,4,3,4,3,3,3],'actual outside costs')
    labels={s:label(s) for s in data}
    charge=lambda ss:group_sum(labels[s] for s in ss)
    # Reconstruct actual local paths, normalized resources, and common source/partner host.
    nodes=list(csv.DictReader((I/'metagraph_order8_20261023/merged_nodes.tsv').open(),delimiter='\t'))
    merged={int(z[k]):int(z['merged_id']) for z in nodes for k in ('class_a','class_b')}
    owner=array('H');owner.frombytes(gzip.decompress((I/'metagraph_order8_20261023/owner_map.u16.gz').read_bytes()))
    if sys.byteorder!='little':owner.byteswap()
    pairs=list(combinations(range(8),2));chords=[(i,j) for i in range(8) for j in range(i+2,8)];full=(1<<21)-1
    def matrix(code):
        a=[[0]*8 for _ in range(8)]
        for k,(i,j) in enumerate(pairs):a[i][j]=(code>>(27-k))&1;a[j][i]=1-a[i][j]
        return a
    def word(a,p):return sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(chords))
    def reflect(m):return sum(((m>>k)&1)<<chords.index((7-j,7-i)) for k,(i,j) in enumerate(chords))
    def partner(a,p):
        pos={x:i for i,x in enumerate(p)}
        return [[0 if i==j else a[i][j] if abs(pos[i]-pos[j])==1 else a[j][i] for j in range(8)] for i in range(8)]
    # Verify original harmonic rows directly from the full pinned multigraph.
    harmonic=[[0]*len(nodes) for _ in potentials];edge_count=0
    with gzip.open(I/'metagraph_order8_20261023/merged_edges.tsv.gz','rt') as stream:
        for z in csv.DictReader(stream,delimiter='\t'):
            a,b=int(z['source']),int(z['target']);w=int(z['blue'])+int(z['black']);edge_count+=1
            for j,f in enumerate(potentials):
                delta=w*(f[a]-f[b]);harmonic[j][a]+=delta;harmonic[j][b]-=delta
    for j,m in enumerate(mods):need(all(v%m==0 for v in harmonic[j]),'every original native Laplacian row harmonic')
    # Independently enumerate the two target canonical classes, not just owner-map aliases.
    target_words={U:set(),V:set()}
    for t in (U,V):
        for canon in set(int(nodes[t][key]) for key in ('canon_a','canon_b')):
            a=matrix(canon)
            for p in permutations(range(8)):
                if all(a[x][y] for x,y in zip(p,p[1:])):target_words[t].add(word(a,p))
    need(not target_words[U]&target_words[V],'distinct literal target canonical classes')
    frames={368:list(range(8)),776:[8,9,4,10,11,12,13,3],3139:[14,15,16,17,18,3,4,19],2125:list(range(20,28))}
    fresh=28
    for s in sorted(set(data)-set(frames)):frames[s]=list(range(fresh,fresh+8));fresh+=8
    matrices={s:matrix(int(nodes[s]['canon_a'])) for s in data};resources={};reserved=set();frozen_clock={}
    for s,z in data.items():
        paths=[r['path'] for r in z['records']]
        actualmenu=set();aa=matrices[s]
        for pp in permutations(range(8)):
            if not all(aa[x][y] for x,y in zip(pp,pp[1:])):continue
            if not all(pp[i]==v for i,v in z['fixed']):continue
            if not all(any(list(pp[i:i+len(block)])==block for i in range(9-len(block))) for block in z['blocks']):continue
            mm=word(aa,pp)
            for tt in (U,V):
                if full^mm in target_words[tt]:actualmenu.add((pp,mm,tt))
        need(actualmenu=={(tuple(r['path']),r['mask'],r['target']) for r in z['records']},'entire guarded native menu exactly two full records')
        frozen_clock[s]=next(i for i in range(8) if paths[0][i]!=paths[1][i])
        resources[s]=[]
        for bit,r in enumerate(z['records']):
            p=r['path'];a=matrices[s];m=r['mask'];t=(U,V)[bit]
            need(sorted(p)==list(range(8)) and all(a[x][y] for x,y in zip(p,p[1:])),'actual Hamiltonian path')
            need(word(a,p)==m and merged[owner[m]]==s and merged[owner[full^m]]==t,'literal normalized source and target classes')
            need(word(partner(a,p),p)==full^m,'literal partner path')
            need(word(list(zip(*a)),p[::-1])==reflect(m),'literal converse path')
            need(all(p[i]==v for i,v in z['fixed']),'all original indexed pins')
            need(all(any(p[i:i+len(block)]==block for i in range(9-len(block))) for block in z['blocks']),'all moving relative blocks')
            orbit={m,reflect(m),full^m,full^reflect(m)}
            need(len(orbit)==4 and not reserved.intersection(orbit),'pairwise disjoint full resource orbits')
            reserved.update(orbit);resources[s].append(sorted(orbit))
            need((X[s]-X[t])%8==r['phase'],'literal ordinary source phase retained')
    need(len(reserved)==656 and fresh==652,'native 652-owner and 656-word supply')
    for s,t in combinations(data,2):
        need(len(set(frames[s])&set(frames[t]))==(2 if {s,t}<={368,776,3139} else 0),'only authenticated shared-owner intersections')
    def glue(amap):
        arcs={}
        for s,a in amap.items():
            em=frames[s]
            for i,j in pairs:
                x,y=em[i],em[j];key=min(x,y),max(x,y);v=a[i][j] if x<y else a[j][i]
                need(key not in arcs or arcs[key]==v,'coherent local host arcs')
                arcs[key]=v
    glue(matrices)
    for bits in product(range(2),repeat=3):
        gp={s:[frames[s][x] for x in data[s]['records'][bit]['path']] for s,bit in zip((368,776,3139),bits)}
        val=(gp[368][1],gp[776][4],gp[3139][7])
        need((val[0]==val[1]==val[2]==4)==(bits==(0,1,1)),'exact external owner4 branch')
        glue({s:partner(matrices[s],data[s]['records'][bit]['path']) for s,bit in zip((368,776,3139),bits)})
    helper_top2=sum(sorted((costs[s] for s in free_helpers),reverse=True)[:2]);helper_top4=sum(sorted((costs[s] for s in free_helpers),reverse=True)[:4])
    core_top14=sum(sorted((costs[s] for s in core),reverse=True)[:14])
    need((helper_top2,helper_top4,core_top14)==(15,29,39),'actual global and helper cost envelopes')
    protected_boundary={str(q):sum(sorted([2]*16+[costs[1923],costs[q]],reverse=True)[:14])+helper_top4 for q in (2087,2835)}
    need(set(protected_boundary.values())=={60},'ten-frozen protected-pair boundary cost60')
    single_boundary={c:sum(sorted([2]*16+[c],reverse=True)[:14])+helper_top4 for c in (3,4)}
    need(single_boundary=={3:58,4:59},'eleven-frozen chosen-good-source boundary costs58and59')
    # A second implementation in C++ enumerates all 317,521,920 base subsets.
    stdin=''.join(f'{s} {labels[s]}\n' for s in cheap+outside)
    with tempfile.TemporaryDirectory(prefix='four-pin-') as td:
        exe=Path(td)/'enumerate';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra',str(P/'enumerate.cpp'),'-o',str(exe)],check=True,capture_output=True)
        raw=subprocess.run([str(exe)],input=stdin,text=True,capture_output=True,check=True).stdout
    rows=[];pairfail=[]
    for line in raw.splitlines():
        a=line.split();v=list(map(int,a[1:]));f=v[:4]
        if a[0]=='R':
            mask=v[4];bad={str(outside[i]):{'grade':v[5+3*i],'charge':v[6+3*i],'holes_across_grades8_9':v[7+3*i]} for i in range(8) if not (mask>>i)&1}
            rows.append({'frozen_cheap':f,'passing':[outside[i] for i in range(8) if (mask>>i)&1],'failure_witnesses':bad})
        elif a[0]=='F':pairfail.append({'frozen_cheap':f,'both_active_pair':v[4:6],'grade':v[6],'charge':v[7],'holes_across_grades8_9_10':v[8]})
        else:raise RuntimeError('unexpected census record')
    need(len(rows)==4845 and len({tuple(r['frozen_cheap']) for r in rows})==4845,'all four-freeze sets exactly once')
    hist=Counter(len(r['passing']) for r in rows)
    need(hist=={8:4522,7:241,6:47,5:17,4:9,3:5,2:4},'independent exact passing histogram')
    banks=[list(p) for p in combinations(outside,2) if all(set(p)&set(r['passing']) for r in rows)]
    need(banks==[[1923,2087],[1923,2835]],'exactly two universal protected choice banks')
    singleton={str(s):next(r['frozen_cheap'] for r in rows if s not in r['passing']) for s in outside}
    need(len(pairfail)==6,'exactly six both-active two-slot failures')
    # Independent Python witness sets, not the C++ Gray traversal or producer DP.
    def supports(ids,grades):
        return {k:{charge(c) for c in combinations(ids,k)} for k in grades}
    chosen={tuple(r['frozen_cheap']) for r in rows if len(r['passing'])==2}|{tuple(z) for z in singleton.values()}|{tuple(z['frozen_cheap']) for z in pairfail}
    direct_controls=[]
    for frozen in sorted(chosen):
        row=next(r for r in rows if tuple(r['frozen_cheap'])==frozen)
        ids=sorted(set(cheap)-set(frozen));a=supports(ids,range(6,11));passed=[]
        for p in outside:
            holes={k:sorted(set(range(2048))-(a[k]|{add(labels[p],x) for x in a[k-1]})) for k in (8,9)}
            if not any(holes.values()):passed.append(p)
            else:
                saved=row['failure_witnesses'][str(p)];need(saved['charge'] in holes[saved['grade']],'direct combinations hostile witness')
        need(passed==row['passing'],'direct combinations agrees with Gray census')
        for fail in pairfail:
            if tuple(fail['frozen_cheap'])!=frozen:continue
            p,q=fail['both_active_pair'];k=fail['grade'];hole=fail['charge']
            coverage=a[k]|{add(labels[p],x) for x in a[k-1]}|{add(labels[q],x) for x in a[k-1]}|{add(add(labels[p],labels[q]),x) for x in a[k-2]}
            need(hole not in coverage,'direct two-active repair failure')
        direct_controls.append({'frozen_cheap':list(frozen),'passing':passed})
    # All free outside subsets and all feasible entry grades, with actual costs.
    envelope=Counter();budget_cases=0
    for p in outside:
        c=costs[p];rc=sorted([2]*16+[c],reverse=True)
        def wr(t):return sum(rc[:max(0,min(17,t))])
        rest=sorted(set(outside)-{p})
        for mask in range(128):
            E=[s for i,s in enumerate(rest) if mask>>i&1];e=len(E);ec=sorted((costs[s] for s in E),reverse=True)
            def we(t):return sum(ec[:max(0,min(e,t))])
            for k in range(8,10+e):
                r0=max(8,k-e)
                need(r0 in (8,9),'canonical receiver grade')
                for r in range(max(0,k-e),min(17,k)+1):
                    if r<8:rn=8;ed=8-r;rd=r+8
                    elif r>9:rn=9;ed=r-9;rd=34-r-9
                    else:rn=r;ed=0;rd=0
                    if ed:
                        need(ed<=e and rd+ed<=18 and wr(rd)+we(ed)<=45+c,'all actual entry costs')
                        envelope['entry_max']=max(envelope['entry_max'],wr(rd)+we(ed));budget_cases+=1
                    if rn>r0:
                        rd=min(2*rn-1,35-2*rn)
                        need(rd+1<=18 and wr(rd)+we(1)<=45+c,'all consecutive-grade costs')
                rd=2*min(r0,17-r0)
                need(rd+2<=18 and wr(rd)+15<=45+c,'all helper-edge costs')
                if e>=2:need(rd+2<=18 and wr(rd)+we(2)<=45+c,'all outside-exchange costs')
                envelope['helper_max_c'+str(c)]=max(envelope['helper_max_c'+str(c)],wr(rd)+15)
    need(envelope['entry_max']==43,'sharp entry envelope over actual outside choices')
    # Actual full endpoint routes for every free-core margin in several hostile patterns.
    helper_pairs=[frozenset(z) for z in combinations(free_helpers,2)];byhq={}
    for hp in helper_pairs:byhq.setdefault(charge(hp),hp)
    need(len(byhq)==972,'actual pinned free-helper pair coverage')
    rng=random.Random(4804);route_results=[];literal_routes=[]
    endpoint_patterns=sorted(chosen)[:5]
    # Include a case where 1923 fails, forcing the adaptive second choice.
    endpoint_patterns.append(tuple(singleton['1923']))
    endpoint_patterns=sorted(set(endpoint_patterns))
    def physical(selected):
        ans={s:[0]*6 for s in list(data)+[U,V]}
        for s in data:
            t=V if s in selected else U
            for j,f in enumerate(potentials):v=2*(f[s]-f[t]);ans[s][j]+=v;ans[t][j]-=v
        return ans
    for ix,frozen in enumerate(endpoint_patterns):
        row=next(r for r in rows if tuple(r['frozen_cheap'])==frozen)
        for bank in banks:
            p=next(s for s in bank if s in row['passing']);R=sorted((set(cheap)-set(frozen))|{p})
            # e=1 realizes four cheap + six outside frozen; e=7 tests maximum frames.
            for e in (0,1,7):
                E=[] if e==0 else sorted(set(bank)-{p}) if e==1 else sorted(set(outside)-{p})
                free=R+E;fixed=set(core)-set(free);fixedV={s for i,s in enumerate(sorted(fixed)) if (i+ix)%2}
                # Full constructive section at both grades via subset enumeration.
                section={8:{},9:{}}
                for k in (8,9):
                    for comb in combinations(R,k):section[k].setdefault(charge(comb),frozenset(comb))
                    need(len(section[k])==2048,'actual full receiver section witnesses')
                for k in range(len(free)+1):
                    a=frozenset(rng.sample(free,k));b=frozenset(rng.sample(free,k));delta=sub(charge(a),charge(b))
                    matched=next(((hp,byhq[add(hq,delta)]) for hq,hp in byhq.items() if add(hq,delta) in byhq),None)
                    need(matched is not None,'two feasible endpoint helper records')
                    ha,hb=matched;qa=add(charge(a),charge(ha));need(qa==add(charge(b),charge(hb)),'endpoint charge equality')
                    start=(a,ha);end=(b,hb)
                    if 8<=k<=9+e:
                        r0=max(8,k-e)
                        def refresh(J,H,r):return (section[r][sub(sub(qa,charge(J)),charge(H))]|frozenset(J),frozenset(H))
                        def norm(state):
                            cset,H=state;J=set(cset)&set(E);r=len(set(cset)&set(R));path=[state]
                            if r<8:
                                for s in sorted(J)[:8-r]:J.remove(s)
                                r=8;path.append(refresh(J,H,r))
                            elif r>9:
                                for s in sorted(set(E)-J)[:r-9]:J.add(s)
                                r=9;path.append(refresh(J,H,r))
                            while r>r0:
                                J.add(next(s for s in E if s not in J));r-=1;path.append(refresh(J,H,r))
                            return path
                        aa=norm(start);bb=norm(end);route=list(aa);J=set(route[-1][0])&set(E);Jb=set(bb[-1][0])&set(E);H=set(ha)
                        for s,t in zip(sorted(J-Jb),sorted(Jb-J)):
                            J.remove(s);J.add(t);route.append(refresh(J,H,r0))
                        for s,t in zip(sorted(H-set(hb)),sorted(set(hb)-H)):
                            H.remove(s);H.add(t);route.append(refresh(J,H,r0))
                        route.append(bb[-1]);route.extend(reversed(bb[:-1]));bound=45+costs[p]
                    else:route=[start,end];bound=sum(sorted((costs[s] for s in free),reverse=True)[:14])+helper_top4
                    need(route[0]==start and route[-1]==end,'literal prescribed endpoints restored')
                    base=None;previous=None;maxsupport=maxcost=0
                    for ci,(C,H) in enumerate(route):
                        selected=set(C)|fixedV|set(H)|{776,3139}
                        need(len(C)==k and len(H)==2 and add(charge(C),charge(H))==qa,'every arrival count and charge')
                        need(not (set(C)&fixed) and all((s in selected)==(s in fixedV) for s in fixed),'every frozen record unchanged')
                        need(368 not in selected and {776,3139}<=selected,'external owner4 preserved')
                        comp=set(free)-set(C);signed=add(sub(charge(free),charge(comp)),charge(H))
                        need(signed==qa,'core-only complement has negative sign and unchanged helper')
                        rowsphys=physical(selected)
                        if base is None:base=rowsphys
                        need(all((rowsphys[s][j]-base[s][j])%mods[j]==0 for s in rowsphys for j in range(6)),'literal all source/target dyadic rows at every arrival')
                        used=[]
                        for s in data:used.extend(resources[s][int(s in selected)])
                        need(len(used)==328 and len(set(used))==328,'capacity-one actual arrival resources')
                        if previous is not None:
                            diff=selected^previous;cost=sum(costs[s] for s in diff);support=len(diff)
                            need(cost<=bound and support<=18,'actual route event support and clocks')
                            maxsupport=max(maxsupport,support);maxcost=max(maxcost,cost)
                        previous=selected
                    route_results.append({'frozen_cheap':list(frozen),'protected_pair':bank,'repair':p,'free_outside':e,'free_grade':k,'events':len(route)-1,'max_support':maxsupport,'max_clocks':maxcost,'budget':bound})
                    if ix==0 and e==1 and k in (7,8,9,10,11):literal_routes.append({'frozen_V':sorted(fixedV),'repair':p,'free_core_grade':k,'route':[{'core_V':sorted(set(C)|fixedV),'helper_V':sorted(set(H)|{776,3139})} for C,H in route]})
    # Reconstruct the separately announced same-target relation failure.
    baseids=sorted(set(cheap)-{176,1425,2071,1457});aa=supports(baseids,range(6,11));lam=add(labels[1841],labels[1457])
    relational={str(k):sorted(set(range(2048))-(aa[k]|{add(lam,x) for x in aa[k-2]})) for k in (8,9,10)}
    need(relational=={'8':[477],'9':[],'10':[1098]},'same-target equality exact grade-shift hostile')
    receipt={'status':'FINITE-EXACT independent audit passed','input_commit':'c09b642a629b7d12b529b88c2697abaee78e73e7','universe':{'cheap':cheap,'outside':outside,'four_freeze_sets':4845,'base_subsets_enumerated':4845*65536,'charge_group':'Z/128 x (Z/2)^4','requested_single_repair_grades':[8,9]},'passing_histogram':dict(sorted(hist.items())),'universal_protected_choice_banks':banks,'singleton_failure_freezes':singleton,'both_active_two_slot_failures':pairfail,'cost_envelopes':dict(envelope),'global_boundary':{'core_top14':core_top14,'helper_top4':helper_top4,'clocks':core_top14+helper_top4,'slots':18,'ten_frozen_protected_pair_clocks':60,'eleven_frozen_chosen_repair_clocks':single_boundary},'budget_cases':budget_cases,'direct_combination_controls':direct_controls,'native_resources':{'slots':82,'owner_host':652,'reserved_words':656,'arrival_words':328},'outside_records':[{'source':s,'label':labels[s],'cost':costs[s],'freeze_clock':frozen_clock[s]} for s in outside],'endpoint_route_summary':{'routes':len(route_results),'events':sum(r['events'] for r in route_results),'max_actual_support':max(r['max_support'] for r in route_results),'max_actual_clocks':max(r['max_clocks'] for r in route_results)},'endpoint_controls':route_results,'literal_routes':literal_routes,'same_target_relation_holes':relational,'checks':dict(CHECKS)}
    receipt=json.loads(json.dumps(receipt))
    detailed={'four_freeze_results':rows}
    for filename,obj in [('receipt.json',receipt),('all_four_freeze_cases.json',detailed)]:
        target=P/filename
        if args.write:target.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
        else:need(json.loads(target.read_text())==obj,'saved exact output matches recomputation')
    print(json.dumps({'status':'PASS','freeze_cases':4845,'passing_histogram':dict(sorted(hist.items())),'protected_pairs':banks,'pair_failures':len(pairfail),'interior_clocks':[48,49],'global_clocks':68,'ten_frozen_clocks':60,'native_slots':82,'endpoint_routes':len(route_results),'checks':sum(CHECKS.values())},sort_keys=True))
if __name__=='__main__':main()
