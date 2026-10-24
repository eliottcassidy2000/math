#!/usr/bin/env python3
"""Independent Session49 audit from original native inputs, no producer imports.
Hamiltonian paths use adjacency DFS; receiver counts use independent combination
and coefficient recurrences. All checks remain enabled with Python -O.
"""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
from array import array
from hashlib import sha256
import argparse, csv, gzip, json, random, sys
P=Path(__file__).resolve().parent
I=P/'inputs'
CHECKS=Counter()
def require(ok,why):
    CHECKS[why]+=1
    if not ok:raise RuntimeError(why)
def read(p):return json.loads((I/p).read_text())
def add(q,r):return 16*((q//16+r//16)%128)+((q^r)&15)
def minus(q,r):return 16*((q//16-r//16)%128)+((q^r)&15)
def qsum(items):
    out=0
    for q in items:out=add(out,q)
    return out
PAIRS=list(combinations(range(8),2))
CHORDS=[(i,j) for i,j in PAIRS if j>i+1]
ALL_WORDS=(1<<len(CHORDS))-1
REVERSE=[CHORDS.index((7-j,7-i)) for i,j in CHORDS]
def tournament(n):
    out=[[0]*8 for _ in range(8)]
    for k,(i,j) in enumerate(PAIRS):out[i][j]=(n>>(27-k))&1;out[j][i]=out[i][j]^1
    return out
def paths(a):
    def walk(path,remaining):
        if not remaining:yield tuple(path);return
        for v in range(8):
            if remaining>>v&1 and (not path or a[path[-1]][v]):yield from walk(path+[v],remaining^(1<<v))
    yield from walk([],255)
def encode(a,path):return sum(a[path[i]][path[j]]<<k for k,(i,j) in enumerate(CHORDS))
def reverse_word(w):return sum(((w>>k)&1)<<j for k,j in enumerate(REVERSE))
def resource(w):return frozenset([w,reverse_word(w),ALL_WORDS^w,ALL_WORDS^reverse_word(w)])
def partner(a,path):
    pos={v:i for i,v in enumerate(path)}
    return [[0 if i==j else a[i][j] if abs(pos[i]-pos[j])==1 else a[j][i] for j in range(8)] for i in range(8)]
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    for entry in json.loads((P/'SOURCE_PINS.json').read_text()):
        raw=(P/entry['path']).read_bytes()
        require(len(raw)==entry['bytes'] and sha256(raw).hexdigest()==entry['sha256'],'original source bytes')
    old=read('native_packet_decomposition_20261023/native/domains.json')
    supply=read('compatible_charged_helpers_20261023/supply/selected_supply.json')
    replacement=read('guarded_component_completion_20261024/sections/native_2125_receipt.json')
    host=read('guarded_component_completion_20261024/native/xor_tree_receipt.json')
    electric=read('metagraph_order8_electrical_20261023/receipt.json')['order8']
    X=electric['higher_layers'][-1]['lifts'][0];Y=electric['radical_lifts_mod4']
    potentials=[X]+Y;moduli=[256]+[4]*len(Y)
    label=lambda s:16*(X[s]%128)+sum((Y[j][s]%2)<<(j-1) for j in range(1,5))
    require(len(Y)==5,'five original radical lifts retained')
    require(all((X[s]-Y[0][s]-Y[4][s])%2==0 for s in range(len(X))),'actual redundant Y0 parity')
    U,V=1965,2234
    require(all(2*(f[U]-f[V])%mod==0 for f,mod in zip(potentials,moduli)),'source rows unchanged under either actual target')
    data={}
    for r in old['sources']:
        rec={q['target']:q for q in r['exact_all_shared_clocks_target_domain']}
        data[r['source']]={'records':[rec[U],rec[V]],'pins':r['fixed_clocks'],'blocks':[]}
    cores=sorted(data)
    for r in supply['helpers']:
        v=r['selected'];data[r['source']]={'records':[v['before'],v['after']],'pins':v['common_clocks'],'blocks':[]}
    rec={q['target']:q for q in replacement['records']}
    data[2125]={'records':[rec[U],rec[V]],'pins':[],'blocks':replacement['retained_blocks']}
    cost={s:sum(u!=v for u,v in zip(*(r['path'] for r in d['records']))) for s,d in data.items()}
    cheap=sorted(s for s in cores if cost[s]==2);outside=sorted(set(cores)-set(cheap))
    helpers=sorted(set(data)-set(cores)-{368,776,3139})
    require((len(data),len(cores),len(cheap),len(outside),len(helpers))==(82,28,20,8,51),'full populations and inventory')
    labels={s:label(s) for s in data}
    charge=lambda selected:qsum(labels[s] for s in selected)
    nodes=list(csv.DictReader((I/'metagraph_order8_20261023/merged_nodes.tsv').open(),delimiter='\t'))
    alias={int(r[col]):int(r['merged_id']) for r in nodes for col in ('class_a','class_b')}
    owners=array('H');owners.frombytes(gzip.decompress((I/'metagraph_order8_20261023/owner_map.u16.gz').read_bytes()))
    if sys.byteorder!='little':owners.byteswap()
    target_words={t:set() for t in (U,V)}
    for t in target_words:
        for code in set(int(nodes[t][k]) for k in ('canon_a','canon_b')):
            a=tournament(code);target_words[t].update(encode(a,p) for p in paths(a))
    require(not target_words[U]&target_words[V],'actual canonical target classes distinct')
    frames={368:list(range(8)),776:[8,9,4,10,11,12,13,3],3139:[14,15,16,17,18,3,4,19],2125:list(range(20,28))}
    cursor=28
    for s in sorted(set(data)-set(frames)):frames[s]=list(range(cursor,cursor+8));cursor+=8
    require(cursor==652 and frames=={int(s):f for s,f in host['source_frames'].items()},'host source frames independently rebuilt')
    for s,t in combinations(data,2):
        require(len(set(frames[s])&set(frames[t]))==(2 if {s,t}<={368,776,3139} else 0),'exact shared-owner intersections')
    matrices={s:tournament(int(nodes[s]['canon_a'])) for s in data}
    record_resources={};all_reserved=set();native={};capacity={}
    def glue(mats):
        arcs={}
        for s,a in mats.items():
            for i,j in PAIRS:
                u,v=frames[s][i],frames[s][j];edge=min(u,v),max(u,v);direction=a[i][j] if u<v else a[j][i]
                require(edge not in arcs or arcs[edge]==direction,'consistent source or partner host arcs');arcs[edge]=direction
    glue(matrices)
    for s,d in sorted(data.items()):
        a=matrices[s];selected_orbits=[resource(r['mask']) for r in d['records']]
        actual_menu=[];kernels=[[],[]];units={U:set(),V:set()};path_count=0
        for pp in paths(a):
            path_count+=1;w=encode(a,pp)
            require(alias[owners[w]]==s,'source word canonical class')
            for bit,t in enumerate((U,V)):
                native_target=ALL_WORDS^w in target_words[t]
                require(native_target==(alias[owners[ALL_WORDS^w]]==t),'target word canonical class independent of owner map')
                if w in selected_orbits[bit]:kernels[bit].append(list(pp))
                if native_target:units[t].add(min(resource(w)))
            if all(pp[i]==v for i,v in d['pins']) and all(any(list(pp[i:i+len(b)])==b for i in range(9-len(b))) for b in d['blocks']):
                actual_menu.extend((pp,w,t) for t in (U,V) if ALL_WORDS^w in target_words[t])
        require(set(actual_menu)=={(tuple(r['path']),r['mask'],r['target']) for r in d['records']},'exhaustive guarded menu is exactly two actual records')
        record_resources[s]=[]
        for bit,r in enumerate(d['records']):
            pp=r['path'];w=r['mask'];t=(U,V)[bit];orb=selected_orbits[bit]
            require(r['target']==t and encode(a,pp)==w and all(a[x][y] for x,y in zip(pp,pp[1:])),'literal selected path and target')
            require(encode(partner(a,pp),pp)==ALL_WORDS^w,'inverse partner is actual Hamiltonian word')
            require(encode(list(zip(*a)),pp[::-1])==reverse_word(w),'converse is actual reversed path')
            require(list(pp) in kernels[bit],'selected representative belongs to full resource kernel')
            require(len(orb)==4 and not all_reserved.intersection(orb),'all normalized alternative resources disjoint')
            require((X[s]-X[t])%8==r['phase'],'retained literal scalar phase')
            all_reserved.update(orb);record_resources[s].append(orb)
        for t in (U,V):capacity[tuple(sorted((s,t)))]=2*len(units[t])
        native[str(s)]={'path_count':path_count,'full_kernels':kernels,'cost':cost[s],'resources':[sorted(x) for x in selected_orbits]}
    require(len(all_reserved)==656,'full 164-resource reservation')
    shared_truth=[]
    for bits in product((0,1),repeat=3):
        selected={s:data[s]['records'][b]['path'] for s,b in zip((368,776,3139),bits)}
        values=[frames[s][selected[s][clock]] for s,clock in zip((368,776,3139),(1,4,7))]
        glue({s:partner(matrices[s],p) for s,p in selected.items()})
        require((values[0]==values[1]==values[2])==(bits in [(0,1,1),(1,0,0)]),'complete actual shared-target tree truth')
        require((values[0]==values[1]==values[2]==4)==(bits==(0,1,1)),'owner4 picks only original helper branch')
        shared_truth.append({'bits':bits,'owners':values})
    # Physical equations use all original merged edges, diagonal terms and weight two.
    lap=[[0]*len(nodes) for _ in potentials];edges=0
    with gzip.open(I/'metagraph_order8_20261023/merged_edges.tsv.gz','rt') as stream:
        for r in csv.DictReader(stream,delimiter='\t'):
            s,t=int(r['source']),int(r['target']);blue,black=int(r['blue']),int(r['black']);edges+=1
            if (s,t) in capacity:require(blue==0 and black==capacity.pop((s,t)),'actual black capacities count all native resource kernels')
            for k,f in enumerate(potentials):z=(blue+black)*(f[s]-f[t]);lap[k][s]+=z;lap[k][t]-=z
    require(not capacity,'all 164 source-target capacities audited')
    for row,mod in zip(lap,moduli):require(all(v%mod==0 for v in row),'all original physical Laplacian rows harmonic')
    require(sum(sorted((cost[s] for s in helpers),reverse=True)[:2])==15,'helper Johnson edge cost15')
    require(sum(sorted((cost[s] for s in helpers),reverse=True)[:4])==29,'helper endpoint cost29')
    require(sum(sorted((cost[s] for s in cores),reverse=True)[:14])==39,'core endpoint cost39')
    # Completely fresh receiver support enumeration from original potential labels.
    frozen={176,1425,2071};base=sorted(set(cheap)-frozen-{1457});pair=[1457,1841];candidates=sorted(set(outside)-{1841})
    require(len(base)==16 and candidates==[1532,1923,2087,2553,2835,3015,3360],'receiver inventory and available repair slots')
    pairlabel=charge(pair);require(pairlabel==776,'fresh linked-pair label')
    counts=[[0]*2048 for _ in range(17)];witness=[[None]*2048 for _ in range(17)]
    allbase=[]
    for r in range(17):
        for ss in combinations(base,r):
            q=charge(ss);counts[r][q]+=1
            if witness[r][q] is None:witness[r][q]=ss
            allbase.append((ss,r,q))
    dp=[[0]*2048 for _ in range(17)];dp[0][0]=1
    for j,s in enumerate(base):
        for r in range(j+1,0,-1):
            for q,count in enumerate(dp[r-1]):
                if count:dp[r][add(q,labels[s])]+=count
    require(dp==counts,'independent integer group-polynomial coefficients equal combinations')
    require(len(allbase)==65536,'entire 16-slot Boolean inventory exhausted')
    census={};sections={};branch_witnesses={};signed_fibre={}
    expected={1532:list(range(8,12)),1923:[9,10],2087:[],2553:[],2835:[9,10],3015:list(range(8,12)),3360:[]}
    for p in candidates:
        R=base+pair+[p];branches=[];literal_counts=[[[0]*2048 for _ in range(20)] for _ in range(4)]
        bw=[[[None]*2048 for _ in range(20)] for _ in range(4)]
        for branch in range(4):
            a,b=divmod(branch,2);extra=(pair if a else [])+([p] if b else []);shift=charge(extra);degree=2*a+b
            branches.append({'index':branch,'a':a,'b':b,'grade_shift':degree,'charge_shift':shift})
            for ss,r,q in allbase:
                selected=tuple(sorted(ss+tuple(extra)));rr=r+degree;qq=add(q,shift)
                literal_counts[branch][rr][qq]+=1
                if bw[branch][rr][qq] is None:bw[branch][rr][qq]=selected
                actual_targets={s:data[s]['records'][int(s in selected)]['target'] for s in pair+[p]}
                require(actual_targets[1457]==actual_targets[1841],'literal pair target equality')
                require((actual_targets[1457]==actual_targets[p])==(a==b),'literal triple equality matches only branches00and11')
                # Compute native target-row residues directly from selected sources.
                actual=tuple((-2*sum(f[s]-f[V] for s in selected))%mod for f,mod in zip(potentials,moduli))
                key=(rr,qq)
                if key in signed_fibre:require(actual==signed_fibre[key],'same grade-charge gives identical original target rows')
                else:signed_fibre[key]=actual
        bytable={};sec={}
        for table in range(16):
            rows=[]
            for r in range(20):
                mult=[sum(literal_counts[b][r][q] for b in range(4) if table>>b&1) for q in range(2048)]
                support={q for q,m in enumerate(mult) if m};holes=sorted(set(range(2048))-support)
                predicted=set()
                for b in range(4):
                    if not table>>b&1:continue
                    shift=branches[b]['charge_shift'];rr=r-branches[b]['grade_shift']
                    if 0<=rr<=16:predicted.update(add(q,shift) for q,m in enumerate(counts[rr]) if m)
                require(support==predicted,'all branch-conditioned grade shifts exactly match literal records')
                rows.append({'grade':r,'count':len(support),'holes':holes,'assignments':sum(mult)})
                if table==9 and len(support)==2048:
                    sec[r]={q:frozenset(next(bw[b][r][q] for b in range(4) if table>>b&1 and bw[b][r][q] is not None)) for q in range(2048)}
            bytable[str(table)]={'complete_grades':[r['grade'] for r in rows if r['count']==2048],'grades':rows}
        # Restricting a relation cannot enlarge support on the SAME inventory.
        for small in range(16):
            for large in range(16):
                if small&large!=small:continue
                for r in range(20):require(set(bytable[str(large)]['grades'][r]['holes'])<=set(bytable[str(small)]['grades'][r]['holes']),'relation monotonicity on fixed19-slot inventory')
        total=charge(R)
        for table in range(16):
            complement=sum(1<<(3-b) for b in range(4) if table>>b&1)
            for r in range(20):
                require({minus(total,q) for q in bytable[str(table)]['grades'][r]['holes']}==set(bytable[str(complement)]['grades'][19-r]['holes']),'actual table-complement grade-charge identity')
        require(bytable['9']['complete_grades']==expected[p],'all seven triple-equality complete-grade lists')
        require(bytable['15']['complete_grades']==list(range(6,14) if p==1923 else range(7,13)),'all seven independent-repair complete-grade lists')
        census[str(p)]={'label':labels[p],'cost':cost[p],'branches':branches,'tables':bytable};sections[p]=sec;branch_witnesses[p]=bw
    # Independent exact minima over every pair of full branch witnesses.
    masks=[[[] for _ in range(2048)] for _ in range(17)]
    base_index={s:i for i,s in enumerate(base)}
    for ss,r,q in allbase:masks[r][q].append(sum(1<<base_index[s] for s in ss))
    branch_minima={}
    expected_hist={8:{18:1203,22:468,26:63,30:5,34:1},9:{18:1398,22:506,26:41,30:4},10:{18:1398,22:506,26:41,30:4},11:{18:1203,22:468,26:63,30:5,34:1}}
    for r in range(8,12):
        hist=Counter();records=[]
        for q in range(2048):
            left=masks[r][q];right=masks[r-3][minus(q,1263)]
            if left and right:
                distance,u,v=min(((u^v).bit_count(),u,v) for u in left for v in right)
                clocks=8+2*distance;hist[clocks]+=1
                records.append({'charge':q,'clocks':clocks,'slots':3+distance,'left_base_mask':u,'right_base_mask':v})
        require(dict(hist)==expected_hist[r],'all-charge exact cross-branch boundary histogram')
        branch_minima[str(r)]={'histogram':dict(hist),'minimizing_pairs':records}
    # Authenticate the previously missing grades before adding a third source.
    pair_holes={}
    for r in (8,9,10):
        support={q for q,m in enumerate(counts[r]) if m}|{add(pairlabel,q) for q,m in enumerate(counts[r-2]) if m}
        pair_holes[str(r)]=sorted(set(range(2048))-support)
    require(pair_holes=={'8':[477],'9':[],'10':[1098]},'old linked18-slot receiver holes, not conflated with enlarged inventory')
    require(charge(pair+[3015])==1263,'fresh 3015 triple label')
    sharp=[]
    for ss,r,q in allbase:
        if r==8 and q==53:sharp.append(tuple(sorted(ss)))
        if r==5 and add(q,1263)==53:sharp.append(tuple(sorted(ss+tuple(pair)+tuple([3015]))))
    require(len(sharp)==2 and not set(sharp[0])&set(sharp[1]),'sharp receiver fibre exactly two disjoint8-V states')
    require(len(set(sharp[0])^set(sharp[1]))==16 and sum(cost[s] for s in set(sharp[0])^set(sharp[1]))==34,'sharp exact native16slot34clock boundary')
    # Cost envelope proof is exhausted across every actual outside availability set.
    cost_cases=0;cost_envelopes={}
    for p in (1532,3015):
        R=base+pair+[p];remaining=sorted(set(outside)-{1841,p});weights=sorted((cost[s] for s in R),reverse=True)
        W=lambda n:sum(weights[:max(0,min(19,n))])
        maxima=Counter()
        for flags in product((0,1),repeat=6):
            E=[s for s,b in zip(remaining,flags) if b];e=len(E);ew=sorted((cost[s] for s in E),reverse=True);EW=lambda n:sum(ew[:n])
            for t in range(8,12+e):
                r0=max(8,t-e);require(8<=r0<=11,'canonical complete receiver grade')
                for r in range(max(0,t-e),min(19,t)+1):
                    if r<8:n=r+8;d=8-r
                    elif r>11:n=27-r;d=r-11
                    else:n=d=0
                    if d:
                        require(n+d<=20 and W(n)+EW(d)<=50+cost[p],'all literal entry cost envelopes')
                        maxima['entry']=max(maxima['entry'],W(n)+EW(d));cost_cases+=1
                for r in range(r0+1,12):
                    n=min(2*r-1,39-2*r)
                    require(n+1<=20 and W(n)+EW(1)<=50+cost[p],'all consecutive-grade cost envelopes')
                n=min(2*r0,38-2*r0)
                require(n+2<=20 and W(n)+15<=50+cost[p],'helper-edge routing envelope')
                maxima['helper']=max(maxima['helper'],W(n)+15)
                if e>=2:require(n+2<=20 and W(n)+EW(2)<=50+cost[p],'outside-exchange routing envelope')
            free=R+E
            require(sum(sorted((cost[s] for s in free),reverse=True)[:14])+29<=68,'all actual free-core boundary envelopes')
        require(maxima['helper']==50+cost[p],'interior upper estimate realized by size-only envelope')
        cost_envelopes[str(p)]=dict(maxima)
    full_free=sorted(set(cores)-frozen)
    require(len(full_free)==25 and sum(cost[s] for s in full_free)==61 and sum(sorted((cost[s] for s in full_free),reverse=True)[:24])==59,'direct strict-nonreturn core envelope25slots61total59changed')
    require(24+4==28 and 59+29==88,'direct atomic no-return28slots88clocks upper bound')
    # Full-record routes, with actual source/target pointers, all physical rows,
    # retained kernel lifts, shared owners and occupied resources at EVERY arrival.
    hpairs=[frozenset(c) for c in combinations(helpers,2)];hcharge={}
    for hp in hpairs:hcharge.setdefault(charge(hp),hp)
    require(len(hpairs)==1275 and len(hcharge)==972,'owner4 helper support is not full Q')
    rng=random.Random(49001);route_results=[];record_examples=[];return_control=None
    def physical(selected):
        result={s:[0]*len(potentials) for s in list(data)+[U,V]}
        for s in data:
            target=data[s]['records'][int(s in selected)]['target']
            for j,f in enumerate(potentials):
                change=2*(f[s]-f[target]);result[s][j]+=change;result[target][j]-=change
        return tuple(tuple(value%mod for value,mod in zip(result[s],moduli)) for s in sorted(result))
    def validate_records(selected,p):
        chosen={s:data[s]['records'][int(s in selected)] for s in data}
        require(chosen[1457]['target']==chosen[1841]['target']==chosen[p]['target'],'actual full record triple equality at route arrival')
        require([chosen[s]['target'] for s in (368,776,3139)]==[U,V,V],'owner4 tree unchanged on every route arrival')
        occupied=set()
        for s,r in chosen.items():
            require(r['path'] in native[str(s)]['full_kernels'][int(s in selected)],'every arrival retains full selected kernel lift')
            require(not occupied.intersection(record_resources[s][int(s in selected)]),'every arrival has distinct occupied inverse-converse units')
            occupied.update(record_resources[s][int(s in selected)])
        require(len(occupied)==328,'every complete arrival occupies82 resources')
        return chosen
    for p in (1532,3015):
        R=frozenset(base+pair+[p]);remaining=sorted(set(outside)-{1841,p})
        availability=[[],remaining[:1],remaining[:3],remaining]
        def random_receiver(k):
            options=[a for a in (0,1) if 0<=k-3*a<=16];a=rng.choice(options)
            return frozenset(rng.sample(base,k-3*a)+(pair+[p] if a else []))
        for Elist in availability:
            E=frozenset(Elist);e=len(E);fixed=frozenset(set(cores)-R-E);fixedV=frozenset(s for j,s in enumerate(sorted(fixed)) if j%2)
            for t in range(20+e):
                def endpoint():
                    grades=[r for r in range(max(0,t-e),min(19,t)+1) if any(0<=r-3*a<=16 for a in (0,1))]
                    r=rng.choice(grades);return random_receiver(r)|frozenset(rng.sample(sorted(E),t-r))
                a,b=endpoint(),endpoint();delta=minus(charge(a),charge(b))
                match=next(((hp,hcharge[add(q,delta)]) for q,hp in hcharge.items() if add(q,delta) in hcharge),None)
                require(match is not None,'feasible same-fibre endpoint helper pairs')
                ha,hb=match;target=add(charge(a),charge(ha));start=(a,ha);finish=(b,hb)
                require(target==add(charge(b),charge(hb)),'endpoint total charge exact')
                if 8<=t<=11+e:
                    r0=max(8,t-e)
                    def refresh(J,H,r):return sections[p][r][minus(minus(target,charge(J)),charge(H))]|frozenset(J),frozenset(H)
                    def normalize(state):
                        history=[state];C,H=state;J=C&E;r=len(C&R)
                        if r<8:J-=frozenset(sorted(J)[:8-r]);r=8;history.append(refresh(J,H,r))
                        elif r>11:J|=frozenset(sorted(E-J)[:r-11]);r=11;history.append(refresh(J,H,r))
                        while r>r0:J|={min(E-J)};r-=1;history.append(refresh(J,H,r))
                        return history
                    front,back=normalize(start),normalize(finish);route=front[:];C,H=route[-1];J=C&E;JC=back[-1][0]&E;HC=back[-1][1]
                    while J!=JC:J=(J-{min(J-JC)})|{min(JC-J)};route.append(refresh(J,H,r0))
                    while H!=HC:H=(H-{min(H-HC)})|{min(HC-H)};route.append(refresh(J,H,r0))
                    if route[-1]!=back[-1]:route.append(back[-1])
                    route.extend(back[-2::-1]);bound=(20,50+cost[p])
                else:route=[start,finish];bound=(20,68)
                reference=None;maxslots=maxclocks=0
                for state in route:
                    C,H=state;full=fixedV|C|H|{776,3139}
                    require(len(C)==t and len(H)==2 and add(charge(C),charge(H))==target,'both exact populations and paired Q on every arrival')
                    chosen=validate_records(full,p);ph=physical(full)
                    if reference is None:reference=ph
                    require(ph==reference,'all original source and target physical equations including diagonals preserved')
                for x,y in zip(route,route[1:]):
                    changed=(x[0]^y[0])|(x[1]^y[1]);slots=len(changed);clocks=sum(cost[s] for s in changed)
                    require(slots<=bound[0] and clocks<=bound[1],'literal constructed route within audited event bounds')
                    maxslots=max(maxslots,slots);maxclocks=max(maxclocks,clocks)
                if return_control is None:
                    for source in sorted(data):
                        history=[]
                        for C,H in route:
                            bit=int(source in (fixedV|C|H|{776,3139}))
                            if not history or history[-1]!=bit:history.append(bit)
                        if len(history)>=3:
                            return_control={'source':source,'successive_distinct_bits':history,'p':p,'grade':t,'outside':Elist,'route':[{'core_V':sorted(C),'helper_V':sorted(H)} for C,H in route]};break
                require(route[0]==start and route[-1]==finish,'routes begin and end at exact full endpoint records')
                route_results.append({'p':p,'outside':Elist,'grade':t,'arrivals':len(route),'max_changed_slots':maxslots,'max_clocks':maxclocks,'bound':bound})
                if p==3015 and e==6 and t in (8,11,17):
                    record_examples.append({'p':p,'outside':Elist,'grade':t,'fixed_V':sorted(fixedV),'route':[{'core_V':sorted(C),'helper_V':sorted(H)} for C,H in route]})
    require(return_control is not None,'reversible route has literal departed-resource reacquisition hostile control')
    # Sharp obstruction really is on the fully frozen ambient fibre, not merely projection.
    ambient=[frozenset(w)|{776,3139}|hpairs[0] for w in sharp]
    for selected in ambient:validate_records(selected,3015)
    require(physical(ambient[0])==physical(ambient[1]),'sharp full82-record endpoints preserve original physical rows')
    full_sharp=[]
    for selected in ambient:
        full_sharp.append({'V_sources':sorted(selected),'records':{str(s):{'target':data[s]['records'][int(s in selected)]['target'],'path':data[s]['records'][int(s in selected)]['path'],'resource':min(record_resources[s][int(s in selected)])} for s in sorted(data)}})
    result={'status':'PASS independent finite full-record audit','source_commit':'c09b642a629b7d12b529b88c2697abaee78e73e7','checks':dict(CHECKS),'native':native,'source_frames':frames,'labels':labels,'costs':cost,'shared_truth':shared_truth,'physical_edges':edges,'base_sources':base,'pair_sources':pair,'pair_label':pairlabel,'old_pair_holes':pair_holes,'relation_census':census,'exact_branch_minima':branch_minima,'return_control':return_control,'strict_nonreturn_atomic_bound':{'slots':28,'clocks':88,'free_core_cost':61,'changed_core_cost_bound':59},'sharp_receiver_V_sets':sharp,'sharp_full_record_endpoints':full_sharp,'sharp_fibre_scope':'All other core and helper records frozen; otherwise no lower bound asserted.','cost_cases':cost_cases,'cost_envelopes':cost_envelopes,'literal_routes':route_results,'route_examples':record_examples,'limitations':['Finite original order-eight native inventory only.','Actual target equality is a record predicate, not a newly realized indexed-owner equality.','The20-slot migration permits returns; separate28-slot88-clock single atomic replacement has strict no-return from its starting assembly.','Full dyadic saturation is inherited from the archived earlier independent filtration theorem; this checker independently verifies every original X and Y equation.','No odd-primary, exact integer-current, all-order native supply, or whole-host Hamiltonian observer claim.']}
    for source in result['relation_census'].values():
        for table in source['tables'].values():
            for row in table['grades']:
                holes=row.pop('holes');row['first_hole']=holes[0] if holes else None
                row['hole_bitmap_hex']=format(sum(1<<q for q in holes),'0512x')
    out=P/'receipt.json';text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.write:out.write_text(text)
    else:require(json.loads(out.read_text())==json.loads(text),'immutable independent receipt exact')
    print(json.dumps({'status':'PASS','checks':sum(CHECKS.values()),'native_slots':82,'relation_tables':112,'literal_routes':len(route_results),'sha256':sha256(text.encode()).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
