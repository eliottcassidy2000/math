"""Independent controls for escape, transport, birth and symmetry certificates."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations
from math import factorial, gcd, isqrt
from pathlib import Path
import argparse
import importlib.util
import json
import carrier_scout as carrier
import birth_scout as birth
from decode17 import decode
from scout import rows,h_count

HERE=Path(__file__).resolve().parent
GATES=0

def need(b,why):
    global GATES
    GATES+=1
    if not b:raise RuntimeError(why)

def path_count(out):
    return sum(all(out[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(len(out))))

def relabel_mask(out,p):
    return sum(((out[p[a]]>>p[b])&1)<<e for e,(a,b) in enumerate(combinations(range(len(out)),2)))

def symmetries_and_escapes():
    receipt=json.loads((HERE/'escape_scout.json').read_text())
    traps={r['source'] for r in receipt['records']};need(len(traps)==768,'complete trap count')
    counts=[h_count(rows(6,x)) for x in range(1<<15)]
    independent=[0]*(1<<15);edges=list(combinations(range(6),2));indices={e:k for k,e in enumerate(edges)}
    for p in permutations(range(6)):
        fixed=value=0
        for a,b in zip(p,p[1:]):
            bit=1<<indices[min(a,b),max(a,b)];fixed|=bit
            if a<b:value|=bit
        free=[e for e in range(15) if not fixed>>e&1]
        for assignment in range(1<<10):
            x=value|sum(((assignment>>j)&1)<<e for j,e in enumerate(free));independent[x]+=1
    need(independent==counts,'complete path-cylinder and endpoint-DP six-vertex censuses agree')
    need(traps=={x for x,h in enumerate(counts) if (h-17)**2==4 and all((counts[x^(1<<e)]-17)**2>4 for e in range(15))},'complete trap predicate independently checked')
    target_scores={tuple(sorted(r.bit_count() for r in rows(6,x))) for x,h in enumerate(counts) if h==17}
    need(target_scores=={(1,1,2,3,4,4)},'all target17 score profiles agree')
    records={r['source']:r for r in receipt['records']}
    covered=set();types=[]
    for mask in (72,140,178,314):
        out=rows(6,mask);orbit={relabel_mask(out,p) for p in permutations(range(6))}
        group=[p for p in permutations(range(6)) if relabel_mask(out,p)==mask]
        need(not covered&orbit,'four ordinary source orbits are disjoint');covered|=orbit
        need(len(group) in (3,5) and len(orbit)*len(group)==factorial(6),'actual group and orbit orders')
        need(path_count(out)==15,'trap count by independent permutation enumeration')
        universal=[i for i,r in enumerate(out) if r.bit_count() in (0,5)]
        need(len(universal)==1,'one universal source or sink')
        core=[i for i in range(6) if i!=universal[0]]
        degrees=sorted(sum(out[i]>>j&1 for j in core) for i in core)
        need(degrees==([2]*5 if len(group)==5 else [1,2,2,2,3]),'five-vertex core classification')
        shell=2 if len(group)==3 else 3
        source_scores=sorted(r.bit_count() for r in out)
        need(sum(abs(a-b) for a,b in zip(source_scores,(1,1,2,3,4,4)))==2*shell,'score transport explains the sharp edit lower bound')
        targets={mask^sum(1<<e for e in es) for es in combinations(range(15),shell)
                 if h_count(rows(6,mask^sum(1<<e for e in es)))==17}
        need(len(targets)==(3 if len(group)==3 else 20),'full shortest shell census')
        targetorbits=[];remaining=set(targets)
        while remaining:
            y=min(remaining);orb={relabel_mask(rows(6,y),p) for p in group}
            need(orb<=targets and len(orb)==len(group),'source symmetry acts freely on target shell')
            remaining-=orb;targetorbits.append(sorted(orb))
        ports=[v for v in range(6) if sum(p[v]==v for p in group)==1]
        need(len(ports)==(3 if len(group)==3 else 5),'one moved vertex destroys residual stabilizer')
        need(universal[0] not in ports,'marking only the already fixed universal vertex is insufficient')
        for y in targets:need(path_count(rows(6,y))==17,'escape targets counted by permutations')
        types.append({'representative':mask,'orbit_size':len(orbit),'automorphism_order':len(group),
                      'core_scores':degrees,'shortest_flip_count':shell,'shortest_target_orbits':targetorbits,
                      'symmetry_breaking_vertex_ports':ports})
    need(covered==traps,'all traps are precisely these four ordinary classes')
    depth=json.loads((HERE/'escape_depth.json').read_text())
    need(depth['traps_reachable_at_peak16']==0,'strict lower barrier independently flooded')
    need(depth['all_source_distance_histogram']=={'0':1440,'1':13520,'2':16080,'3':1728},'target17 covering radius is three')
    for row in depth['trap_types']:
        need(row['key'][1:]==[36,3],'every trap has optimal barrier36 and three-step barrier path')
        for kind in ('shortest','minimax'):
            path=row[kind];masks=path['masks']
            need(all((a^b).bit_count()==1 for a,b in zip(masks,masks[1:])),'every intermediate step is a real edge flip')
            need([path_count(rows(6,x)) for x in masks]==path['counts'],'independent whole escape trajectory')
        decoded=decode(row['source'],36)
        need(len(decoded['masks'])==4 and decoded['counts'][-1]==17 and max((h-17)**2 for h in decoded['counts'])==36,'public constrained decoder reproduces the sharp escape')
    return types

def literal_carrier_audit():
    certificate=json.loads((HERE/'carrier_scout.json').read_text())
    cases=0
    for row in certificate:
        n=row['n'];free,ff,expand=carrier.chart(n);size=1<<len(free)
        graphs=[rows(n,expand(x)) for x in range(size)]
        for direction,i,j in row['relocation_labels']:
            p=list(range(n));v=p.pop(i);p.insert(j,v)
            table=carrier.permutation_map(n,p,free,ff)
            for x,out in enumerate(graphs):
                legal=all(out[a]>>b&1 for a,b in zip(p,p[1:]))
                expected=sum(((out[p[a]]>>p[b])&1)<<e for e,(a,b) in enumerate(free)) if legal else -1
                need(table[x]==expected,'compiled transition agrees with literal arcs and path guard');cases+=1
        need(row['legality_only']['classes']==size,'legality alone reconstructs the full tested cube')
        need(row['C_relocations']['classes']==size,'count-observer closure retains all tested orientations')
    last=certificate[-1]
    need(last['C_one_mark']['classes']==10455 and last['C_one_mark_converse']['classes']==16149,'restricted exact minima')
    return {'literal_transition_cases':cases,'n7_C_only':3008,'n7_C_and_R':10455,
            'n7_C_R_converse':16149,'n7_all_relocations':32768,
            'all_order_injectivity':'OPEN; exact through seven'}

def birth_audit():
    E,P=birth.native.bank(25);EE,PP=birth.native.bank(26);new=[e for e in EE if e not in E]
    old=[birth.native.encode(p,E) for p in P];fresh=[birth.native.encode(p,E+new) for p in PP]
    A=birth.active(old,32);B=birth.active(fresh,34);mask=(1<<32)-1
    need(len(A)==4864 and len(B)==11232,'literal source counts')
    dead={x&mask for x in B if (x&mask) not in A}
    need(len(dead)==9184,'all old-dead parents recovered')
    # Independent cylinder intersection formula for second moments, not completion counting.
    for a in range(4):
        first=second=0
        for m,v in fresh:
            if ((a<<32)^v)&m&~mask:continue
            first+=1<<(32-(m&mask).bit_count())
            for mm,vv in fresh:
                if ((a<<32)^vv)&mm&~mask or (v^vv)&m&mm:continue
                second+=1<<(32-((m|mm)&mask).bit_count())
        r=int(a in (1,2))
        need(first==512+5120*r and second==512+6240*r,'first and second moments by path-pair cylinders')
        direct=Counter(h for x,h in B.items() if x>>32==a)
        need(direct==Counter({1:512+4096*r,2:464*r,3:32*r}), 'complete positive output distribution')
        need(sum((x&mask) not in A for x in B if x>>32==a)==4592*r,'opposite-port resurrection law')
    certificate=json.loads((HERE/'birth_scout.json').read_text())
    decoders=[]
    for row in certificate['dead_source_decoders']:
        x=row['source25_mask'];y=row['source26_mask'];path=row['square_sum_path']
        need(sum(x&m==v for m,v in old)==0,'actual source has zero guarded old paths')
        need(sum(y&m==v for m,v in fresh)==1,'actual target has one guarded new path')
        need(sorted(path)==list(range(1,27)) and all(isqrt(a+b)**2==a+b for a,b in zip(path,path[1:])),'native square-sum Hamilton path')
        need(y&mask==x,'old orientation bits preserved exactly')
        need(F(row['initial_energy'])==1-F(4144,(1<<32)-4864),'exact conditional old-dead energy')
        for step in row['trace']:
            vals=[None if v=='None' else F(v) for v in step['conditional_energies']]
            pops=step['conditional_populations'];before=F(step['old_energy'])
            need(sum(v*m for v,m in zip(vals,pops) if m)/sum(pops)==before,'conditional energies use actual population weights')
            need(vals[step['chosen']]<=before,'decoder energy never increases')
        decoders.append({k:v for k,v in row.items() if k!='trace'})
    return {'two_port_disagreement_parameter':'r=P(bit10 != bit23); independent of old uniform source',
            'mean':'(1+10r)/2^23','second_moment':'(16+195r)/2^27',
            'old_dead_activation_probability':'287r/2^28','dead_source_decoders':decoders}

def hostile_control():
    file=HERE.parent/'regular_parent_symmetry_20261029'/'verify.py'
    spec=importlib.util.spec_from_file_location('symmetry_source',file);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    inherited=module.hostile_checks();sources=[]
    for flip in (False,True):
        out=[0]*12
        for u,v in combinations(range(12),2):
            b,s=divmod(u,4);c,t=divmod(v,4)
            arc=(c-b)%3==1 if b!=c else s==0 or (s>0 and (t-s)%3==1)
            if flip and b==c==2:arc=not arc
            out[u if arc else v]|=1<<(v if arc else u)
        sources.append(out)
    # Three ordered C3 components and three singleton components, on the same named vertices.
    blocks=[(0,),(1,2,3),(4,),(5,6,7),(8,),(9,10,11)];target=[0]*12
    for i,block in enumerate(blocks):
        if len(block)==3:
            for u,v in zip(block,block[1:]+block[:1]):target[u]|=1<<v
        for later in blocks[i+1:]:
            for u in block:
                for v in later:target[u]|=1<<v
    need(h_count(target)==27,'actual common-vertex target count27')
    groups=[module.iso_maps(out,out) for out in sources]
    def preserves(p):return all(((target[u]>>v)&1)==((target[p[u]]>>p[v])&1) for u,v in combinations(range(12),2))
    need(all(preserves(p) for p in groups[1]),'all27 source automorphisms preserve target')
    need(not all(preserves(p) for p in groups[0]),'81 source cannot preserve this target')
    need(27%len(groups[0])!=0 and 27%len(groups[1])==0,'all-target divisibility obstruction separates hostile pair')
    orbit={relabel_mask(target,p) for p in groups[0]}
    need(len(orbit)==3,'sharp three-object equivariant support for the81 source')
    return {'inherited':inherited,'H27_target_out_rows':target,
            'equivariant_map_on_27_orbit':True,'equivariant_map_on_81_orbit':False,
            'minimum_equivariant_target_supports':[3,1]}

def arithmetic_controls():
    cycle=[4,8,3,7,2,6,10,5,9,4]
    step=lambda x:2+(x-2+4)%9
    need(all(step(a)==b for a,b in zip(cycle,cycle[1:])),'exact requested cycle as translation modulo nine')
    cycletournament=[sum(1<<j for j in range(9) if (j-i)%9 in (1,2,3,4)) for i in range(9)]
    cycleH=h_count(cycletournament)
    need(cycleH%9==0 and path_count(cycletournament)==cycleH,'cyclic tournament path count independently retains the nine-phase symmetry')
    square=[[9,2,7],[4,6,8],[5,10,3]]
    image=[[step(x) for x in row] for row in square]
    lines=square+[list(col) for col in zip(*square)]+[[square[i][i] for i in range(3)],[square[i][2-i] for i in range(3)]]
    need(all(sum(row)==18 for row in lines) and len({sum(row) for row in image})>1,'nine-cycle does not preserve magic equations')
    # Exact genus formula for X_0(13), using the two elliptic-point counts.
    p=13;e2=sum((x*x+1)%p==0 for x in range(p));e3=sum((x*x+x+1)%p==0 for x in range(p))
    genus=1+F(p+1,12)-F(e2,4)-F(e3,3)-1
    need(genus==0 and (e2,e3)==(2,2),'genus-zero weight-two obstruction')
    # Direct five-vertex isomorphism/converse census.
    full=(1<<10)-1;canonical={}
    for x in range(1024):
        if x not in canonical:
            orb={relabel_mask(rows(5,x),p) for p in permutations(range(5))};key=min(orb)
            canonical.update({y:key for y in orb})
    classes=set(canonical.values());selfdual=[x for x in classes if canonical[x^full]==x]
    need((len(classes),len(selfdual))==(12,8),'valid pasted five-vertex count')
    return {'cycle':cycle,'magic_array':square,'cycle_image':image,'image_row_sums':list(map(sum,image)),
            'cyclic_tournament_out_rows':cycletournament,'cyclic_tournament_H':cycleH,
            'genus_X0_13':str(genus),'five_vertex_classes':12,'self_converse_classes':8,'merged_classes':10}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE/'verification.json');args=parser.parse_args()
    result={'trap_orbits':symmetries_and_escapes(),'carrier':literal_carrier_audit(),
            'birth':birth_audit(),'hostile81_27':hostile_control(),'arithmetic':arithmetic_controls()}
    result['gates']=GATES
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'gates':GATES,'trap_orbits':[(z['orbit_size'],z['automorphism_order']) for z in result['trap_orbits']],
                      'carrier':result['carrier'],'birth_decoders':len(result['birth']['dead_source_decoders'])}))

if __name__=='__main__':main()
