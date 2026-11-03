"""Exact source-preserving chart and multi-edge response experiments."""
from collections import Counter, deque
from itertools import combinations, permutations
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
R=(0,5,6,3,4,1,2)


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def rows(n, mask, fixed_path=False):
    out=[0]*n
    edges=[(i,j) for i,j in combinations(range(n),2) if not fixed_path or j>i+1]
    if fixed_path:
        for i in range(n-1):out[i]|=1<<(i+1)
    for bit,(i,j) in enumerate(edges):
        if mask>>bit&1:out[i]|=1<<j
        else:out[j]|=1<<i
    return out


def ham(out):
    n=len(out);top=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for used in range(1,1<<n):
        for last,count in enumerate(dp[used]):
            if count:
                remaining=out[last]&~used
                while remaining:
                    bit=remaining&-remaining;remaining-=bit
                    dp[used|bit][bit.bit_length()-1]+=count
    return sum(dp[top])


def path_cube_counts(n, fixed_path=False):
    edges=[(i,j) for i,j in combinations(range(n),2) if not fixed_path or j>i+1]
    index={e:k for k,e in enumerate(edges)}
    top=(1<<len(edges))-1;counts=[0]*(top+1)
    for path in permutations(range(n)):
        used=value=0;valid=True
        for a,b in zip(path,path[1:]):
            if fixed_path and abs(a-b)==1:
                if a>b:valid=False;break
                continue
            k=index[tuple(sorted((a,b)))];used|=1<<k
            if a<b:value|=1<<k
        if not valid:continue
        free=top^used;sub=free
        while True:
            counts[value|sub]+=1
            if sub==0:break
            sub=(sub-1)&free
    return counts


def remark(out, path):
    if not all(out[a]>>b&1 for a,b in zip(path,path[1:])):return -1
    return sum(((out[path[a]]>>path[b])&1)<<k
               for k,(a,b) in enumerate((i,j) for i,j in combinations(range(len(out)),2) if j>i+1))


def canonical_labels(signatures):
    bank={};result=[]
    for key in signatures:
        if key not in bank:bank[key]=len(bank)
        result.append(bank[key])
    return result


def refine(counts, transitions):
    labels=canonical_labels(counts);sizes=[len(set(labels))]
    while True:
        nxt=canonical_labels((labels[x],*(labels[t[x]] if t[x]>=0 else -1 for t in transitions))
                             for x in range(len(counts)))
        sizes.append(len(set(nxt)))
        if nxt==labels:return labels,sizes
        labels=nxt


def inverse_table(counts, transitions, target):
    distance=[-1]*len(counts);choice=[None]*len(counts)
    todo=deque()
    for x,value in enumerate(counts):
        if value==target:distance[x]=0;todo.append(x)
    while todo:
        x=todo.popleft()
        for j,t in enumerate(transitions):
            y=t[x]
            if y>=0 and distance[y]<0:
                need(t[y]==x,'reverse BFS uses actual involutions')
                distance[y]=distance[x]+1;choice[y]=j;todo.append(y)
    for start,d in enumerate(distance):
        if d<0:continue
        x=start
        for _ in range(d):
            j=choice[x];y=transitions[j][x]
            need(y>=0 and distance[y]==distance[x]-1,'Boolean inverse descent')
            x=y
        need(counts[x]==target,'Boolean inverse reaches actual target')
    farthest=distance.index(max(distance));x=farthest;trace=[]
    while distance[x]>0:
        j=choice[x];y=transitions[j][x]
        trace.append({'command_index':j,'source':x,'target':y,'counts':[counts[x],counts[y]]});x=y
    return {'target':target,'reachable':sum(d>=0 for d in distance),
            'distance_histogram':dict(sorted(Counter(distance).items())),
            'maximum_distance':max(distance),'farthest_witness':trace}


def carrier():
    n=7;counts=path_cube_counts(n,True);top=len(counts)-1
    graphs=[rows(n,x,True) for x in range(top+1)]
    C=[x^top for x in range(top+1)]
    rem=[remark(out,R) for out in graphs]
    D={x for x,y in enumerate(rem) if y>=0}
    need(len(D)==4096 and all(rem[rem[x]]==x for x in D),'guarded involutive selector')
    need(all(rem[x]!=x for x in D),'nonidentity involution cannot fix an actual tournament source')
    need(not(D&{C[x] for x in D}),'complement reverses selector guard')
    keys=[]
    for x in range(top+1):
        if x in D:keys.append((1,counts[x],counts[C[x]],counts[C[rem[x]]]))
        elif C[x] in D:keys.append((2,counts[x],counts[C[x]],counts[C[rem[C[x]]]]))
        else:keys.append((0,counts[x],counts[C[x]]))
    formula=canonical_labels(keys)
    exact,rounds=refine(counts,[C,rem])
    need(formula==exact,'closed formula is the exact minimal partial observer')
    base=canonical_labels((counts[x],counts[C[x]]) for x in range(top+1))
    need(len(set(base))==3008,'inherited complement-only observer')
    need(rem[24511]==31 and [counts[x] for x in (24511,31,C[24511],C[31])]==[17,17,71,33],
         'literal same-object marking hostile')
    for x in (0,1,31,24511,top,140,314):
        need(ham(graphs[x])==counts[x],'independent endpoint DP controls')
    # Two marked blocks can now follow one another. Choose an explicitly fixed
    # second involution, independently of the observation partition.
    second=(0,3,4,1,2,5,6)
    rem2=[remark(out,second) for out in graphs]
    exact2,rounds2=refine(counts,[C,rem,rem2])
    need(all(y<0 or rem2[y]==x for x,y in enumerate(rem2)),'second partial involution')
    # Find a same first-observer fibre whose second-selector future separates it.
    groups={};hostile=None
    for x,cls in enumerate(exact):
        if cls not in groups:groups[cls]=x
        else:
            y=groups[cls]
            if exact2[x]!=exact2[y]:hostile=(y,x);break
    def distinguish(x,y,transitions):
        todo=deque([(x,y,())]);seen={(x,y)}
        while todo:
            a,b,word=todo.popleft()
            if counts[a]!=counts[b]:return word
            for j,t in enumerate(transitions):
                u,v=t[a],t[b]
                if (u<0)!=(v<0):return word+(j,)
                if u>=0 and (u,v) not in seen:
                    seen.add((u,v));todo.append((u,v,word+(j,)))
        return None
    need(hostile is not None,'second operation really refines')
    separating=distinguish(*hostile,[C,rem,rem2])
    need(separating is not None,'finite distinguishing word')
    operation_cube={};partitions={}
    for subset in range(8):
        commands=[t for j,t in enumerate((C,rem,rem2)) if subset>>j&1]
        p,sizes=refine(counts,commands)
        partitions[subset]=p
        operation_cube[str(subset)]={'commands':[('C','R','S')[j] for j in range(3) if subset>>j&1],
                                    'states':len(set(p)),'refinement_rounds':sizes}
    common=canonical_labels(zip(partitions[1],partitions[2]))
    need(len(set(common))<len(set(partitions[3])),'union of commands needs more than marginal observers')
    out={'n':7,'tiles':len(counts),'selector':R,'domain':len(D),'base_states':len(set(base)),
         'one_selector_states':len(set(exact)),'one_selector_rounds':rounds,
         'phase_state_counts':dict(Counter(key[0] for key in set(keys))),
         'phase_source_counts':dict(Counter(key[0] for key in keys)),
         'marking_hostile':{'tiles':[24511,31],'base_states':[keys[24511],keys[31]],'R':R},
         'second_selector':second,'second_domain':sum(y>=0 for y in rem2),
         'two_selector_states':len(set(exact2)),'two_selector_rounds':rounds2,
         'second_selector_hostile':{'tiles':hostile,'first_state':keys[hostile[0]],
                                    'distinguishing_word':[('C','R','S')[j] for j in separating]},
         'operation_cube':operation_cube,'common_C_and_R_partition':len(set(common)),
         'C_R_target_inverse':inverse_table(counts,[C,rem],17),
         'C_R_S_target_inverse':inverse_table(counts,[C,rem,rem2],17),
         'H_histogram':dict(sorted(Counter(counts).items()))}
    return out


def trap_packets():
    counts=path_cube_counts(6)
    minima=[x for x,h in enumerate(counts) if h!=17 and all(abs(counts[x^(1<<i)]-17)>abs(h-17) for i in range(15))]
    need(len(minima)==768,'inherited strict trap count')
    distances=Counter()
    for x in minima:
        two=any(counts[x^(1<<i)^(1<<j)]==17 for i,j in combinations(range(15),2))
        distances[2 if two else 3]+=1
    need(distances==Counter({2:480,3:288}),'all two-edge tests replayed')
    prior=json.loads((HERE.parent/'prime_inverse_ports_20261102/escape_audit.json').read_text())
    profiles=[]
    for record in prior['orbit_representatives']:
        source=record['source'];packet=record['barrier_walk'];bits=packet['changes']
        values=[counts[source^sum(1<<bits[i] for i in range(3) if u>>i&1)] for u in range(8)]
        energies=[(v-17)**2 for v in values]
        need(values[0]==15 and values[7]==17,'known legal packet endpoints')
        for u in range(8):
            x=source^sum(1<<bits[i] for i in range(3) if u>>i&1)
            graph=rows(6,x)
            need(ham(graph)==values[u],'all cube vertices independently counted')
            P=packet['shared_path'];need(all(graph[a]>>b&1 for a,b in zip(P,P[1:])),
                                        'whole cube preserves the common Hamilton path')
        coefficient=energies.copy()
        for i in range(3):
            for u in range(8):
                if u>>i&1:coefficient[u]-=coefficient[u^(1<<i)]
        even=sum(energies[u] for u in range(8) if u.bit_count()%2==0)
        odd=sum(energies[u] for u in range(8) if u.bit_count()%2==1)
        need(odd-even==coefficient[7],'third mixed difference')
        hits=[sum(values[u]==17 for u in range(8) if u.bit_count()%2==parity) for parity in (0,1)]
        need(hits==[0,1],'same pair marginals but different actual target support')
        for idx in combinations(range(3),2):
            marginals=[]
            for parity in (0,1):
                marginals.append(Counter(tuple(u>>i&1 for i in idx) for u in range(8) if u.bit_count()%2==parity))
            need(marginals[0]==marginals[1],'same complete pair marginals')
        profiles.append({'source':source,'orbit_size':record['orbit_size'],'changes':bits,
                         'shared_path':packet['shared_path'],'counts':values,'energies':energies,
                         'energy_mobius':coefficient,'even_energy_numerator_over4':even,
                         'odd_energy_numerator_over4':odd,'target_hit_numerators_over4':hits,
                         'independent_common_bias_energy_coefficients':[
                             sum(coefficient[u] for u in range(8) if u.bit_count()==degree)
                             for degree in range(4)]})
    values=profiles[0]['counts'];energies=profiles[0]['energies']
    law0=[2,0,0,1,0,1,1,0];law1=[1,1,1,0,1,0,0,1]
    need(sum(law0)==sum(law1)==5,'positive fifth-weight source laws')
    for idx in combinations(range(3),2):
        def marginal(law):
            out=Counter()
            for u,w in enumerate(law):out[tuple(u>>i&1 for i in idx)]+=w
            return out
        need(marginal(law0)==marginal(law1),'target-mean laws retain equal pair marginals')
    need([sum(w*h for w,h in zip(law,values)) for law in (law0,law1)]==[85,85],
         'both means equal exact target17')
    need([sum(w*e for w,e in zip(law,energies)) for law in (law0,law1)]==[144,136],
         'target-mean energy separation')
    return {'traps':768,'two_or_three_distance_histogram':dict(distances),'packet_cubes':profiles,
            'same_target_mean_laws':{'source':72,'denominator':5,'weights':[law0,law1],
                                    'mean':17,'energy_numerators':[144,136],
                                    'target_hit_numerators':[0,1]},
            'scope':'independent requested two-edge replay; three-edge packets inherited'}


if __name__=='__main__':
    result={'carrier':carrier(),'traps':trap_packets()}
    (HERE/'experiment.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'one_selector':result['carrier']['one_selector_states'],
                      'two_selectors':result['carrier']['two_selector_states'],
                      'rounds':result['carrier']['two_selector_rounds'],
                      'trap_cubic_differences':[r['energy_mobius'][7] for r in result['traps']['packet_cubes']]}))
