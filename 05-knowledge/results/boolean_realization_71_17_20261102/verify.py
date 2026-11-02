"""Exact Boolean realization, source coupling and marked-fibre certificates.

Standard library only. Checks remain active under python -O.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations, product
from pathlib import Path
import argparse
import importlib.util
import json
from scout import rows, h_count

HERE=Path(__file__).resolve().parent
GATES=0

def need(ok,why):
    global GATES
    GATES+=1
    if not ok: raise RuntimeError(why)

def paths(out):
    return [p for p in permutations(range(len(out)))
            if all(out[a]>>b&1 for a,b in zip(p,p[1:]))]

def chart(n):
    edges=list(combinations(range(n),2))
    free=[i for i,(a,b) in enumerate(edges) if b>a+1]
    fixed=sum(1<<i for i,(a,b) in enumerate(edges) if b==a+1)
    def expand(x):return fixed+sum((x>>i&1)<<e for i,e in enumerate(free))
    def compress(out,p):
        return sum(((out[p[a]]>>p[b])&1)<<j for j,e in enumerate(free) for a,b in [edges[e]])
    return free,expand,compress

def converse(out):
    allbits=(1<<len(out))-1
    return [allbits^r^(1<<v) for v,r in enumerate(out)]

def chart_key(out):
    _,_,compress=chart(len(out))
    ordinary=min(compress(out,p) for p in paths(out))
    opp=converse(out)
    return min(ordinary,min(compress(opp,p) for p in paths(opp)))

def full_canonical(out):
    n=len(out);edges=list(combinations(range(n),2));top=(1<<len(edges))-1
    values=[sum(((out[p[a]]>>p[b])&1)<<e for e,(a,b) in enumerate(edges))
            for p in permutations(range(n))]
    return min(min(values),min(top^v for v in values))

def is_strong(out):
    for graph in (out,converse(out)):
        seen=1
        while True:
            nxt=seen
            for i,r in enumerate(graph):
                if seen>>i&1:nxt|=r
            if nxt==seen:break
            seen=nxt
        if seen!=(1<<len(out))-1:return False
    return True

def expectation(values,probs):
    result=F(0)
    for mask,v in enumerate(values):
        weight=F(1)
        for e,p in enumerate(probs):weight*=p if mask>>e&1 else 1-p
        result+=weight*v
    return result

def polynomial_reduction():
    n=4;d=6;values=[h_count(rows(n,x)) for x in range(1<<d)]
    coeff=values.copy()
    for e in range(d):
        for x in range(1<<d):
            if x>>e&1:coeff[x]-=coeff[x^(1<<e)]
    coeff[0]-=3
    square=[0]*(1<<d)
    for a,ca in enumerate(coeff):
        for b,cb in enumerate(coeff):square[a|b]+=ca*cb
    examples=[]
    for probs in ([F(1,2)]*d,[F(i+1,d+2) for i in range(d)]):
        value=sum(c*product_value(probs,e) for e,c in enumerate(square))
        direct=expectation([(v-3)**2 for v in values],probs)
        mean=expectation(values,probs)
        need(value==direct,'Boolean OR-monomial square agrees with law')
        need(direct>(mean-3)**2,'naive square loses nonzero variance')
        examples.append({'mean':str(mean),'Boolean_energy':str(value),'naive_square':str((mean-3)**2)})
    return {'n':n,'coefficient_count':sum(c!=0 for c in square),'examples':examples}

def product_value(probs,mask):
    value=F(1)
    for e,p in enumerate(probs):
        if mask>>e&1:value*=p
    return value

def one_edge_checks(censuses):
    cases=[]
    for n,h,x,y,e in ((5,7,4,260,8),(6,21,20,28,3)):
        a,b=censuses[n][x],censuses[n][y]
        need((a,b)==(h-2,h+2) and x^y==1<<e,'sharp one-edge witness')
        px,py=paths(rows(n,x)),paths(rows(n,y))
        need((len(px),len(py))==(a,b),'independent permutation counts')
        common=sorted(set(px)&set(py))
        need(bool(common),'common Hamiltonian path preserves chart guard')
        free,expand,compress=chart(n)
        u,v=compress(rows(n,x),common[0]),compress(rows(n,y),common[0])
        need((u^v).bit_count()==1,'one-edge pair normalizes to one chart edge')
        need(h_count(rows(n,expand(u)))==a and h_count(rows(n,expand(v)))==b,'normalized counts')
        mean=F(a+b,2);second=F(a*a+b*b,2)
        need(mean==h and second-mean*mean==4,'sharp variance four')
        for p in (F(0),F(1,7),F(1,2),F(5,6),F(1)):
            need((1-p)*(a-h)**2+p*(b-h)**2==4,'energy identically four along gap edge')
        # A genuinely interior regular point, not only a fractional cube edge.
        d=n*(n-1)//2;eps=F(1,1024);bins=[[Counter(),Counter()] for _ in range(2)]
        for mask,value in enumerate(censuses[n]):
            bit=(mask>>e)&1;distance=((mask^x)&~(1<<e)).bit_count()
            bins[bit][0][distance]+=value;bins[bit][1][distance]+=value*value
        moments=[[sum(c*eps**k*(1-eps)**(d-1-k) for k,c in b.items()) for b in endpoint] for endpoint in bins]
        A,B=moments[0][0],moments[1][0]
        need(A<h<B,'full-interior interpolation brackets the forbidden target')
        p=(h-A)/(B-A);second=(1-p)*moments[0][1]+p*moments[1][1]
        need(0<p<1 and (1-p)*A+p*B==h and second-h*h>4,'interior mean is exact with positive realization defect')
        probs=[(1-eps if x>>j&1 else eps) for j in range(d)];probs[e]=p
        edge_indices={edge:j for j,edge in enumerate(combinations(range(n),2))}
        weighted_paths=F(0)
        for path in permutations(range(n)):
            weight=F(1)
            for u,v in zip(path,path[1:]):
                q=probs[edge_indices[min(u,v),max(u,v)]];weight*=q if u<v else 1-q
            weighted_paths+=weight
        need(weighted_paths==h,'interior mean independently verified by weighted path products')
        cases.append({'n':n,'target':h,'masks':[x,y],'edge':list(combinations(range(n),2))[e],
                      'counts':[a,b],'mean_at_half':h,'derivative':b-a,'energy':4,
                      'common_path':common[0],'fixed_path_tiles':[compress(rows(n,x),common[0]),compress(rows(n,y),common[0])],
                      'interior':{'epsilon':str(eps),'free_edge_probability':str(p),'mean':h,
                                  'derivative':str(B-A),'energy':str(second-h*h)}})
    return cases

def rounding(n,base,target):
    free,expand,_=chart(n)
    # Work on a three-coordinate face of the fixed-path chart, preserving its guard.
    coords=list(range(3))
    masks=[base^sum(((x>>j)&1)<<e for j,e in enumerate(coords)) for x in range(8)]
    counts=[h_count(rows(n,expand(x))) for x in masks]
    need(counts[0]==target,'rounding starts near an independently checked solution')
    maxerror=max((v-target)**2 for v in counts)
    eps=F(1,8*(maxerror+1))
    probs=[eps]*3
    energies=[(v-target)**2 for v in counts]
    initial=expectation(energies,probs)
    need(initial<4,'certified parity-gap threshold')
    trace=[]
    for e in range(3):
        before=expectation(energies,probs)
        options=[]
        for bit in (0,1):
            q=probs.copy();q[e]=F(bit)
            options.append(expectation(energies,q))
        selected=min(range(2),key=lambda b:options[b]);probs[e]=F(selected)
        need(options[selected]<=before,'conditional-expectation clamping is nonincreasing')
        trace.append({'coordinate':coords[e],'before':str(before),'options':list(map(str,options)),
                      'selected':selected})
    index=sum(int(p)<<e for e,p in enumerate(probs))
    need(counts[index]==target,'rounded Boolean tournament meets target')
    out=rows(n,expand(masks[index]))
    need(len(paths(out))==target,'rounding output independently counted')
    return {'n':n,'target':target,'base_tile':base,'face_counts':counts,'epsilon':str(eps),
            'initial_energy':str(initial),'trace':trace,'output_tile':masks[index],'out_rows':out}

def chart_checks():
    n=7;free,expand,compress=chart(n);top=(1<<len(free))-1
    counts=[h_count(rows(n,expand(x))) for x in range(top+1)]
    floor=min(h for x,h in enumerate(counts) if is_strong(rows(n,expand(x))))
    need(floor==25,'all fixed-path seven-vertex strong floor')
    need(7 not in counts and 21 not in counts,'seven-vertex omissions')
    hits=[x for x,h in enumerate(counts) if h==17 and counts[x^top]==71]
    need(len(hits)==12,'twelve literal 17/71 complement lines')
    fibres=defaultdict(list)
    for x,h in enumerate(counts):
        if h==17:fibres[chart_key(rows(n,expand(x)))].append(x)
    inventory=[];hostile=None
    for key,tiles in sorted(fibres.items()):
        partners=Counter(counts[x^top] for x in tiles)
        inventory.append({'merged_chart_key':key,'mass':len(tiles),'partner_histogram':dict(sorted(partners.items()))})
        xs=[x for x in tiles if counts[x^top]==71]
        ys=[x for x in tiles if counts[x^top]!=71]
        if hostile is None and xs and ys:
            x,y=xs[0],ys[0];out=rows(n,expand(x));opp=converse(out)
            transport=next(((p,False) for p in paths(out) if compress(out,p)==y),None)
            if transport is None:transport=next((p,True) for p in paths(opp) if compress(opp,p)==y)
            need(chart_key(rows(n,expand(x)))==chart_key(rows(n,expand(y))),'same actual merged node')
            hostile={'tiles':[x,y],'counts':[17,17],'complement_tiles':[x^top,y^top],
                     'complement_counts':[counts[x^top],counts[y^top]],
                     'transport_path':transport[0],'uses_converse':transport[1]}
    need(hostile is not None,'actual same-node complement non-descent')
    first=hits[0]
    for x in (first,first^top):need(len(paths(rows(n,expand(x))))==counts[x],'17/71 paths independent')
    return {'n':n,'universe':top+1,'strong_minimum':floor,'histogram':dict(sorted(Counter(counts).items())),
            'count_and_complement_observer_states':len(set((h,counts[x^top]) for x,h in enumerate(counts))),
            '17_to_71_lines':hits,'H17_partner_values':sorted(set(counts[x^top] for x,h in enumerate(counts) if h==17)),
            'H17_fibres':inventory,'same_node_non_descent':hostile,
            'witness':{'tiles':[first,first^top],'out_rows':[rows(n,expand(first)),rows(n,expand(first^top))]},
            'rounding':[rounding(n,next(x for x,h in enumerate(counts) if h==target),target) for target in (17,71)]}

def fake_fibre_check():
    free,expand,compress=chart(5)
    q=[full_canonical(rows(5,expand(x))) for x in range(64)]
    need((q[4],q[26],q[9])==(12,41,12),'inherited cube hostile coordinate convention')
    fake=q.copy();fake[4],fake[26]=fake[26],fake[4]
    need(Counter(q)==Counter(fake),'false partition has exactly correct masses')
    path=(0,4,1,2,3);out=rows(5,expand(4))
    need(all(out[a]>>b&1 for a,b in zip(path,path[1:])),'legal path transporter')
    need(compress(out,path)==9,'transport normalizes tiling four to nine')
    need(q[4]==q[9] and fake[4]!=fake[9],'Boolean partition fails actual path transport')
    return {'swapped_tiles':[4,26],'path':path,'transported_tile':9,
            'true_labels':[q[4],q[9]],'false_labels':[fake[4],fake[9]],
            'masses':dict(sorted(Counter(q).items()))}

def source_checks():
    weights=[1,-3,3];values=[1,3,5]
    need(sum(weights)==1 and sum(w*v for w,v in zip(weights,values))==7,'signed mean hostile')
    need(sum(w*(v-7)**2 for w,v in zip(weights,values))==0,'signed zero-energy hostile')
    # Actual ordered tournament joins: all left vertices dominate all right vertices.
    joins={}
    for a,b in product((4,260),repeat=2):
        A,B=rows(5,a),rows(5,b)
        out=[r|(((1<<5)-1)<<5) for r in A]+[r<<5 for r in B]
        value=h_count(out)
        need(value==h_count(A)*h_count(B),'actual ordered join multiplicativity')
        joins[(a,b)]=value
    independent=F(sum(joins.values()),4);shared=F(joins[(4,4)]+joins[(260,260)],2)
    need((independent,shared)==(49,53),'shared-source correlation is retained')
    clocks=[]
    for word in ((1,7),(7,1)):
        P=3**len(word);Q=2**sum(word);B=0;prefix=0
        for i,a in enumerate(word):B+=3**(len(word)-1-i)*2**prefix;prefix+=a
        c=((Q-B)*pow(P,-1,2*Q))%(2*Q);d=(P*c+B)//Q
        trajectory=[c];valuations=[]
        for a in word:
            v=3*trajectory[-1]+1;k=0
            while v%2==0:v//=2;k+=1
            valuations.append(k);trajectory.append(v)
        need(tuple(valuations)==word and trajectory[-1]==d,'exact ordered valuation clock')
        need(B==Q*d-P*c,'source-target affine identity')
        clocks.append({'valuation_word':word,'P':P,'Q':Q,'B':B,'source':c,'target':d,'trajectory':trajectory})
    need([x['source'] for x in clocks]==[483,469] and all(x['target']==17 for x in clocks),'recovered 17/71 source collision')
    return {'signed_hostile':{'weights':weights,'H_values':values,'mean':7,'energy':0},
            'join_values':{str(k):v for k,v in joins.items()},'independent_mean':str(independent),
            'shared_mean':str(shared),'clock_words':clocks}

def local_trap(counts,target):
    # Complete connected sublevel components at energy <=4; neutral steps are allowed.
    remaining={x for x,h in enumerate(counts) if (h-target)**2<=4}
    traps=[];components=0;dimension=(len(counts)-1).bit_length()
    while remaining:
        seed=min(remaining);remaining.remove(seed);component=[seed]
        for x in component:
            for e in range(dimension):
                y=x^(1<<e)
                if y in remaining:remaining.remove(y);component.append(y)
        components+=1
        if all(counts[x]!=target for x in component):traps.append(component)
    seed=72;h=counts[seed]
    neighbors=[counts[seed^(1<<e)] for e in range(dimension)]
    need(h==15 and len(paths(rows(6,seed)))==15,'strict local minimum count independently checked')
    need(all((v-target)**2>4 for v in neighbors),'every edge strictly increases target error')
    return {'target':target,'energy_cap':4,'components':components,'trapped_components':len(traps),
            'trap_sizes':dict(sorted(Counter(map(len,traps)).items())),
            'first_trap':traps[0] if traps else [],
            'strict_local_witness':{'mask':seed,'H':h,'out_rows':rows(6,seed),'neighbor_counts':neighbors},
            'conclusion':'finite sublevel test only; no global algorithm inferred'}

def hostile_81_27():
    path=HERE.parent/'regular_parent_symmetry_20261029'/'verify.py'
    spec=importlib.util.spec_from_file_location('inherited_hostile',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    result=module.hostile_checks()
    need([x['group_order'] for x in result['tournaments']]==[81,27],'retained hostile symmetry groups')
    return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE/'verification.json')
    args=parser.parse_args();censuses={};report=[]
    for n in range(3,7):
        counts=[h_count(rows(n,x)) for x in range(1<<(n*(n-1)//2))]
        censuses[n]=counts
        need(all(h%2==1 and h not in (7,21) for h in counts),'labelled parity and forbidden controls')
        if n<=5:
            for x,h in enumerate(counts):need(len(paths(rows(n,x)))==h,'full permutation-vs-DP census through five')
        report.append({'n':n,'universe':len(counts),'histogram':dict(sorted(Counter(counts).items()))})
    result={'labelled_censuses':report,'Boolean_reduction':polynomial_reduction(),
            'sharp_gap_witnesses':one_edge_checks(censuses),'chart_seven':chart_checks(),
            'false_Boolean_partition':fake_fibre_check(),'sources':source_checks(),
            'H17_local_energy_probe':local_trap(censuses[6],17),'hostile_81_27':hostile_81_27()}
    result['gates']=GATES
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'gates':GATES,'chart_lines_17_71':len(result['chart_seven']['17_to_71_lines']),
                      'same_node_non_descent':result['chart_seven']['same_node_non_descent'],
                      'H17_fibres':result['chart_seven']['H17_fibres'],
                      'local_energy_probe':result['H17_local_energy_probe']}))

if __name__=='__main__':main()
