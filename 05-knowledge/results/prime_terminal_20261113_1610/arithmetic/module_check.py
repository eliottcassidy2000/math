"""Complete small-core response fibres and padded-module clock controls."""
from collections import Counter
from itertools import combinations,product
from math import comb,lcm
from pathlib import Path
import json
import sys
import check as core

HERE=Path(__file__).resolve().parent

def subset_counts(rows):
    core.need(core.valid(rows),'literal tournament')
    n=len(rows);h=[1]+[0]*((1<<n)-1);dp=[[0]*n for _ in h]
    for mask in range(1,1<<n):
        for v in range(n):
            if mask>>v&1:
                rest=mask^(1<<v)
                dp[mask][v]=1 if not rest else sum(dp[rest][u] for u in range(n) if rows[u]>>v&1)
        h[mask]=sum(dp[mask])
    return h

def bins(rows,J=2):
    h=subset_counts(rows);previous=[1]+[0]*(len(h)-1);result=[0]
    for j in range(1,J+1):
        current=[0]*len(h)
        for mask in range(len(h)):
            sub=mask
            while True:
                current[mask]+=h[sub]*previous[mask^sub]
                if not sub:break
                sub=(sub-1)&mask
        result.append(current[-1]);previous=current
    return result

def run_to_bin(a):
    return [0]+[sum((-1)**(j-i)*comb(j,i)*a[j] for j in range(i,len(a))) for i in range(1,len(a))]

def main():
    profiles={11:Counter(),13:Counter()};representatives={};checked=permutations_checked=0;head_max={}
    for n in range(1,6):
        pairs=list(combinations(range(n),2));largest=0
        for bits in range(1<<len(pairs)):
            rows=[0]*n
            for index,(a,b) in enumerate(pairs):
                if bits>>index&1:rows[a]|=1<<b
                else:rows[b]|=1<<a
            h=subset_counts(rows);H=h[-1];largest=max(largest,H);checked+=1
            if H in profiles:
                P2=sum(h[s]*h[(len(h)-1)^s] for s in range(len(h)))
                core.need(core.strong(rows),'prime core in the bounded head')
                core.need(core.permutation_hp(rows)==H,'independent complete matching-head permutation count')
                permutations_checked+=1
                profiles[H][P2]+=1;representatives.setdefault(H,rows)
        head_max[n]=largest
    core.need(set(profiles[11])=={100} and set(profiles[13])=={104},'full small-core profile fibres')
    core.need(max(head_max[n] for n in range(1,5))==5 and 13**3<5**5,'certified strong-order cutoff inputs')
    modules={11:tuple(representatives[11]),13:tuple(representatives[13])}
    sources=[((2,20,9,35,45,7),0),((2,12,9,49,135,215,159,15),4),(core.apex('10010110'),4)]
    frame_a=((0,29,15),(0,89,45),(0,91,45))
    physical=[]
    for H,base in modules.items():
        base_bins=bins(base,4)
        for t in range(3):
            module=core.join(base,core.transitive(t)) if t else base
            P=bins(module,4)
            core.need(P[1:]==[base_bins[j]*j**t for j in range(1,5)],'ordinal bin clock')
            actual=[]
            for (frame,root),a in zip(sources,frame_a):
                rows=core.inflate(frame,root,module)
                value=core.hp(rows)
                c=run_to_bin(a)
                core.need(value==sum(c[j]*P[j] for j in range(1,len(c))),'physical run/bin expansion')
                actual.append(value)
            core.need(actual[1]-3*actual[0]==2*H and actual[2]-3*actual[0]==4*H,'scaled constant digits')
            if H==13:
                core.need(all(v%13==0 for v in actual[1:]),'all-size obstruction finite controls')
            physical.append({'module_H':H,'padding':t,'module_rows':module,'P1_to_P4':P[1:],'values':actual})
    positive=next(r for r in physical if r['module_H']==11 and r['padding']==1)
    core.need(positive['values']==[2989,8989,9011] and core.factor(9011)==[9011],'literal nonunit-module prime output')
    module=tuple(positive['module_rows'])
    h91=core.apex('10010110');exterior=core.join((2,4,1),h91)
    core.need(core.hp(exterior)==273,'actual smaller exterior')
    canonical=core.join(module,exterior)
    raw=core.inflate(*sources[2],module)
    order=list(range(4,4+len(module)))+[v for v in range(len(raw)) if v not in range(4,4+len(module))]
    target_core=tuple(sum(((raw[v]>>w)&1)<<j for j,w in enumerate(order)) for v in order)
    target=core.join(target_core,core.transitive(len(canonical)-len(target_core)))
    core.need(core.hp(module)*core.hp(exterior)==3003 and core.hp(target_core)==9011,'canonical unit-digit reframe counts')
    core.need(len(canonical)==len(target)==18 and all((canonical[i]&((1<<len(module))-1))==target[i]&((1<<len(module))-1)
              for i in range(len(module))),'same named module and vertex set')
    for graph in (canonical,target):
        core.need(core.valid(graph) and all(len({(graph[v]>>w)&1 for w in range(len(module))})==1
                                           for v in range(len(module),len(graph))),'module remains uniform')
    flips=[(i,j) for i in range(18) for j in range(i+1,18) if ((canonical[i]^target[i])>>j)&1]
    core.need(all(j>=len(module) for i,j in flips),'no old module arc changed')
    # A strict superset of physically possible bin coefficient differences:
    # all coefficients in[-2,2] at bases1,2,3,4. No physicality is inferred.
    formal_checks=0;max_hits=0
    P=bins(modules[11],4)
    for coefficients in product(range(-2,3),repeat=4):
        for digit in (2,4):
            hits=[t for t in range(21) if sum(coefficients[j-1]*P[j]*j**t for j in range(1,5))==digit]
            core.need(len(hits)<=3,'finite Chebyshev zero-count control')
            max_hits=max(max_hits,len(hits));formal_checks+=1
    isolated=(-18,1)
    isolated_hits=[t for t in range(21) if sum(isolated[j-1]*P[j]*j**t for j in (1,2))==2]
    core.need(isolated_hits==[1],'one-time unit digit can occur in the signed response lattice')
    modulus=lcm(P[2],P[3]//3,P[4])
    q=modulus*((2*pow(modulus,-1,11))%11)
    sharp=((2-q)//11,3*q//P[2],-3*q//P[3],q//P[4])
    core.need(all(denom and numerator%denom==0 for numerator,denom in
                  ((2-q,11),(3*q,P[2]),(-3*q,P[3]),(q,P[4]))),'integral sharp coefficients')
    sharp_hits=[t for t in range(21) if sum(sharp[j-1]*P[j]*j**t for j in range(1,5))==2]
    core.need(sharp_hits==[0,1,2],'sharp three-hit signed response example')
    # The unit-module boundary is positive, not covered by the nonunit no-go.
    unit=[]
    for t in range(3):
        module=core.transitive(t+1)
        values=[core.hp(core.inflate(frame,root,module)) for frame,root in sources]
        core.need(values[1]-3*values[0]==2 and values[2]-3*values[0]==4,'unit-module positive boundary')
        unit.append({'t':t,'values':values})
    report={'status':'PASS','labelled_tournaments_orders1_to5':checked,
            'head_H_maxima':head_max,'profiles':{str(h):dict(c) for h,c in profiles.items()},
            'representative_rows':modules,'matching_head_independent_permutation_checks':permutations_checked,
            'physical_padding_controls':physical,'unit_module_controls':unit,
            'canonical_9011_reframe':{'module_rows':module,'exterior_rows':exterior,'source_rows':canonical,'target_rows':target,
                                      'source_count':3003,'target_count':9011,'source_exterior_count':273,
                                      'same_order':18,'changed_arcs':flips,'scope':'Batch exterior rewire; no bounded uniform edit claim'},
            'formal_clock_universe':'coefficients[-2,2]^4, digits2/4, padding0..20; not assumed physical',
            'formal_clock_checks':formal_checks,'formal_maximum_hits':max_hits,
            'formal_positive_controls':{'isolated_coefficients':isolated,'isolated_hits':isolated_hits,
                                        'sharp_coefficients':sharp,'sharp_hits':sharp_hits,
                                        'scope':'Integer bin-response combinations, not asserted actual frame differences'},
            'cited_cutoff':{'paper':'Arthur H. Busch, A Note on the Number of Hamiltonian Paths in Strong Tournaments, EJC13(2006)N3',
                            'url':'https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/',
                            'statement':'H(C)^3>=5^(|C|-1) for strong tournaments C; 13^3<5^5'}}
    suffix='.optimized' if sys.flags.optimize else ''
    (HERE/('module_receipt'+suffix+'.json')).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k in ('status','labelled_tournaments_orders1_to5','profiles',
                                                          'matching_head_independent_permutation_checks','formal_clock_checks','formal_maximum_hits')}))


if __name__=='__main__':main()
