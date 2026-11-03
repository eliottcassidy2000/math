"""Reproduce the seam theorem, finite atom atlas and bounded-value certificates.

No removable asserts. Source words are retained with their literal charts.
The optional full 594-matrix run is exploratory; the 65-symbol certificate
already gives the complete positive odd prefix through30000 except7,21.
"""
from itertools import permutations
from pathlib import Path
from math import gcd, isqrt
import argparse
import hashlib
import json
import random
import transfer as t
import extension as e
import search as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SELECTED=('s14','s85','s248','s261','s262','s268','s330','s453')
NUM=(1,0,-1,-6,-76,-72,1,10,-7)
DEN=(1,-1,-1,-8,-95,-87,-64,-157,-629,-202,-27,-227,-585)


def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p


def add(a,b,sign=1):
    c=[0]*max(len(a),len(b))
    for j,x in enumerate(a):c[j]+=x
    for j,x in enumerate(b):c[j]+=sign*x
    return trim(c)


def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)


def determinant(a):
    n=len(a);answer=0
    for p in permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        product=sign
        for i in range(n):product*=a[i][p[i]]
        answer+=product
    return answer


def bareiss(a):
    a=[list(row) for row in a];n=len(a);previous=1;sign=1
    for k in range(n-1):
        if a[k][k]==0:
            pivot=next((i for i in range(k+1,n) if a[i][k]),None)
            if pivot is None:return 0
            a[k],a[pivot]=a[pivot],a[k];sign=-sign
        value=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=value*a[i][j]-a[i][k]*a[k][j]
                t.need(numerator%previous==0,'fraction-free determinant')
                a[i][j]=numerator//previous
            a[i][k]=0
        previous=value
    return sign*a[-1][-1]


def operator(word):
    M=[[int(i==j) for j in range(3)] for i in range(3)]
    for letter in reversed(word):
        P=t.POSITIVE[letter]
        M=[[sum(M[i][k]*P[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    return M


def gcd_mod(a,b,p):
    a=trim([x%p for x in a]);b=trim([x%p for x in b])
    while b!=[0]:
        while len(a)>=len(b) and a!=[0]:
            d=len(a)-len(b);c=a[-1]*pow(b[-1],-1,p)%p
            a=trim([(a[j]-c*b[j-d] if d<=j<d+len(b) else a[j])%p for j in range(len(a))])
        a,b=b,a
    return [(x*pow(a[-1],-1,p))%p for x in a]


def algebra():
    t.verify_projection()
    h=[t.fast_vector(('a',)*n)[0] for n in range(5)]
    hankel=[[h[i+j] for j in range(3)] for i in range(3)]
    t.need(determinant(hankel)==-4,'rational minimum dimension three')
    observation=[[1,0,0],[1,2,0],[1,2,2]]
    for k,P in t.POSITIVE.items():
        t.need(P[0][0]%2==1 and P[1][0]%2==P[2][0]%2==0,('parity lattice',k))
    # Polynomial adjugate/determinant identity, independent of a series fit.
    R=[[[int(i==j)] for j in range(3)] for i in range(3)]
    for name,P in t.POSITIVE.items():
        for i in range(3):
            for j in range(3):R[i][j]=add(R[i][j],[0]*t.COST[name]+[P[i][j]],-1)
    D=[0]
    for p in permutations(range(3)):
        sign=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
        term=[1]
        for i in range(3):term=mul(term,R[i][p[i]])
        D=add(D,term,sign)
    N=add(mul(R[1][1],R[2][2]),mul(R[1][2],R[2][1]),-1)
    t.need(N==list(NUM) and D==list(DEN),'exact rational mass function')
    t.need(D[-1]%7 and gcd_mod(N,D,7)==[1],'coprime mass denominator, minimal order12')
    mass=[(1,0,0)]
    for n in range(1,101):
        v=[0,0,0]
        for name,P in t.POSITIVE.items():
            if t.COST[name]<=n:
                child=s.transition(mass[n-t.COST[name]],P)
                for j in range(3):v[j]+=child[j]
        mass.append(tuple(v))
    for n in range(101):
        coeff=sum(D[j]*mass[n-j][0] for j in range(min(n,len(D)-1)+1))
        t.need(coeff==(N[n] if n<len(N) else 0),('mass recurrence',n))
    periods={}
    for name,P in t.POSITIVE.items():
        tr=sum(P[i][i] for i in range(3))
        e2=sum(P[i][i]*P[j][j]-P[i][j]*P[j][i] for i in range(3) for j in range(i+1,3))
        det=determinant(P)
        seq=[t.fast_vector((name,)*k)[0] for k in range(9)]
        for k in range(6):t.need(seq[k+3]==tr*seq[k+2]-e2*seq[k+1]+det*seq[k],('period law',name))
        periods[name]={'characteristic':[1,-tr,e2,-det],'first_values':seq,'determinant':det}
    basis=[(),('a',),('b',),('q0',),('q1',),('q2',),('q4',),('a','b'),('a','q1')]
    full=[sum(operator(w),[]) for w in basis]
    contexts=[(i,j) for i in range(3) for j in range(3)]
    two_hankel=[[t.fast_vector(('a',)*i+w+('a',)*j)[0] for i,j in contexts] for w in basis]
    t.need(abs(bareiss(full))==128 and abs(bareiss(two_hankel))==8192,'two-sided rank nine')
    t.need(t.fast_vector(('b','a'))==t.fast_vector(('q5',))==(5,4,5),'one-sided twin')
    t.need(t.fast_vector(('b','a','a'))[0]==9 and t.fast_vector(('q5','a'))[0]==13,'right-attachment separator')
    return {'hankel':hankel,'hankel_det':-4,'observation':observation,
            'observation_smith':[1,2,2],'mass_numerator':N,'mass_denominator':D,
            'mass_first21':[v[0] for v in mass[:21]],'periods':periods,
            'two_sided':{'basis_words':basis,'matrix_span_det':bareiss(full),
                         'hankel_contexts':contexts,'hankel':two_hankel,'hankel_det':bareiss(two_hankel)}}


def cone_check(matrices,K):
    for name,P in matrices.items():
        t.need(all(x>=0 for row in P for x in row) and P[0][0]>=1,('positive monotonicity',name))
        t.need(P[0][0]%2==1 and P[1][0]%2==P[2][0]%2==0,('odd observer',name))
        for v in [(1,0,0),(1,K,0),(1,0,K),(1,K,K)]:
            h,x,y=s.transition(v,P)
            t.need(0<=x<=K*h and 0<=y<=K*h,('invariant cone ray',name,K,v))


def literal(word,atoms):
    rows=()
    for name in word:
        A=atoms[name];k,d=len(rows),len(A)
        rows=tuple(row|(((1<<d)-1)<<k) for row in rows)+tuple(row<<k for row in A)
    return t.complement(rows)


def atom_atlas():
    atlas={name:{'source_rows':rows,'matrix':t.POSITIVE[name]} for name,rows in t.ATOMS.items()}
    five=six=0
    for n in (5,6):
        for mask in range(1<<((n-1)*(n-2)//2)):
            rows=t.chart(n,mask)
            if not e.strong(rows):continue
            name=('r' if n==5 else 's')+str(mask)
            atlas[name]={'source_rows':rows,'matrix':e.compile_atom(rows)}
            if n==5:five+=1
            else:six+=1
    t.need((five,six)==(50,903),'complete marked strong atom inventory')
    return atlas


def direct_checks(atlas):
    atoms={k:tuple(v['source_rows']) for k,v in atlas.items()}
    mats={k:tuple(map(tuple,v['matrix'])) for k,v in atlas.items()}
    checked=0
    for n in range(1,11):
        for word in t.words(n):
            actual=t.subset_count(literal(word,atoms))
            t.need(t.count(word)==actual==t.fast_vector(word)[0],('three independent paths',word))
            checked+=1
    # All 960 atoms with each of the four two-sided singleton contexts.
    atom_cases=0
    for name in atlas:
        for left in ((),('a',)):
            for right in ((),('a',)):
                word=left+(name,)+right
                v=(1,0,0)
                for k in reversed(word):v=s.transition(v,mats[k])
                t.need(v[0]==t.subset_count(literal(word,atoms)),('literal atom context',word))
                atom_cases+=1
    # Literal permutations are an independent small-order count route.
    permutation_checks=0
    for n in range(1,7):
        for word in t.words(n):
            R=literal(word,atoms)
            count=sum(all(R[p[i]]>>p[i+1]&1 for i in range(n-1)) for p in permutations(range(n)))
            t.need(count==t.fast_vector(word)[0],('permutation count',word))
            permutation_checks+=1
    controls={}
    for word in [('a','q1'),('a','q4'),('a','q2'),('a','q5'),('r13','a'),('q0','r13')]:
        R=literal(word,atoms);h=t.subset_count(R)
        controls[' '.join(word)]={'H':h,'rows':R,'strong':e.strong(R)}
    t.need(controls['r13 a']['H']==37 and controls['q0 r13']['H']==133,'single-atom repairs')
    t.need(t.subset_count(atoms['r13'])==11,'uncomplemented r13 source count')
    patches={47:('a','s14'),83:('a','s268'),103:('a','s330'),
             491:('a','a','s248','a'),709:('a','s261','a','a','a'),
             833:('a','s85','b'),3227:('a','r35','s262'),4619:('a','s453','r29'),
             55:('a','q1','a','a'),111:('q2','a','a','b')}
    for expected,word in patches.items():
        R=literal(word,atoms)
        t.need(t.subset_count(R)==expected and e.strong(R),('displayed literal witness',expected,word))
        controls[' '.join(word)]={'H':expected,'rows':R,'strong':True}
    return {'golden_words_through10':checked,'all_atom_singleton_contexts':atom_cases,
            'literal_permutation_words_through6':permutation_checks,'controls':controls}


def independent_closure(bound,matrices,K,capture=None):
    # LIFO, no vertex costs or parents: independent of Dijkstra and its ties.
    seen={(1,0,0)};stack=[(1,0,0)]
    items=tuple(reversed(sorted(matrices.items())))
    while stack:
        h,x,y=stack.pop()
        for _,P in items:
            nh=P[0][0]*h+P[1][0]*x+P[2][0]*y
            if nh>bound:continue
            state=(nh,P[0][1]*h+P[1][1]*x+P[2][1]*y,P[0][2]*h+P[1][2]*x+P[2][2]*y)
            if state not in seen:seen.add(state);stack.append(state)
    digest=hashlib.sha256()
    if capture is not None:capture.extend(sorted(seen))
    for h,x,y in sorted(seen):digest.update(f'{h},{x},{y}\n'.encode('ascii'))
    return len(seen),digest.hexdigest()


def lanes(d):
    present={r['H'] for r in d['values']};bound=d['bound']
    p=s.primes(bound);p7=s.primes(bound//7)
    d['prime_lane']={'eligible':len([q for q in p if q!=7]),'missing':[q for q in p if q!=7 and q not in present]}
    d['seven_prime_lane']={'eligible':len([q for q in p7 if q!=3]),'missing_p':[q for q in p7 if q!=3 and 7*q not in present]}
    return d


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    al=algebra(); atlas=atom_atlas(); direct=direct_checks(atlas)
    cone_check({k:v['matrix'] for k,v in atlas.items()},16)
    five={k:v['matrix'] for k,v in atlas.items() if k.startswith('r')}
    fcost={k:5 for k in five}
    selected={k:atlas[k]['matrix'] for k in SELECTED}
    matrices=t.POSITIVE|five|selected
    costs=t.COST|fcost|{k:6 for k in selected}
    families=[('golden',t.POSITIVE,t.COST,4,100000),
              ('golden_plus_r13',t.POSITIVE|{'r13':five['r13']},t.COST|{'r13':5},8,100000),
              ('all_strong_through_five',t.POSITIVE|five,t.COST|fcost,8,10000),
              ('selected_six',matrices,costs,16,30000)]
    coverages={}
    summaries={}
    for name,M,cost,K,bound in families:
        cone_check(M,K)
        d=lanes(s.explore(bound,M,cost,K))
        size,digest=independent_closure(bound,M,K)
        t.need((size,digest)==(d['states'],d['state_sha256']),('independent closure',name))
        prime_targets=set(s.primes(bound))|{7*p for p in s.primes(bound//7) if p!=3}
        physical_atoms={k:tuple(atlas[k]['source_rows']) for k in M}
        strong_checks=0
        for item in d['values']:
            if item['H'] not in prime_targets:continue
            t.need(e.strong(literal(item['word'],physical_atoms)),('strong lane realization',name,item['H']))
            strong_checks+=1
        d['strong_lane_witness_checks']=strong_checks
        coverages[name]=d
        summaries[name]={k:v for k,v in d.items() if k not in ('values','missing_odds','prime_lane','seven_prime_lane')}
        summaries[name].update(attained=len(d['values']),
            prime_hits=d['prime_lane']['eligible']-len(d['prime_lane']['missing']),
            prime_total=d['prime_lane']['eligible'],
            seven_prime_hits=d['seven_prime_lane']['eligible']-len(d['seven_prime_lane']['missing_p']),
            seven_prime_total=d['seven_prime_lane']['eligible'],first_missing_odds=d['missing_odds'][:12])
    t.need(coverages['selected_six']['missing_odds']==[7,21],'complete65-symbol prefix')
    t.need(coverages['golden']['prime_lane']['missing'][0]==37,'all-size golden prime obstruction')
    t.need(coverages['golden']['seven_prime_lane']['missing_p'][0]==19,'all-size golden7p obstruction')
    # The small closed witnesses are separately convenient for hand inspection.
    exclusions={}
    for B in (37,133):
        d=s.explore(B);closed=[]
        independent_closure(B,t.POSITIVE,4,closed)
        closed_set=set(closed)
        for v in closed:
            for P in t.POSITIVE.values():
                w=s.transition(v,P)
                t.need(w[0]>B or w in closed_set,('finite closure certificate',B,v,w))
        d['closed_states']=closed
        exclusions[str(B)]=d
    receipt={'algebra':al,'direct':direct,'coverage_summaries':summaries,
             'selected_six_atoms':SELECTED,'atom_counts':{'five':50,'six':903},
             'distinct_atom_matrices':{'five':len(set(five.values())),
               'six':len({v['matrix'] for k,v in atlas.items() if k.startswith('s')})}}
    sources=[
        '05-knowledge/results/fibonacci_parity_precision_20261103/COMPLEMENTS.md',
        '05-knowledge/results/charged_recurrence_memory_20261028/COMPLEMENT.md',
        '05-knowledge/results/overnight_rank_return_20261029/DEFECT_TRANSFER.md',
        '05-knowledge/results/square_birth_windows_20261103/GROWTH.md',
        '05-knowledge/results/scc_basins_20261020d/arithmetic.md',
        '05-knowledge/results/paired_recursions_20261103/README.md',
        '05-knowledge/results/session150_hamiltonian_spectrum_20261103/transfer/TRANSFER.md',
        '05-knowledge/results/session150_hamiltonian_spectrum_20261103/prime_induction/RESULT.md',
        '05-knowledge/results/prime_lane_ports_five_groups_20261102/BOUNDARY_BUDGET.md',
        '01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md',
    ]
    receipt['source_sha256_lf']={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').encode()).hexdigest() for p in sources}
    outputs={'receipt.json':receipt,'atom_atlas.json':atlas,'coverage.json':coverages,'exclusions.json':exclusions}
    for name,obj in outputs.items():
        text=json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n'
        if args.write:(HERE/name).write_text(text,encoding='utf-8',newline='\n')
        else:t.need(json.loads((HERE/name).read_text(encoding='utf-8'))==json.loads(text),('saved receipt',name))
    print(json.dumps({'status':'PASS','direct':{k:v for k,v in direct.items() if k!='controls'},
                      'families':summaries},sort_keys=True))


if __name__=='__main__':main()
