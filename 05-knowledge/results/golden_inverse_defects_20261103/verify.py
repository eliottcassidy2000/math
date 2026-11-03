"""Exact inverse, attachment obstructions, and transported offset families."""
from collections import Counter
from fractions import Fraction as F
from itertools import product,combinations
from pathlib import Path
import argparse,hashlib,json,math
import inverse as i

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
C=((-1,0,1),(1,-1,1),(1,1,-1))
D=((0,F(1,2),F(1,2)),(1,0,1),(1,F(1,2),F(1,2)))
SEEDS=(('a','a','r5'),('a','a','q1','a','a'),('a','a','a','a','a','q0'))
TAU=(2,-1,0)
FIRST=(15,17,19,23,25,27,29,33,37,41)
SECOND=(2,3,5,7,11,13,17,19,23,27)


def coefficients(P):
    Q=i.mm(i.mm(D,P),C)
    i.need(all(F(x).denominator==1 for row in Q for x in row),'integer defect matrix')
    return tuple(tuple(int(x) for x in row) for row in Q)


def matrix_certificates():
    Q={c:coefficients(P) for c,P in i.BANK.items()}
    identity=tuple(tuple(int(r==s) for s in range(3)) for r in range(3))
    for c,P in i.BANK.items():
        det,adj=i.INVERSES[c]
        i.need(det and i.mm(P,adj)==tuple(tuple(det*x for x in row) for row in identity),('adjugate',c))
        i.need(all(P[r][s]>=i.BANK['a'][r][s] for r in range(3) for s in range(3)),('singleton domination',c))
    # The complete 960-chart table is closed under reverse-converse.
    phases={};allQ={}
    lookup={tuple(a['source_rows']):c for c,a in i.ATLAS.items()}
    for c,a in i.ATLAS.items():
        P=a['matrix'];q=coefficients(P);allQ[c]=q
        i.need(min(i.defect(P[0]))>=0,('initial cone',c))
        if c!='a':
            i.need(min(x for row in q for x in row)>=1,('positive non-singleton defect action',c))
            i.need(P[0][1]+P[0][2]>=4,('non-singleton initial rank',c))
            i.need(all(P[r][s]>=i.BANK['b'][r][s] for r in range(3) for s in range(3)),('triangle domination',c))
        rows=a['source_rows'];n=len(rows)
        opposite=tuple(sum(1<<s for s in range(n) if s!=r and rows[n-1-s]>>(n-1-r)&1) for r in range(n))
        i.need(opposite in lookup,('reverse-converse closure',c))
        phases[c]=lookup[opposite]
    i.need(Q['a']==((0,0,1),(1,0,1),(0,1,1)),'defect shift')
    i.need(i.step(TAU,i.BANK['q2'])==TAU and i.step(TAU,i.BANK['b'])==tuple(-x for x in TAU),'sign transport')
    Pa2=i.mm(i.BANK['a'],i.BANK['a'])
    i.need(min(map(sum,Pa2))==3,'two-step length contraction')
    E=tuple(tuple(i.BANK['b'][r][s]-2*Pa2[r][s]+int(r==s) for s in range(3)) for r in range(3))
    i.need(min(x for row in E for x in row)>=0,'arbitrary-position insertion bound')
    mixed={}
    for c,q in Q.items():
        if c=='a':continue
        functional=tuple(row[0]-2*row[1]-row[2] for row in q)
        if min(functional)<=0<=max(functional):mixed[c]=functional
    i.need(mixed=={'r12':(9,0,-3),'r58':(-3,2,3)},'remaining mixed singleton guard')
    return {'defect_matrices':Q,'determinants':{c:d for c,(d,_) in i.INVERSES.items()},
            'golden_reverse_converse':{c:phases[c] for c in i.GOLDEN},
            'golden_tau_images':{c:i.step(TAU,i.BANK[c]) for c in i.GOLDEN},
            'full_atlas_checked':len(allQ),'insertion_difference':E,'mixed_singleton_guards':mixed}


def source_words(v):
    if v==i.ROOT:yield ();return
    for c,u in i.predecessors(v):
        if i.polynomial(u):
            for w in source_words(u):yield (c,)+w


def exact_fibres():
    states=i.forward_states(133)
    forward=i.forward_polynomials(states)
    for v,poly in forward.items():i.need(i.polynomial(v)==poly,('two routes to source fibre',v))
    candidates=0;nonmembers=0
    for h in range(1,22,2):
        for x in range(16*h+1):
            for y in range(max(x,h-x),min(x+h,16*h)+1):
                v=(h,x,y)
                if not i.valid(v):continue
                candidates+=1
                found=i.witness(v)
                i.need((found is not None)==(v in states),('complete cone universe',v))
                if found is None:nonmembers+=1
                else:i.need(i.state(found)==v,('decoded witness',v))
    # Independent literal Hamilton DP and word-by-word coefficient counting.
    actual={};literal_checks=0;endpoint_checks=0
    for n in range(1,10):
        for word in i.words(n):
            v=i.state(word)
            i.need(i.hamilton(i.literal(word))==v[0],('literal word',word))
            if n<=8:
                H,E,Cover=i.endpoint_stats(i.literal(word))
                i.need(i.defect(v)==(Cover-H,E,H-E),('physical defect coordinates',word))
                endpoint_checks+=1
            literal_checks+=1
            row=actual.setdefault(v,Counter());row[n]+=1
    for v,row in actual.items():
        exact=dict(i.polynomial(v))
        for n in range(1,10):i.need(exact.get(n,0)==row.get(n,0),('coefficient comparison',v,n))
    hist={}
    for h in sorted(set(FIRST+SECOND)):
        terms=Counter()
        for v,poly in forward.items():
            if v[0]==h:
                for n,count in poly:terms[n]+=count
        hist[h]=dict(sorted(terms.items()))
    return {'states_through133':len(states),'cone_candidates_through21':candidates,
            'cone_nonmembers':nonmembers,'literal_words_through9':literal_checks,'endpoint_checks':endpoint_checks,
            'count_fibre_polynomials':hist},states,forward


def attachment_checks(states):
    checks=0
    for v in states:
        if v[0]<3:continue
        hits=[]
        for offset in (2,4):
            predicted=i.ternary_left(v,offset)
            actual,_=i.extension_search(v,3*v[0]+offset)
            i.need((predicted is not None)==(actual is not None),('all-length left selector',v,offset))
            if predicted is not None:
                q=v
                for c in reversed(predicted):q=i.step(q,i.BANK[c])
                i.need(q[0]==3*v[0]+offset,('selector realization',v,offset))
                hits.append(offset)
            checks+=1
        i.need(len(hits)<2,('incompatible two offsets',v))
    low=i.low_defect_words()
    for v in set(low.values()):
        i.need(sum(dict(i.polynomial(v)).values())==sum(z==v for z in low.values()),('all low-defect words',v))
    h25=[]
    for v in sorted(states):
        if v[0]!=25:continue
        for w in source_words(v):
            table=[[i.state(('a',)*l+w+('a',)*r)[0] for r in range(3)] for l in range(3)]
            i.need(not any(x in (77,79) for row in table for x in row),('all H25 two-sided wrappers',w))
            h25.append({'word':w,'state':v,'table':table})
    i.need(len(h25)==11,'all physical H25 source words')
    q=[i.state(('a',)*n)[0] for n in range(15)]
    i.need(all(q[n+1]<3*q[n]+2<3*q[n]+4<q[n+2] for n in range(6,13)),'singleton insertion hostile seeds')
    return {'left_saturation_checks':checks,'low_defect_words':[{'word':w,'state':v,'d':i.defect(v)[0]} for w,v in sorted(low.items())],
            'low_defect_counts':dict(Counter(i.defect(v)[0] for v in low.values())),
            'low_defect_states':len(set(low.values())),'largest_low_defect_H':max(v[0] for v in low.values()),
            'H25_all_source_words':h25,'singleton_counts':q}


def offset_checks():
    vectors=[i.state(w) for w in SEEDS]
    for j in (1,2):i.need(vectors[j]==tuple(3*x+j*y for x,y in zip(vectors[0],TAU)),('seed relation',j))
    count=0;literal_count=0
    for length in range(9):
        for prefix in product(('b','q2'),repeat=length):
            sign=(-1)**prefix.count('b');h=i.state(prefix+SEEDS[0])[0]
            for j in (1,2):
                target=i.state(prefix+SEEDS[j])[0]
                i.need(target==3*h+2*j*sign,('all-phase algebra replay',prefix,j))
                count+=1
                if length<=1:
                    i.need(i.hamilton(i.literal(prefix+SEEDS[j]))==target,('literal offset',prefix,j))
                    literal_count+=1
    seq=[i.state(('q2',)*k+SEEDS[0])[0] for k in range(15)]
    for k in range(13):i.need(seq[k+2]==18*seq[k+1]-5*seq[k]-72,('quadratic clock',k))
    # Both output seeds fail the next signed-line ternary update.
    missing=[]
    for seed in vectors[1:]:
        for j in (1,2):
            v=tuple(3*x+j*y for x,y in zip(seed,TAU))
            i.need(not i.polynomial(v),('nonclosure of proposed affine carry',v))
            missing.append(v)
    return {'seed_words':SEEDS,'seed_vectors':vectors,'mixed_word_checks':count,
            'literal_offset_checks':literal_count,'q2_clock':seq,
            'q2_generating_numerator':[31,-172,69],'q2_generating_denominator':[1,-19,23,-5],
            'unrealizable_next_carries':missing}


def phase_clock_checks():
    # Integer coordinates (Y,X2,g)=(y,2x+h-y,h-y), with X2=g mod2.
    J=((0,1,1),(0,2,0),(1,-1,-1))
    q=((11,18,0),(4,7,0),(0,0,1));b=((5,8,0),(2,3,0),(0,0,-1))
    i.need(i.mm(i.BANK['q2'],J)==i.mm(J,q),'q2 split transport')
    i.need(i.mm(i.BANK['b'],J)==i.mm(J,b),'triangle split transport')
    orders={}
    for name,p,expected in [('b',17,68),('q2',19,171)]:
        P=i.BANK[name];M=((1,0,0),(0,1,0),(0,0,1));identity=M
        for exponent in range(1,expected+1):
            M=tuple(tuple(x%p for x in row) for row in i.mm(M,P))
            if M==identity:break
        i.need(exponent==expected and M==identity,('complete modular clock',name,p))
        orders[name]={'modulus':p,'order':exponent}
    Uq=((11,18),(4,7));Ub=((5,8),(2,3))
    even=[(0,0)]*81;odd=[(0,0)]*81;E=[0]*81;O=[0]*81
    even[0]=(25,34);E[0]=1
    def rowmul(v,P):return tuple(sum(v[j]*P[j][k] for j in range(2)) for k in range(2))
    def plus(a,b):return tuple(x+y for x,y in zip(a,b))
    for n in range(1,81):
        if n>=4:
            even[n]=plus(even[n],rowmul(even[n-4],Uq));odd[n]=plus(odd[n],rowmul(odd[n-4],Uq))
            E[n]+=E[n-4];O[n]+=O[n-4]
        if n>=3:
            even[n]=plus(even[n],rowmul(odd[n-3],Ub));odd[n]=plus(odd[n],rowmul(even[n-3],Ub))
            E[n]+=O[n-3];O[n]+=E[n-3]
    bruteE=Counter();bruteO=Counter();massE=Counter();massO=Counter();count=0
    stack=[((),0)]
    while stack:
        w,n=stack.pop();parity=w.count('b')%2
        (bruteO if parity else bruteE)[n]+=1
        (massO if parity else massE)[n]+=i.state(w+SEEDS[0])[0];count+=1
        for c,cost in [('b',3),('q2',4)]:
            if n+cost<=28:stack.append((w+(c,),n+cost))
    for n in range(29):
        i.need((E[n],O[n],even[n][0]+6*E[n],odd[n][0]-6*O[n])==
               (bruteE[n],bruteO[n],massE[n],massO[n]),('parity generating recurrence',n))
    for n in range(81):
        coefficient=E[n]-(2*E[n-4] if n>=4 else 0)-(E[n-6] if n>=6 else 0)+(E[n-8] if n>=8 else 0)
        i.need(coefficient==int(n==0)-int(n==4),('even-prefix rational function',n))
    return {'integer_change':J,'q2_matrix':q,'triangle_matrix':b,'modular_orders':orders,
            'prefix_enumerations':count,'even_prefix_counts_first29':E[:29],
            'even_mass_first29':[even[n][0]+6*E[n] for n in range(29)]}


def ear_recovery():
    rows=[0]*5
    for bit,(a,b) in enumerate(combinations(range(5),2)):
        if 8>>bit&1:rows[a]|=1<<b
        else:rows[b]|=1<<a
    i.need(i.hamilton(rows)==9,'inherited H9 parent')
    counts=Counter();by_mask={}
    for mask in range(1,31):
        child=tuple(row|((1<<5) if not mask>>v&1 else 0) for v,row in enumerate(rows))+(mask,)
        H=i.hamilton(child);counts[H]+=1;by_mask[mask]=H
    i.need(tuple(sorted(counts))==FIRST,'inherited exact ear image')
    return {'source_code':8,'source_rows':rows,'parent_H':9,'ear_multiplicities':dict(sorted(counts.items())),
            'half_defect_polynomial':{(h-9)//2:c for h,c in sorted(counts.items())},'cut_responses':by_mask}


def correction_carrier_checks():
    basis=(TAU,i.step(TAU,i.BANK['q4']),i.step(TAU,i.BANK['q5']))
    i.need(i.adjugate(basis)[0]==-2,'full correction lattice basis')
    i.need(all(P[1][0]%2==P[2][0]%2==0 for P in i.BANK.values()),'correction lattice closure')
    observers=[]
    for m in range(1,19):
        responses={(2*u%m,2*(u+x)%m,2*(u+x+y)%m)
                   for u,x,y in product(range(m),repeat=3)}
        expected=(m//math.gcd(m,2))**3
        i.need(len(responses)==expected,('correction observer quotient',m))
        observers.append((m,expected))
    return {'basis':basis,'determinant':-2,'modular_observer_sizes':observers}


def transported_lane_scan():
    bound=1000000;seen={};stack=[(i.state(SEEDS[0]),())];values={};even=0
    while stack:
        v,w=stack.pop()
        if v in seen:continue
        seen[v]=w;parity=w.count('b')%2
        i.need(v[0]-v[2]==6*(-1)**parity,('phase recovered from full state',w))
        if not parity:
            even+=1
            for j in (1,2):
                n=3*v[0]+2*j
                if n<=bound:
                    target=w+SEEDS[j]
                    i.need(i.state(target)[0]==n,('lane source word',n,target))
                    values[n]={'prefix':w,'seed':j,'source_word':target,'source_H':v[0]}
        for c in ('b','q2'):
            u=i.step(v,i.BANK[c])
            i.need(u[0]>=5*v[0],('uniform prefix growth',c,v))
            if 3*u[0]+2<=bound:stack.append((u,(c,)+w))
    def prime(n):return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))
    primes=[n for n in sorted(values) if prime(n)]
    seven=[n for n in sorted(values) if n%7==0 and prime(n//7)]
    i.need(primes==[97,4729,21839,21841,81197,81199,83597] and seven==[1253],'complete bounded prefix lanes')
    return {'output_bound':bound,'prefix_states':len(seen),'even_prefix_states':even,
            'distinct_outputs':len(values),'prime_outputs':primes,'seven_prime_outputs':seven,
            'all_output_witnesses':values}


def seven_field_clock():
    # F49=F7[theta], theta^2=4theta+2. Map (Y,X2) to4X2+Ytheta.
    def mul(a,b):
        return ((a[0]*b[0]+2*a[1]*b[1])%7,
                (a[0]*b[1]+a[1]*b[0]+4*a[1]*b[1])%7)
    def power(a,n):
        answer=(1,0)
        while n:
            if n&1:answer=mul(answer,a)
            a=mul(a,a);n//=2
        return answer
    theta=(0,1);alpha=(3,4)
    powers=[power(theta,k) for k in range(48)]
    i.need(len(set(powers))==48 and power(theta,48)==(1,0),'primitive F49 clock')
    i.need(power(theta,7)==(4,6) and power(theta,8)==(5,0),'F49 Frobenius and norm')
    i.need(powers[38]==alpha,'actual initial field phase')
    hits={2:[],4:[]};hist=Counter();vectors=[]
    v=i.state(SEEDS[0])
    for k in range(48):
        h,x,y=v;field=(4*(2*x+h-y)%7,y%7)
        i.need(field==mul(alpha,powers[k]),('literal count clock to field',k))
        i.need(power(field,7)==mul(alpha,powers[(7*k+36)%48]),('affine Frobenius phase',k))
        hist[h%7]+=1;vectors.append(tuple(t%7 for t in v))
        for offset in (2,4):
            if (3*h+offset)%7==0:hits[offset].append(k)
        v=i.step(v,i.BANK['q2'])
    i.need(len(set(vectors))==48 and tuple(t%7 for t in v)==vectors[0],'exact initial orbit period')
    i.need(hist==Counter({0:7,1:7,2:7,3:7,4:7,5:7,6:6}),'F49 linear fibres')
    i.need(hits=={2:[1,4,5,19,30,40,47],4:[6,16,23,25,28,29,43]},'all divisibility clocks')
    return {'period':48,'seed_field_element':alpha,'seed_exponent':38,'frobenius_shift':36,
            'H_residue_counts':dict(sorted(hist.items())),'seven_divisibility_positions':hits}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    matrices=matrix_certificates();fibres,states,forward=exact_fibres()
    attachments=attachment_checks(states);offsets=offset_checks();ear=ear_recovery();phase=phase_clock_checks()
    correction=correction_carrier_checks();lanes=transported_lane_scan();field=seven_field_clock()
    scan=i.forward_states(10000)
    both=[];one=0
    for h,x,y in scan:
        hit=[j for j in (1,2) if (3*h+2*j,3*x-j,3*y) in scan]
        if hit:one+=1
        if len(hit)==2:both.append((h,x,y))
    i.need(sorted(both)==[(31,14,25),(417,341,411),(1575,1222,1569),(3159,2679,3321)],'finite seed discovery replay')
    sources=[str((i.PREVIOUS/name).relative_to(ROOT)).replace('\\','/') for name in ('THEORY.md','atom_atlas.json')]
    sources += ['01-canon/theorems/THM-4118-ear-response-lattice-and-stateful-unit-component-intervals.md',
        '05-knowledge/results/clocks_carriers_genus_20261102/GENUS_BRIDGE.md',
        '05-knowledge/results/session151_spectrum_metagraph_20261103/strategy/TERNARY_COMPILER.md',
        '05-knowledge/results/golden_hamilton_transfer_20261103/RECURRENCE.md',
        '05-knowledge/results/golden_hamilton_transfer_20261103/INCOMING.md',
        '05-knowledge/results/session153_gadgets_20261103/cuts/README.md']
    receipt={'matrices':matrices,'fibres':fibres,'attachments':attachments,'offsets':offsets,'ear':ear,'phase':phase,
        'correction_carrier':correction,'transported_lanes':lanes,'seven_field_clock':field,
        'seed_scan':{'bound':10000,'states':len(scan),'one_or_more':one,'both':sorted(both)},
        'source_sha256_lf':{p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').encode()).hexdigest() for p in sources}}
    path=HERE/'receipt.json';text=json.dumps(receipt,sort_keys=True,separators=(',',':'))+'\n'
    if args.write:path.write_text(text,encoding='utf-8',newline='\n')
    else:i.need(json.loads(path.read_text(encoding='utf-8'))==json.loads(text),'receipt reproduction')
    print(json.dumps({'status':'PASS','fibres':{k:v for k,v in fibres.items() if k!='count_fibre_polynomials'},
        'left_checks':attachments['left_saturation_checks'],'low_defect_words':len(attachments['low_defect_words']),
        'H25_words':len(attachments['H25_all_source_words']),'offset_checks':offsets['mixed_word_checks'],
        'seed_scan':receipt['seed_scan']},sort_keys=True))


if __name__=='__main__':main()
