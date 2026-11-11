"""Reuse the complete inherited class bank; construct exact multiplicity tables.

Default: save full orders1..8 table and all-order fixed-H formulas for H<=375.
Query: python3 -B build.py --n 100 --h 39
The query counts labelled tournaments on [n], not converse-merged classes.
"""
from collections import defaultdict,Counter
from fractions import Fraction
from functools import lru_cache
from math import comb,factorial
from pathlib import Path
import argparse,csv,hashlib,io,json,subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
SOURCE='05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/'
GLOBAL_SOURCE='05-knowledge/results/session173_modular_clocks_20261103/counts527/inventory_375.json'


def load():
    pins={}
    def read(name):
        data=subprocess.check_output(['git','show','HEAD:'+SOURCE+name],cwd=ROOT)
        pins[name]=hashlib.sha256(data).hexdigest()
        return data
    old=json.loads(read('extended_census.json'))
    records=[]
    for n in(1,2):records.append(dict(n=n,H=1,aut=1,multiplicity=1,strong=n==1))
    for x in old['records']:
        records.append(dict(n=x['n'],H=x['H'],aut=x['aut'],
                            multiplicity=x['ordinary_multiplicity'],strong=x['strong']))
    for x in csv.DictReader(io.StringIO(read('n8_classes.csv').decode())):
        records.append(dict(n=8,H=int(x['H']),aut=int(x['aut']),multiplicity=1,strong=bool(int(x['strong']))))
    return old,records,pins


def atom_bank(records):
    atoms=Counter()
    for x in records:
        if x['strong']:atoms[x['n'],x['H'],x['aut']]+=x['multiplicity']
    return atoms


def complete_atom_bank():
    raw=subprocess.check_output(['git','show','HEAD:'+GLOBAL_SOURCE],cwd=ROOT)
    inventory=json.loads(raw);atoms=Counter({(1,1,1):1})
    for n,rows in inventory['orders'].items():
        for code,h,self_converse,aut in rows:
            if h<=375:atoms[int(n),h,aut]+=1
    return atoms,hashlib.sha256(raw).hexdigest()


def core_compiler(atoms):
    by_H=defaultdict(list)
    for (n,h,aut),multiplicity in atoms.items():
        if h>1:by_H[h].append((n,aut,multiplicity))
    @lru_cache(None)
    def cores(h):
        if h==1:return {(0,0):(Fraction(1),Fraction(1))}
        result=defaultdict(lambda:[Fraction(0),Fraction(0)])
        for q,rows in by_H.items():
            if h%q:continue
            for (r,m),(ordinary,mass) in cores(h//q).items():
                for n,aut,count in rows:
                    result[r+1,m+n][0]+=count*ordinary
                    result[r+1,m+n][1]+=Fraction(count,aut)*mass
        return {key:tuple(value) for key,value in result.items()}
    return cores


def exact(n,h,cores):
    ordinary=Fraction(0);mass=Fraction(0)
    for (r,m),(a,b) in cores(h).items():
        if n>=m:
            placements=comb(n-m+r,r)
            ordinary+=a*placements;mass+=b*placements
    labelled=factorial(n)*mass
    if ordinary.denominator!=1 or labelled.denominator!=1:raise RuntimeError(('noninteger',n,h))
    return int(ordinary),int(labelled)


def through(n,h,cores):
    ordinary=Fraction(0)
    for (r,m),(a,b) in cores(h).items():
        if n>=m:ordinary+=a*comb(n-m+r+1,r+1)
    # The empty core would count the order-zero tournament as well.
    if h==1:ordinary-=1
    labelled=sum(exact(j,h,cores)[1] for j in range(1,n+1))
    if ordinary.denominator!=1:raise RuntimeError(('noninteger cumulative',n,h))
    return int(ordinary),labelled


def polynomial(terms,channel):
    result=[Fraction(0)]
    for (r,m),weights in terms.items():
        poly=[weights[channel]/factorial(r)]
        for j in range(r):
            new=[Fraction(0)]*(len(poly)+1)
            for k,c in enumerate(poly):new[k]+=c*(-m+r-j);new[k+1]+=c
            poly=new
        result.extend([Fraction(0)]*max(0,len(poly)-len(result)))
        for k,c in enumerate(poly):result[k]+=c
    while len(result)>1 and result[-1]==0:result.pop()
    return list(map(str,result))


def prime_columns(bins):
    source=GLOBAL_SOURCE
    raw=subprocess.check_output(['git','show','HEAD:'+source],cwd=ROOT)
    inventory=json.loads(raw)
    primes=[p for p in range(2,100) if all(p%d for d in range(2,int(p**.5)+1))]
    rows=[];all_terms={}
    for p in primes:
        count=Counter();mass=defaultdict(Fraction)
        for k,classes in inventory['orders'].items():
            for code,h,self_converse,aut in classes:
                if h==p:count[int(k)]+=1;mass[int(k)]+=Fraction(1,aut)
        for n in range(1,9):
            ordinary=sum(c*max(0,n-k+1) for k,c in count.items())
            labelled=factorial(n)*sum(c*max(0,n-k+1) for k,c in mass.items())
            if (ordinary,labelled)!=tuple(bins[n,p][:2]):raise RuntimeError(('prime inventory overlap',n,p))
        slope=sum(count.values());intercept=sum(c*(1-k) for k,c in count.items())
        mass_slope=sum(mass.values());mass_intercept=sum(c*(1-k) for k,c in mass.items())
        rows.append(dict(prime=p,attainable=bool(count),minimum_order=min(count,default=''),
                         global_strong_classes=slope,strong_by_order=';'.join(f'{k}:{c}' for k,c in sorted(count.items())),
                         ordinary_linear_slope=slope,ordinary_linear_intercept=intercept,
                         labelled_over_factorial_slope=str(mass_slope),labelled_over_factorial_intercept=str(mass_intercept),
                         linear_from_order=max((k-1 for k in count),default=1)))
        all_terms[p]=[dict(order=k,classes=count[k],automorphism_mass=str(mass[k])) for k in sorted(count)]
    with (HERE/'prime_columns_lt100.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    result=dict(status='FINITE-EXACT globally complete prime columns below100 from THM-4744 inventory',
                source=source,sha256=hashlib.sha256(raw).hexdigest(),prime_count=len(primes),
                attained_prime_count=sum(r['attainable'] for r in rows),all_attained_by_order=max(r['minimum_order'] for r in rows if r['attainable']),
                exact_formula='N_n(p)=sum_k classes[k]*max(0,n-k+1); L_n(p)/n!=sum_k automorphism_mass[k]*max(0,n-k+1)',
                terms=all_terms,ordinary_sum_slope=sum(r['ordinary_linear_slope'] for r in rows),
                ordinary_sum_intercept=sum(r['ordinary_linear_intercept'] for r in rows))
    (HERE/'prime_columns_lt100.json').write_text(json.dumps(result,indent=2)+'\n')
    return {k:v for k,v in result.items() if k!='terms'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--n',type=int);parser.add_argument('--h',type=int)
    args=parser.parse_args()
    old,records,pins=load();atoms=atom_bank(records)
    complete_atoms,complete_pin=complete_atom_bank();cores=core_compiler(complete_atoms)
    if args.n is not None or args.h is not None:
        if args.n is None or args.h is None or args.n<1 or not 1<=args.h<=375:
            parser.error('supply n>=1 and 1<=h<=375; this all-order bank is complete only in that H range')
        ordinary,labelled=exact(args.n,args.h,cores) if args.h&1 else(0,0)
        ordinary_through,labelled_through=through(args.n,args.h,cores) if args.h&1 else(0,0)
        print(json.dumps(dict(n=args.n,H=args.h,ordinary_isomorphism_classes=ordinary,
                              labelled_tournaments_on_n_vertices=labelled,
                              up_to_n_ordinary_classes=ordinary_through,
                              up_to_n_labelled_tournaments=labelled_through),indent=2));return
    bins=defaultdict(lambda:[0,0,0,0])
    for x in records:
        n,h,count,aut=x['n'],x['H'],x['multiplicity'],x['aut']
        lab=count*factorial(n)//aut
        if factorial(n)%aut:raise RuntimeError(('orbit denominator',x))
        b=bins[n,h];b[0]+=count;b[1]+=lab
        if x['strong']:b[2]+=count;b[3]+=lab
    # Exact first-SCC recurrence in both counting measures.
    unlab=[defaultdict(int) for _ in range(9)];lab=[defaultdict(int) for _ in range(9)]
    unlab[0][1]=lab[0][1]=1
    atom_u=Counter();atom_l=Counter()
    for (k,h,aut),count in atoms.items():
        atom_u[k,h]+=count;atom_l[k,h]+=count*factorial(k)//aut
    for n in range(1,9):
        for (k,q),count in atom_u.items():
            if k>n:continue
            for h,value in unlab[n-k].items():unlab[n][q*h]+=count*value
            for h,value in lab[n-k].items():lab[n][q*h]+=comb(n,k)*atom_l[k,q]*value
        if dict(unlab[n])!={h:v[0] for (nn,h),v in bins.items() if nn==n}:raise RuntimeError(('ordinary SCC',n))
        if dict(lab[n])!={h:v[1] for (nn,h),v in bins.items() if nn==n}:raise RuntimeError(('labelled SCC',n))
        if sum(lab[n].values())!=2**(n*(n-1)//2):raise RuntimeError(('labelled total',n))
        if sum(h*c for h,c in lab[n].items())!=factorial(n)*2**((n-1)*(n-2)//2):raise RuntimeError(('path moment',n))
        if n>=3:
            for h,b in old['counts_by_order'][str(n)].items():
                if bins[n,int(h)][:2]!=[b['ordinary'],b['labelled']]:raise RuntimeError(('source extraction',n,h))
    max_H=max(h for n,h in bins);cumulative=defaultdict(lambda:[0,0]);table=[]
    for n in range(1,9):
        for h in range(1,max_H+1,2):
            b=bins[n,h];cumulative[h][0]+=b[0];cumulative[h][1]+=b[1]
            table.append(dict(n=n,H=h,ordinary_classes=b[0],labelled_tournaments=b[1],
                              strong_ordinary_classes=b[2],strong_labelled_tournaments=b[3],
                              up_to_n_ordinary_classes=cumulative[h][0],up_to_n_labelled_tournaments=cumulative[h][1]))
    with (HERE/'spectrum_n1_n8.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
    saved_atoms=[dict(n=n,H=h,automorphism_order=aut,ordinary_classes=count)
                 for (n,h,aut),count in sorted(atoms.items())]
    (HERE/'strong_atoms_n1_n8.json').write_text(json.dumps(saved_atoms,indent=2)+'\n')
    low={}
    for h in range(1,376,2):
        terms=cores(h)
        for n in range(1,9):
            if exact(n,h,cores)!=tuple(bins[n,h][:2]):raise RuntimeError(('core words',n,h))
            direct=tuple(sum(bins[j,h][k] for j in range(1,n+1)) for k in (0,1))
            if through(n,h,cores)!=direct:raise RuntimeError(('cumulative core words',n,h))
        low[h]=dict(core_words=[dict(strong_components=r,total_size=m,ordinary_coefficient=str(a),automorphism_mass=str(b))
                               for (r,m),(a,b) in sorted(terms.items())],
                    ordinary_polynomial_ascending=polynomial(terms,0),labelled_over_factorial_polynomial_ascending=polynomial(terms,1),
                    exact_polynomial_from=max((m-r for r,m in terms),default=0),
                    degree=max((r for r,m in terms),default=None))
    (HERE/'all_order_H_le73.json').write_text(json.dumps({h:v for h,v in low.items() if h<=73},indent=2)+'\n')
    (HERE/'all_order_H_le375.json').write_text(json.dumps(low,indent=2)+'\n')
    summaries=[]
    for n in range(1,9):
        present=sorted(h for (nn,h),b in bins.items() if nn==n and b[0])
        maximum=max(present)
        missing=[h for h in range(1,maximum+1,2) if not bins[n,h][0]]
        modes_u=[h for h in present if bins[n,h][0]==max(bins[n,t][0] for t in present)]
        modes_l=[h for h in present if bins[n,h][1]==max(bins[n,t][1] for t in present)]
        summaries.append(dict(n=n,ordinary_classes=sum(bins[n,h][0] for h in present),
                              labelled_tournaments=sum(bins[n,h][1] for h in present),
                              strong_classes=sum(bins[n,h][2] for h in present),
                              distinct_H=len(present),maximum_H=maximum,missing_odd_up_to_max=missing,
                              ordinary_modes=modes_u,labelled_modes=modes_l))
    prime_report=prime_columns(bins)
    report=dict(status='FINITE-EXACT inherited complete orders1..8, with PROVED labelled weighted extension; all-n fixed-H formulas complete for H<=375',
                no_new_tournament_enumeration=True,source_sha256=pins,table_rows=len(table),summary=summaries,
                prime_columns=prime_report,
                complete_strong_source=GLOBAL_SOURCE,complete_strong_sha256=complete_pin,
                scope_H_completeness='THM-4744 inventory through H375 and order12; cited strong minimum at order13 is625, excluding all larger orders',
                sample_queries={str(h):dict(zip(('ordinary','labelled'),exact(20,h,cores))) for h in(3,9,15,39,63)},
                checks=['exact first-SCC convolution in ordinary and labelled channels for n1..8',
                        'labelled total2^(n choose2) and first Hamiltonian-path moment',
                        'all original per-order ordinary/labelled bins',
                        'all fixed-H<=375 core-word formulas against every order1..8 bin',
                        'cumulative core-word formula against all direct order sums through8',
                        'global prime inventory below100 agrees with every order1..8 ordinary and labelled bin'])
    (HERE/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(table_rows=len(table),summary=summaries,
                         sample_polynomials={h:low[h] for h in(9,15,39,63)}),indent=2))


if __name__=='__main__':main()
