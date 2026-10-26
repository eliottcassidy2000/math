"""Mixed-precision Newton quotient of exact integral clone-excess moments.

The integer polynomial evaluation image, not the larger arbitrary constant-tail
module, is the carrier.  Sparse keys retain semantic support at zero residues.
"""
from importlib.util import module_from_spec,spec_from_file_location
from math import comb,gcd,prod
from pathlib import Path
from random import Random
import argparse
import json


def need(ok,message):
    if not ok: raise RuntimeError(message)


def triangular(j): return j*(j+1)//2


def precision_list(s,K):
    if not s: return [K]
    return [K-triangular(j) for j in range(s) if triangular(j)<K]


def zeros(s,K): return [0]*len(precision_list(s,K))


def add(out,term,s,K,factor=1):
    for j,p in enumerate(precision_list(s,K)):
        out[j]=(out[j]+factor*term[j])%(1<<p)
    return out


def multiply_x(c,s,K):
    need(s>0,'nonterminal polynomial multiplication')
    # X q_j = q_(j+1) + 2^(j+1) q_j.
    return [((c[j]<<(j+1))+(c[j-1] if j else 0))%(1<<p)
            for j,p in enumerate(precision_list(s,K))]


def monomial(exponent,s,K):
    need(exponent>=0 and s>0,'nonnegative exponent and nonterminal horizon')
    out=zeros(s,K); out[0]=1
    for _ in range(exponent): out=multiply_x(out,s,K)
    return out


def evaluate(c,s,a,K):
    need(1<=a<=s,'available positive source evaluation')
    modulus=1<<K; point=pow(2,a,modulus); q=1; total=0
    for j,cj in enumerate(c):
        total=(total+cj*q)%modulus
        q=q*(point-pow(2,j+1,modulus))%modulus
    return total


def transform(c,s,rest,a,delta,j0,factor,K):
    """P'(X)=factor*2^j0*X^delta*P(2^a X), at destination nodes."""
    need(delta>=0 and j0>=0 and 0<=a<=s,'integral transform data')
    need(a+rest<=s,'decreasing continuation horizon')
    modulus=1<<K
    if not s:
        need(rest==a==0,'already terminal identity')
        return [(factor*c[0]*(1<<j0))%modulus]
    if not rest:
        need(a>0,'a positive completing batch reads a positive source node')
        return [(factor*evaluate(c,s,a,K)*(1<<j0))%modulus]
    out=zeros(rest,K); basis=zeros(rest,K); basis[0]=1
    for j,cj in enumerate(c):
        add(out,basis,rest,K,cj)
        # q_(j+1)(2^a X) = (2^a X - 2^(j+1))*q_j(2^a X).
        xbasis=multiply_x(basis,rest,K)
        basis=[((x<<a)-(y<<(j+1)))%(1<<p)
               for x,y,p in zip(xbasis,basis,precision_list(rest,K))]
    for _ in range(delta): out=multiply_x(out,rest,K)
    return [(factor*x*(1<<j0))%(1<<p)
            for x,p in zip(out,precision_list(rest,K))]


def bound(D):
    return max(((x+2*i-1)//(2*i) for i,x in enumerate(D) if i),default=0)


def sparse_to_newton(s,D,bank,K):
    if not s: return [sum(bank.values())%(1<<K)]
    b=bound(D); out=zeros(s,K)
    for m,x in bank.items():
        need(m>=b,'integral polynomial exponent')
        add(out,monomial(m-b,s,K),s,K,x)
    return out


def odd_factorials(n,mod):
    out=[1]
    for j in range(1,n+1): out.append(out[-1]*(j//(j&-j))%mod)
    return out


def coefficient(odd,u,ell,k,e,mod):
    beta=(-1)**e*comb(k+e-1,e) if k else int(e==0)
    denominator=odd[u]*pow(ell,k,mod)*odd[k]%mod
    return beta*odd[u+ell*(k+e)]*pow(denominator,-1,mod)%mod


def next_defect(ell,rest,D,k):
    return tuple(D[(d-1)//2]+k*(d-gcd(ell,d))
                 for d in range(1,min(ell-2,rest)+1,2))


def audit_stage(n,ell,states,reference,K):
    banks={}; mod=1<<K
    # Independent original f-response weights, not Newton or moment transfer.
    from math import factorial
    for (s,response),weight in reference.items():
        if not s:
            x=weight*(factorial(n)>>(n-n.bit_count()))
            need(x.denominator==1,'terminal reference integer')
            if int(x)%mod: banks[(0,())]={0:int(x)%mod}
            continue
        m=response[0]
        D=tuple(d*m-f for d,f in zip(range(1,min(ell,s)+1,2),response))
        u=n-s; x=weight*(factorial(u)>>(u-u.bit_count()))
        need(x.denominator==1,'source reference integer')
        banks.setdefault((s,D),{})[m]=int(x)%mod
    need(set(banks)==set(states),'exact reference structural key equality')
    for key,(kind,data) in states.items():
        s,D=key
        want=banks[key] if kind=='sparse' else sparse_to_newton(s,D,banks[key],K)
        need(data==want,('exact stagewise Newton bank',n,ell,key))
    return len(states)


def newton_target(n,K,reference=None):
    if type(n) is not int or n<0 or type(K) is not int or K<1:
        raise ValueError('integer n>=0 and integer K>=1 required')
    if not n:
        return 1,{'keys':1,'cells':1,'residue_bits':K,'rows':[],
                  'conversions':0,'zero_sparse_contributions':0,'audited_states':0}
    mod=1<<K; odd=odd_factorials(n,mod); top=n if n%2 else n-1
    states={(n,(0,)*((top+1)//2)):('sparse',{0:1})}; rows=[]
    stats={'conversions':0,'zero_sparse_contributions':0,'audited_states':0}
    for ell in range(top,0,-2):
        if reference is not None:
            stats['audited_states']+=audit_stage(n,ell,states,reference[ell],K)
        rows.append({'ell':ell,'keys':len(states),
                     'cells':sum(len(data) for kind,data in states.values()),
                     'residue_bits':sum(len(data)*K if kind=='sparse' else sum(precision_list(key[0],K))
                                        for key,(kind,data) in states.items())})
        nxt={}
        def add_dense(key,term):
            rest,Dp=key
            if key not in nxt: nxt[key]=('dense',zeros(rest,K))
            kind,data=nxt[key]
            if kind=='sparse':
                data=sparse_to_newton(rest,Dp,data,K); nxt[key]=('dense',data)
                stats['conversions']+=1
            add(data,term,rest,K)
        def add_sparse(key,m,x):
            rest,Dp=key
            if rest and not x: stats['zero_sparse_contributions']+=1
            if key not in nxt: nxt[key]=('sparse',{}) if rest else ('dense',[0])
            kind,data=nxt[key]
            if kind=='sparse':
                # Semantic support survives every modular zero.
                data[m]=(data.get(m,0)+x)%mod
                if len(data)>len(precision_list(rest,K)):
                    nxt[key]=('dense',sparse_to_newton(rest,Dp,data,K))
                    stats['conversions']+=1
            else:
                add(data,sparse_to_newton(rest,Dp,{m:x},K),rest,K)
        for (s,D),(kind,data) in states.items():
            selected=D[(ell-1)//2] if ell<=s else 0
            for k in range(s//ell+1):
                for e in (range((s-ell*k)//ell+1) if k else (0,)):
                    rest=s-ell*(k+e)
                    if ell==1 and rest: continue
                    Dp=next_defect(ell,rest,D,k); key=(rest,Dp)
                    factor=coefficient(odd,n-s,ell,k,e,mod)
                    base=k*(ell-1)//2+ell*k*(k-1)//2-(k-k.bit_count())
                    if kind=='sparse':
                        for m,x in data.items():
                            shift=base+k*(ell*m-selected)
                            need(shift>=0,'sparse nonnegative shift')
                            add_sparse(key,m+k,0 if shift>=K else (factor*x<<shift)%mod)
                    else:
                        b,bp=bound(D),bound(Dp)
                        term=transform(data,s,rest,ell*k,k+b-bp,
                                       base+k*(ell*b-selected),factor,K)
                        add_dense(key,term)
        # No semantic zero is dropped before the terminal boundary.
        states={key:entry for key,entry in nxt.items() if key[0] or entry[1][0]}
    value=states.get((0,()),('dense',[0]))[1][0]
    stats.update({'keys':sum(r['keys'] for r in rows),'cells':sum(r['cells'] for r in rows),
                  'residue_bits':sum(r['residue_bits'] for r in rows),'rows':rows})
    return value*pow(odd[n],-1,mod)%mod,stats


def load(path,name):
    spec=spec_from_file_location(name,path); module=module_from_spec(spec);spec.loader.exec_module(module)
    return module


def controls(repo):
    root=repo/'05-knowledge/results/port_residual_integral_20261025/counting'
    exact=load(root/'check.py','exact_reference'); old=load(root/'modular.py','old_modular')
    comparisons=0; audited_states=0; samples=[]; zero_seen=False
    literal=load(repo/'05-knowledge/results/ports_regular_actions_hn_20261025/counting/check_excess.py','literal_core')
    literal_counts=[]
    for n in range(45):
        count,_=exact.integral_target(n,True)
        reference,reference_count=exact.rational_stages(n,True)
        need(reference_count==count,'independent f-response exact count')
        if n<=5:
            census=literal.literal_pairfree_classes(n)
            need(census==count,'literal pair-free class count')
            literal_counts.append(census)
        for K in range(1,10):
            got,stats=newton_target(n,K,reference); audited_states+=stats['audited_states']; prior,before=old.modular_target(n,K,True)
            need(got==prior==count%(1<<K),('core residue',n,K))
            if n:
                need(all(a['keys']==b['keys'] and a['cells']<=b['cells']
                         for a,b in zip(stats['rows'],before['rows'])),('stagewise coordinate domination',n,K))
            zero_seen|=bool(stats['zero_sparse_contributions']); comparisons+=1
            if n in (18,26,40,44) and K in (2,4,9):
                samples.append({'n':n,'K':K,'old_cells':before['cells'],'new_cells':stats['cells'],
                                'old_allocated_residue_bits':K*before['cells'],
                                'new_allocated_residue_bits':stats['residue_bits']})
    need(zero_seen,'zero sparse residues retained in actual runs')
    rng=Random(20261026); squares=0; kernel_cases=0
    for K in range(1,11):
        for s in range(1,16):
            # Any integer polynomial, including signed coefficients.
            bank={j:rng.randrange(-100,101) for j in range(9)}
            c=sparse_to_newton(s,(0,),bank,K)
            for r in range(1,s+1):
                need(evaluate(c,s,r,K)==sum(x*pow(2,r*m,1<<K) for m,x in bank.items())%(1<<K),
                     ('all-horizon polynomial evaluation',s,K,r))
            for rest in range(s+1):
                for a in range(s-rest+1):
                    if rest==0 and a==0: continue
                    delta=rng.randrange(5);j0=rng.randrange(4);factor=rng.randrange(-10,11)
                    out=transform(c,s,rest,a,delta,j0,factor,K)
                    dest={m+delta:factor*x*(1<<(j0+a*m)) for m,x in bank.items()}
                    target=sparse_to_newton(rest,(0,),dest,K)
                    need(out==target,('exact mixed transfer square',s,K,rest,a)); squares+=1
            # Every omitted/mixed coefficient modulus evaluates to zero.
            for j in range(min(s,12)):
                p=max(0,K-triangular(j))
                for r in range(1,s+1):
                    q=prod((1<<r)-(1<<i) for i in range(1,j+1))
                    need(((1<<p)*q)%(1<<K)==0,'Newton kernel generator')
                kernel_cases+=1
    # Kernel converse/minimality: exhaustive small mixed-coordinate representatives.
    from itertools import product
    unique_cases=0
    for K in range(1,5):
        for s in range(1,5):
            precisions=precision_list(s,K); outputs=set()
            for c in product(*(range(1<<p) for p in precisions)):
                image=tuple(evaluate(c,s,r,K) for r in range(1,s+1))
                need(image not in outputs,'distinct mixed residues have distinct observations')
                outputs.add(image)
            need(len(outputs)==1<<sum(precisions),'exact quotient cardinality');unique_cases+=1
    # Bounded DEGREE h is distinct from h arbitrary sparse exponent slots.
    bounded_degree_cases=0
    for K in range(1,4):
        for h in range(1,5):
            for s in range(1,5):
                images=set()
                for coeffs in product(range(1<<K),repeat=h):
                    images.add(tuple(sum(c*pow(2,r*j,1<<K) for j,c in enumerate(coeffs))%(1<<K)
                                     for r in range(1,s+1)))
                p=[K-triangular(j) for j in range(min(h,s)) if triangular(j)<K]
                need(len(images)==1<<sum(p),'bounded-degree exact image invariant factors')
                bounded_degree_cases+=1
    # Sparse h=1 is NOT always module rank 1: the only exponent can be too large.
    for K in range(1,7):
        c=sparse_to_newton(8,(0,),{K:1},K)
        need(all(evaluate(c,8,r,K)==0 for r in range(1,9)),
             'arbitrary sparse exponent gaps can lower rank further')
    # Division by 2 or X is outside the kernel-preserving operation class.
    need(sparse_to_newton(3,(0,),{0:4},2)==[0,0] and
         evaluate(sparse_to_newton(3,(0,),{0:2},2),3,1,2)==2,
         'exact coefficient division by 2 does not descend')
    need(sparse_to_newton(3,(0,),{2:1},2)==[0,0] and
         evaluate(sparse_to_newton(3,(0,),{1:1},2),3,1,2)==2,
         'negative gauge exponent does not descend')
    # Full constant-tail module is strictly larger than actual polynomial image.
    # K=2, Y_1=0,Y_2=1 cannot be a polynomial bank: P(4)-P(2) is even.
    need(all((evaluate(c,2,2,2)-evaluate(c,2,1,2))%2==0
             for c in product(range(4),range(2))),'invalid constant-tail input excluded')
    for n,K in [(-1,2),(3,0),(True,2),(3,1.5)]:
        try:newton_target(n,K)
        except ValueError:pass
        else:raise RuntimeError('invalid modular inputs accepted')
    return {'status':'PASS','core_orders':[0,44],'absolute_precision_bits':[1,9],
            'complete_count_comparisons':comparisons,'exact_intermediate_state_checks':audited_states,
            'literal_pairfree_class_counts_orders_0_to_5':literal_counts,'arbitrary_integer_polynomial_transfer_squares':squares,
            'kernel_generator_checks':kernel_cases,'exhaustive_quotient_cardinality_cases':unique_cases,
            'bounded_degree_image_cardinality_cases':bounded_degree_cases,
            'samples':samples,'precision_example_K64':precision_list(100,64),
            'scope':'Exact clone-excess absolute residues; image-lattice Newton compression, not an ordinary scalar-state reduction or normalized-precision theorem.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path)
    args=parser.parse_args()
    if args.repo is None:
        for candidate in [Path(__file__).resolve().parent.parent/'session64-repo',*Path(__file__).resolve().parents]:
            if (candidate/'05-knowledge/results/port_residual_integral_20261025/counting/check.py').is_file():
                args.repo=candidate;break
    if args.repo is None:raise FileNotFoundError('supply --repo')
    print(json.dumps(controls(args.repo),indent=2,sort_keys=True))
