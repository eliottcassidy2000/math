"""Independent native Newton audit; explicit gates, no producer controls/audit reuse.

Imports the production target only. Fresh f-response, unmerged odd partitions,
exact interpolation, and direct all-node observations are written here.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import factorial, gcd, prod
from pathlib import Path
from random import Random
import argparse, json, sys
sys.dont_write_bytecode = True


def need(ok, msg):
    if not ok: raise RuntimeError(msg)


def load(path, name):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec); spec.loader.exec_module(module)
    return module


def oddpart(n):
    while n and not n % 2: n //= 2
    return n


def partitions(total, lower=1):
    if not total:
        yield ()
        return
    for length in range(lower, total+1, 2):
        for tail in partitions(total-length, length):
            yield (length,) + tail


def partition_series(N):
    """Unmerged Burnside sum, repeated geometric factors (no beta/binomial)."""
    ordinary = [Fraction(0) for _ in range(N+1)]
    core = [Fraction(0) for _ in range(N+1)]
    cases = 0
    for actual in range(N+1):
        for part in partitions(actual):
            mult = Counter(part)
            pair_orbits = sum((d-1)//2 for d in part)
            pair_orbits += sum(gcd(part[i], part[j]) for i in range(len(part)) for j in range(i))
            centralizer = prod(d**k*factorial(k) for d,k in mult.items())
            weight = Fraction(2**pair_orbits, centralizer)
            ordinary[actual] += weight
            # One factor (1+z^d)^(-1) per actual cycle.
            excess = [1] + [0]*(N-actual)
            for d in part:
                nxt = [0]*len(excess)
                for old, value in enumerate(excess):
                    if not value: continue
                    for e in range((N-actual-old)//d+1):
                        nxt[old+d*e] += value*(-1)**e
                excess = nxt
            for e,c in enumerate(excess): core[actual+e] += weight*c
            cases += 1
    need(all(x.denominator==1 for x in ordinary+core),'integer unmerged partition totals')
    return [int(x) for x in ordinary], [int(x) for x in core], cases


def f_stages(n, core):
    """Original gcd-response carrier, no defects, gauge, moments or Newton."""
    if not n: return {}, 1
    ell = n if n%2 else n-1
    states = {(n, (0,)*((ell+1)//2)): Fraction(1)}
    snapshots = {}
    for d in range(ell,0,-2):
        snapshots[d] = states
        nxt = defaultdict(Fraction)
        for (s,f), weight in states.items():
            selected = f[(d-1)//2] if d<=s else 0
            for k in range(s//d+1):
                e_values = range((s-d*k)//d+1) if core and k else (0,)
                # Expand k geometric factors directly to avoid production beta.
                coefficients = [1]+[0]*((s-d*k)//d)
                for _ in range(k):
                    coefficients = [sum(coefficients[j]*(-1)**(i-j) for j in range(i+1))
                                    for i in range(len(coefficients))]
                exponent = k*(d-1)//2+d*k*(k-1)//2+k*selected
                for e in e_values:
                    rest = s-d*(k+e)
                    if d==1 and rest: continue
                    fn = tuple(f[(t-1)//2]+k*gcd(t,d)
                               for t in range(1,min(d-2,rest)+1,2))
                    nxt[(rest,fn)] += weight*Fraction(coefficients[e]*2**exponent, d**k*factorial(k))
        states = {key:x for key,x in nxt.items() if x}
    final = states.get((0,()),Fraction(0))
    need(final.denominator==1,('integer fresh f count',n,core))
    return snapshots,int(final)


def b_of(D):
    return max([0]+[(D[i]+2*i-1)//(2*i) for i in range(1,len(D))])


def moment(bank,b,r): return sum(x*2**(r*(m-b)) for m,x in bank.items())


def exact_newton(bank,b,s,K):
    """Integer divided differences at 2,4,...; not multiply_x conversion."""
    coefficients=[]
    for j in range(s):
        tau=j*(j+1)//2
        if tau>=K: break
        node=2**(j+1)
        numerator=moment(bank,b,j+1)
        for i,c in enumerate(coefficients):
            numerator -= c*prod(node-2**t for t in range(1,i+1))
        denominator=prod(node-2**t for t in range(1,j+1))
        q,r=divmod(numerator,denominator)
        need(r==0,'integral divided difference')
        coefficients.append(q)
    return [c%(2**(K-j*(j+1)//2)) for j,c in enumerate(coefficients)]


def mixed_value(c,r,K):
    node=2**r
    return sum(x*prod(node-2**i for i in range(1,j+1)) for j,x in enumerate(c))%(2**K)


def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True)
    p.add_argument('--producer',type=Path,required=True);args=p.parse_args()
    source_hash=sha256(args.producer.read_bytes()).hexdigest()
    need(source_hash=='0342e4682a78de8505fde96b71f7d34c9038dbde3ace163793a3375c7c93ab6a','pinned production source')
    native=load(args.producer,'native_newton')
    exact=load(args.repo/'05-knowledge/results/port_residual_integral_20261025/counting/check.py','native_exact')
    old=load(args.repo/'05-knowledge/results/port_residual_integral_20261025/counting/modular.py','native_old')
    N=48
    ordinary,core,partition_cases=partition_series(N)
    stats=Counter(); precision=(1,2,3,4,6,9,12,16,24)
    samples=[]
    def independent_stage(n,ell,states,reference,K):
        expected={};mod=2**K
        for (s,f),w in reference.items():
            x=w*oddpart(factorial(n-s))
            need(x.denominator==1,'fresh f integer prefix scaling')
            if not s:
                if x%mod: expected[(0,())]={0:int(x)}
                continue
            m=f[0];D=tuple(d*m-value for d,value in zip(range(1,min(ell,s)+1,2),f))
            bank=expected.setdefault((s,D),{})
            need(m not in bank,'fresh f-to-bank injection');bank[m]=int(x)
        need(set(expected)==set(states),('native structural state keys',n,ell,K))
        for (s,D),(kind,data) in states.items():
            bank=expected[(s,D)]
            if not s:
                need(kind=='dense' and data==[sum(bank.values())%mod],'separate native terminal scalar')
                stats['terminal_stage_states']+=1;continue
            b=b_of(D);need(all(m>=b for m in bank),'nonnegative source exponents')
            t=sum(j*(j+1)//2<K for j in range(s))
            need(len(data)==min(len(bank),t),('exact adaptive coordinate count',n,ell,K,s,D))
            if kind=='sparse':
                need(data=={m:x%mod for m,x in bank.items()},'native sparse coefficients and semantic support')
                stats['zero_sparse_keys']+=sum(x==0 for x in data.values())
            else:
                need(data==exact_newton(bank,b,s,K),('native independent divided differences',n,ell,K,s,D))
                stats['zero_dense_vectors']+=int(not any(data))
            for r in range(1,s+1):
                got=moment(data,b,r)%mod if kind=='sparse' else mixed_value(data,r,K)
                need(got==moment(bank,b,r)%mod,'native all-horizon moment observation')
                stats['all_node_observations']+=1
            stats['native_intermediate_states']+=1
        return len(states)
    native.audit_stage=independent_stage
    for n in range(N+1):
        snapshots,count=f_stages(n,True)
        need(count==core[n],('independent core f versus partitions',n))
        exact_count,_=exact.integral_target(n,True)
        need(exact_count==count,('inherited exact core engine',n))
        for K in precision:
            got,metrics=native.newton_target(n,K,snapshots)
            need(got==count%(2**K),('native core count',n,K))
            stats['native_complete_core_counts']+=1
            if n in (2,12,24,44,48) and K in (4,9,24):
                samples.append({'n':n,'K':K,'residue':got,'cells':metrics['cells'],'bits':metrics['residue_bits']})
        _,ordinary_f=f_stages(n,False)
        exact_ord,_=exact.integral_target(n,False)
        need(ordinary_f==ordinary[n]==exact_ord,('ordinary is distinct and unrestricted',n))
        stats['unrestricted_ordinary_exact_counts']+=1
    # Actual legal kernel transfer checks, including positive terminal reads.
    rng=Random(6464)
    for n in range(1,25):
        snapshots,_=f_stages(n,True)
        for ell,states in snapshots.items():
            # Every distinct reachable (s,D), not every redundant f response.
            keys={(s,tuple(d*f[0]-v for d,v in zip(range(1,min(ell,s)+1,2),f)))
                  for (s,f) in states if s}
            for s,D in keys:
                b=b_of(D);selected=D[(ell-1)//2] if ell<=s else 0
                for k in range(s//ell+1):
                    for e in (range((s-ell*k)//ell+1) if k else (0,)):
                        rest=s-ell*(k+e)
                        if ell==1 and rest:continue
                        Dp=tuple(D[(d-1)//2]+k*(d-gcd(ell,d)) for d in range(1,min(ell-2,rest)+1,2))
                        delta=k+b-b_of(Dp);a=ell*k
                        j0=k*(ell-1)//2+ell*k*(k-1)//2-(k-k.bit_count())+k*(ell*b-selected)
                        need(delta>=0 and j0>=0 and a+rest<=s,'actual legal shift and horizon')
                        numerator=oddpart(factorial(n-rest))
                        denominator=oddpart(factorial(n-s))*ell**k*oddpart(factorial(k))
                        factor,remainder=divmod(numerator,denominator)
                        need(not remainder,'actual exact odd quotient')
                        # beta computed by geometric power coefficient.
                        bb=[1]+[0]*e
                        for _ in range(k):bb=[sum(bb[j]*(-1)**(i-j) for j in range(i+1)) for i in range(e+1)]
                        factor*=bb[e]
                        K=1+rng.randrange(12);mod=2**K
                        need(native.coefficient(native.odd_factorials(n,mod),n-s,ell,k,e,mod)==factor%mod,'native signed coefficient')
                        bank={d:rng.randrange(-9,10) for d in (0,1,3,8)}
                        c=exact_newton(bank,0,s,K)
                        out=native.transform(c,s,rest,a,delta,j0,factor,K)
                        for r in (range(1,rest+1) if rest else (0,)):
                            want=factor*2**(j0+r*delta)*moment(bank,0,a+r)%mod
                            got=mixed_value(out,r,K) if rest else out[0]
                            need(want==got,'actual legal signed transition arithmetic')
                            stats['legal_transition_node_checks']+=1
                        # Add an invisible kernel polynomial to one random Newton basis coordinate.
                        j=rng.randrange(s+3);scale=2**max(K-j*(j+1)//2,0)
                        for r in (range(1,rest+1) if rest else (0,)):
                            need(a+r>0,'terminal reads positive node')
                            q=prod(2**(a+r)-2**i for i in range(1,j+1))
                            need(factor*2**(j0+r*delta)*scale*q%mod==0,'actual legal kernel stability')
                            stats['legal_kernel_node_checks']+=1
                        stats['actual_legal_batches']+=1
                        stats['signed_excess_batches']+=int(e>0)
                        stats['terminal_batches']+=int(rest==0)
    # Sparse support rank and invalid-operation boundaries, independently observed.
    sparse_images={}
    for exponents in ((0,), (1,), (4,), (0,4), (1,4), (0,1,4), (0,2,5)):
        K=4;S=4;mod=2**K
        images={tuple(sum(c*pow(2,r*d,mod) for c,d in zip(coef,exponents))%mod for r in range(1,S+1))
                for coef in product(range(mod),repeat=len(exponents))}
        sparse_images[str(exponents)]=len(images)
    need(sparse_images['(4,)']==1 and sparse_images['(0,)']==16,'sparse exponent hostile')
    need(all(moment({0:4},0,r)%4==0 for r in range(1,4)) and moment({0:2},0,1)%4==2,'division by two hostile')
    need(all(moment({2:1},0,r)%4==0 for r in range(1,4)) and moment({1:1},0,1)%4==2,'division by X hostile')
    # Concrete unrestricted ordinary count, independent partition check at 50.
    _,o50=f_stages(50,False);ordinary50=sum(Fraction(2**(sum((d-1)//2 for d in part)+sum(gcd(part[i],part[j]) for i in range(len(part)) for j in range(i))),prod(d**k*factorial(k) for d,k in Counter(part).items())) for part in partitions(50))
    need(o50==ordinary50 and o50%(2**25)==0,'ordinary 50 exact baseline')
    abs2,_=old.modular_target(50,2,False);wide,_=old.modular_target(50,27,False)
    need(abs2==0 and o50//2**25%4==wide//2**25==2,'absolute versus normalized precision hostile')
    # Pair-free/core is not ordinary; normalizing numerator does not change target semantics.
    need(core[2]==0 and ordinary[2]==1 and core[5]==5 and ordinary[5]==12,'core/ordinary semantic distinction')
    print(json.dumps({'status':'PASS','production_sha256':source_hash,'core_orders':[0,N],'precision_bits':precision,
       'partition_terms':partition_cases,'native_controls':dict(sorted(stats.items())),
       'core_counts_0_to_16':core[:17],'ordinary_counts_0_to_16':ordinary[:17],
       'sparse_image_sizes_K4_S4':sparse_images,'normalized_order50_hostile':{'absolute_mod4':abs2,'normalized_mod4':2,'required_absolute_bits':27},
       'samples':samples,'scope':'Independent native signed clone-excess validation; no runtime, normalized-Newton, literal-clone quotient, or ordinary-state reduction claim.'},sort_keys=True,indent=2))

if __name__=='__main__':main()
