"""Exact content-aware counting, dyadic FFT loss, and signed Fano crossing checks."""
from collections import Counter
from fractions import Fraction
from functools import reduce
from itertools import combinations, product
from math import comb, factorial, gcd, lcm
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
GATES = 0


def need(value, reason):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(reason)


def v2(n):
    return (n & -n).bit_length()-1


def effective_period(scale, moduli):
    return lcm(*(m//gcd(m, scale) for m in moduli))


def merge_crt(a, mod_a, b, mod_b):
    g = gcd(mod_a, mod_b)
    need((b-a) % g == 0, "compatible generalized CRT")
    q = mod_b//g
    t = 0 if q == 1 else ((b-a)//g)*pow(mod_a//g, -1, q) % q
    mod = mod_a*q
    return (a+mod_a*t) % mod, mod


def unscale_residues(scale, pairs):
    a, mod = 0, 1
    for m, r in pairs:
        d = scale % m
        g = gcd(m, d)
        need(r % g == 0, "scaled residue belongs to the multiplier image")
        q = m//g
        b = 0 if q == 1 else (r//g)*pow(d//g, -1, q) % q
        a, mod = merge_crt(a, mod, b, q)
    return a, mod


def content_checks():
    cases = 0
    for scale in range(1, 49):
        for cap in range(2, 20):
            moduli = list(range(1, cap+1))
            whole = lcm(*moduli)
            period = effective_period(scale, moduli)
            need(period == whole//gcd(whole, scale), "exact content quotient of lcm")
            for c in [0,1,2,7,period-1,period,period+1]:
                residues = [(m, scale*c % m) for m in moduli]
                actual, actual_period = unscale_residues(scale, residues)
                need((actual,actual_period) == (c % period,period), "constructive complete fibre")
                need(all(scale*c % m == scale*(c+period) % m for m in moduli),
                     "period is an actual collision")
                cases += 1
            if period <= 256:
                signatures = {tuple(scale*c % m for m in moduli) for c in range(period)}
                need(len(signatures) == period, "independent full-period injectivity")
    odd_cases=[]
    for n in range(0, 65):
        M=n+3
        bank=list(range(3,2*M,2))
        Q=lcm(*bank)
        B=comb(2*M,M)
        need(v2(B)==M.bit_count(), "central-binomial two-valuation")
        oddB=B//(1<<v2(B))
        need(Q%oddB==0, "odd central binomial divides odd lcm")
        need(Q>2**n, "linear-size odd bank exceeds original count height")
        for E in [0,1,7,31,1000]:
            need(effective_period(1<<E,bank)==Q, "dyadic scale invisible to odd bank")
        odd_cases.append({'bits':n,'M':M,'moduli':len(bank),'largest_modulus':bank[-1],
                          'effective_period':Q})
    saturation=[]
    for cap in [5,9,15,31]:
        scale=lcm(*range(1,cap+1))
        need(effective_period(scale,range(1,cap+1))==1, "full-content hostile erases bounded bank")
        saturation.append({'cap':cap,'scale':scale,'zero_and_one_indistinguishable':True})
    # The external engine's documented generic grid multiplier, not just 2^E.
    grid=[]
    for g in range(1,6):
        rect=comb(g,2)**2
        nonrect=g**4-rect
        for lam in range(1,16):
            t=(lam+2)//3
            scale=(t+1)**nonrect * 2**(t*rect)
            odd=scale >> v2(scale)
            need(odd==((t+1) >> v2(t+1))**nonrect, "actual generic grid odd content")
            grid.append({'G':g,'lambda':lam,'t':t,'rectangles':rect,
                         'nonrectangles':nonrect,'scale':scale,'odd_part':odd})
    return {'recovered_fibres':cases,'odd_banks':odd_cases,'saturation_hostiles':saturation,
            'generic_grid_multipliers':grid}


def partitions(n, low=1):
    if n == 0:
        yield ()
    for a in range(low,n+1,2):
        for rest in partitions(n-a,a):
            yield (a,)+rest


def orbit_terms(n, colors):
    for lam in partitions(n):
        multiplicities=Counter(lam)
        z=1
        for a,c in multiplicities.items(): z*=a**c*factorial(c)
        energy=sum((a-1)//2 for a in lam)+sum(gcd(a,b)
                   for i,a in enumerate(lam) for b in lam[i+1:])
        yield z,energy,len(lam)


def orbit_checks():
    bank=[]
    known=[1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168]
    for n in range(0,21):
        for ports in range(3):
            Q=2**ports; group=factorial(n)
            terms=list(orbit_terms(n,Q))
            direct=sum((Fraction(2**e*Q**c,z) for z,e,c in terms),Fraction())
            need(direct.denominator==1, "exact rational Burnside integrality")
            count=direct.numerator
            if ports==0 and n<len(known):
                need(count==known[n], "inherited A000568 values, no new terms")
            # Height certificate uses labelled orientations and individually named ports.
            bound=2**(n*(n-1)//2+ports*n)
            whole=1; moduli=[]
            for m in range(1,1000):
                whole=lcm(whole,m); moduli.append(m)
                if whole//gcd(whole,group)>bound: break
            else: raise RuntimeError("modulus budget insufficient")
            pairs=[]
            for m in moduli:
                num=sum((group//z)*pow(2,e,m)*pow(Q,c,m) for z,e,c in terms)%m
                pairs.append((m,num))
            decoded, period=unscale_residues(group,pairs)
            need(decoded==count and period>bound, "full orbit count from scaled modular numerator")
            bank.append({'n':n,'named_ports':ports,'count':count,'last_modulus':moduli[-1],
                         'residue_period':period,'labelled_upper_bound':bound})
    # Native canonical hostile; import immutable helper, not a guessed group order.
    src=HERE.parent/'regular_parent_symmetry_20261029'/'verify.py'
    spec=importlib.util.spec_from_file_location('inherited_regular_parent',src)
    old=importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
    hostile=old.hostile_checks()
    receipts=[]
    for row in hostile['tournaments']:
        group=row['group_order']; raw=row['regular_binary_masks']; free=raw//group
        need(raw%group==0, "free-mask orbits divide by actual full group")
        residue=raw%(3*group)
        need(residue//group==free%3, "one trit after exact group-order precision")
        receipts.append({'group':group,'regular_masks':raw,'free_orbits':free,
                         'raw_mod_3':raw%3,'safe_modulus':3*group,
                         'recovered_orbits_mod_3':residue//group})
    need([(r['group'],r['free_orbits']) for r in receipts]==[(81,20),(27,64)],
         "81/27 full source discrimination")
    return {'orbit_counts':bank,'hostile_81_27':hostile,'division_receipts':receipts}


def dft(x, root, mod):
    T=len(x)
    return [sum(a*pow(root,i*j,mod) for i,a in enumerate(x))%mod for j in range(T)]


def conv(x,y,mod):
    T=len(x)
    return [sum(x[j]*y[(i-j)%T] for j in range(T))%mod for i in range(T)]


def eps_coordinates(bits):
    # X=1+epsilon; Lucas parity gives binom(j,i) odd iff i is a bit-subset of j.
    return [sum(bits[j] for j in range(i,len(bits)) if i & j == i)%2
            for i in range(len(bits))]


def truncated_product(x,y):
    return [sum(x[j]*y[i-j] for j in range(i+1))%2 for i in range(len(x))]


def fft_checks():
    result=[]
    for p in [3,5,7,11,13]:
        for a in range(1,5):
            m=2**a; T=2*m; N=p**m+1; odd=N//2
            need(N%4==2 and gcd(T,odd)==1, "exact 2-primary boundary")
            root=p%N
            need(pow(root,m,N)==N-1, "primitive dyadic half-turn")
            need(all(sum(pow(root,i*j,N) for i in range(T))%N==
                     (T%N if j==0 else 0) for j in range(T)), "principal root sums")
            # Actual kernel generators and independent inverse on odd component.
            generators=[]
            for j in range(T-1):
                x=[0]*T; x[j]=odd; x[-1]=odd
                need(not any(dft(x,root,N)), "explicit missing parity generator")
                generators.append(x)
            for seed in range(8):
                x=[(seed*17+j*j+3*j+seed*j)%N for j in range(T)]
                y=[(seed+j*7+j*j*j)%N for j in range(T)]
                fx=dft(x,root,N)
                inv=dft(fx,pow(root,-1,odd),odd)
                inv=[v*pow(T,-1,odd)%odd for v in inv]
                need(inv==[v%odd for v in x], "odd component inverse DFT")
                parities=[v%2 for v in x[:-1]]
                parities.append((fx[0]-sum(parities))%2)
                rebuilt=[v+odd*((bit-v)%2) for v,bit in zip(inv,parities)]
                need(rebuilt==x, "optimal parity-sidecar reconstruction")
                need(dft(conv(x,y,N),root,N)==[u*v%N for u,v in zip(fx,dft(y,root,N))],
                     "convolution homomorphism including singular part")
                bx=[v%2 for v in x]; by=[v%2 for v in y]
                need(eps_coordinates(conv(bx,by,2))==
                     truncated_product(eps_coordinates(bx),eps_coordinates(by)),
                     "retained parity algebra is F2[epsilon]/epsilon^T")
            # Shift's faithful parity clock; (1+epsilon) has order T.
            e=[0]*T; e[1]=1
            need(eps_coordinates(e)==[1,1]+[0]*(T-2), "cyclic shift equals one plus epsilon")
            clock=[1,1]+[0]*(T-2)
            powers=[[1]+[0]*(T-1)]
            for _ in range(T):
                powers.append(truncated_product(powers[-1],clock))
            need(powers[T]==powers[0] and powers[T//2]!=powers[0],
                 "actual nilpotent shift has full dyadic order")
            hidden=[1,1]+[0]*(T-2)
            shifts={tuple(hidden[-j:]+hidden[:-j]) if j else tuple(hidden) for j in range(T)}
            need(len(shifts)==T, "hidden parity word has full shift orbit")
            need(all(sum(w)%2==0 for w in shifts), "entire hidden orbit has zero transform")
            need(all(not any(dft([odd*b for b in w],root,N)) for w in shifts),
                 "full lifted hidden shift orbit has zero Fourier output")
            need(eps_coordinates(eps_coordinates(hidden))==hidden, "Pascal change of coordinates invertible")
            square=eps_coordinates(conv(hidden,hidden,2))
            eps=eps_coordinates(hidden)
            need(square==[eps[i//2] if i%2==0 else 0 for i in range(T)],
                 "actual characteristic-two Frobenius is doubling epsilon degree")
            result.append({'p':p,'a':a,'m':m,'length':T,'modulus':N,'odd_component':odd,
                           'missing_bits':T-1,'kernel_order':2**(T-1),
                           'faithful_shift_order':T,'frobenius_image_dimension':T//2})
    # Complete kernel fibre enumeration, separately from the inverse tests.
    exhaustive=[]
    for T in [4,8,16]:
        p=3; m=T//2; N=p**m+1; odd=N//2
        images=Counter()
        for bits in product(range(2),repeat=T):
            # All p^i are odd, so these binary-lift images are constant vectors.
            output=tuple(odd*(sum(bits)%2) for _ in range(T))
            images[output]+=1
        need(len(images)==2 and images[(0,)*T]==2**(T-1), "complete parity kernel census")
        exhaustive.append({'length':T,'parity_inputs':2**T,'image_sizes':sorted(images.values())})
    full_images=Counter(tuple(dft(x,3,10)) for x in product(range(10),repeat=4))
    need(len(full_images)==1250 and set(full_images.values())=={8},
         "independent complete Z/10 length-four transform census")
    # Positive control in the paper's odd Fermat-type modulus: inversion is legal.
    positive=[]
    for a in range(1,5):
        m=2**a; T=2*m; N=2**m+1
        need(gcd(T,N)==1, "base-two modulus permits division by transform length")
        x=list(range(T)); fx=dft(x,2,N)
        rebuilt=[v*pow(T,-1,N)%N for v in dft(fx,pow(2,-1,N),N)]
        need(rebuilt==[v%N for v in x], "base-two full inverse")
        positive.append({'m':m,'T':T,'modulus':N})
    return {'odd_base_transforms':result,'parity_censuses':exhaustive,
            'full_Z10_length4':{'inputs':10000,'images':len(full_images),'fibre_size':8},
            'base_two_controls':positive}


def oct_sign(x,y):
    p=[x>>j&1 for j in range(3)]; q=[y>>j&1 for j in range(3)]
    phi=sum(p[i]*q[j] for i in range(3) for j in range(i,3))
    phi+=q[0]*p[1]*p[2]+p[0]*q[1]*p[2]+p[0]*p[1]*q[2]
    return (-1)**(phi%2)


def norm_polynomial_energy(sign,q):
    coefficients=Counter()
    for z in range(q):
        for x,v in product(range(q),repeat=2):
            y=x^z; u=v^z
            key=(tuple(sorted((x,v))),tuple(sorted((y,u))))
            coefficients[key]+=sign(x,y)*sign(v,u)
    for x,y in product(range(q),repeat=2):
        coefficients[((x,x),(y,y))]-=1
    need(all(not c or (a[0]!=a[1] and b[0]!=b[1])
             for (a,b),c in coefficients.items()), "norm defect contains only distinct-coordinate monomials")
    return sum(c*c for c in coefficients.values())


def general_crossing_energy_checks():
    rows=[]
    for q in [1,2,4,8,16]:
        for seed in range(6):
            def sign(x,y):
                return (-1)**(((x*7+y*11+x*y*3+seed)*(seed+1)^(x>>(seed%3)))&1)
            energy=0
            for x,y,u,v in product(range(q),repeat=4):
                supported=x^y==u^v
                defect=(sign(x,y)*sign(v,u)+sign(x,u)*sign(v,y) if supported else 0)
                defect-=2*int(x==v and y==u)
                energy+=defect*defect
            need(energy==norm_polynomial_energy(sign,q),
                 "all-size signed-XOR crossing norm equals norm-composition mean-square defect")
            rows.append({'dimension':q,'seed':seed,'defect_energy':energy})
    return rows


def crossing_checks():
    lines=sorted({tuple(sorted((x,y,x^y))) for x in range(1,8)
                  for y in range(x+1,8)})
    index={frozenset(line):i for i,line in enumerate(lines)}
    need(len(lines)==7, "actual XOR Fano incidence")
    normals=[next(a for a in range(1,8) if all((a&x).bit_count()%2==0 for x in line)) for line in lines]
    need(sorted(normals)==list(range(1,8)), "all seven distinct nonzero line normals")
    gauge_masks=set()
    for f in range(128):
        mask=sum((sum(f>>(x-1)&1 for x in line)%2)<<i for i,line in enumerate(lines))
        gauge_masks.add(mask)
    need(len(gauge_masks)==16, "inherited four-dimensional Fano gauge image")
    good=[]; bad=[]; base_neg=0; base_nonzero=0
    first=None
    for mask in range(128):
        syndrome=0
        for i,a in enumerate(normals):
            if mask>>i&1: syndrome^=a
        need((syndrome==0)==(mask in gauge_masks), "explicit Hamming syndrome equals gauge image")
        def s(x,y):
            if not x or not y or x==y: return oct_sign(x,y)
            return oct_sign(x,y)*(-1 if mask>>index[frozenset((x,y,x^y))]&1 else 1)
        failures=0; defect_energy=0
        for x,y,u,v in product(range(8),repeat=4):
            r=s(x,y)*s(v,u) if x^y==v^u else 0
            c=s(x,u)*s(v,y) if x^u==v^y else 0
            target=2*int(x==v and y==u)
            if mask==0:
                base_nonzero+=bool(r or c); base_neg+=int(r<0)+int(c<0)
                need(r+c==target, "external crossing identity replayed literally")
                if first is None and target==0 and r==-c and r:
                    first={'x':x,'y':y,'u':u,'v':v,'R':r,'C':c,'absolute_sum':abs(r)+abs(c)}
            if r+c!=target: failures+=1
            defect_energy+=(r+c-target)**2
        need((failures==0)==(mask in gauge_masks), "crossing-valid iff inherited Fano gauge code")
        need(defect_energy==norm_polynomial_energy(s,8), "four ordered ports per norm rectangle")
        need(defect_energy==(0 if failures==0 else 384), "exact crossing and norm energy gap")
        (good if failures==0 else bad).append({'mask':mask,'syndrome':syndrome,'failures':failures,
                                             'defect_energy':defect_energy})
    need(len(good)==16 and len(bad)==112, "inherited good/defective Fano split")
    need(all(row['failures']==96 for row in bad), "every defective table has exactly 96 failed ports")
    # All future named line flips: each present syndrome has a different acceptance row.
    sigma=lambda mask: reduce(int.__xor__,
                         (normals[i] for i in range(7) if mask>>i&1),0)
    futures={mask:tuple(sigma(mask^a)==0 for a in range(128)) for mask in range(128)}
    need(len(set(futures.values()))==8, "smallest future-legality carrier has exactly eight states")
    need(all((futures[x]==futures[y])==(sigma(x)==sigma(y)) for x in range(128) for y in range(128)),
         "future-legality equivalence is exactly Hamming syndrome")
    branch_signatures=set()
    for mask in range(128):
        readings=[]
        for i,line in enumerate(lines):
            x,y=line[:2]
            # The named R branch at u=0, v=x XOR y reads the multiplication sign.
            readings.append(oct_sign(x,y)*(-1 if mask>>i&1 else 1))
        branch_signatures.add(tuple(readings))
    need(len(branch_signatures)==128, "actual named branch coefficients retain all seven line bits")
    need(all(sigma(mask^127)==sigma(mask) for mask in range(128)), "converse closes on syndrome")
    # Every linear point re-marking acts contragrediently on the syndrome.
    linear_maps=[]; base_orbit=set(); base_stabilizer=0; syndrome_maps=set()
    for cols in product(range(1,8),repeat=3):
        A=tuple(reduce(int.__xor__,(cols[j] for j in range(3) if x>>j&1),0)
                for x in range(8))
        if len(set(A))!=8: continue
        line_perm=[index[frozenset(A[x] for x in line)] for line in lines]
        dual=[0]*8
        for i,normal in enumerate(normals): dual[normal]=normals[line_perm[i]]
        need(all(dual[x^y]==dual[x]^dual[y] for x,y in product(range(8),repeat=2)),
             "contragredient action on line normals is linear")
        offset=0
        for i,line in enumerate(lines):
            x,y=line[:2]
            bit=int(oct_sign(x,y)!=oct_sign(A[x],A[y]))
            need(all(int(oct_sign(x,y)!=oct_sign(A[x],A[y]))==bit
                     for x,y in combinations(line,2)), "one orientation bit per transported line")
            offset |= bit << line_perm[i]
        need(sigma(offset)==0, "reference sign cocycle stays in the Hamming code")
        for mask in range(128):
            mapped=offset ^ sum(((mask>>i)&1)<<line_perm[i] for i in range(7))
            need(sigma(mapped)==dual[sigma(mask)], "native Fano re-marking closes on syndrome")
        linear_maps.append(A); base_orbit.add(offset); base_stabilizer+=int(offset==0)
        syndrome_maps.add(tuple(dual))
    need(len(linear_maps)==len(syndrome_maps)==168, "faithful GL(3,2) syndrome action")
    affine_maps={tuple(s^t for s in A) for A in syndrome_maps for t in range(8)}
    need(len(affine_maps)==1344, "future flips and re-markings induce AGL(3,2)")
    need(len(base_orbit)*base_stabilizer==168, "native Fano base orbit-stabilizer")
    need(len(base_orbit | {m^127 for m in base_orbit})==16,
         "linear re-markings plus converse reach the entire good code")
    return {'xor_lines':lines,'line_normals':normals,'good_masks':good,'defective_masks':bad,
            'negative_base_terms':base_neg,'supported_quadruples':base_nonzero,
            'absolute_value_hostile':first,'future_legality_states':len(set(futures.values())),
            'named_branch_states':len(branch_signatures),
            'remarking':{'linear_group_order':len(linear_maps),'affine_syndrome_group_order':len(affine_maps),
                         'base_orbit_size':len(base_orbit),'base_stabilizer_order':base_stabilizer},
            'coarse_acceptance_hostile':{'masks':[1,2],'future_flip':1,
                                        'result_syndromes':[sigma(0),sigma(3)]}}


def divisor_search_checks():
    result=[]
    for N in range(2,501):
        divisors=[d for d in range(2,N+1) if N%d==0]
        bits=N.bit_length(); prefix=0; path=[]
        for remaining in range(bits-1,-1,-1):
            trial=2*prefix
            c=sum(d>>remaining==trial for d in divisors)
            scale=2**(bits+3)
            bank=list(range(3,2*(bits+3),2))
            zero=all(scale*c%m==0 for m in bank)
            need(zero==(c==0), "scaled modular zero gives correct prefix branch")
            prefix=trial+int(zero)
            path.append(prefix)
        need(prefix==divisors[0], "native divisor-prefix reconstruction")
        if N in [25,27,57,91,169,343,499]:
            result.append({'N':N,'smallest_factor':prefix,'prefix_path':path})
    return result


def recovered_trace_checks():
    rows=[]
    for a,mod in [(5,13),(2,11),(5,29)]:
        def mul(x,y):
            r,s=x; t,u=y
            return ((r*t+s*u)%mod,(r*u+s*t+a*s*u)%mod)
        def trace(x): return (2*x[0]+a*x[1])%mod
        power=(1,0)
        for _ in range(12): power=mul(power,(0,1))
        sampled=set(); primitive=set()
        for x in product(range(mod),repeat=2):
            z=x; history=[]
            for _ in range(4):
                history.append(trace(z)); z=mul(z,power)
            sampled.add(tuple(history))
            primitive.add((trace(x),trace(mul(x,(0,1)))))
        delta=gcd(2,a)
        need(len(primitive)==mod**2//(gcd(mod,delta)*gcd(mod,(a*a+4)//delta)),
             "inherited primitive trace Smith factors replayed")
        if mod==13:
            need(power==(1,0) and len(sampled)==13 and len(primitive)==169,
                 "12/13 clock collapse and exact primitive recovery")
        if mod==11:
            need(power==(10,0) and len(sampled)==11 and len(primitive)==121,
                 "twelfth clock alternates at eleven")
        if mod==29:
            need(len(primitive)==29, "ramified trace cannot be made faithful by primitive clock")
        rows.append({'a':a,'modulus':mod,'source_states':mod**2,'u12':power,
                     'four_sampled_readings_image':len(sampled),'primitive_pair_image':len(primitive)})
    return rows


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=HERE/'certificate.json')
    args=parser.parse_args()
    data={'content':content_checks(),'orbit':orbit_checks(),'fft':fft_checks(),
          'crossing':crossing_checks(),'prefix':divisor_search_checks(),
          'recovered_trace_clocks':recovered_trace_checks(),
          'general_crossing_energy':general_crossing_energy_checks()}
    data['gates']=GATES
    raw=json.dumps(data,indent=2)+'\n'
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps({'gates':GATES,'sha256':hashlib.sha256(raw.encode()).hexdigest()}))


if __name__=='__main__':
    main()
