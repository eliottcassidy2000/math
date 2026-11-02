"""Exact checks; run with and without -O. No external packages or assert gates."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import factorial, gcd, lcm
from random import Random
import argparse
import hashlib
import json
from pathlib import Path


def need(condition, label):
    if not condition:
        raise ArithmeticError(label)


def stability_constants():
    n, c, eta, a = 500, F(24999, 100000), F(1, 62), F(1, 15)
    checks = {
        'exp12_taylor_slack': sum((F(12**j, factorial(j)) for j in range(11)), F()) - 25000,
        'small_variance_slack': c / 3 - F(8, 499),
        'retained_mass_endpoint': c * F(n-1, 8) - F(4*n, n-1) - F(n, 124) - (F(n, 64)-F(5,16)),
        'retained_mass_derivative_lower': c / 8 - F(1,124) - F(1,64),
        'epsilon_integer_endpoint': 4**30-5500*3**30,
        'baseline_absorption_endpoint': c-F(32*n,(n-1)**2)-4*eta-F(4,77)-a,
        'joint_budget_slack': a-4*eta,
        'matching_size_slack': F(n//2)-4/eta,
        'variance_upper_bound_slack': F(n,7)-F(n,8)-F(4*n,n-1),
    }
    for key, value in checks.items():
        need(value > 0, key)
    need(checks['retained_mass_endpoint'] == F(654055969,12375200000), 'mass proof constant')
    need(checks['baseline_absorption_endpoint'] == F(464029788839,178309616100000), 'absorption proof constant')
    need(checks['joint_budget_slack'] == F(1,465), 'joint proof constant')
    # Independent algebra: the last contradiction for every finite endpoint
    # before the coefficient a-R*eta^2 changes sign. The proof handles all n.
    for R in range(250, 1000):
        left = eta * (a - 4*eta)
        right_coefficient = a-R*eta**2
        if right_coefficient > 0:
            need(F(4,R)*right_coefficient < left, 'last scalar contradiction')
        else:
            need(left > 0, 'negative right coefficient')
    return {k: str(v) for k,v in checks.items()}


def stability_276():
    n, R = 276, 138
    c, eta, a, delta = F(24999,100000), F(2,69), F(7,60), F(1,59)
    Q=R+F(33,2)+F(8,R)
    checks={
        'small_variance_slack': c/3-F(8,275),
        'retained_floor_slack': c*F(n-1,2)-F(4*n,n-1)-2*eta*n-14,
        'retained_derivative_lower':c/2-2*eta,
        'tail_at_one_removed':delta-Q*F(3,4)**(R-1),
        'tail_at_max_removed':delta-Q*F(3,4)**15/F(R-15),
        'absorption_slack':c-4*eta-delta-a,
        'joint_positive_slack':a-4*eta,
    }
    for key,value in checks.items():
        need(value>0, 'n276 '+key)
    need(checks['retained_floor_slack']==F(31599,88000),'n276 floor constant')
    need(checks['tail_at_max_removed']==F(18763479311,119479547723776),'n276 tail constant')
    need(checks['absorption_slack']==F(58643,135700000),'n276 absorption constant')
    need(checks['joint_positive_slack']==F(1,1380),'n276 shared constant')
    need(R*eta==4,'n276 endpoint size equality')
    # Check the whole finite L interval at the first two parities. Uniform
    # extension is proved by convexity and monotonicity in the text.
    for RR in [138,139,150,250]:
        QQ=RR+F(33,2)+F(8,RR)
        for L in range(1,RR-15+1):
            need(QQ*F(3,4)**(RR-L)/L < delta, 'n276 full endpoint interval')
        coeff=a-RR*eta**2
        if coeff>0:
            need(F(4,RR)*coeff<=eta*(a-4*eta),'n276 final non-strict bound')
        else:
            need(eta*(a-4*eta)>0,'n276 negative-coefficient case')
    # An explicit boundary rejects an unproved extension of this certificate.
    hostile=F(1,4)-F(32,137)-(137+F(33,2)+F(8,137))*F(3,4)**15/F(122)
    need(hostile<0,'n274 scalar certificate hostile')
    checks['n274_scalar_hostile']=hostile
    return {key:str(value) for key,value in checks.items()}


def graph_check(n, weights):
    edges = list(combinations(range(n), 2))
    positive = [(w, e) for w,e in zip(weights,edges) if w]
    positive.sort(key=lambda pair: (-pair[0], pair[1]))
    owners, selected, groups = {}, [], []
    discarded = []
    for w,(u,v) in positive:
        if u not in owners and v not in owners:
            j = len(selected)
            selected.append(w)
            groups.append(0)
            owners[u] = owners[v] = j
        else:
            j = min(owners[x] for x in (u,v) if x in owners)
            need(selected[j] >= w, 'greedy domination')
            groups[j] += w
            discarded.append((w,(u,v)))
    # W computed two independent ways: explicit incident pairs and degrees.
    W = sum(w*z for i,(w,e) in enumerate(positive) for z,f in positive[i+1:] if set(e)&set(f))
    deg = [0]*n
    for w,(u,v) in positive:
        deg[u] += w
        deg[v] += w
    need(2*W == sum(d*d for d in deg)-2*sum(w*w for w,e in positive), 'wedge degree identity')
    D = sum(w for w,e in discarded)
    ddeg = [0]*n
    for w,(u,v) in discarded:
        ddeg[u] += w
        ddeg[v] += w
    need(2*W >= sum(d*d for d in ddeg), 'discarded degree cancellation')
    need(n*W >= 2*D*D, 'sharp greedy bound')
    need(4*W >= sum(max(4*w*S,S*S+2*w*S) for w,S in zip(selected,groups)), 'disjoint group allocation')
    L = sum(w<=1 for w in selected)
    Dl = sum(S for w,S in zip(selected,groups) if w<=1)
    Dh = D-Dl
    if L:
        need(4*L*(W-Dh)>=Dl*Dl, 'joint heavy light budget')
    else:
        need(W>=Dh, 'all heavy budget')
    need(W>=D-L, 'completed square joint budget')


def graphs():
    total = 0
    universes = []
    for n, alphabet in [(5,range(3)), (6,range(2))]:
        count = 0
        for weights in product(alphabet, repeat=n*(n-1)//2):
            graph_check(n, weights)
            count += 1
        universes.append({'n': n, 'integer_weight_alphabet': list(alphabet), 'graphs': count})
        total += count
    rng = Random(20261102)
    for n in [7,11,25,64]:
        for _ in range(12):
            weights = [rng.randrange(8) if rng.randrange(3)==0 else 0 for _ in range(n*(n-1)//2)]
            graph_check(n, weights)
            total += 1
    for r in range(2,31):
        n=2*r
        W=r*r*(r-1)
        D=r*(r-1)
        need(F(n*W,D*D)==F(2*r,r-1), 'sharpness family')
    return {'universes':universes, 'seeded_stress_graphs':48, 'total':total,
            'sharpness_family':'K_(r,r), 2<=r<=30; symbolic ratio 2r/(r-1)'}


def reachable(program):
    return [j for j,(a,b) in enumerate(program) if all(b % program[i][1] for i in range(j))]


def phi(program, m):
    return lcm(*(b*m//gcd(a,m) for j,(a,b) in enumerate(program) if j in reachable(program)))


def compile_precision(program, horizon, m):
    """Guarded-word proof compiler, not literal trajectory simulation."""
    if horizon == 0:
        return m, []
    modulus = 1
    frontier = [((), F(1), 1, frozenset())]
    audit = []
    for t in range(horizon):
        children = []
        for word, ratio, D, forbidden in frontier:
            thresholds = [ratio.denominator*b//gcd(ratio.numerator,b) for a,b in program]
            for j,(a,b) in enumerate(program):
                next_D = lcm(D, thresholds[j])
                next_forbidden = forbidden | frozenset(thresholds[:j])
                if any(next_D % f == 0 for f in next_forbidden):
                    continue
                next_ratio = ratio*F(a,b)
                next_word = word+(j,)
                need(next_D % next_ratio.denominator == 0, 'prefix integrality')
                modulus = lcm(modulus, next_D)
                if t+1 == horizon:
                    modulus = lcm(modulus, next_ratio.denominator*m//gcd(next_ratio.numerator,m))
                row = (next_word,next_ratio,next_D,next_forbidden)
                children.append(row)
                audit.append(row)
        frontier = children
    return modulus, audit


def observe(program, N, horizon, m):
    """Independent literal first-applicable-fraction execution."""
    word = []
    for _ in range(horizon):
        for j,(a,b) in enumerate(program):
            if N % b == 0:
                word.append(j)
                N = a*(N//b)
                break
        else:
            return tuple(word), 'HALT'
    return tuple(word), N % m


def least_cyclic_period(values):
    # KMP prefix function: does not consult divisibility, guards, or the compiler.
    pi = [0]*len(values)
    for i in range(1,len(values)):
        j=pi[i-1]
        while j and values[i] != values[j]:
            j=pi[j-1]
        if values[i] == values[j]:
            j += 1
        pi[i]=j
    p=len(values)-pi[-1]
    return p if len(values)%p==0 else len(values)


def fractran():
    fractions = [(a,b) for a in range(1,5) for b in range(1,5) if gcd(a,b)==1]
    programs = [()] + [(f,) for f in fractions] + list(product(fractions,repeat=2))
    cases=states=maximum=word_count=0
    for program in programs:
        for r in range(1,4):
            for m in range(1,5):
                M,audit = compile_precision(program,r,m)
                sufficient=m
                for _ in range(r):
                    sufficient=phi(program,sufficient)
                need(sufficient%M==0, 'exact divides iterated compiler')
                if r==1:
                    need(M==phi(program,m), 'one step formula')
                for word,ratio,D,forbidden in audit:
                    trace=observe(program,D,len(word),1)
                    need(trace==(word,0), 'realization witness')
                values=[observe(program,N,r,m) for N in range(1,2*M+1)]
                need(values[:M]==values[M:], 'literal periodicity')
                need(least_cyclic_period(values[:M])==M, 'literal minimal period')
                cases+=1
                states+=2*M
                maximum=max(maximum,M)
                word_count+=len(audit)
    # Complete guarded-word comparison on an infinite-family parameter box;
    # the proof, not this finite box, supplies the universal quantifiers.
    family_cases=0
    for a in range(2,9):
        for b in range(2,9):
            if gcd(a,b)>1:
                continue
            program=((a,b),(b,a))
            for r in range(1,7):
                for m in range(1,17):
                    M,rows=compile_precision(program,r,m)
                    expected=a*b**r*m if r%2 else b**r*lcm(a,m)
                    need(M==expected, 'reciprocal parity formula')
                    need(sum(len(w)==r for w,ratio,D,forbidden in rows)==r+1, 'reciprocal live words')
                    family_cases+=1
    need(observe(((1,4),(3,2)),2,1,100)==((1,),3), 'squareful hostile state 2')
    need(observe(((1,4),(3,2)),8,1,100)==((0,),2), 'squareful hostile state 8')
    need(2%6==8%6, 'hostile same primorial digit')
    need(compile_precision(((2,3),(3,2)),2,1)[0]==18, 'precision shadowing hostile')
    need(phi(((2,3),(3,2)),phi(((2,3),(3,2)),1))==36, 'naive precision hostile')
    need(observe(((2,3),(3,2)),4,2,1)==((1,0),0), 'unrealizable repeated second rule')
    for p in [2,3,5,11]:
        for r in range(1,9):
            for m in range(1,13):
                need(compile_precision(((1,p),),r,m)[0]==m*p**r, 'unbounded horizon family')
    carrier_states=0
    for a,b in [(2,3),(3,5),(5,2)]:
        def state(N,h):
            if not h:
                return 'END'
            t=0
            while t<h and N%b==0:
                N//=b
                t+=1
            if t:
                return t
            return 'ALT2' if N%a==0 else 'HALT'
        for h in range(1,7):
            observations={}
            for N in range(1,a*b**h+1):
                s=state(N,h)
                value=observe(((a,b),(b,a)),N,h,1)
                if s in observations:
                    need(observations[s]==value,'minimal carrier sufficient')
                observations[s]=value
                if h>=2 and s!='HALT':
                    nxt=a*N//b if N%b==0 else b*N//a
                    predicted=1 if s=='ALT2' else ('ALT2' if s==1 else s-1)
                    need(state(nxt,h-1)==predicted,'clocked carrier transition')
                carrier_states+=1
            need(len(observations)==h+2,'minimal carrier size')
            need(len(set(observations.values()))==h+2,'minimal carrier necessity')
    return {'fraction_alphabet':fractions,'programs':len(programs),'horizons':[1,2,3],
            'target_moduli':[1,2,3,4], 'literal_cases':cases,'literal_states':states,
            'maximum_minimal_modulus':maximum,'realized_word_checks':word_count,
            'reciprocal_family_checks':family_cases,
            'minimal_clocked_carrier_integer_states':carrier_states,
            'hostiles':{'same_mod_6_different_rules':[2,8], 'minimal_two_step':18,'iterated_one_step':36}}


def determinant(A):
    B=[list(map(F,row)) for row in A]
    ans=F(1)
    for j in range(len(B)):
        pivot=next((i for i in range(j,len(B)) if B[i][j]),None)
        if pivot is None:
            return F(0)
        if pivot!=j:
            B[j],B[pivot]=B[pivot],B[j]
            ans=-ans
        d=B[j][j]
        ans*=d
        for i in range(j+1,len(B)):
            q=B[i][j]/d
            for k in range(j+1,len(B)):
                B[i][k]-=q*B[j][k]
    return ans


def native(A):
    n=len(A)
    dp={}
    for i in range(n):
        dp[(1<<i,i)]=F(1)
    for S in range(1,1<<n):
        for v in range(n):
            if (S,v) not in dp:
                continue
            for w in range(n):
                if not S>>w&1:
                    key=(S|1<<w,w)
                    dp[key]=dp.get(key,F())+dp[S,v]*A[v][w]
    H=sum(dp.get(((1<<n)-1,i),F()) for i in range(n))
    u=determinant([[A[i][j]+(i==j) for j in range(n)] for i in range(n)])
    v=(-1)**n*determinant(A)
    V=sum(A[i][j]*(1-A[i][j]) for i,j in combinations(range(n),2))
    return H,u,v,V


def pad(A):
    n=len(A)
    B=[[F() for _ in range(n+2)] for _ in range(n+2)]
    for i in range(n):
        for j in range(n):
            B[i][j]=A[i][j]
        B[i][n]=B[i][n+1]=F(1)
    B[n][n+1]=B[n+1][n]=F(1,2)
    return B


def padding():
    sources=[[[F(0) if i==j else F(1,2) for j in range(3)] for i in range(3)],
             [[F(0),F(1),F(0)],[F(0),F(0),F(1)],[F(1),F(0),F(0)]]]
    rows=[]
    for A in sources:
        initial=native(A)
        for r in range(1,3):
            A=pad(A)
            H,u,v,V=native(A)
            H0,u0,v0,V0=initial
            need((H,u,abs(v),V)==(H0,u0*F(3,4)**r,abs(v0)*F(1,4)**r,V0+F(r,4)), 'literal ordinal padding')
            rows.append({'r':r,'H':str(H),'u':str(u),'abs_v':str(abs(v)),'V':str(V)})
    return rows


def pair_reserve():
    cases=0
    for p,a,b in product((F(i,8) for i in range(9)),repeat=3):
        tau=p*a*(1-b)+(1-p)*(1-a)*b
        alpha=min(p,1-p)
        qsum=a*(1-a)+b*(1-b)
        need(tau>=alpha*(qsum+(a-b)**2), 'joint profile reserve')
        if p==F(1,2):
            need(tau==alpha*(qsum+(a-b)**2), 'fair-edge identity')
        need(alpha>=p*(1-p), 'variance coefficient improvement')
        cases+=1
    # Same marginal variances, different cyclic probability.
    p=F(1,2)
    need(p*0*(1-0)+(1-p)*(1-0)*0==0, 'aligned port hostile')
    need(p*1*(1-0)+(1-p)*(1-1)*0==F(1,2), 'discordant port hostile')
    return {'triangle_grid_cases':cases,'grid_denominator':8,
            'same_variances_different_cycles':['0','1/2']}


def square_sum():
    counts={}
    visits=0
    for n in range(1,26):
        adj=[0]*n
        squares={i*i for i in range(2,8)}
        for u,v in combinations(range(n),2):
            if u+v+2 in squares:
                adj[u]|=1<<v
                adj[v]|=1<<u
        leaves=[i for i in range(n) if adj[i].bit_count()==1]
        if n>1 and (any(a==0 for a in adj) or len(leaves)>2):
            counts[str(n)]=0
            continue
        full=(1<<n)-1
        @lru_cache(None)
        def extend(S,last):
            if S==full:
                return 1
            candidates=adj[last]&~S
            total=0
            while candidates:
                bit=candidates&-candidates
                candidates-=bit
                total+=extend(S|bit,bit.bit_length()-1)
            return total
        if leaves:
            u=leaves[0]
            undirected=extend(1<<u,u)
        else:
            directed=sum(extend(1<<u,u) for u in range(n))
            need(n==1 or directed%2==0, 'square-sum reversal parity')
            undirected=directed if n==1 else directed//2
        counts[str(n)]=undirected
        visits+=extend.cache_info().misses
    need({n for n in range(2,26) if counts[str(n)]}=={15,16,17,23,25}, 'small onset set')
    need([counts[str(n)] for n in [15,16,17,23,25]]==[1,1,1,3,10], 'small path counts')
    return {'counts_up_to_reversal':counts,'memo_states':visits,
            'scope':'finite n=1..25 only; no all-size construction claimed'}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    report={'status':'PASS','stability_276':stability_276(),
            'intermediate_stability_500':stability_constants(),
            'graphs':graphs(),'fractran':fractran(),'native_padding':padding(),
            'pair_reserve':pair_reserve(),'square_sum':square_sum()}
    payload=json.dumps(report,sort_keys=True,separators=(',',':')).encode()
    report['semantic_sha256']=hashlib.sha256(payload).hexdigest()
    report['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':report['status'],'graphs':report['graphs']['total'],
                      'literal_fractran_cases':report['fractran']['literal_cases'],
                      'semantic_sha256':report['semantic_sha256']}))


if __name__=='__main__':
    main()
