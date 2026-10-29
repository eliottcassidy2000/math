#!/usr/bin/env python3
"""Exact A000568 recreation and composition audit; Python standard library only.

No tournament-class list is used for the counts. OEIS is an audit, not an input
to construction. Prime recovery/recomposition is explicitly a Burnside round trip.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import factorial, gcd, lcm
import argparse
import json


# Audit reference only: https://oeis.org/A000568/b000568.txt, retrieved 2026-10-28.
# Construction below never consults this tuple; comparison occurs afterward.
OEIS_A000568_0_50 = (
    1,
    1,
    1,
    2,
    4,
    12,
    56,
    456,
    6880,
    191536,
    9733056,
    903753248,
    154108311168,
    48542114686912,
    28401423719122304,
    31021002160355166848,
    63530415842308265100288,
    244912778438520759443245824,
    1783398846284777975419600287232,
    24605641171260376770598003978281472,
    645022068557873570931850526424042500096,
    32207364031661175384456332260036660040346624,
    3070169883150468336193188889176239554269865953280,
    559879382429394075397997876821117309031348506639435776,
    195692027657521876084316842660833482785173437775365039898624,
    131326696677895002131450257709457767457170027052967027982788816896,
    169484335125246268100514597385576342667201246238506672327765919863947264,
    421255599848131447082003884098323929861369544621589389269735693986231100612608,
    2019284625667208265086928694043799677058780746074756618649807453554008410636526845952,
    18691296182213712407784892577100643237772159079535345610331272616359410643727554822061146112,
    334493774260141796028606267674709437232608940215918926763659414050175507824571200950884097540096000,
    11585805009377407483760329507477613962387207099806115603395226601918143978183445560553337846836297396682752,
    777510212704809602577105489198198692851891665391100208885324555434113885980852261704841726000960370177169043226624,
    101193361693125457990849046332017257500954450121854613136962529367013247492336823109239997244028630456633654846271092097024,
    25566010532016016687406724976631466697641013194578052960760500860844908999003850753699781560435594001532021694030270428757225373696,
    12549163328480034384022586604627435683222335285625319361365805379963449319586848339828196059799465957306927089740943745163632094948026875904,
    11977388019552054692083922151958065860796846899026839549253416480609097326545504675055283701368974668002855110699910350933703101581087772125367369728,
    22245401009936554231489011970059411118940446748866376454574020373414881769002409820397140027794167192183948888599628431109193420826617865995050529774636892160,
    80457490378175048878236241093915935987362717294991651875299540184497475761612685289771920455116187482060432330331817434121649131892824379319217954379657337242063994880,
    567076578285122471704586780782720308902114040542520919186644143060785152994318019686632497789949430912316036731894249260999702300870338609414965732474947326740618008554287661056,
    7793841145798991288339541567431082530560199377689930717403316856994498241370267414973048275931891882666312041777200765060752406653234738641596384861093820021561209315748075967788164841472,
    209010218654756435461084627123436598880395661738635597818413624103552400467500965951449254650876548284378862411200483352322220803409319721629824988334580040609936151806821055449455814143969084833792,
    10943293606424234728894381377588824215781834358634661636384420748138320373316940088367934824116446827210456074517628129515093533683892507110724627355895755814537618864508356606665235598809366626681079560404992,
    1119281727109786401619969084598953639932635647019577789486358300929322778125436049834173527372065702578768334343342346438379022244832696522465596152241689667369911087992777361948021957036903259214109391266883646940774400,
    223756958857166158963681037861726180036702556083780548679770723997088438981299856082958914061037863408379154922270825551856972794406643078714186069382334929988917527072961645849037502172291091334184458027644196374615900824193204224,
    87474978865511190799790879819946602468651463196700303025331822925608712240755318590862808175942822666915371488638435605019149846325035884452279153919715904824713139598609118229549714329948250882252270659436382961312997790007671905617912528896,
    66907656627540457170501802895638731118154547896184908665803219707779217672450123352587217995764427523603237210632321128283699630153779980478664552078038130261217473130129635131994761104111595479162606596972730174987406977495342158864248759129286134202368,
    100174633462774035901042871788797518465544988805495671189667445723448496499822925385704689467399853127956227466515273357659586240806376614658646206117923114350626056227448060939476281661431267903042052614164351205856445781131568634152583606075485729733889874412437504,
    293715131468050239237393586298190886514034052337748256266893125857197421465509729077400622506985110333384537639252042267380141886525721900488681470438795099330564673130908450859567033860130949293669445559871075429367677508771221785675250135041356372820059642983984004209534042112,
    1687213465092585843361181736380139541941055912056738893805579966474686655501866254673931726160227237600108252660206909117152060507311356946815848464425089091314431812637433889290136314482038634434254654575046115063583638291749189326981428796797199387785877941192685629316977695634589218766848,
    18996334831713632411474577006236471641541836655076071576812389618726770007568256674117763406499904656118507200178560110821023466608873904455425696482070096778094878929660913666378608972920617911174855467425281513498391999800181323362375354207355488730927747317479653572915278014559169658331645189344460800,
)


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def partitions(n, low=1, odd=True):
    if not n:
        yield ()
    for d in range(low, n + 1, 2 if odd else 1):
        for rest in partitions(n-d, d, odd):
            yield (d,) + rest


def centralizer(part):
    answer = 1
    for d, count in Counter(part).items():
        answer *= d**count * factorial(count)
    return answer


def gcd_exponent(part):
    return sum((d-1)//2 for d in part) + sum(gcd(d,e) for d,e in combinations(part,2))


def phi(n):
    return sum(gcd(k,n) == 1 for k in range(1,n+1))


def cyclic_burnside_exponent(part):
    """Independent pair-orbit count by averaging powers of the permutation."""
    order = lcm(*part) if part else 1
    numerator = 0
    for d in range(1, order+1):
        if order % d == 0:
            fixed = sum(k for k in part if d % k == 0)
            numerator += phi(order//d) * fixed * (fixed-1)//2
    require(numerator % order == 0, "pair-orbit average must be integral")
    return numerator//order


def direct_fixed_orientations(part):
    """Literally follow ordered pairs. No odd-cycle or gcd formula is used."""
    perm, start = [], 0
    for d in part:
        perm.extend(list(range(start+1,start+d))+[start])
        start += d
    unseen = set(combinations(range(start),2))
    orbits = 0
    while unseen:
        a,b = origin = min(unseen)
        while True:
            unseen.discard(tuple(sorted((a,b))))
            a,b = perm[a],perm[b]
            if (a,b) == origin[::-1]:
                return 0
            if (a,b) == origin:
                break
        orbits += 1
    return 2**orbits


def add(a,b,scale=1):
    result = a.copy()
    for part, value in b.items():
        result[part] = result.get(part,F(0)) + scale*value
        if not result[part]:
            del result[part]
    return result


def multiply(a,b,bound):
    result = {}
    for p,v in a.items():
        for q,w in b.items():
            if sum(p)+sum(q) <= bound:
                part = tuple(sorted(p+q))
                result[part] = result.get(part,F(0)) + v*w
    return {p:v for p,v in result.items() if v}


def prime_index(bound):
    B = {p:F(2**gcd_exponent(p),centralizer(p))
         for n in range(bound+1) for p in partitions(n)}
    A = add(B,{():F(1)},-1)
    # Invert multiplicatively by homogeneous recursion, independently of P.
    inverse = {():F(1)}
    for n in range(1,bound+1):
        term = multiply(A,inverse,n)
        inverse.update({p:-v for p,v in term.items() if sum(p)==n})
    residual = add(add({():F(1)},inverse,-1),{(1,):F(1)},-1)
    shifted = {d:{tuple(d*i for i in p):v for p,v in A.items() if d*sum(p)<=bound}
               for d in range(1,bound+1,2)}

    @lru_cache(None)
    def substituted_monomial(part):
        if not part:
            return {():F(1)}
        return multiply(substituted_monomial(part[:-1]),shifted[part[-1]],bound)

    P = {}
    for n in range(3,bound+1):
        homogeneous = {p:v for p,v in residual.items() if sum(p)==n}
        P.update(homogeneous)
        for part,value in homogeneous.items():
            residual = add(residual,substituted_monomial(part),-value)
    require(not residual,"full prime plethystic identity failed")
    require(all(v>0 for v in P.values()),"negative prime coefficient")
    for part in P:
        for i,d in enumerate(part):
            require(d==1 or len(part)==1 or any(gcd(d,e)>1 for j,e in enumerate(part) if j!=i),
                    "forbidden coprime-isolated prime cycle")
    return P


def scalar_product(a,b,bound):
    out = [F(0)]*(bound+1)
    for i,v in enumerate(a):
        for j,w in enumerate(b[:bound+1-i]):
            out[i+j] += v*w
    return out


def recompose(P,bound):
    """Compute t from P alone, with no tournament count input."""
    t, strong = [F(0)]*(bound+1), [F(0)]*(bound+1)
    t[0] = 1
    for n in range(1,bound+1):
        strong[n] = F(n==1)
        for part,v in P.items():
            if sum(part)>n:
                continue
            term = [v]+[F(0)]*n
            for d in part:
                child = [F(0)]*(n+1)
                for j in range(1,n//d+1):
                    child[d*j] = t[j]
                term = scalar_product(term,child,n)
            strong[n] += term[n]
        t[n] = sum(strong[j]*t[n-j] for j in range(1,n+1))
    require(all(v.denominator==1 for v in t),"composition gave fractional class count")
    return list(map(int,t))


def palette_separation_audit(P,bound):
    """Invert the universal-color observation matrix without large TT sizes."""
    @lru_cache(None)
    def coefficient(part,capacities):
        if not part:
            return int(not any(capacities))
        d=part[0]
        total=0
        for capacity,multiplicity in Counter(capacities).items():
            if capacity>=d:
                rest=list(capacities)
                rest.remove(capacity)
                rest.append(capacity-d)
                total+=multiplicity*coefficient(part[1:],tuple(sorted(rest)))
        return total
    dimension=0
    for n in range(1,bound+1):
        parts=sorted(partitions(n),key=lambda p:(-len(p),p))
        matrix={(mu,lam):coefficient(lam,mu) for mu in parts for lam in parts}
        recovered={}
        for mu in parts:
            require(all(not matrix[mu,lam] for lam in parts
                        if len(lam)<len(mu) or (len(lam)==len(mu) and lam!=mu)),
                    "palette observation matrix not triangular")
            require(matrix[mu,mu]>0,"palette matrix diagonal vanished")
            observed=sum(matrix[mu,lam]*P.get(lam,0) for lam in parts)
            recovered[mu]=F(observed-sum(matrix[mu,lam]*v for lam,v in recovered.items()))/matrix[mu,mu]
        require(all(recovered[p]==P.get(p,0) for p in parts),"palette recovery failed")
        dimension+=len(parts)
    return {"through":bound,"independent_coordinates":dimension,
            "check":"triangular observation matrices invert actual aggregate prime data"}


def hostile_pair():
    # A=1 -> C3; reverse only the first A in the second tournament.
    def tournament(reverse_first):
        rows = [[0]*12 for _ in range(12)]
        for u,v in combinations(range(12),2):
            i,a,j,b = u//4,u%4,v//4,v%4
            if i!=j:
                arc = (j-i)%3==1
            else:
                arc = a==0 or (a>0 and b>0 and (b-a)%3==1)
                if reverse_first and i==0:
                    arc = not arc
            rows[u][v],rows[v][u] = int(arc),int(not arc)
        return rows

    def automorphisms(rows):
        n=len(rows)
        degree=list(map(sum,rows))
        mapping={}
        used=set()
        found=[]
        def visit():
            if len(mapping)==n:
                found.append(tuple(mapping[i] for i in range(n)))
                return
            choices=[]
            for u in range(n):
                if u not in mapping:
                    possible=[v for v in range(n) if v not in used and degree[u]==degree[v]
                              and all(rows[u][a]==rows[v][b] for a,b in mapping.items())]
                    choices.append((len(possible),u,possible))
            _,u,possible=min(choices)
            for v in possible:
                mapping[u]=v
                used.add(v)
                visit()
                used.remove(v)
                del mapping[u]
        visit()
        return found

    def index(perms):
        counts=Counter()
        for p in perms:
            seen=set()
            lengths=[]
            for v in range(len(p)):
                if v not in seen:
                    w,length=v,0
                    while w not in seen:
                        seen.add(w)
                        length+=1
                        w=p[w]
                    lengths.append(length)
            counts[tuple(sorted(lengths))]+=1
        return {p:F(c,len(perms)) for p,c in counts.items()}

    first,second=tournament(False),tournament(True)
    triples=lambda rows:{t for t in combinations(range(12),3)
                         if all(sum(rows[v][w] for w in t)==1 for v in t)}
    require(triples(first)==triples(second),"hostile cyclic triples not equal")
    actions=[automorphisms(rows) for rows in (first,second)]
    require(list(map(len,actions))==[81,27],"hostile symmetry orders changed")
    Z=[index(action) for action in actions]
    require(Z[0]!=Z[1],"individual indices must distinguish hostile parents")
    return {"cyclic_triples":len(triples(first)),"automorphism_orders":list(map(len,actions)),
            "cycle_indices":[{','.join(map(str,p)):str(v) for p,v in sorted(z.items())} for z in Z]}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--count",type=int,default=50)
    parser.add_argument("--prime",type=int,default=20)
    args=parser.parse_args()
    require(3<=args.prime<=args.count,"require 3 <= prime bound <= count bound")
    counts=[]
    odd_states=0
    exponent_checks=0
    for n in range(args.count+1):
        total=F(0)
        for part in partitions(n):
            e=gcd_exponent(part)
            require(e==cyclic_burnside_exponent(part),"independent cyclic Burnside discrepancy")
            total+=F(2**e,centralizer(part))
            odd_states+=1
            exponent_checks+=1
        require(total.denominator==1,"nonintegral unlabelled count")
        counts.append(int(total))
    direct_checks=0
    for n in range(min(14,args.count)+1):
        for part in partitions(n,odd=False):
            expected=0 if any(d%2==0 for d in part) else 2**gcd_exponent(part)
            require(direct_fixed_orientations(part)==expected,"literal pair-orbit audit failed")
            direct_checks+=1
    # Audit only after construction. Never feed known counts to the recurrence.
    covered=min(args.count,len(OEIS_A000568_0_50)-1)
    require(all(counts[n]==OEIS_A000568_0_50[n] for n in range(covered+1)),"OEIS reference disagreement")
    P=prime_index(args.prime)
    recomposed=recompose(P,args.prime)
    require(recomposed==counts[:args.prime+1],"canonical composition round trip failed")
    prime_counts=[sum(v for p,v in P.items() if sum(p)==n) for n in range(args.prime+1)]
    require(all(v.denominator==1 for v in prime_counts),"nonintegral prime count")
    prior=[0,0,0,1,0,3,15,197,4008,136655,7906161,796876429,142980435791,
           46453867277461,27687212017612715,30571668704134867685]
    require(prime_counts[:min(16,args.prime+1)]==prior[:min(16,args.prime+1)],"Session79 prime mismatch")
    report={"status":"PASS","source_pin":"5ed1835235b36a57b8cd11eb705aa259b6bc0f08",
            "counts_through":args.count,"A000568":counts,"OEIS_verified_through":covered,
            "odd_partition_states":odd_states,"independent_cyclic_exponent_checks":exponent_checks,
            "literal_pair_orbit_checks":direct_checks,"literal_checks_through":min(14,args.count),
            "prime_recovery_and_scalar_recomposition_through":args.prime,
            "prime_counts":list(map(int,prime_counts)),"prime_index_monomials":len(P),
            "aggregate_prime_cycle_index":{','.join(map(str,p)):str(v) for p,v in sorted(P.items(),key=lambda t:(sum(t[0]),t[0]))},
            "hostile_81_27":hostile_pair(),
            "universal_palette_separation":palette_separation_audit(P,min(12,args.prime)),
            "limits":["Prime recovery/recomposition is a Burnside round trip, not independent discovery.",
                      "No prime representatives, individual group atlas, or metagraph edges are reconstructed.",
                      "Explicit literal orientation audit is per permutation cycle type, not tournament-class census."]}
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
