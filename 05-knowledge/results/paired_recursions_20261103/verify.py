"""Exact recurrences, converse-compatible pairings, and scoped hostile controls.
Run from any directory: python -X utf8 verify.py [--output receipt.json].
No assertion statements are used for validation.
"""
from pathlib import Path
from fractions import Fraction
from collections import Counter, defaultdict
from itertools import combinations, permutations
from math import gcd, factorial, prod
import argparse, importlib.util, json, sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def load_module(name,rel):
    spec=importlib.util.spec_from_file_location(name,ROOT/rel)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
old=load_module("precision_inherited","05-knowledge/results/port_residual_integral_20261025/precision/check.py")
ports=load_module("ports_inherited","05-knowledge/results/primegame_port_spectrum_20261103/port_scout.py")
def need(x,msg):
    if not x: raise RuntimeError(msg)
def v2(x): return (x & -x).bit_length()-1
def phi(n): return sum(gcd(i,n)==1 for i in range(1,n+1))
def channel(n,r=1,s=0):
    value=Fraction(0)
    for parts in old.all_odd_partitions(n):
        e,z,_=old.term(parts)
        value+=Fraction(1 << (r*e+s*len(parts)),z)
    need(value.denominator==1,("integrality",n,r,s))
    return int(value)
def coprime_pair_sum(n,bits=4):
    mod=1<<bits
    return sum(pow(a*(n-a),-1,mod) for a in range(1,(n+1)//2,2)
               if gcd(a,n)==1)%mod
def distinct_prime_factors(n):
    ans=[];p=2
    while p*p<=n:
        if n%p==0:
            ans.append(p)
            while n%p==0:n//=p
        p=3 if p==2 else p+2
    if n>1:ans.append(n)
    return ans
def phi_factored(n):
    ans=n
    for p in distinct_prime_factors(n):ans=ans//p*(p-1)
    return ans
def crt_pair_sum(n,bits=4):
    need(n>=4 and n%2==0 and bits>=1,'CRT domain')
    modulus=1<<bits;limit=n//2-1;divisors=[(1,1)]
    for p in distinct_prime_factors(n):
        if p!=2:divisors += [(d*p,-mu) for d,mu in divisors[:]]
    answer=0
    for residue in range(1,modulus,2):
        count=0
        for divisor,mu in divisors:
            first=residue*pow(divisor,-1,modulus)%modulus
            bound=limit//divisor
            count+=mu*max(0,1+(bound-first)//modulus)
        answer+=count*pow(residue*(n-residue),-1,modulus)
    return answer%modulus
def odd_four_bits(n):
    if n==3:return 1
    need(n>=5 and n%2,("odd input",n))
    ans=pow(n,-1,16)+2*pow(n-2,-1,16)
    triple_count=0
    for a in range(1,n,2):
        for b in range(a+2,n,2):
            c=n-a-b
            if c<=b:break
            if gcd(a,b)==gcd(a,c)==gcd(b,c)==1:
                ans+=4*pow(a*b*c,-1,16);triple_count+=1
    ans+=8*(n>=7 and n%3!=0)
    return ans%16
def even_four_bits(n):
    need(n>=4 and n%2==0,("even input",n))
    ans=coprime_pair_sum(n,4)
    if n>=12 and n%6==0:ans+=4*coprime_pair_sum(n//3,2)
    if n>=6:ans+=8*(phi(n-2)//2)
    ans+=6*(n==4)+2*(n==6)+8*(n==10)
    return ans%16
def even_four_bits_crt(n):
    need(n>=4 and n%2==0,'even CRT input')
    ans=crt_pair_sum(n,4)
    if n>=12 and n%6==0:ans+=4*crt_pair_sum(n//3,2)
    if n>=6:ans+=8*(phi_factored(n-2)//2)
    ans+=6*(n==4)+2*(n==6)+8*(n==10)
    return ans%16
def anti_pair_orbit_count(parts):
    # Independent direct permutation action, including endpoint swap sign.
    n=sum(parts);p=[];base=0
    for k in parts:
        p.extend(base+(j+1)%k for j in range(k));base+=k
    unseen=set(combinations(range(n),2));cycles=0
    while unseen:
        start=min(unseen);i,j=start;length=0
        while True:
            unseen.discard(tuple(sorted((i,j))))
            i,j=p[i],p[j];length+=1
            if {i,j}==set(start):break
        swap=(i,j)!=start
        if swap != bool(length%2):return 0
        cycles+=1
    return 1<<cycles
def all_parts(n,least=1):
    if not n:yield ()
    for a in range(least,n+1):
        for tail in all_parts(n-a,a):yield (a,)+tail
def direct_converse_count(n):
    ans=Fraction(0)
    for parts in all_parts(n):
        z=prod(a**m*factorial(m) for a,m in Counter(parts).items())
        ans+=Fraction(anti_pair_orbit_count(parts),z)
    need(ans.denominator==1,("anti-integrality",n));return int(ans)
def scc(a):
    n=len(a);reach=[r|(1<<i) for i,r in enumerate(a)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    return all(r==(1<<n)-1 for r in reach)
def full_class_census(n):
    ps=list(permutations(range(n))); edges=list(combinations(range(n),2))
    code_to_class={};entries=[]
    for mask in range(1<<len(edges)):
        if mask in code_to_class:continue
        a=ports.rows(n,mask);cid=len(entries)
        for p in ps:
            code=sum(((a[p[i]]>>p[j])&1)<<k for k,(i,j) in enumerate(edges))
            code_to_class[code]=cid
        entries.append({"mask":mask,"H":ports.h(a),"strong":scc(a)})
    for cid,e in enumerate(entries):
        e["converse"]=code_to_class[((1<<len(edges))-1)^e["mask"]]
        e["self_converse"]=e["converse"]==cid
    return entries,code_to_class
def construct_free_action(entries,d):
    # No metagraph-locality or H preservation is imposed.
    block=1<<d;fixed=[i for i,e in enumerate(entries) if e["converse"]==i]
    pairs=[(i,e["converse"]) for i,e in enumerate(entries) if i<e["converse"]]
    need(len(fixed)%block==0 and 2*len(pairs)%block==0,"action divisibility")
    blocks=[fixed[i:i+block] for i in range(0,len(fixed),block)]
    if d:
        for j in range(0,len(pairs),block//2):
            blocks.append([x for pair in pairs[j:j+block//2] for x in pair])
    else:blocks.extend([[x] for pair in pairs for x in pair])
    generators=[]
    for bit in range(d):
        g=list(range(len(entries)))
        for b in blocks:
            for i,x in enumerate(b):g[x]=b[i^(1<<bit)]
        generators.append(g)
    R=[e["converse"] for e in entries]
    for i in range(len(entries)):
        orbit=set()
        for word in range(block):
            j=i
            for bit,g in enumerate(generators):
                if word>>bit&1:j=g[j]
            orbit.add(j)
        need(len(orbit)==block,"action not free")
        for g in generators:
            need(g[g[i]]==i and g[R[i]]==R[g[i]],"not commuting involutions")
            for h in generators:need(g[h[i]]==h[g[i]],"generators do not commute")
    return generators
def weighted_recurrence(layers):
    # In the coefficient ring [h][k]=[hk].
    A=[Counter({1:1})];S=[Counter()];C=[Counter({1:1})];U=[Counter()]
    for entries in layers:
        A.append(Counter(e["H"] for e in entries))
        S.append(Counter(e["H"] for e in entries if e["strong"]))
        C.append(Counter(e["H"] for e in entries if e["self_converse"]))
        U.append(Counter(e["H"] for e in entries if e["strong"] and e["self_converse"]))
    for n in range(1,len(A)):
        predicted=Counter()
        for k in range(1,n+1):
            for h,ct in S[k].items():
                for t,dt in A[n-k].items():predicted[h*t]+=ct*dt
        need(predicted==A[n],("SCC recurrence",n))
        predicted=Counter()
        for m in range(n//2+1):
            center=U[n-2*m] if n>2*m else Counter({1:1})
            for h,ct in A[m].items():
                for t,dt in center.items():predicted[h*h*t]+=ct*dt
        need(predicted==C[n],("converse SCC recurrence",n))
    return A,S,C,U
def native_matching_probe(code_map):
    n=5;edges=list(combinations(range(n),2))
    match=[edges.index((0,1)),edges.index((2,3))]
    other=[i for i in range(len(edges)) if i not in match]
    blocks=[];frame_witness={}
    for f in range(1<<len(other)):
        base=sum(((f>>j)&1)<<k for j,k in enumerate(other))
        b=frozenset(code_map[base|sum(((v>>j)&1)<<k for j,k in enumerate(match))] for v in range(4))
        if len(b)==4:
            blocks.append(b);frame_witness.setdefault(b,base)
    unique=sorted(set(blocks),key=lambda b:tuple(sorted(b)))
    tilings=sum(len(a|b|c)==12 for a,b,c in combinations(unique,3))
    minimum_cover=next(k for k in range(1,len(unique)+1)
                       if any(len(set().union(*bs))==12 for bs in combinations(unique,k)))
    cover=next(bs for bs in combinations(unique,minimum_cover) if len(set().union(*bs))==12)
    ids_to_mask={}
    for mask,cid in code_map.items():ids_to_mask[cid]=min(mask,ids_to_mask.get(cid,mask))
    return {"full_frames":len(blocks),"unique_four_class_blocks":len(unique),
            "disjoint_three_block_tilings":tilings,"minimum_full_block_cover":minimum_cover,
            "four_frame_witness":[{"base_mask":frame_witness[b],
              "class_masks":sorted(ids_to_mask[i] for i in b)} for b in cover]}
def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output",default="receipt.json")
    args=parser.parse_args();checks=0
    exact=[channel(n) for n in range(41)]
    known=[1,1,1,2,4,12,56,456,6880,191536,9733056,903753248]
    need(exact[:len(known)]==known,"external OEIS prefix")
    for n in range(3,201):
        expected=(odd_four_bits(n) if n%2 else even_four_bits(n))
        inherited,_=old.normalized_residue(n,4)
        need(expected==inherited,("four bits",n,expected,inherited));checks+=1
        if n<=40:need(exact[n]//(1<<(n//2))%16==expected,("independent full sum",n))
    crt_cases=list(range(4,1001,2))+[20402,338,5510,30030,1000000,1000002]
    for n in crt_cases:
        need(even_four_bits_crt(n)==even_four_bits(n),('CRT formula',n))
        for bits in (1,2,3,4,5):
            need(crt_pair_sum(n,bits)==coprime_pair_sum(n,bits),('CRT weights',n,bits))
    counts=[]
    for n in range(41):
        f=channel(n//2,2,n%2);a=exact[n]
        if n<=12:need(f==direct_converse_count(n),("twisted direct",n))
        d=min(v2(a),v2(f))
        if n>=3 and n%2:need(d==(n-1)//2,("odd axes",n))
        if n%4==2:need(d==n//2-1,("2 mod 4 axes",n))
        counts.append({"n":n,"a":a,"self_converse":f,"max_free_commuting_rank":d})
    layers=[];fine=[];native=None
    for n in range(1,7):
        es,cmap=full_class_census(n);layers.append(es)
        need(len(es)==exact[n],"class census")
        need(sum(e["self_converse"] for e in es)==counts[n]["self_converse"],"self census")
        need(sum(e["H"]==1 for e in es)==1,"unique H1")
        gs=construct_free_action(es,counts[n]["max_free_commuting_rank"])
        fixed=Counter(e["H"] for e in es if e["self_converse"])
        # H-preserving involution: converse on chiral pairs, pair fixed classes by H.
        J=[e["converse"] for e in es]
        for h in fixed:
            ids=[i for i,e in enumerate(es) if e["H"]==h and e["self_converse"]]
            for i in range(0,len(ids)-1,2):J[ids[i]],J[ids[i+1]]=ids[i+1],ids[i]
        need(sum(J[i]==i for i in range(len(es)))==sum(v%2 for v in fixed.values()),"minimum defect")
        for i,e in enumerate(es):
            need(J[J[i]]==i and es[J[i]]["H"]==e["H"],"H pairing")
            need(J[e["converse"]]==es[J[i]]["converse"],"H pairing commutes")
        fine.append({"n":n,"minimum_H_preserving_fixed_classes":sum(v%2 for v in fixed.values()),
                     "odd_self_H_fibres":[h for h,v in sorted(fixed.items()) if v%2],
                     "entries":es,"generators":gs,"H_preserving_pairing":J})
        if n==5:native=native_matching_probe(cmap)
    A,S,C,U=weighted_recurrence(layers)
    # Independent literal Hamilton-path counts at n <= 5.
    literal=0
    for n,es in enumerate(layers[:5],1):
        for e in es:
            rows=ports.rows(n,e["mask"])
            brute=sum(all(rows[p[i]]>>p[i+1]&1 for i in range(n-1)) for p in permutations(range(n)))
            need(brute==e["H"],"literal path enumeration");literal+=1
    need(native["full_frames"]==136 and native["disjoint_three_block_tilings"]==0,"native matching control")
    fails=[]
    for n in range(4,77,2):
        residue,_=old.normalized_residue(n,12)
        need(residue!=0,("valuation unresolved",n))
        actual=v2(residue);guess=v2(phi(n))-1
        if actual!=guess:fails.append({"n":n,"relative_valuation":actual,"totient_guess":guess})
    need(fails[0]["n"]==40,"smallest false stronger valuation")
    for row in fails:
        actual=channel(row['n'])
        need(v2(actual)==row['n']//2+row['relative_valuation'],'independent valuation hostile')
    sys.path.insert(0,str(ROOT/'05-knowledge/results/primegame_port_spectrum_20261103'))
    hostile_module=load_module('growth_hostile','05-knowledge/results/primegame_port_spectrum_20261103/verify.py')
    hostile=hostile_module.hostile()
    result={"status":"PASS","scope":{"exact_Burnside":40,"independent_pair_orbits":12,
             "four_binary_digit_law":200,"all_labelled_class_orbits":6,
             "literal_path_order_check_classes":literal},
            "counts":counts,"four_bit_checks":checks,"CRT_cases":crt_cases,
            "class_pairings":fine,"hostile_81_27":hostile,
            "native_matching_n5":native,"stronger_totient_valuation_failures_through76":fails,
            "provenance":"Inherited partition/precision and port code imported by exact repository path; source hashes in manifest."}
    (HERE/args.output).write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"PASS","checks":checks,"native":native,"failures":fails,
                      "H_defects":[(x["n"],x["minimum_H_preserving_fixed_classes"]) for x in fine]}))
if __name__=="__main__":main()

