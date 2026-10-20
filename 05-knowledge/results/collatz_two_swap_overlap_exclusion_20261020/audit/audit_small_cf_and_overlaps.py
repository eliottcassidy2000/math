#!/usr/bin/env python3
"""Independent exact-log and actual-edge audit of the small overlap residual.

No producer implementations are imported. Python assertions are not used.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import json

checks = 0

def need(test, why):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(why)


def log_interval(t, terms=80):
    # log((1+t)/(1-t)); positive tail bounded geometrically.
    power=t
    lower=F(0)
    for j in range(terms):
        lower += 2*power/F(2*j+1)
        power *= t*t
    return lower, lower+2*power/F((2*terms+1)) / (1-t*t)


def cf_prefix(lower, upper):
    coeffs=[]
    pairs=[]
    a0,a1,b0,b1=0,1,1,0
    while True:
        a=lower.numerator//lower.denominator
        need(a == upper.numerator//upper.denominator, 'ambiguous exact CF digit')
        coeffs.append(a)
        a0,a1=a1,a*a1+a0
        b0,b1=b1,a*b1+b0
        pairs.append((a1,b1))
        if b1>=4096:
            break
        need(lower>a, 'cannot invert zero interval endpoint')
        lower,upper=1/(upper-a),1/(lower-a)
    return coeffs,pairs


def carry(bits):
    # Reverse scan independent of producer's forward affine composition.
    value=0
    multiplier=1
    for bit in reversed(bits):
        value=2*value + (multiplier if bit else 0)
        if bit:
            multiplier*=3
    return value


def cycle_data(values, sign):
    need(len(set(values)) == len(values), 'known cycle is not primitive')
    for i,n in enumerate(values):
        next_n = n//2 if n%2==0 else (3*n+sign)//2
        need(next_n == values[(i+1)%len(values)], 'known cycle edge wrong')
    return {'p':sum(v%2 for v in values), 'q':len(values), 'minimum':min(values), 'sign':sign}


def actual_paths(p,q):
    original=[((i+1)*p)//q-(i*p)//q for i in range(q)]
    labels=[i if bit else None for i,bit in enumerate(original)]
    for edge in range(q):
        j=(edge+1)%q
        if original[edge] == original[j]:
            continue
        middle=labels[:]
        first_label = middle[edge] if original[edge] else middle[j]
        direction=1 if original[edge] else -1
        middle[edge],middle[j]=middle[j],middle[edge]
        for next_edge in sorted({(edge-1)%q,(edge+1)%q}):
            k=(next_edge+1)%q
            if (middle[next_edge] is None) == (middle[k] is None):
                continue
            second_label=middle[next_edge] if middle[next_edge] is not None else middle[k]
            if second_label == first_label:
                continue
            second_direction=1 if middle[next_edge] is not None else -1
            need(direction == second_direction, 'distinct-one overlap directions differ')
            final=middle[:]
            final[next_edge],final[k]=final[k],final[next_edge]
            bits=[int(x is not None) for x in final]
            need(sum(bits)==p and len(bits)==q, 'shape changed')
            yield direction,edge,next_edge,bits,first_label


def main():
    l2,u2=log_interval(F(1,3))
    # Use log3 = log(3/2)+log2, instead of producer's log3 series.
    l32,u32=log_interval(F(1,5))
    l3,u3=l2+l32,u2+u32
    need(l2>F(2,3), 'log2 lower bound')
    lower,upper=l3/u2,u3/l2
    coeffs,pairs=cf_prefix(lower,upper)
    need(pairs[-1] == (24727,15601), 'first out-of-range denominator')
    candidates=[(p,q) for q,p in pairs if 2<=p<4096 and p<q<2*p]
    need(candidates == [(2,3),(5,8),(12,19),(41,65),(53,84),(306,485),(665,1054)], 'candidate list')
    cycles=[cycle_data([1,2],1),cycle_data([1],-1),cycle_data([5,7,10],-1),
            cycle_data([17,25,37,55,82,41,61,91,136,68,34],-1)]
    extras=sorted({(c['p'],c['q']) for c in cycles if c['p']>=2 and c['p']<c['q']<2*c['p']} - set(candidates))
    need(extras==[(7,11)], 'known-cycle shape coverage')
    need((4095**2)<(1<<24) == 4096**2, 'Legendre strict threshold')
    rows=[]
    hits=[]
    digest=sha256()
    all_paths=0
    crossing_paths=0
    producer_compatible_records=[]
    for p,q in sorted(candidates+extras):
        r=q-p
        need(gcd(p,q)==1,'primitive shape')
        D=(1<<q)-3**p
        directions={1:0,-1:0}
        min_nonintegral_denom=None
        for direction,e1,e2,bits,first_label in actual_paths(p,q):
            directions[direction]+=1
            all_paths+=1
            crossing_paths += (e1==q-1 or e2==q-1)
            B=carry(bits)
            content=gcd(B,abs(D))
            denominator=abs(D)//content
            first_residue=(-((first_label*p)//q+1)*r)%p
            u=first_residue if direction==1 else first_residue-r
            family='plus' if direction==1 else 'minus'
            producer_compatible_records.append({'p':p,'q':q,'r':r,'family':family,
                'u':u,'edges':[e1,e2],'B':B,'D':D,'content':content,
                'word':''.join(map(str,bits))})
            # Positional carry independently checks the reverse scan.
            after=0
            positional=0
            for i in range(q-1,-1,-1):
                if bits[i]:
                    positional+=(1<<i)*3**after
                    after+=1
            need(B==positional,'independent carry mismatch')
            if denominator==1:
                n=B//D
                orbit=[]
                for bit in bits:
                    need(n%2==bit,'integral output has wrong itinerary')
                    orbit.append(n)
                    n=n//2 if n%2==0 else (3*n+1)//2
                need(n==B//D,'integer orbit failed closure')
                hits.append({'p':p,'q':q,'direction':direction,'edges':[e1,e2],
                             'word':''.join(map(str,bits)),'B':B,'D':D,
                             'fixed_point':B//D,'orbit':orbit})
            elif min_nonintegral_denom is None or denominator<min_nonintegral_denom:
                min_nonintegral_denom=denominator
            encoded=[p,q,direction,e1,e2,''.join(map(str,bits)),str(B),str(D),str(content)]
            digest.update((json.dumps(encoded,separators=(',',':'))+'\n').encode())
        expected=min(r,p-r)
        need(directions == {1:expected,-1:expected},'actual-edge versus residue interval count')
        rows.append({'p':p,'q':q,'r':r,'D_sign':1 if D>0 else -1,
                     'plus_paths':directions[1],'minus_paths':directions[-1],
                     'minimum_nonintegral_denominator':min_nonintegral_denom})
    need(all_paths==906,'complete candidate path count')
    need(sorted((h['p'],h['q'],h['direction'],h['fixed_point']) for h in hits)
         == [(2,3,-1,-5),(2,3,1,-7)],'negative exceptions changed')
    producer_compatible_records.sort(key=lambda z:(z['p'],0 if z['family']=='plus' else 1,z['u']))
    producer_digest=sha256()
    for row in producer_compatible_records:
        producer_digest.update(json.dumps(row,sort_keys=True,separators=(',',':')).encode()+b'\n')
    producer=Path(__file__).resolve().parent.parent
    receipt_path=producer/'candidates.receipt.json'
    receipt=json.loads(receipt_path.read_text())
    need(receipt_path.read_bytes()==(producer/'candidates.receipt.optimized.json').read_bytes(), 'producer normal/-O mismatch')
    need(sha256((producer/'verify_small_convergent_overlaps.py').read_bytes()).hexdigest()==receipt['script_sha256'], 'producer source hash')
    need(producer_digest.hexdigest()==receipt['all_case_trace_sha256'], 'independent exact full case trace mismatch')
    need(receipt['total_paths']==all_paths,'producer path count mismatch')
    result={'status':'FINITE-EXACT independent signed small-CF and actual-edge audit',
            'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'checks':checks,'log_decomposition':'log3=log2+log(3/2), 80 terms each',
            'continued_fraction_prefix':coeffs,'convergents_through_first_excess':pairs,
            'convergent_shapes_p_q':candidates,'extra_known_cycle_shapes_p_q':extras,
            'known_cycles':cycles,'rows':rows,'path_count':all_paths,
            'original_cut_crossing_paths':crossing_paths,'integral_outputs':hits,
            'exact_trace_sha256':digest.hexdigest(),
            'producer_compatible_trace_sha256':producer_digest.hexdigest(),
            'producer_script_sha256':receipt['script_sha256']}
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
