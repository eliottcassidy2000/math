#!/usr/bin/env python3
"""Exact observer spectra and native continuation probabilities. Standard library only."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, gcd
from pathlib import Path
import hashlib
import json
import sys

CHECKS = 0

def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)

@lru_cache(None)
def transforms(n):
    ee = tuple(combinations(range(n), 2))
    pos = {e: k for k, e in enumerate(ee)}
    return tuple((tuple(pos[tuple(sorted((p[a], p[b])))] for a, b in ee),
                  tuple(int(p[a] > p[b]) for a, b in ee)) for p in permutations(range(n)))

def images(n, mask, affine=False):
    return [sum((((mask >> k) & 1) ^ (eps[k] if affine else 0)) << image[k]
                for k in range(len(image))) for image, eps in transforms(n)]

@lru_cache(None)
def canon(n, mask, affine=False):
    return min(images(n, mask, affine))

def actual_arc(n, mask, a, b):
    ee = tuple(combinations(range(n), 2))
    bit = (mask >> ee.index(tuple(sorted((a, b))))) & 1
    return (a < b) != bool(bit)

def native_records(n):
    ce = tuple(combinations(range(n-1), 2))
    ne = tuple(combinations(range(n), 2))
    records = []
    for code in range(1 << len(ce)):
        native = sum(((code >> k) & 1) << ne.index((a, b+1)) for k, (a, b) in enumerate(ce))
        need(all(actual_arc(n, native, i, i+1) for i in range(n-1)), 'Actual forced path')
        records.append((code, native, canon(n, native), canon(n-1, code), canon(n, native, True)))
    return records

def incidence(records, source):
    labels = sorted({tuple(row[i] for i in source) for row in records})
    targets = sorted({row[-1] for row in records})
    w = [[F(0) for _ in targets] for _ in labels]
    for row in records:
        w[labels.index(tuple(row[i] for i in source))][targets.index(row[-1])] += 1
    return labels, targets, w

def transpose(a):
    return [list(row) for row in zip(*a)]

def mm(a, b):
    bt = transpose(b)
    return [[sum((x*y for x, y in zip(row, col)), F(0)) for col in bt] for row in a]

def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]

def coupling(w):
    r = [sum(row) for row in w]
    c = [sum(col) for col in transpose(w)]
    need(all(r) and all(c), 'Strictly positive incidence marginals')
    k = [[sum(w[i][t]*w[j][t]/(r[i]*c[t]) for t in range(len(c)))
          for j in range(len(r))] for i in range(len(r))]
    need(all(sum(row) == 1 for row in k), 'Stochastic coupling')
    need(all(r[i]*k[i][j] == r[j]*k[j][i] for i in range(len(r)) for j in range(len(r))),
         'Weighted reversibility')
    return r, c, k

def rank(a):
    a = [row[:] for row in a]
    row = 0
    for col in range(len(a[0])):
        pivot = next((j for j in range(row, len(a)) if a[j][col]), None)
        if pivot is None:
            continue
        a[pivot], a[row] = a[row], a[pivot]
        scale = a[row][col]
        a[row] = [x/scale for x in a[row]]
        for j in range(len(a)):
            if j != row and a[j][col]:
                scale = a[j][col]
                a[j] = [x-scale*y for x, y in zip(a[j], a[row])]
        row += 1
        if row == len(a):
            break
    return row

def charpoly(a):
    """Newton trace identities, returned in descending powers."""
    power = identity(len(a))
    traces = [F(0)]
    coeff = [F(1)]
    for k in range(1, len(a)+1):
        power = mm(power, a)
        traces.append(sum(power[i][i] for i in range(len(a))))
        coeff.append(-sum(coeff[k-i]*traces[i] for i in range(1, k+1))/k)
    return coeff

def apply(a, v):
    return [sum(x*y for x, y in zip(row, v)) for row in a]

def expected_char(factor, zeros):
    f = [F(x, factor[0]) for x in factor]
    out = [F(0)]*(len(f)+1)
    for i, x in enumerate(f):
        out[i] += x
        out[i+1] -= x
    return out+[F(0)]*zeros

FULL = [34479695250,-101537044560,124233482025,-81343863297,30498542418,-6453737759,695751836,-27925620]
COMP = [3377602800,-8397030240,8429380020,-4360639860,1225690398,-180039440,12012273,-267391]
PAIR = [562933800,-914443020,553227084,-153488811,19184700,-874496]

# Exact polynomial arithmetic uses increasing powers.
def trim(a):
    a = list(a)
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a

def divrem(a, b, modulus=None):
    a, b = trim(a), trim(b)
    if b == [0]:
        raise ZeroDivisionError
    out = [0]*max(1, len(a)-len(b)+1)
    while a != [0] and len(a) >= len(b):
        i = len(a)-len(b)
        c = a[-1]*pow(b[-1], -1, modulus) % modulus if modulus else F(a[-1], b[-1])
        out[i] = c
        for j in range(len(b)):
            a[i+j] -= c*b[j]
            if modulus:
                a[i+j] %= modulus
        a = trim(a)
    return trim(out), a

def poly_gcd(a, b, p):
    while trim(b) != [0]:
        a, b = b, divrem(a, b, p)[1]
    scale = pow(a[-1], -1, p)
    return [(v*scale) % p for v in a]

def mulmod(a, b, f, p):
    c = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] = (c[i+j]+x*y) % p
    return divrem(c, f, p)[1]

def powmod(a, k, f, p):
    out = [1]
    while k:
        if k & 1:
            out = mulmod(out, a, f, p)
        a = mulmod(a, a, f, p)
        k //= 2
    return out

def irreducible_certificate(descending, p):
    degree = len(descending)-1
    need(degree in (5, 7), 'Rabin specialization uses a prime degree')
    lead = descending[0] % p
    need(lead != 0, 'Degree retained modulo prime')
    f = [(c*pow(lead,-1,p)) % p for c in reversed(descending)]
    xp = powmod([0,1], p, f, p)
    minus_x = xp[:]+[0]*max(0, 2-len(xp))
    minus_x[1] = (minus_x[1]-1) % p
    g = poly_gcd(f, trim(minus_x), p)
    need(g == [1], 'No degree-one finite-field factor')
    z = [0,1]
    frobenius = []
    for _ in range(degree):
        z = powmod(z, p, f, p)
        frobenius.append(z)
    need(z == [0,1], 'Every irreducible factor degree divides prime degree')
    return {'prime':p, 'monic_mod_p_ascending':f, 'gcd_with_xp_minus_x':g,
            'frobenius_remainders':frobenius}

def sturm(descending):
    f = [F(x) for x in reversed(descending)]
    chain = [f, [i*f[i] for i in range(1,len(f))]]
    while chain[-1] != [0]:
        rem = divrem(chain[-2], chain[-1])[1]
        if rem == [0]:
            break
        chain.append([-x for x in rem])
    return chain

def evaluate(poly, x):
    out = F(0)
    for a in reversed(poly):
        out = out*x+a
    return out

def variations(chain, x):
    signs = [1 if y > 0 else -1 for p in chain if (y:=evaluate(p,x))]
    return sum(a != b for a,b in zip(signs,signs[1:]))

def root_certificate(f, lower, upper):
    chain = sturm(f)
    v = [variations(chain,x) for x in (F(0),lower,upper,F(1))]
    need(v[0]-v[3] == len(f)-1, 'All nonconstant roots lie in (0,1)')
    need(v[1]-v[2] == 1 and v[2]-v[3] == 0, 'Unique largest root in exact bracket')
    need(all(evaluate(factor, x) != 0 for factor in (list(reversed(f)),) for x in (lower,upper)),
         'Isolation endpoints are not roots')
    return {'largest_root_interval':[str(lower),str(upper)],'sturm_variations_0_low_high_1':v,
            'sturm_chain_sha256':hashlib.sha256(repr(chain).encode()).hexdigest()}

def serialization(a):
    return [[str(x) for x in row] for row in a]

def trace_defect(w, labels, coordinate):
    r,c,k = coupling(w)
    blocks = sorted({label[coordinate] for label in labels})
    ww = [[sum(w[i][t] for i,label in enumerate(labels) if label[coordinate] == block)
           for t in range(len(c))] for block in blocks]
    rr,_,kk = coupling(ww)
    variance = F(0)
    for i,label in enumerate(labels):
        b = blocks.index(label[coordinate])
        variance += sum(r[i]*(w[i][t]/r[i]-ww[b][t]/rr[b])**2/c[t] for t in range(len(c)))
    direct = sum(k[i][i] for i in range(len(r)))-sum(kk[i][i] for i in range(len(rr)))
    need(variance == direct, 'Exact trace defect equals profile variance')
    return variance

def check_converse(n, records, source):
    labels,targets,w = incidence(records,source)
    totalmask = (1 << comb(n,2))-1
    merged = []
    for row in records:
        row = list(row)
        row[-1] = min(row[-1],canon(n,row[1]^totalmask,True))
        merged.append(tuple(row))
    _,_,wm = incidence(merged,source)
    need(coupling(w)[2] == coupling(wm)[2], 'Native converse merging leaves source operator unchanged')
    for j,t in enumerate(targets):
        opposite = canon(n,t^totalmask,True)
        need([row[j] for row in w] == [row[targets.index(opposite)] for row in w],
             'Converse target columns exactly equal')


def spectral_controls():
    out = {}
    for n in (4,5):
        all_records = native_records(n)
        positive = all_records[1:]
        # Independent literal Hamiltonian orders and stabilizer normalization.
        target_counts = Counter(row[-1] for row in all_records)
        ht = []
        for target,mass in sorted(target_counts.items()):
            h = sum(all(actual_arc(n,target,a,b) for a,b in zip(p,p[1:])) for p in permutations(range(n)))
            aut = images(n,target,True).count(target)
            need(F(h,aut) == mass, 'Actual H/Aut native column mass')
            ht.append([target,h,aut,mass])
        cases = {}
        for name,source in (('full',(2,)),('compressed',(3,)),('refinement',(2,3))):
            labels,targets,w = incidence(positive,source)
            r,c,k = coupling(w)
            cp = charpoly(k)
            # Independently swap source and target, allowing dimension zero differences.
            _,_,dual = coupling(transpose(w))
            dcp = charpoly(dual)
            nonzero_cp = cp[:]
            while nonzero_cp[-1] == 0:
                nonzero_cp.pop()
            nonzero_dual = dcp[:]
            while nonzero_dual[-1] == 0:
                nonzero_dual.pop()
            need(nonzero_cp == nonzero_dual, 'Independent source/target characteristic agreement')
            cases[name] = {'labels':labels,'targets':targets,'incidence':serialization(w),
                           'masses':[str(x) for x in r],'rank':rank(w),
                           'trace':str(sum(k[i][i] for i in range(len(k)))),
                           'charpoly_descending':[str(x) for x in cp]}
            check_converse(n,positive,source)
            if n == 4 and name in ('full','compressed'):
                slow = [4]+[-3]*(len(r)-1)
                need(apply(k,slow) == [F(8,15)*x for x in slow], 'Explicit 8/15 eigenspace')
                zeros = ([0,1,-2,0],[0,0,1,-1]) if name=='full' else ([0,1,-3],)
                for v in zeros:
                    need(apply(k,v) == [0]*len(r), 'Explicit zero eigenspace')
                need(cp == expected_char([15,-8],len(r)-2), 'Complete n4 spectrum')
            if n == 5 and name in ('full','compressed'):
                factor = FULL if name=='full' else COMP
                need(cp == expected_char(factor,len(r)-8), 'Exact native degree-seven factor')
        labels,_,w = incidence(positive,(2,3))
        defects = [trace_defect(w,labels,i) for i in (0,1)]
        need(defects == ([F(0),F(0)] if n==4 else [F(599548,675675),F(8662,6435)]),
             'Full/compressed refinement defects')
        # Complete cube, compressed complement pairs; target converse already tested.
        paired=[]
        for row in all_records:
            row=list(row)
            row[3]=min(row[3],canon(n-1,row[0]^((1<<comb(n-1,2))-1)))
            paired.append(tuple(row))
        labels,targets,w = incidence(paired,(3,))
        r,c,k=coupling(w)
        cp=charpoly(k)
        need(cp == (expected_char([15,-7],0) if n==4 else expected_char(PAIR,0)),
             'Complement-paired compressed characteristic polynomial')
        check_converse(n,paired,(3,))
        if n==4:
            need(r==[2,6] and apply(k,[3,-1])==[F(7,15)*3,-F(7,15)], 'Paired mass and eigenspace')
        cases['paired_compressed']={'labels':labels,'targets':targets,'incidence':serialization(w),
                                    'masses':[str(x) for x in r],'charpoly_descending':[str(x) for x in cp]}
        out[str(n)]={'word_records':all_records,'target_H_Aut_mass':ht,'cases':cases,
                     'refinement_defects':[str(x) for x in defects]}
        if n==5:
            bycode={row[0]:row for row in positive}
            need(bycode[3][2]==bycode[33][2] and bycode[3][3]!=bycode[33][3], 'Full does not refine compressed')
            need(bycode[3][3]==bycode[9][3] and bycode[3][2]!=bycode[9][2], 'Compressed does not refine full')
            out[str(n)]['incomparability_witnesses']=[[bycode[3],bycode[33]],[bycode[3],bycode[9]]]
    certs={}
    for name,f,p,lo,hi in (('full',FULL,31,635789,635790),('compressed',COMP,23,642470,642471),
                          ('paired_compressed',PAIR,19,573377,573378)):
        certs[name]={'irreducibility':irreducible_certificate(f,p),
                     'root_isolation':root_certificate(f,F(lo,10**6),F(hi,10**6))}
    need(F(635790,10**6)<F(642470,10**6), 'Compressed native gap strictly smaller')
    out['algebraic_certificates']=certs
    return out


def native_host(t):
    def arc(x,y):
        i,u=x; j,v=y
        if i==j:
            return (v-u)%3==1
        if (j-i)%t <= (t-1)//2:
            return (v-u)%3 in (1,2)
        return (u-v)%3 not in (1,2)
    vertices=list(product(range(t),range(3)))
    need(all(sum(arc(x,y) for y in vertices if y!=x)==(3*t-1)//2 for x in vertices),
         'Literal host regularity')
    return arc


def transfer_controls():
    correction=(2,1,-1,-2,-1,1)
    counts=[1,0,0]
    j=[0,1]
    for r in range(2,83):
        j.append(j[-1]+2*j[-2])
    rows=[]
    for r in range(81):
        predicted=[(2**r+correction[(r-2*d)%6])//3 for d in range(3)]
        need(all((2**r+correction[(r-2*d)%6])%3==0 for d in range(3)), 'Exact division by three')
        need(counts==predicted, 'Integer transfer agrees with period-six formula')
        need(counts==[sum(comb(r,k) for k in range(d,r+1,3)) for d in range(3)],
             'Weight-refined binomial transfer identity')
        need(j[r]==(2**r-(-1)**r)//3, 'Jacobsthal recurrence and closed form')
        if r:
            need(gcd(j[r],j[r+1])==1, 'Consecutive Jacobsthal terms coprime')
            need(j[r+1]*j[r-1]-j[r]**2==(-1)**r*2**(r-1), 'Determinant-two Cassini identity')
        if r>=3:
            previous=[(2**(r-1)+correction[(r-1-2*d)%6])//3 for d in range(3)]
            observed=[]
            for d in range(3):
                a,b=previous[d],previous[(d-1)%3]
                need(a+b==counts[d] and a>0 and b>0, 'Actual conditional splitting counts')
                pair=tuple(sorted((F(a,a+b),F(b,a+b))))
                expected=(F(1,2),F(1,2)) if a==b else (F(j[r]-1,2*j[r]),F(j[r]+1,2*j[r]))
                need(pair==expected, 'Exact native Jacobsthal splitting fractions')
                need(F(7,15) not in pair and F(8,15) not in pair, 'No accidental spectral-pair probability')
                if a!=b:
                    p,q=pair[0].numerator,pair[0].denominator
                    need(q==j[r] and abs(2*p-q)==1, 'True Farey adjacency to half')
                    need(abs(pair[0].numerator*pair[1].denominator-pair[1].numerator*pair[0].denominator)==j[r],
                         'Complementary pair is not Farey adjacent')
                    k=(j[r]-1)//2
                    need(F(k,2*k+1)==pair[0] and (k-1)==(j[r]-3)//2, 'Exact Berggren parabolic address depth')
                observed.append([str(x) for x in pair])
            if r<=12:
                g=F(1,3)
                heap=0
                for _ in range((j[r]-3)//2):
                    g=(1-g)/(3-4*g)
                    heap=3*heap+1
                need(g==F(j[r]-1,2*j[r]), 'Independent Berggren branch replay')
                need(heap==(3**((j[r]-3)//2)-1)//2, 'Independent ternary heap replay')
                rows.append({'r':r,'A':counts[:],'J':j[r],'conditional_pairs':observed,
                             'berggren_depth':(j[r]-3)//2})
        counts=[counts[d]+counts[(d-1)%3] for d in range(3)]
    need((2**2+correction[2])%9==(2**6+correction[0])%9==3, 'Lossy mod9 numerators collide')
    need(((2**2+correction[2])//3)%9==1 and ((2**6+correction[0])//3)%9==4,
         'Counts do not collide modulo nine')
    for r in range(81):
        for d in range(3):
            numerator=2**r+correction[(r-2*d)%6]
            need((3*(numerator%8))%8==(numerator//3)%8, 'Division by three modulo eight')
            need(((numerator%27)//3)%9==(numerator//3)%9, 'Modulo twenty-seven repairs modulo-nine division')
    actual=[]
    for t in (3,5,7):
        arc=native_host(t); r=t-1
        cnt=Counter(); first=Counter()
        for alpha in range(3):
            for bits in product((0,1),repeat=r):
                phases=[alpha]
                for bit in bits:
                    phases.append((phases[-1]+bit)%3)
                path=[(i,(a+u)%3) for i,a in enumerate(phases) for u in range(3)]
                need(len(set(path))==3*t and all(arc(a,b) for a,b in zip(path,path[1:])),
                     'Actual retained owner path and every boundary arc')
                delta=(phases[-1]-alpha)%3
                cnt[alpha,delta]+=1
                first[alpha,delta,bits[0]]+=1
        for alpha in range(3):
            for d in range(3):
                a=(2**r+correction[(r-2*d)%6])//3
                need(cnt[alpha,d]==a, 'Literal pinned native count')
                for bit in (0,1):
                    child=(2**(r-1)+correction[(r-1-2*(d-bit))%6])//3
                    need(first[alpha,d,bit]==child, 'Literal native conditional numerator')
        actual.append({'cells':t,'order':3*t,'accepted_with_common_phase':sum(cnt.values()),
                       'fixed_start_endpoint_counts':[cnt[0,d] for d in range(3)]})
    # Independent polynomial determinant by all six permutations.
    # Entries are constant 1 on diagonal and z on cyclic successor.
    determinant=Counter()
    for p in permutations(range(3)):
        if all(p[i] in (i,(i+1)%3) for i in range(3)):
            exponent=sum(p[i]!=i for i in range(3))
            sign=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
            determinant[exponent]+=sign
    need(dict(determinant)=={0:1,3:1}, 'Native polynomial transfer determinant 1+z^3')
    return {'formula_range':[0,80],'small_rows':rows,'literal_hosts':actual,
            'transfer_determinant':dict(determinant),'mod9_hostile':{'numerators':[3,66],'counts':[1,22]}}


def main():
    out={'status':'PASS','spectra':spectral_controls(),'native_transfer':transfer_controls(),'checks':CHECKS}
    target=Path(__file__).with_name('receipt.json')
    text=json.dumps(out,indent=2,sort_keys=True)+'\n'
    if '--write' in sys.argv:
        target.write_text(text,encoding='utf-8')
    else:
        need(target.read_text(encoding='utf-8')==text,'Frozen exact receipt')
    print(json.dumps({'status':'PASS','checks':CHECKS,'normalization':'exact integers/Fractions only',
                      'full_degree':7,'compressed_degree':7,'paired_degree':5},sort_keys=True))

if __name__=='__main__':
    main()
