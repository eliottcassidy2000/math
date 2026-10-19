"""Source-owned 15/35 guards and retained ROOT words for the 125 residuals.

Default verifies supplied certificates, with no orbit discovery. --discover
creates the finite receipt for exactly the baseline-residual j in 1..1024.
All checks are active under -O. Native words use U(n)=(3n+1)/2^v2(3n+1).
No universal Collatz assertion or random source replacement is used.
"""
import argparse
import base64
from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path
import json
import zlib

ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'05-knowledge/results/guarded_children_105_20261018g.json'
CHECKS=0


def check(ok,why):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(why)


def v2(n):
    if n<=0: raise ValueError('positive argument required')
    return (n & -n).bit_length()-1


def baseline_accepts(j):
    # Independently replay the inherited sufficient grammar through20 T ticks.
    x=(3*pow(125,j)+1)//4
    a,b,D=3,1,4
    for depth in range(1,19):
        if x&1:
            x=(3*x+1)//2
            a,b=3*a,3*b+D
        else: x//=2
        D*=2
        if 2*(125*a+b)<=125*D: return True
    return False


def discover_word(source,cap=100000):
    x=source; word=[]
    for _ in range(cap):
        if x==1: return bytes(word)
        y=3*x+1; a=v2(y)
        check(a<256,'byte encoding bound')
        word.append(a); x=y>>a
    raise RuntimeError(('UNRESOLVED cap',source.bit_length()))


def word_digest(word): return sha256(word).hexdigest()


def carrier(word):
    p,b,A=1,0,0
    for a in word:
        b=3*b+(1<<A); p*=3; A+=a
    return p,b,A


def first_half(source,word):
    x=source
    for index,a in enumerate(word,1):
        x=(3*x+1)>>a
        if 2*x<=source: return index,x
    raise RuntimeError('grounded nontrivial source has no half child')


def threshold(p,b,A,minimum):
    c=1
    while (minimum << (A+c)) < 2*(p*minimum+b): c+=1
    return c


def template(j,word):
    n=125**j
    length,h=first_half(n,word)
    prefix=word[:length-1]
    p0,b0,A=carrier(prefix)
    p,b=3*p0,3*b0+(1<<A)
    c0=threshold(p,b,A,n)
    K=A+c0-2
    denominator=(1 << (A+c0))-2*p
    check(denominator>0,'strict asymptotic half slope')
    cutoff=(2*b+denominator-1)//denominator
    check(K>=19 and 0<j<(1<<K),'positive least exponent / inherited residual retained')
    check(word[length-1]>=c0,'actual valuation meets coarse gate')
    check(n>=cutoff,'family minimum pays carry')
    q,bb,AA=1,0,0
    for a in prefix:
        bb=3*bb+(1<<AA); q*=3; AA+=a
        check(2*q>=(1<<AA),'every proper odd checkpoint stays above half')
    return dict(j=j,prefix=list(prefix),prefix_cost=A,P=str(p),B=str(b),
                c_min=c0,guard_bits=K,source_cutoff=str(cutoff),
                actual_c=word[length-1],child=str(h),child_rank=len(word)-length,
                child_cost=sum(word[length:]),child_word_sha256=word_digest(word[length:]),
                first_odd_half_index=length,
                half_shortcut_clock=A+word[length-1],
                half_ordinary_clock=A+word[length-1]+length,
                source_residues=[n%15,n%35],child_residues=[h%15,h%35],
                child_mod105=h%105)


def source_word_check(j,word):
    n=125**j; x=n
    check(len(word)>0,'nonempty source word')
    for a in word:
        check(x>1 and a>0,'strict first hit / positive valuation')
        y=3*x+1
        check(v2(y)==a,'literal exact valuation')
        x=y>>a
    check(x==1,'literal ROOT terminal')
    # A separate reverse path does not assume intermediate forward states.
    x=1
    for a in reversed(word):
        numerator=(x<<a)-1
        check(numerator%3==0,'reverse integrality')
        x=numerator//3
        check(x>1 and x&1,'reverse first-hit and oddness')
    check(x==n,'reverse recovers original source')


def exact_last_valuation(j,rec):
    p,b,A=int(rec['P']),int(rec['B']),rec['prefix_cost']
    bits=A+rec['c_min']+16
    while True:
        modulus=1<<bits
        numerator=(p*pow(125,j,modulus)+b)%modulus
        if numerator: return v2(numerator)-A
        bits+=16
        if bits>A+rec['c_min']+256:
            raise RuntimeError('valuation probe cap, not an inferred valuation')


def modular_guard_check(j,rec):
    A,K=rec['prefix_cost'],rec['guard_bits']
    check(j% (1<<K)==rec['j'],'source-specific exponent guard')
    c=exact_last_valuation(j,rec)
    check(c>=rec['c_min'],'variable-depth completion')
    # Preserve enough denominator precision for every division and both odd moduli.
    modulus=105*(1 << (A+c+1))
    x=pow(125,j,modulus)
    for a in rec['prefix']:
        y=3*x+1
        check(v2(y)==a,'native prefix in quotient retaining denominator')
        x=y>>a
    y=3*x+1
    check(v2(y)==c,'actual terminal valuation')
    h=(y>>c)%105
    p,b=int(rec['P']),int(rec['B'])
    pair=[]
    for m in (15,35):
        phase=(A+c)%12  # sufficient for both, since ord_105(2)=12.
        value=((p*pow(125,j,m)+b)*pow(pow(2,phase,m),-1,m))%m
        pair.append(value)
    check((pair[0]-pair[1])%5==0,'shared phase compatibility')
    decoded=(pair[0]+15*((pair[1]-pair[0])%7))%105
    check(decoded==h,'mod15/mod35 decoder keeps actual child')
    return dict(j=str(j),c=c,child_mod105=h,pair=pair)


def discover():
    residual=[j for j in range(1,1025) if not baseline_accepts(j)]
    check(len(residual)==63,'frozen baseline residual count')
    rows=[]
    for j in residual:
        word=discover_word(125**j)
        rec=template(j,word)
        rec.update(word_rank=len(word),word_cost=sum(word),word_sha256=word_digest(word),
                   word_zlib_base85=base64.b85encode(zlib.compress(word,9)).decode('ascii'))
        rows.append(rec)
    output=dict(status='FINITE-EXACT supplied ROOT words; PROVED guarded half-descent families',
                input_revision='621a39554e38',baseline='digit_125_nested_105_20261018q.md P2, depth18',
                universe='Positive integer exponents 1..1024 rejected by the inherited sufficient atlas',
                cap_odd_steps=100000,rows=rows)
    DEST.write_text(json.dumps(output,indent=2,sort_keys=True)+'\n',encoding='utf-8')


def verify():
    data=json.loads(DEST.read_text(encoding='utf-8'))
    residual=[j for j in range(1,1025) if not baseline_accepts(j)]
    check([r['j'] for r in data['rows']]==residual,'all and only baseline residuals')
    totals=[0,0]; controls=[]; broadened=0
    for rec in data['rows']:
        j=rec['j']
        word=zlib.decompress(base64.b85decode(rec['word_zlib_base85']))
        check(word_digest(word)==rec['word_sha256'],'decoded certificate hash')
        check((len(word),sum(word))==(rec['word_rank'],rec['word_cost']),'word metadata')
        source_word_check(j,word)
        expected=template(j,word)
        check(all(rec[k]==v for k,v in expected.items()),'actual source / guard / child binding')
        child_word=word[rec['first_odd_half_index']:]
        h=int(rec['child'])
        check(word_digest(child_word)==rec['child_word_sha256'],'grounded child suffix')
        child=h
        for a in child_word:
            check(child>1 and v2(3*child+1)==a,'child exact suffix')
            child=(3*child+1)>>a
        check(child==1,'grounded child terminal')
        # Decimal lift has a fully authenticated fixed prefix to the same125^j.
        decimal=(2048*1000**j-23)//27
        check(27*decimal==2048*1000**j-23,'decimal source membership')
        x=decimal
        for a in (1,2,8+3*j):
            check(x>1 and v2(3*x+1)==a,'decimal native prefix')
            x=(3*x+1)>>a
        check(x==125**j,'decimal lifts this exact grounded source')
        for lift in (0,1,2,7,10**12):
            ctrl=modular_guard_check(j+(lift<<rec['guard_bits']),rec)
            if j in (10,74,712): controls.append(ctrl)
        # A neighboring cell fails the retained valuation threshold.
        failed=j+(1<<(rec['guard_bits']-1))
        modulus=1<<(rec['prefix_cost']+rec['c_min'])
        check((int(rec['P'])*pow(125,failed,modulus)+int(rec['B']))%modulus !=0,
              'one-bit-short guard rejects its other half')
        totals[0]+=len(word); totals[1]+=sum(word)
        broadened+=rec['actual_c']>rec['c_min']
    for i,a in enumerate(data['rows']):
        for b in data['rows'][i+1:]:
            k=min(a['guard_bits'],b['guard_bits'])
            check((a['j']-b['j'])%(1<<k)!=0,'new guards disjoint')
    mass=sum((Fraction(1,1<<r['guard_bits']) for r in data['rows']),Fraction())
    check(totals==[575652,1152705],'independent residual-word totals')
    report=dict(active_checks=CHECKS,residual_count=len(residual),
                residual_exponents=residual,total_odd_rank=totals[0],total_halving_cost=totals[1],
                guard_precision_range=[min(r['guard_bits'] for r in data['rows']),
                                       max(r['guard_bits'] for r in data['rows'])],
                exact_added_exponent_density=str(mass),
                threshold_lower_than_observed_valuation=broadened,
                control_lifts=controls)
    print(json.dumps(report,indent=2))
    return report


def main():
    global CHECKS
    parser=argparse.ArgumentParser()
    parser.add_argument('--discover',action='store_true')
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    if args.discover: discover()
    CHECKS=0
    report=verify()
    if args.report:
        args.report.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')


if __name__=='__main__': main()
