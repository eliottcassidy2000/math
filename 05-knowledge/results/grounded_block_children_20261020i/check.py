"""Exact stable reset blocks and source-bound grounded child factories.

All finite words are supplied in inputs.json. Production does not discover
orbits or ROOT suffixes. Symbolic final valuations are never materialized by
the modular reader. Default saves the receipt; --verify is read only.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product
from hashlib import sha256
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def integer(x, lo=0):
    if type(x) is not int or x < lo:
        raise ValueError('exact integer outside domain')
    return x


def odd(x):
    integer(x, 1)
    if not x & 1:
        raise ValueError('positive odd integer required')
    return x


def word(w):
    if type(w) is not tuple:
        raise ValueError('literal immutable valuation word required')
    for a in w:
        integer(a, 1)
    return w


def v2(x):
    integer(x, 1)
    return (x & -x).bit_length() - 1


def state(w):
    word(w)
    lam, pos, debt = F(1), F(0), F(0)
    P, Q, B = 1, 1, 0
    for a in w:
        mu = F(3, 2**a)
        lam, pos = mu*lam, mu*pos + (a == 1)
        debt = mu*debt + (1-F(4, 2**a) if a >= 3 else 0)
        P, Q, B = 3*P, Q*2**a, 3*B+Q
    return lam, pos, debt, P, Q, B


def replay(n, w, root=False):
    odd(n); word(w)
    if type(root) is not bool:
        raise ValueError('literal ROOT flag required')
    x, xs = n, [n]
    for a in w:
        if x == 1:
            raise ValueError('ROOT padding rejected')
        z = 3*x+1
        if v2(z) != a:
            raise ValueError('wrong native source or valuation')
        x = z >> a
        xs.append(x)
    if root and x != 1:
        raise ValueError('supplied word has no ROOT endpoint')
    return tuple(xs)


def reverse(w, h):
    word(w); odd(h)
    x = h
    for a in reversed(w):
        z = (x << a)-1
        if z % 3:
            raise ValueError('nonintegral inverse seam')
        x = odd(z//3)
    return x


def stable(w):
    lam, _, debt, *_ = state(w)
    return lam+debt <= 1


def cone(w):
    lam, _, debt, *_ = state(w)
    return lam >= 1 and debt <= lam-1


def cut(w):
    word(w)
    return next((k for k in range(len(w)-1, -1, -1)
                 if cone(w[:k]) and stable(w[k:])), None)


@dataclass(frozen=True)
class Root:
    source: int
    word: tuple


def audit_root(r):
    if type(r) is not Root:
        raise ValueError('typed ROOT certificate required')
    odd(r.source); word(r.word)
    replay(r.source, r.word, root=True)
    if reverse(r.word, 1) != r.source:
        raise RuntimeError('independent inverse ROOT identity failed')
    return r


@dataclass(frozen=True)
class Paid:
    source: int
    word: tuple
    cut: int


def audit_paid(r):
    if type(r) is not Paid:
        raise ValueError('typed paid receipt required')
    odd(r.source); word(r.word); integer(r.cut)
    if not 0 <= r.cut < len(r.word):
        raise ValueError('empty terminal block')
    if not cone(r.word[:r.cut]) or not stable(r.word[r.cut:]):
        raise ValueError('cone or stable block premise failed')
    h = replay(r.source, r.word)[-1]
    if h >= r.source:
        raise ValueError('no actual original-source descent')
    if not state(r.word)[2] < r.source-h:
        raise RuntimeError('proved positive-envelope payment failed')
    return h


def discharge(r, child):
    h = audit_paid(r); audit_root(child)
    if child.source != h:
        raise ValueError('child ROOT certificate is from a different source')
    result = r.word+child.word
    replay(r.source, result, root=True)
    if reverse(result, 1) != r.source:
        raise RuntimeError('complete inverse source recovery failed')
    return result


def coarse(r, t):
    """Variable-final-valuation paid ray; no automatic child ROOT claim."""
    h = audit_paid(r); integer(t)
    if r.word[-1] < 2:
        raise ValueError('variable final letter must be at least two')
    _, _, _, P, Q, _ = state(r.word)
    a = r.word[-1]+v2(h+P*t)
    result = Paid(r.source+Q*t, r.word[:-1]+(a,), r.cut)
    endpoint = audit_paid(result)
    if endpoint != (h+P*t) >> v2(h+P*t):
        raise RuntimeError('coarse endpoint identity failed')
    return result


def cap_stable(w, cap):
    integer(cap, 1)
    lam, _, debt, *_ = state(w)
    return lam+debt/cap <= 1


def cap_cut(w, cap):
    word(w); integer(cap, 1)
    return next((k for k in range(len(w)-1, -1, -1)
                 if cone(w[:k]) and cap_stable(w[k:], cap)), None)


@dataclass(frozen=True)
class CapPaid:
    source: int
    word: tuple
    cut: int
    cap: int


def audit_cap(r):
    if type(r) is not CapPaid:
        raise ValueError('typed cap-payment certificate required')
    odd(r.source); word(r.word); integer(r.cut); integer(r.cap, 1)
    if not 0 <= r.cut < len(r.word):
        raise ValueError('nonempty cap terminal required')
    if not cone(r.word[:r.cut]) or not cap_stable(r.word[r.cut:], r.cap):
        raise ValueError('cap cone or terminal premise failed')
    h = replay(r.source, r.word)[-1]
    if r.source-h < r.cap+1:
        raise ValueError('original-source gap does not pay cap')
    if not state(r.word)[2] < r.cap+1 <= r.source-h:
        raise RuntimeError('proved cap payment inequality failed')
    return h


def cap_discharge(r, child):
    h = audit_cap(r); audit_root(child)
    if h != child.source:
        raise ValueError('wrong cap child source')
    whole = r.word+child.word
    replay(r.source, whole, root=True)
    if reverse(whole, 1) != r.source:
        raise RuntimeError('inverse cap ROOT source failed')
    return whole


def cap_coarse(r, t):
    h = audit_cap(r); integer(t)
    if r.word[-1] < 2:
        raise ValueError('variable cap final valuation at least two')
    _, _, _, P, Q, _ = state(r.word)
    extra = v2(h+P*t)
    result = CapPaid(r.source+Q*t, r.word[:-1]+(r.word[-1]+extra,), r.cut, r.cap)
    endpoint = audit_cap(result)
    if endpoint != (h+P*t) >> extra:
        raise RuntimeError('cap coarse native endpoint failed')
    return result


def phase2(target, precision):
    """Least exponent modulo 2*3^(precision-1), for a literal 3-adic unit."""
    integer(target, 1); integer(precision, 1)
    if target % 3 == 0:
        raise ValueError('target is not a ternary unit')
    a, period = (0 if target % 3 == 1 else 1), 2
    for k in range(1, precision):
        modulus = 3**(k+1)
        found = [a+i*period for i in range(3)
                 if pow(2, a+i*period, modulus) == target % modulus]
        if len(found) != 1:
            raise RuntimeError('primitive-root phase lift failed')
        a = found[0]
        period *= 3
    return a, period


def greater_shift(x, e, y):
    """Exact x*2^e > y without expanding a huge exponent."""
    integer(x, 1); integer(e); integer(y)
    left_bits = x.bit_length()+e
    if left_bits != y.bit_length():
        return left_bits > y.bit_length()
    return (x << e) > y


@dataclass(frozen=True)
class Factory:
    head: tuple
    cut: int
    child: Root
    minimum: int
    first: int
    period: int
    P: int
    A: int
    C: int


def factory(head, k, child, minimum=2):
    word(head); integer(k); integer(minimum, 2); audit_root(child)
    if child.source <= 1 or child.source % 3 == 0:
        raise ValueError('this factory retains a non-ROOT ternary-unit child')
    if k > len(head) or not cone(head[:k]):
        raise ValueError('bad cone prefix')
    l, _, d, *_ = state(head[k:])
    if l+d > F(4, 3):
        raise ValueError('last-letter completion cannot stabilize this block')
    _, _, _, oldP, Q, B = state(head)
    P, A, C = 3*oldP, sum(head), 3*B+Q
    target = (C*pow(child.source*Q, -1, P)) % P
    phase, period = phase2(target, len(head)+1)
    lower = minimum
    # This loop depends on the finite head coefficients, not the giant phase.
    while not greater_shift(child.source, A+lower, C+P*child.source):
        lower += 1
    first = phase + max(0, (lower-phase+period-1)//period)*period
    return Factory(head, k, child, minimum, first, period, P, A, C)


def audit_factory(f):
    if type(f) is not Factory:
        raise ValueError('typed grounded factory required')
    word(f.head); integer(f.cut); audit_root(f.child)
    for field in ('minimum', 'first', 'period', 'P', 'A', 'C'):
        integer(getattr(f, field), 0 if field == 'A' else 1)
    expected = factory(f.head, f.cut, f.child, f.minimum)
    if f != expected:
        raise ValueError('factory source/head/phase metadata changed')
    return f


def member(f, index):
    audit_factory(f); integer(index)
    a = f.first+f.period*index
    return {'index': index, 'final_valuation': a,
            'head': f.head, 'child': f.child.source,
            'source': {'numerator_coefficient': f.child.source,
                       'power_of_two': f.A+a, 'subtract': f.C, 'divide': f.P},
            'ROOT_odd_steps': len(f.head)+1+len(f.child.word),
            'ROOT_valuation_cost': f.A+a+sum(f.child.word),
            'ROOT_ordinary_steps': f.A+a+sum(f.child.word)+len(f.head)+1+len(f.child.word)}


def source_mod(f, index, modulus):
    audit_factory(f); integer(index); integer(modulus, 1)
    a = f.first+f.period*index
    numerator = f.child.source*pow(2, f.A+a, f.P*modulus)-f.C
    if numerator % f.P:
        raise RuntimeError('denominator-preserving reader is not integral')
    return (numerator//f.P) % modulus


def materialize(f, index, max_bits=4096):
    audit_factory(f); integer(index); integer(max_bits, 1)
    a = f.first+f.period*index
    if f.child.source.bit_length()+f.A+a > max_bits:
        raise ValueError('symbolic source exceeds the declared materialization cap')
    n = ((f.child.source << (f.A+a))-f.C)//f.P
    receipt = Paid(n, f.head+(a,), f.cut)
    discharge(receipt, f.child)
    return receipt


def reject(fn, message):
    try:
        fn()
    except ValueError:
        need(True, message)
    else:
        need(False, message)


def frac(x):
    return [x.numerator, x.denominator]


def verify():
    global CHECKS
    CHECKS = 0
    input_path = HERE/'inputs.json'
    inp = json.loads(input_path.read_text(encoding='utf-8'))
    oldpath = HERE.parent/'debt_annulus_20261020h/receipt.json'
    need(sha256(oldpath.read_bytes().replace(b'\r\n',b'\n')).hexdigest()==inp['old_receipt_LF_SHA256'],
         'immutable old pending universe')
    old = json.loads(oldpath.read_text(encoding='utf-8'))
    roots = {}
    for r in inp['roots']:
        if type(r) is not dict or type(r['word']) is not list:
            raise ValueError('literal supplied ROOT row required')
        cert = audit_root(Root(r['source'], tuple(r['word'])))
        need(cert.source not in roots, 'no repeated child source')
        roots[cert.source] = cert
    rows, unresolved, used, records = [], [], set(), []
    for row in inp['rows']:
        if type(row) is not dict or type(row['word']) is not list:
            raise ValueError('literal supplied first-exit row required')
        n = odd(row['source']); w = word(tuple(row['word'])); h = odd(row['child'])
        xs = replay(n,w)
        need(xs[-1]==h and h<n and all(x>=n for x in xs[:-1]), 'actual first exit and original source')
        k = cut(w)
        if row['cut'] is not None:
            integer(row['cut'])
        need(k==row['cut'], 'selected exact stable-block cut')
        if k is None:
            unresolved.append(n); continue
        p = Paid(n,w,k); child=roots[h]; used.add(h)
        whole = discharge(p,child); records.append(p)
        st=state(w[k:])
        need(st[0]+st[2]<=1, 'stable terminal exact inequality')
        f = factory(w[:-1],k,child,w[-1])
        need(materialize(f,0).source==n and f.first==w[-1], 'factory retains its actual seed')
        for t in (0,1,2,3,7,16):
            c=coarse(p,t)
            need(audit_paid(c)<c.source, 'whole native coarse-ray payment')
        for index in (0,1,2,5,10**12):
            info=member(f,index)
            need(info['ROOT_ordinary_steps']==info['ROOT_valuation_cost']+info['ROOT_odd_steps'], 'literal clock distinction')
            for modulus in (1,3,15,35,105,2**(sum(w)+1),3**(len(w)+2)):
                nmod=source_mod(f,index,modulus);a=info['final_valuation']
                need((f.P*nmod+f.C-f.child.source*pow(2,f.A+a,modulus))%modulus==0, 'symbolic source equation in full query')
            need(source_mod(f,index,105)%15==source_mod(f,index,15) and
                 source_mod(f,index,105)%35==source_mod(f,index,35), 'compatible source registers')
        need(source_mod(f,1,2**sum(w))==n%2**sum(w), 'grounded inflation stays in the paid coarse ray')
        rows.append({'source':n,'child':h,'cut':k,'terminal':w[k:], 'stable_sum':frac(st[0]+st[2]),
                     'ROOT_rank':len(whole),'ROOT_cost':sum(whole),'factory_period':f.period,
                     'registers_index0_5':[source_mod(f,j,105) for j in range(6)]})
    need([r['source'] for r in inp['rows']]==old['pending_sources'], 'complete exact inherited 175 source universe')
    need(len(rows)==146 and len(unresolved)==29 and min(unresolved)==927, 'new controller finite coverage')
    need(used==set(roots) and len(roots)==95, 'all95 supplied child certificates used')
    # Separate source-gap cap extension; the146 unit factories stay unchanged.
    cap_roots = dict(roots)
    for row in inp['cap_roots']:
        if type(row) is not dict or type(row['word']) is not list:
            raise ValueError('literal cap child ROOT row required')
        child = audit_root(Root(row['source'], tuple(row['word'])))
        need(child.source not in cap_roots, 'cap suffix is new or reuses retained old suffix')
        cap_roots[child.source] = child
    cap_rows = []
    for row in inp['rows']:
        if row['source'] not in unresolved:
            continue
        n, w, h = row['source'], tuple(row['word']), row['child']
        k = cap_cut(w, 3)
        need(k is not None, 'all29 residual words admit cap3')
        cp = CapPaid(n,w,k,3); whole = cap_discharge(cp,cap_roots[h])
        for t in (0,1,2,3,7,16):
            lifted = cap_coarse(cp,t)
            need(audit_cap(lifted) <= h+3**len(w)*t, 'cap coarse actual endpoint')
        lam, _, debt, *_ = state(w[k:])
        cap_rows.append({'source':n,'child':h,'cut':k,'cap':3,'terminal':w[k:],
                         'stable_sum':frac(lam+debt/3),'actual_gap':n-h,
                         'whole_debt':frac(state(w)[2]),'ROOT_rank':len(whole),'ROOT_cost':sum(whole)})
    need(len(cap_rows)==29 and min(r['actual_gap'] for r in cap_rows)==50,
         'all29 actual original gaps pay cap3')
    need(len(cap_roots)==113, 'eighteen additional strict child suffixes')
    hostile=next(r for r in cap_rows if r['source']==4575)
    need(F(*hostile['whole_debt'])>2 and hostile['actual_gap']==110,
         'fixed-threshold-two obstruction bypassed by retained source gap')
    consecutive=sum(all(a>=2 for a in p.word[p.cut:]) for p in records)
    # Classification by existence of any consecutive-reset cut, independently of selected cut.
    consecutive=0
    for p in records:
        consecutive += any(cone(p.word[:k]) and all(a>=2 for a in p.word[k:]) for k in range(len(p.word)))
    need(consecutive==117, '117 simple and29 genuinely mixed stable suffixes')

    # Independent suffix formulas, exact interval iff, last-letter completion,
    # and finite composition closure; no fitted finite statement is promoted.
    short=[]; unstable_witnesses=[]
    for length in range(0,6):
        for w in product(range(1,5),repeat=length):
            l,p,d,P,Q,B=state(w)
            expected=sum((1-F(4,2**a))*F(3**(length-i-1),2**sum(w[i+1:]))
                         for i,a in enumerate(w) if a>=3)
            need(d==expected and l==F(P,Q), 'independent suffix debt and slope')
            for a in (2,3,4,8):
                need(stable(w+(a,))==(l+d<=F(4,3)), 'sharp variable-final stable completion iff')
                for cap in (1,3,7):
                    value=state(w+(a,))[0]+state(w+(a,))[2]/cap
                    expected=F(3,4)*(l+d/cap) if a==2 else F(1,cap)+(3*l+(3*d-4)/cap)/2**a
                    need(value==expected, 'exact cap final-letter formula')
                    if value<=1:
                        need(cap_stable(w+(a+1,),cap), 'raising actual final valuation preserves cap')
            if l+d>1:
                x=max(F(0),(1-d)/l);x=(1+x)/2
                need(0<=x<1 and l*x+d>1, 'literal rational unstable-block hostile')
                if len(unstable_witnesses)<4:unstable_witnesses.append([w,frac(x),frac(l*x+d)])
            else:
                for x in (F(0),F(1,3),F(999,1000)):
                    need(0<=l*x+d<1, 'stable positive interval control')
                if length<=3:short.append(w)
    for u in short:
        for v in short:
            need(stable(u+v), 'stable-block monoid closure')
    need(stable((4,3)) and not stable((1,4)), 'positive block and expanding interval hostile')
    need(stable((1,2,3)) and not stable((2,3,1)), 'same letters, different ordered stable predicate')
    for r in range(0,13):
        l,_,d,*_=state((4,3)*r)
        need(l==F(9,128)**r and d==F(100,119)*(1-F(9,128)**r), 'exact repeated-block recursion')
    for precision in range(1,7):
        modulus=3**precision;period=2*3**(precision-1)
        powers={pow(2,a,modulus):a for a in range(period)}
        need(len(powers)==period, 'independent complete ternary-unit power table')
        for target,a in powers.items():
            need(phase2(target,precision)==(a,period), 'all small exponent phases')

    # All-r architecture: exact heads differ with r; final phase remains marked.
    seed=next(p for p in records if p.source==31);u=seed.word[:-2]
    need(len(u)==33 and sum(u)==49 and seed.cut==33 and seed.word[-2:]==(4,3), 'literal31 architecture')
    need(all(state(seed.word[:j])[0]>1 for j in range(1,len(seed.word))), '31 proper odd prefixes grow on its entire coarse ray')
    h=roots[23]; architecture=[]
    for r in (1,2,3,4,8,16):
        f=factory(u+(4,3)*(r-1)+(4,),len(u),h,3)
        need(f.P==3**(33+2*r) and f.A==46+7*r and f.period==2*3**(32+2*r), 'all-r exact clock and denominator')
        table=[source_mod(f,j,105) for j in range(6)]
        need(len(set(table))==6, 'six distinct marked31 architecture phases')
        for j in range(12):
            need(source_mod(f,j,105)==table[j%6], 'period-six source register controls')
        for j in (0,1,10**20):
            a=f.first+f.period*j
            need(pow(2,a,f.P)==(f.C*pow(23*(1<<f.A),-1,f.P))%f.P, 'all-r source divisibility phase')
            need(greater_shift(23,f.A+a,f.C+23*f.P), 'symbolic original-source strict descent')
        architecture.append({'r':r,'head':f.head,'cut':f.cut,'first_valuation':f.first,'period':f.period,
                             'P':f.P,'head_cost':f.A,'C':f.C,'child':23,'ROOT_suffix':h.word,
                             'registers':table,'index0':member(f,0)})
    f1=factory(seed.word[:-1],seed.cut,h,3)
    need(materialize(f1,0)==seed, '31 recovered literally, not just modulo a guard')
    need(len(discharge(seed,h))==39 and sum(discharge(seed,h))==67, '31 strict ROOT clocks39/67/106')
    reject(lambda:materialize(f1,1,4096), 'giant source protected by explicit bit cap')

    # Source and type hostiles; no unsupported ROOT or source substitution.
    reject(lambda:audit_paid(Paid(True,seed.word,seed.cut)), 'Boolean source rejected')
    reject(lambda:audit_paid(Paid(31,list(seed.word),seed.cut)), 'mutable word rejected')
    reject(lambda:audit_paid(Paid(31,seed.word,True)), 'Boolean cut rejected')
    reject(lambda:discharge(seed,roots[19]), 'wrong actual child rejected')
    reject(lambda:audit_root(Root(23,h.word+(2,))), 'ROOT padding rejected')
    reject(lambda:coarse(seed,-1), 'negative progression parameter rejected')
    reject(lambda:phase2(3,2), 'nonunit ternary target rejected')
    reject(lambda:source_mod(f1,0,False), 'Boolean query modulus rejected')
    reject(lambda:source_mod(replace(f1,first=float(f1.first)),0,105), 'float phase metadata rejected')
    reject(lambda:source_mod(replace(f1,head=f1.head+(2,)),0,105), 'changed factory head rejected')
    reject(lambda:factory((1,),0,h), 'expanding final-prefix cannot be stabilized by last valuation')
    reject(lambda:audit_cap(CapPaid(4575,tuple(next(r['word'] for r in inp['rows'] if r['source']==4575)),True,3)), 'Boolean cap cut rejected')
    reject(lambda:cap_stable((4,3),3.0), 'floating cap rejected')
    reject(lambda:cap_discharge(cp,roots[23]), 'cap wrong child rejected')
    reject(lambda:audit_cap(CapPaid(3,(1,4),0,3)), 'source gap smaller than cap budget rejected')
    mersenne=[p.source for p in records if (p.source+1)&p.source==0]
    need(mersenne==[31,127,2047], 'exact new Mersenne controller anchors')
    for p in records:
        if p.source not in mersenne:continue
        K=sum(p.word)
        need(p.source+1 < 2**K and p.source != 2**K-1, 'coarse Mersenne guard has only the anchor')
    return {'status':'PROVED block monoid and exact grounded factories; FINITE-EXACT175 inherited pending anchors',
            'checks':CHECKS,'old_receipt_LF_SHA256':inp['old_receipt_LF_SHA256'],
            'input_LF_SHA256':sha256(input_path.read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
            'new_paid_grounded':rows,'unit_cap_pending':unresolved,'cap3_paid_grounded':cap_rows,'pending':[],
            'counts':{'old_pending':175,'unit_cap_new':146,'cap3_additional':29,'new':175,'remaining':0,
                'consecutive_suffix':117,'additional_mixed_blocks':29,'unit_cap_supplied_children':95,
                'additional_supplied_children':18,'supplied_children':113,'inherited_plus_new_controller':5000},
            'mersenne_controller_anchors':mersenne,'all_r_31_architecture':architecture,
            'unstable_interval_hostiles':unstable_witnesses,
            'scope':'coarse rays are paid; only marked fixed-child subfamilies are grounded; no universal ROOT or first-exit claim for repeated architecture'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    payload=verify();text=json.dumps(payload,sort_keys=True,indent=2)+'\n';out=HERE/'receipt.json'
    if args.verify:
        if out.read_text(encoding='utf-8').replace('\r\n','\n')!=text:
            raise RuntimeError('saved receipt differs')
    else:
        out.write_text(text,encoding='utf-8',newline='\n')
    print(json.dumps({'checks':CHECKS,'counts':payload['counts'],'LF_SHA256':sha256(text.encode()).hexdigest()}))


if __name__=='__main__':
    main()
