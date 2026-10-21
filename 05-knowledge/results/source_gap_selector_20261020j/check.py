"""Positive complement register and source-owned gap selectors.

All actual words are supplied. No production first-exit or ROOT discovery.
Default writes receipt.json; --verify independently recomputes and compares.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent/'grounded_block_children_20261020i'
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def integer(x, minimum=0):
    if type(x) is not int or x < minimum:
        raise ValueError('literal integer outside domain')
    return x


def rational(x, minimum=F(0)):
    if type(x) not in (int, F) or x < minimum:
        raise ValueError('exact rational outside domain')
    return F(x)


def odd(n):
    integer(n, 1)
    if not n & 1:
        raise ValueError('positive odd source required')
    return n


def word(w):
    if type(w) is not tuple:
        raise ValueError('immutable literal valuation word required')
    for a in w:
        integer(a, 1)
    return w


def v2(n):
    integer(n, 1)
    return (n & -n).bit_length()-1


def coefficients(w):
    word(w)
    lam, pos, debt = F(1), F(0), F(0)
    P, Q, B = 1, 1, 0
    inverse_nonones = F(0)
    for a in w:
        if a >= 2:
            inverse_nonones += 1/lam
        mu = F(3, 2**a)
        lam, pos, debt = mu*lam, mu*pos+(a == 1), mu*debt+max(F(0), 1-F(4, 2**a))
        P, Q, B = 3*P, 2**a*Q, 3*B+Q
    return lam, pos, debt, P, Q, B, inverse_nonones


def replay(n, w):
    odd(n); word(w)
    xs = [n]
    for a in w:
        if xs[-1] == 1:
            raise ValueError('strict ROOT stopping excludes padding')
        z = 3*xs[-1]+1
        if v2(z) != a:
            raise ValueError('word is not native at the retained source')
        xs.append(z >> a)
    return tuple(xs)


def inverse(w, h):
    word(w); odd(h)
    for a in reversed(w):
        z = (h << a)-1
        if z % 3:
            raise ValueError('nonintegral inverse source')
        h = odd(z//3)
    return h


@dataclass(frozen=True)
class Register:
    source: int
    word: tuple
    endpoint: int
    slope: F
    debt: F
    complement: F
    normalized: F


def register(n, w):
    xs = replay(n,w)
    lam, pos, debt, P, Q, B, reset_sum = coefficients(w)
    h = xs[-1]
    K = F(h+1)-debt
    need(K == lam*(n+1+F(2,3)*reset_sum) > 0, 'positive complement and normalized register')
    need(h-1+debt == lam*(n-1)+pos == 2*h-K, 'exact positive envelope')
    need(Q*h == P*n+B and inverse(w,h)==n, 'independent complete affine and inverse source')
    return Register(n,w,h,lam,debt,K,K/lam)


def audit(r):
    if type(r) is not Register:
        raise ValueError('source register required')
    odd(r.source); word(r.word); odd(r.endpoint)
    for field in ('slope','debt','complement','normalized'):
        rational(getattr(r,field))
    expected = register(r.source,r.word)
    if r != expected:
        raise ValueError('source, history or complement register changed')
    return r


def half_paid(r):
    audit(r)
    if r.source <= 1 or not r.word or 2*r.endpoint+1 > r.source:
        raise ValueError('no original-source half descent')
    need(r.debt < r.endpoint+1 <= r.source-r.endpoint, 'any-history half-descent pays')
    return r.source-r.endpoint-r.debt


def repair(n, head, tail):
    """Three native non-one letters after any descent, or earlier strict ROOT."""
    odd(n); word(head); word(tail)
    if n <= 1 or not head:
        raise ValueError('nontrivial original source and descent head required')
    first = register(n,head)
    if first.endpoint >= n:
        raise ValueError('head has no actual original-source descent')
    if len(tail)>3 or any(a<2 for a in tail):
        raise ValueError('repair must use at most three non-one letters')
    result = register(n,head+tail)
    if len(tail)!=3 and result.endpoint!=1:
        raise ValueError('short repair must terminate at first ROOT')
    need(result.debt < n-result.endpoint, 'same-source bounded repair payment')
    return result


def last_threshold(r):
    """Exact pre-step threshold, not a claim that the actual next valuation meets it."""
    audit(r)
    if r.source<=1 or r.endpoint==1:
        raise ValueError('non-ROOT head at nontrivial original source required')
    E = r.endpoint-1+r.debt
    a = 2
    while F(3,2**a)*E >= r.source-1:
        a += 1
    return a


def cap_requirement(v):
    word(v)
    if not v:
        raise ValueError('nonempty terminal word required')
    lam, _, debt, *_ = coefficients(v)
    if lam>=1:
        return None
    return max(F(1),debt/(1-lam))


def cap_selector(n, w):
    """Exact architecture test with the actual gap retained; never proves arrival."""
    r = register(n,w)
    g = n-r.endpoint
    if g<2:
        return None
    options=[]
    for k in range(len(w)):
        lam,_,debt,*_=coefficients(w[:k])
        if lam<1 or debt>lam-1:
            continue
        c=cap_requirement(w[k:])
        if c is not None and c<=g-1:
            options.append((c,k))
    return min(options) if options else None


def affine_barrier(alpha, beta):
    alpha,beta=rational(alpha),rational(beta)
    if alpha+beta==0:
        raise ValueError('empty nonnegative-debt domain is excluded')
    return alpha>=beta>=1


def barrier_hostile(alpha,beta):
    """Native state transition, synthetic incoming debt: not a claimed source history."""
    alpha,beta=rational(alpha),rational(beta)
    if affine_barrier(alpha,beta):
        raise ValueError('the barrier is invariant')
    if alpha<beta:
        x,a=3,1
    else:
        a=4
        while beta-1+(alpha-3*beta+4)/2**a>=0:
            a+=2
        x=(2**a-1)//3
    y=(3*x+1)//2**a
    mu=F(3,2**a);delta=max(F(0),1-F(4,2**a))
    upper=alpha*x+beta
    slack=alpha*y+beta-mu*upper-delta
    eps=min(upper/2,-slack/(2*mu))
    d=upper-eps
    need(upper>0 and slack<0 and v2(3*x+1)==a, 'actual native barrier-hostile edge')
    need(0<=d<upper and mu*d+delta>alpha*y+beta, 'strict relaxed incoming-debt violation')
    return {'x':x,'a':a,'y':y,'incoming_debt':d,'outgoing_debt':mu*d+delta}


def rejects(fn):
    try:fn()
    except ValueError:return True
    return False


def serial(x):
    if type(x) is F:return [x.numerator,x.denominator]
    if type(x) is Register:return {k:serial(getattr(x,k)) for k in x.__dataclass_fields__}
    if type(x) is dict:return {str(k):serial(v) for k,v in x.items()}
    if type(x) in (tuple,list):return [serial(v) for v in x]
    return x


def lfhash(p):
    return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def verify():
    global CHECKS
    CHECKS=0
    inp=json.loads((HERE/'inputs.json').read_text(encoding='utf-8'))
    old=json.loads((PRIOR/'inputs.json').read_text(encoding='utf-8'))
    expected=json.loads((PRIOR/'receipt.json').read_text(encoding='utf-8'))
    need(expected['counts']['new']==175, 'declared inherited175 controls only')
    roots={r['source']:tuple(r['word']) for r in old['roots']+old['cap_roots']}
    finite=[]; repairs=[];half_count=0
    for row in old['rows']:
        n=row['source'];w=tuple(row['word']);r=register(n,w)
        need(r.endpoint==row['child'] and all(x>=n for x in replay(n,w)[:-1]), 'retained actual first-exit control')
        half=2*r.endpoint+1<=n
        if half:
            half_paid(r);half_count+=1
        cap=cap_selector(n,w)
        need(cap is not None, 'rational architecture selector covers old175 supplied words')
        suffix=roots[r.endpoint];whole=w+suffix
        root=register(n,whole)
        need(root.endpoint==1 and root.debt<2<=n-1, 'all supplied strict ROOT histories pay independently of cone')
        repair_tail=[]
        for a in suffix[:3]:
            if a<2:break
            repair_tail.append(a)
        if len(repair_tail)==3 or replay(r.endpoint,tuple(repair_tail))[-1]==1:
            repaired=repair(n,w,tuple(repair_tail))
            repairs.append({'source':n,'head_endpoint':r.endpoint,'tail':repair_tail,'paid_endpoint':repaired.endpoint})
        finite.append({'source':n,'endpoint':r.endpoint,'half_gate':half,'least_cap':cap[0],'cut':cap[1]})
    need(len(finite)==175 and half_count==124, 'declared half-gate count124 of175')
    special=[]; root_over_one=None
    for row in inp['rows']:
        n=odd(row['source']);whole=word(tuple(row['ROOT_word']))
        xs=replay(n,whole);need(xs[-1]==1,'supplied special ROOT endpoint')
        last_norm=F(n+1)
        for t in range(len(whole)+1):
            r=register(n,whole[:t])
            need(r.normalized>=last_norm and r.complement>0,'all-prefix positive monotone register')
            if t:
                prev=register(n,whole[:t-1]);a=whole[t-1]
                wanted=prev.normalized+(F(2,3)/prev.slope if a>=2 else 0)
                need(r.normalized==wanted,'exact normalized refinement increment')
            last_norm=r.normalized
        end=register(n,whole)
        first=next(i for i,x in enumerate(xs) if i and x<n)
        exit_record=register(n,whole[:first])
        special.append({'source':n,'ROOT':end,'first_exit':exit_record,
                        'first_exit_gap':n-exit_record.endpoint,'first_exit_cap':cap_selector(n,whole[:first])})
        if n==55:
            bad=register(n,whole[:39]);need(bad.endpoint==53 and bad.debt>=2,'actual later descent not automatically paid')
            need(last_threshold(bad)==2 and whole[39]>=2,'next actual reset meets the pre-step threshold')
            fixed=register(n,whole[:40]);need(fixed.debt<n-fixed.endpoint,'literal next reset repairs55 same-source debt')
        if n==123027327:
            need(end.debt>1 and end.debt<2,'D<h hostile at ROOT; D<h+1 survives')
            root_over_one=end
    # Independent formal word check; no arbitrary word is labelled a native source.
    for length in range(6):
        for w in product(range(1,5),repeat=length):
            lam,pos,d,P,Q,B,J=coefficients(w)
            independent=sum(F(2**sum(w[:i]),3**i) for i,a in enumerate(w) if a>=2)
            need(J==independent and lam==F(P,Q),'independent inverse-slope sum')
            for n in (3,19):
                x=F(P*n+B,Q)
                K=x+1-d
                need(K==lam*(n+1+F(2,3)*J)>0,'formal affine complement identity')
            if length:
                c=cap_requirement(w)
                for trial in (F(1),F(3,2),F(3),F(7)):
                    need((lam+d/trial<=1)==(c is not None and c<=trial),'sharp rational cap iff')
    # Full classification of nonnegative affine state-only half-planes.
    hostiles=[]
    for alpha,beta in product((F(i,4) for i in range(9)),repeat=2):
        if alpha+beta==0:continue
        if affine_barrier(alpha,beta):
            for a in range(1,21):
                mu=F(3,2**a);delta=max(F(0),1-F(4,2**a))
                need(alpha/F(2**a)+beta*(1-mu)-delta>=0,'affine barrier local slack')
        else:
            hostile=barrier_hostile(alpha,beta)
            if len(hostiles)<5:hostiles.append({'alpha':alpha,'beta':beta,**hostile})
    # Three is sharp given only D<x+1, not a claim about reachable source histories.
    x,n,d=129,131,F(129);E=x-1+d
    need(replay(x,(2,2,2))==(129,97,73,55),'literal native tail for relaxed sharpness')
    need(d<x+1 and F(3,4)**2*E>=n-1 and F(3,4)**3*E<n-1,
         'two-letter failure and three-letter repair in the relaxed invariant state space')
    base=register(55,tuple(inp['rows'][0]['ROOT_word'])[:39])
    need(rejects(lambda:audit(replace(base,source=True))), 'Boolean source rejected')
    need(rejects(lambda:audit(replace(base,normalized=float(base.normalized)))), 'float register rejected')
    need(rejects(lambda:audit(replace(base,complement=base.complement+1))), 'changed complement rejected')
    need(rejects(lambda:register(3,(1,4,2))), 'ROOT padding rejected')
    need(rejects(lambda:half_paid(register(31,(1,)))), 'half gate cannot infer future descent')
    need(rejects(lambda:repair(55,base.word,(1,2,2))), 'wrong actual/non-one tail rejected')
    need(rejects(lambda:affine_barrier(0,0)), 'vacuous empty debt region rejected')
    need(rejects(lambda:cap_requirement(())), 'empty terminal cap rejected')
    return {'status':'PROVED exact complement, state-barrier classification and conditional repair; arrival OPEN',
            'checks':CHECKS,'input_hash':lfhash(HERE/'inputs.json'),'inherited_inputs_hash':lfhash(PRIOR/'inputs.json'),
            'inherited_receipt_hash':lfhash(PRIOR/'receipt.json'),'finite175':finite,'half_gate_count':half_count,
            'immediate_three_letter_repairs':repairs,'special_sources':special,'ROOT_debt_above_one_hostile':root_over_one,
            'affine_barrier_hostiles':hostiles,'scope':'no universal arrival, no new-source inflation counted as fixed-source progress; no signed ROOT-energy implication'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=verify();out=json.dumps(serial(result),sort_keys=True,indent=2)+'\n';p=HERE/'receipt.json'
    if args.verify:
        if p.read_text(encoding='utf-8').replace('\r\n','\n')!=out:
            raise ValueError('saved receipt differs')
    else:p.write_text(out,encoding='utf-8',newline='\n')
    print(json.dumps({'checks':CHECKS,'LF_SHA256':sha256(out.encode()).hexdigest(),
                      'half_gate':result['half_gate_count'],'immediate_repair_count':len(result['immediate_three_letter_repairs'])}))


if __name__=='__main__':main()
