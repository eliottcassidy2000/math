"""Supplied Mersenne histories and order-forgetting source-gap burst guards.

Production replays inputs only. All checks remain active with python -O.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
CHECKS=0


def need(ok, label):
    global CHECKS
    CHECKS+=1
    if not ok: raise ValueError(label)


def v2(x): return (x & -x).bit_length()-1


def validate_row(row):
    for field,lo,hi in (('exponent',14,80),('mod15',0,14),('mod35',0,34)):
        if type(row[field]) is not int or not lo<=row[field]<=hi:
            raise ValueError('literal integer source register: '+field)
    if type(row['source']) is not int or row['source']!=(1<<row['exponent'])-1:
        raise ValueError('literal retained Mersenne source')
    if row['mod15']!=row['source']%15 or row['mod35']!=row['source']%35:
        raise ValueError('changed modular source guard')
    return row


def rejects(fn):
    try:fn()
    except ValueError:return True
    return False


def affine(w):
    mu,P=F(1),F(0)
    for a in w:
        t=F(3,2**a)
        mu,P=t*mu,t*P+(a==1)
    return mu,P


def histogram(w):
    need(all(type(a) is int and a>=1 for a in w), 'positive valuation letters')
    mu=F(3**len(w),2**sum(w))
    r=w.count(1)
    return mu,2*(F(3,2)**r-1)


def replay(n,w):
    need(type(n) is int and n>1 and n%2==1, 'positive odd original source')
    x=n;D=F(0);E=F(n-1);xs=[n];debts=[D]
    for a in w:
        need(x!=1,'strict first-ROOT stop')
        z=3*x+1
        need(type(a) is int and a==v2(z), 'actual source valuation')
        t=F(3,2**a);x=z>>a
        D=t*D+max(F(0),1-F(4,2**a))
        E=t*E+(a==1)
        need(E==x-1+D and 0<=D<x+1,'positive-complement recurrence')
        xs.append(x);debts.append(D)
    return xs,debts


def inverse(h,w):
    for a in reversed(w):
        z=(h<<a)-1
        need(z>0 and z%3==0,'inverse integral positive source')
        h=z//3
        need(h%2==1,'inverse odd source')
    return h


def run():
    for length in range(7):
        for w in product(range(1,5), repeat=length):
            mu,P=affine(w);bound,Q=histogram(w)
            need(mu==bound and P<=Q,'all-order emission envelope')
            sorted_word=tuple(a for a in w if a>=2)+(1,)*w.count(1)
            need(affine(sorted_word)==(mu,Q),'sharpness with all ones last')
    need(F(1)/(1-F(27,32))==F(32,5),'three-letter common floor')
    need(F(27,32)**10<F(4,19)<F(27,32)**9,'ten-block uniform relaxed repair threshold')
    need(F(32,5)+F(27,32)**10*(14-F(32,5))<8,'worst original source n9 pays within30 guarded letters')
    for w in product(range(1,5),repeat=3):
        if w.count(1)<=1:
            mu,P=affine(w)
            need(mu<=F(27,32) and P<=1,'every guarded triple contracts above common floor')
    need(histogram((1,2))==(F(9,8),F(1)), 'one reset per pair is not contraction')
    need(affine((1,2,2))[1]==F(9,16)<1,'histogram deliberately forgets favourable order')
    inp=json.loads((HERE/'inputs.json').read_text(encoding='utf-8'))
    sample=next(row for row in inp['rows'] if row['exponent']==16)
    need(rejects(lambda:validate_row(dict(sample,mod15=False))),'Boolean register alias rejected')
    need(rejects(lambda:validate_row(dict(sample,mod35=float(sample['mod35'])))),'floating register alias rejected')
    need(rejects(lambda:validate_row(dict(sample,source=sample['source']+2))),'changed original source rejected')
    rows=[];counts=Counter()
    for row in inp['rows']:
        validate_row(row)
        e,n,w=row['exponent'],row['source'],tuple(row['ROOT_word'])
        need(n==(1<<e)-1,'literal retained Mersenne source')
        need(row['mod15']==n%15 and row['mod35']==n%35,'both retained modular guards')
        need(row['mod15']%5==row['mod35']%5,'compatible registers')
        need(inverse(1,w)==n,'independent ROOT inverse reconstruction')
        xs,debts=replay(n,w)
        need(xs[-1]==1,'supplied grounded ROOT child')
        j=next(i for i,x in enumerate(xs) if i and x<n);h=xs[j]
        need(all(x>=n for x in xs[:j]),'retained first-exit cut')
        gate=0 if 2*h+1<=n else None
        if gate is None:
            for k in range(1,len(w)-j+1):
                mu,P=histogram(w[j:j+k])
                if 2*h*mu+P<=n-1:
                    gate=k;break
        need(gate is not None,'declared finite bank covered by half or burst guard')
        stop=j+gate;child=xs[stop]
        need(debts[stop]<n-child,'strict original-source payment')
        need(inverse(1,w[stop:])==child,'retained child is grounded by its own suffix')
        old_tail=[]
        for a in w[j:j+3]:
            if a<2:break
            old_tail.append(a)
        old_repair=len(old_tail)==3 or xs[j+len(old_tail)]==1
        category='half' if gate==0 else 'burst'
        counts[category]+=1
        if gate and not old_repair:counts['burst_without_three_reset_gate']+=1
        r=w[j:stop].count(1);A=sum(w[j:stop])
        need((row['mod15'],row['mod35'])==(n%15,n%35),'source unchanged at child')
        rows.append({'exponent':e,'source':n,'mod15':n%15,'mod35':n%35,
                     'first_exit_step':j,'first_exit':h,'repair_length':gate,
                     'repair_ones':r,'repair_valuation_cost':A,'paid_child':child,
                     'child_ROOT_word':list(w[stop:]),'head_steps':stop,
                     'head_valuation_clock':sum(w[:stop]),'ROOT_steps':len(w),
                     'ROOT_valuation_clock':sum(w),'old_three_reset_gate':old_repair,
                     'payment_margin':str(n-child-debts[stop])})
    need([r['exponent'] for r in rows]==list(range(14,81)),'complete declared67-source universe')
    return {'status':'PROVED scoped histogram and common-floor lemmas; FINITE-EXACT supplied67-source certificates',
            'checks':CHECKS,'counts':dict(counts),'rows':rows,
            'checker_sha256_lf':sha256(Path(__file__).read_text(encoding='utf-8').encode()).hexdigest(),
            'inputs_sha256_lf':sha256((HERE/'inputs.json').read_text(encoding='utf-8').encode()).hexdigest()}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    result=run();path=HERE/'receipt.json'
    if a.verify:need(result==json.loads(path.read_text(encoding='utf-8')),'saved receipt matches exact replay')
    else:path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
