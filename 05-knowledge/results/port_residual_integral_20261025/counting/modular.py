"""Fixed absolute binary precision for the integral adaptive recurrence.

Sparse keys record semantic support even when their stored residue is zero.
The production calculation uses only residues and odd modular inverses.
"""
from importlib.util import spec_from_file_location, module_from_spec
from math import comb, gcd
from pathlib import Path
from random import Random
import json


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def bound(defect):
    return max(((x + 2*i - 1)//(2*i) for i,x in enumerate(defect) if i), default=0)


def intrinsic(ell,k):
    return k*(ell-1)//2+ell*k*(k-1)//2


def beta(k,e):
    return (-1)**e*comb(k+e-1,e) if k else int(e==0)


def odd_factorials(n,modulus):
    out=[1]
    for j in range(1,n+1):
        out.append(out[-1]*(j//(j&-j))%modulus)
    return out


def shifted(value,exponent,K):
    need(exponent>=0,'nonnegative shift')
    return 0 if exponent>=K else (value<<exponent)% (1<<K)


def to_dense(s,defect,support,K):
    if not s:
        return [sum(support.values())%(1<<K)]
    b=bound(defect)
    need(all(m>=b for m in support),'valid sparse lower bound')
    return [0]+[sum(shifted(x,r*(m-b),K) for m,x in support.items())%(1<<K)
                for r in range(1,min(s,K)+1)]


def read_channel(data,s,a,K):
    need(0<=a<=s and (a>0 or s==0),'affordable source coordinate')
    return data[0] if not s else data[min(a,K)]


def integral_coeff(odd,u,ell,k,e,modulus):
    up=u+ell*(k+e)
    denominator=odd[u]*pow(ell,k,modulus)*odd[k]%modulus
    need(denominator%2==1,'all residual denominators are odd units')
    return beta(k,e)*odd[up]*pow(denominator,-1,modulus)%modulus


def next_defect(ell,rest,defect,k):
    return tuple(defect[(d-1)//2]+k*(d-gcd(ell,d))
                 for d in range(1,min(ell-2,rest)+1,2))


def dense_batch(n,ell,s,defect,data,k,e,K,odd):
    rest=s-ell*(k+e); dn=next_defect(ell,rest,defect,k)
    b,bp=bound(defect),bound(dn)
    selected=defect[(ell-1)//2] if ell<=s else 0
    j0=intrinsic(ell,k)-(k-k.bit_count())+k*(ell*b-selected)
    delta=k+b-bp
    need(j0>=0 and delta>=0,'integral constant-tail transfer')
    factor=integral_coeff(odd,n-s,ell,k,e,1<<K)
    out=[0]*(min(rest,K)+1 if rest else 1)
    for r in (range(1,min(rest,K)+1) if rest else (0,)):
        out[r]=shifted(factor*read_channel(data,s,ell*k+r,K),j0+r*delta,K)
    return (rest,dn),out


def modular_target(n,K,core=True):
    if type(n) is not int or n<0 or type(K) is not int or K<1:
        raise ValueError('n>=0 and K>=1 must be integers')
    modulus=1<<K
    if n==0:
        return 1,{'keys':1,'cells':1,'conversions':0,'zero_sparse_contributions':0,'rows':[]}
    odd=odd_factorials(n,modulus)
    top=n if n%2 else n-1
    states={(n,(0,)*((top+1)//2)):('sparse',{0:1}) if core else ('ordinary',1)}
    stats={'conversions':0,'zero_sparse_contributions':0}; rows=[]
    for ell in range(top,0,-2):
        rows.append({'length':ell,'keys':len(states),
                     'cells':sum(len(data) if kind=='sparse' else max(1,min(key[0],K)) if kind=='dense' else 1
                                 for key,(kind,data) in states.items())})
        nxt={}

        def add_sparse(key,m,x):
            rest,dn=key
            if rest and not x:
                stats['zero_sparse_contributions']+=1
            if key not in nxt:
                nxt[key]=('sparse',{}) if rest else ('dense',[0])
            kind,data=nxt[key]
            if kind=='sparse':
                # Do not drop zero residues: keys are the exact semantic support.
                data[m]=(data.get(m,0)+x)%modulus
                if len(data)>min(rest,K):
                    nxt[key]=('dense',to_dense(rest,dn,data,K))
                    stats['conversions']+=1
            elif rest:
                bp=bound(dn)
                need(m>=bp,'destination sparse lower bound')
                for r in range(1,min(rest,K)+1):
                    data[r]=(data[r]+shifted(x,r*(m-bp),K))%modulus
            else:
                data[0]=(data[0]+x)%modulus

        def add_dense(key,vector):
            rest,dn=key
            if key not in nxt:
                nxt[key]=('dense',[0]*(min(rest,K)+1 if rest else 1))
            kind,data=nxt[key]
            if kind=='sparse':
                data=to_dense(rest,dn,data,K); nxt[key]=('dense',data)
                stats['conversions']+=1
            for r in (range(1,min(rest,K)+1) if rest else (0,)):
                data[r]=(data[r]+vector[r])%modulus

        for (s,defect),(kind,data) in states.items():
            selected=defect[(ell-1)//2] if ell<=s else 0
            for k in range(s//ell+1):
                for e in (range((s-ell*k)//ell+1) if core and k else (0,)):
                    rest=s-ell*(k+e)
                    if ell==1 and rest:
                        continue
                    dn=next_defect(ell,rest,defect,k); key=(rest,dn)
                    factor=integral_coeff(odd,n-s,ell,k,e,modulus)
                    base=intrinsic(ell,k)-(k-k.bit_count())
                    if kind=='sparse':
                        for m,x in data.items():
                            add_sparse(key,m+k,shifted(factor*x,base+k*(ell*m-selected),K))
                    elif kind=='dense':
                        checked,vector=dense_batch(n,ell,s,defect,data,k,e,K,odd)
                        need(checked==key,'modular key agreement'); add_dense(key,vector)
                    else:
                        b,bp=bound(defect),bound(dn)
                        shift=base+k*(ell*b-selected)+rest*(k+b-bp)
                        need(ell*k+rest==s,'ordinary one-channel source')
                        value=shifted(factor*data,shift,K)
                        nxt[key]=('ordinary',(nxt.get(key,('ordinary',0))[1]+value)%modulus)
        # Nonterminal zero residues are retained, in both sparse and dense form.
        states={key:entry for key,entry in nxt.items()
                if key[0] or (entry[1][0] if entry[0]=='dense' else entry[1])}
    need(set(states)<={(0,())},'terminal modular key')
    kind,data=states.get((0,()),('ordinary',0))
    value=data[0] if kind=='dense' else data
    result=value*pow(odd[n],-1,modulus)%modulus
    stats.update({'keys':sum(r['keys'] for r in rows),'cells':sum(r['cells'] for r in rows),'rows':rows})
    return result,stats


def tail_controls(exact):
    rng=Random(20261025); checks=commutations=0
    for K in range(1,7):
        for ell in (3,5,7):
            for _ in range(20):
                s=ell+rng.randrange(12); defect=(0,)+tuple(rng.randrange(15) for _ in range((ell-1)//2))
                b=bound(defect); u=(ell+2)*(b+5); n=u+s
                odd=odd_factorials(n,1<<K)
                data=[0]+[rng.randrange(1<<K) for _ in range(min(s,K))]
                expanded=[0]+[data[min(r,K)] for r in range(1,s+1)]
                support={b+j:rng.randrange(1<<K) for j in range(5)}
                bank=to_dense(s,defect,support,K)
                for k in range(s//ell+1):
                    for e in (range((s-ell*k)//ell+1) if k else (0,)):
                        key,out=dense_batch(n,ell,s,defect,data,k,e,K,odd)
                        fullkey,full=exact.dense_contribution(n,ell,s,defect,expanded,k,e)
                        need(key==fullkey,'independent full key')
                        rest,dn=key
                        for r in (range(1,rest+1) if rest else (0,)):
                            need(full[r]%(1<<K)==read_channel(out,rest,r,K),'arbitrary constant-tail closure')
                        _,result=dense_batch(n,ell,s,defect,bank,k,e,K,odd)
                        dest={}
                        for m,x in support.items():
                            shift=intrinsic(ell,k)-(k-k.bit_count())+k*(ell*m-defect[-1])
                            dest[m+k]=shifted(integral_coeff(odd,u,ell,k,e,1<<K)*x,shift,K)
                        need(result==to_dense(rest,dn,dest,K),'modular sparse/dense commutation')
                        checks+=1; commutations+=1
    # A real skip transition preserves an invalid product-lattice tail (0,1).
    fullkey,full=exact.dense_contribution(2,3,2,(0,),[0,0,1],0,0)
    falsekey,false=dense_batch(2,3,2,(0,),[0,0],0,0,1,odd_factorials(2,2))
    need(fullkey==falsekey and full[2]%2==1 and read_channel(false,2,2,1)==0,
         'arbitrary product lattice does not admit false tail compression')
    need(to_dense(2,(0,),{0:1},1)==[0,1],'actual nonzero tail is retained')
    _,completion=dense_batch(2,1,2,(0,),[0,1],2,0,1,odd_factorials(2,2))
    need(completion==[1],'a high-index actual tail has a nonzero individual continuation')
    zero_support={0:0,1:0}
    need(len(zero_support)==2 and to_dense(3,(0,),zero_support,2)==[0,0,0],
         'zero residues do not remove semantic support')
    return {'arbitrary_constant_tail_transfers':checks,'sparse_dense_commutations':commutations}


def main():
    p=Path(__file__).with_name('check.py')
    spec=spec_from_file_location('exact_integer_reference',p)
    exact=module_from_spec(spec); spec.loader.exec_module(exact)
    comparisons=0; samples=[]; zero_seen=False
    for n in range(41):
        for core in (False,True):
            value,before=exact.integral_target(n,core)
            for K in range(1,7):
                residue,stats=modular_target(n,K,core)
                need(residue==value%(1<<K),('absolute residue',n,K,core))
                if n:
                    need(all(a['length']==b['length'] and a['cells']<=b['cells']
                             for a,b in zip(stats['rows'],before['rows'])),('modular storage domination',n,K,core))
                zero_seen |= bool(stats['zero_sparse_contributions'])
                comparisons+=1
                if n in (18,26,40) and core and K in (1,2,6):
                    samples.append({'n':n,'K':K,'exact_cells':before['cells'],
                                    'modular_cells':stats['cells'],'residue':residue})
    value44,_=exact.integral_target(44,True)
    for K in range(1,7):
        residue,_=modular_target(44,K,True)
        need(residue==value44%(1<<K),('order44 dense boundary',K)); comparisons+=1
    need(zero_seen,'complete runs exercise zero-residue sparse support')
    absolute2,_=modular_target(50,2,False)
    absolute27,_=modular_target(50,27,False)
    need(absolute2==0 and absolute27%(1<<25)==0 and absolute27>>25==2,
         'normalized precision requires baseline digits')
    controls=tail_controls(exact)
    for args in ((-1,2),(3,0),(True,2),(3,1.5)):
        try:
            modular_target(*args)
        except ValueError:
            pass
        else:
            raise RuntimeError('invalid modular inputs accepted')
    print(json.dumps({'status':'PASS','ordinary_and_core_orders':[0,40],
                      'absolute_precision_bits':[1,6],'exact_count_comparisons':comparisons,
                      'additional_core_order':44,'controls':controls,'samples':samples,
                      'normalized_hostile':{'n':50,'baseline_bits':25,'absolute_2_bits':absolute2,
                                            'absolute_27_bits':absolute27,'normalized_2_bits':absolute27>>25},
                      'scope':'constant-tail module and absolute residues; no positivity or minimum-rank claim'},indent=2))


if __name__=='__main__':
    main()
