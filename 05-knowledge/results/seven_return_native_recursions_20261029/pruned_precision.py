"""Exact precision recursion with a bound on unfinished multiplicity refunds."""
from math import gcd
from baseline_precision import base_valuation,cycle_floor,term

def count(n,bits,r=1,s=0,mode='completion'):
    if any(type(x) is not int for x in (n,bits,r,s)) or n<0 or bits<1 or r<1 or s<0 or mode not in ('none','collision','completion','sum'):
        raise ValueError('invalid precision request')
    stats=dict(nodes=0,pruned_nodes=0,completed=0,retained=0)
    if n<=2:return 1,dict(stats,base=base_valuation(n,r,s),cycle_cap=n)
    odd=n%2;ks=[1 if odd else 2];k=ks[0]+2
    while k<=n and cycle_floor(k,r,s,bool(odd))<bits:ks.append(k);k+=2
    twos=[j-j.bit_count() for j in range(ks[-1]+1)]
    caps=[]
    for size in range(ks[-1]+1):
        row=[-1]
        for v in range(1,size+1):row.append(max(row[-1],twos[size-v]+twos[v]))
        caps.append(row)
    modulus=1<<bits;answer=0
    for k in ks:
        T=r*(k-1)**2//2+s*(k-1) if odd else r*k*(k-2)//2+s*(k-2)
        def rec(total,left,parts,run,frozen,delta):
            nonlocal answer
            stats['nodes']+=1
            refund=(frozen+twos[run+left]) if parts else twos[left]
            if mode=='sum' and parts and total>left*parts[-1]:
                limit=min(left,(total-left*parts[-1])//2)
                if limit<1:raise ArithmeticError('nonpositive odd-part excess capacity')
                refund=frozen+caps[run+left][limit]
            ceiling=refund if mode in ('completion','sum') else twos[k]
            if mode!='none' and T+r*delta-ceiling>=bits:
                stats['pruned_nodes']+=1;return
            if left==0:
                if total:raise ArithmeticError('unfinished partition')
                stats['completed']+=1
                exponent,denominator=term(parts,r,s,bool(odd))
                if exponent!=T+r*delta-frozen-twos[run]:raise ArithmeticError('incremental valuation mismatch')
                if exponent<0:raise ArithmeticError('negative normalized exponent')
                if exponent<bits:
                    answer+=(1<<exponent)*pow(denominator,-1,modulus)
                    stats['retained']+=1
                return
            lower=parts[-1] if parts else 1
            choices=([total] if total>=lower and total%2 else []) if left==1 else range(lower,total//left+1,2)
            for a in choices:
                dd=delta+sum(gcd(a,x)-1 for x in parts)
                if parts and a==parts[-1]:rr=run+1;ff=frozen
                else:rr=1;ff=frozen+twos[run]
                rec(total-a,left-1,parts+(a,),rr,ff,dd)
        rec(n,k,(),0,0,0)
    return answer%modulus,dict(stats,base=base_valuation(n,r,s),cycle_cap=ks[-1])

if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser();p.add_argument('n',type=int);p.add_argument('bits',type=int)
    a=p.parse_args();value,stats=count(a.n,a.bits)
    print(json.dumps(dict(residue=value,stats=stats),indent=2))
