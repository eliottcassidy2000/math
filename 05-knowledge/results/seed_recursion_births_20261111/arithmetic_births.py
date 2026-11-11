"""Exact node-birth profiles; independent pair enumeration checks formulas."""
from pathlib import Path
from math import isqrt,comb
import json

HERE=Path(__file__).resolve().parent


def require(ok,message):
    if not ok:raise RuntimeError(message)


def kernel(n):
    out=1;p=2
    while p*p<=n:
        parity=0
        while n%p==0:n//=p;parity^=1
        if parity:out*=p
        p+=1
    return out*n


def square(n):return isqrt(n)**2==n


def main():
    limit=4096
    tau=[0]*(limit+1)
    for d in range(1,limit+1):
        for n in range(d,limit+1,d):tau[n]+=1
    records=[];D=0;weak_total=strict_total=0;classes={};edges=0
    for n in range(1,limit+1):
        D+=tau[n]
        if n==1:weak=strict=0
        else:
            weak=(tau[n]-2+square(n))//2
            strict=(tau[n]-2-square(n))//2
        weak_total+=weak;strict_total+=strict
        require(weak_total==(D-2*n+isqrt(n))//2,'cumulative weak factors')
        require(strict_total==(D-2*n-isqrt(n)+2)//2,'cumulative strict factors')
        d=kernel(n);a=isqrt(n//d)
        require(classes.get(d,0)==a-1,'square-product older neighbours')
        classes[d]=a;edges+=a-1
        require(edges==sum(comb(v,2) for v in classes.values()),'clique edge total')
        records.append(dict(n=n,sum_strict=(n-1)//2,sum_weak=n//2,
                            product_strict=strict,product_weak=weak,
                            product_strict_atom=n>=2 and strict==0,
                            square_sum_birth_degree=isqrt(2*n-1)-isqrt(n),
                            square_product_birth_degree=a-1,squarefree_kernel=d,
                            square_product_components=len(classes),square_product_edges=edges))
    pair_controls=0
    for n in range(1,129):
        r=records[n-1]
        require(r['sum_strict']==sum(a+b==n for a in range(1,n) for b in range(a+1,n)),'direct sums')
        require(r['product_strict']==sum(a*b==n for a in range(2,n) for b in range(a+1,n)),'direct proper factors')
        require(r['square_sum_birth_degree']==sum(square(n+x) for x in range(1,n)),'square-sum birth')
        require(r['square_product_birth_degree']==sum(square(n*x) for x in range(1,n)),'square-product birth')
        pair_controls+=4
    # Independent causal closure from {2,3}, using distinct parents.
    additive={2,3};multiplicative={2,3}
    for n in range(4,limit+1):
        if any(a in additive and n-a in additive for a in range(1,(n+1)//2)):
            additive.add(n)
        if any(n%a==0 and a in multiplicative and n//a in multiplicative for a in range(2,isqrt(n)+1) if a*a<n):
            multiplicative.add(n)
    smooth={2,3}
    a=2
    while a<=limit:
        b=3
        while a*b<=limit:smooth.add(a*b);b*=3
        a*=2
    require(additive==set(range(1,limit+1))-{1,4,6},'strict additive closure')
    require(multiplicative==smooth,'strict multiplicative closure')
    report=dict(limit=limit,recorded_prefix=512,direct_pair_controls=pair_controls,records=records[:512],at_limit=records[-1],
                strict_sum_closure_holes=[1,4,6],strict_product_closure_size=len(smooth),
                square_sum_components_at_14=1)
    # Check the connectivity boundary without relying on its stored value.
    comp=[]
    for n in range(1,33):
        unseen=set(range(1,n+1));parts=[]
        while unseen:
            todo=[unseen.pop()];part=set(todo)
            while todo:
                v=todo.pop();new={w for w in unseen if square(v+w)}
                unseen-=new;part|=new;todo.extend(new)
            parts.append(sorted(part))
        comp.append(dict(n=n,components=parts))
    require(len(comp[12]['components'])==2 and len(comp[13]['components'])==1,'13/14 connectivity boundary')
    report['square_sum_component_controls']=comp
    (HERE/'arithmetic_births.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(limit=limit,pair_controls=pair_controls,at_4096=records[-1]),indent=2))


if __name__=='__main__':main()
