"""All positive chain lengths on a streamed multilinear profile at an exact target."""
import argparse
from itertools import permutations, product
from math import gcd
import sys


def order_coeff(f, perm):
    return [f[sum((mask >> i & 1) << perm[i] for i in range(3))] for mask in range(8)]


def all_widths(f, target):
    found = set()
    for perm in permutations(range(3)):
        d,b,c,bc,a,ab,ac,abc = order_coeff(f,perm)
        for m in range(1,target.bit_length()+1):
            x = 1 << m
            if (target - d) % x:
                break
            intercept, linear = d + b*x, c+bc*x
            slope, mixed = a+ab*x, ac+abc*x
            if intercept + (linear+slope)*x + mixed*x*x > target:
                break
            for k in range(m,target.bit_length()+1):
                y = 1 << k
                if (target-intercept) % y:
                    break
                den = slope + mixed*y
                if den <= 0:
                    raise ValueError("Strong port must have positive slope")
                residual = target - intercept - linear*y
                if residual < den*y:
                    break
                q, rem = divmod(residual, den)
                if rem or q & (q-1):
                    continue
                widths = [0,0,0]
                for i,w in zip(perm,(m,k,q.bit_length()-1)): widths[i]=w
                found.add(tuple(widths))
    return found


def multi_widths(f, target, minimums=None):
    K = len(f).bit_length()-1
    if minimums is None:
        minimums = [1]*K
    widths = [0]*K
    out = set()
    def visit(coeff, ids, last, minimum):
        delta=target-coeff[0]
        if delta % gcd(*coeff[1:]):
            return
        if len(ids)==1:
            power,rem=divmod(delta,coeff[1])
            if not rem and power>0 and power & (power-1)==0:
                w=power.bit_length()-1
                if w>=minimums[ids[0]] and (w>minimum or (w==minimum and ids[0]>last)):
                    widths[ids[0]]=w
                    out.add(tuple(widths))
            return
        # The minimum remaining width divides target - constant in the 2-adic sense.
        degree = [0]*(len(ids)+1)
        for mask,value in enumerate(coeff):
            degree[mask.bit_count()]+=value
        pairs=[]
        for i in range(len(ids)):
            low=[c for mask,c in enumerate(coeff) if not mask>>i&1]
            high=[c for mask,c in enumerate(coeff) if mask>>i&1]
            pairs.append((low,high))
        for m in range(minimum,target.bit_length()+1):
            x=1<<m
            if delta % x:
                break
            if all(m>=minimums[idx] for idx in ids):
                if sum(c << (m*j) for j,c in enumerate(degree))>target:
                    break
            for i,idx in enumerate(ids):
                if m<minimums[idx] or (m==minimum and idx<=last):
                    continue
                widths[idx]=m
                low,high=pairs[i]
                reduced=[a+b*x for a,b in zip(low,high)]
                visit(reduced, ids[:i]+ids[i+1:], idx, m)
    visit(f, list(range(K)), -1, 1)
    return out


def compositions(total, size):
    if size == 1:
        if total >= 1: yield (total,)
        return
    for x in range(1,total-size+2):
        for tail in compositions(total-x,size-1):
            yield (x,)+tail


def solve(f, target):
    k=len(f).bit_length()-1
    # Indistinguishable incidence columns in the monomial support enter
    # solely via their product. Add their widths, then recover all splits.
    groups={}
    support=[mask for mask,c in enumerate(f) if c]
    for i in range(k):
        signature=tuple((mask>>i)&1 for mask in support)
        groups.setdefault(signature,[]).append(i)
    parts=list(groups.values())
    u=len(parts)
    reduced=[0]*(1<<u)
    for mask in support:
        tag=sum(((mask>>group[0])&1)<<i for i,group in enumerate(parts))
        reduced[tag]=f[mask]
    if u==3 and all(len(g)==1 for g in parts): compact=all_widths(reduced,target)
    else: compact=multi_widths(reduced,target,[len(g) for g in parts])
    answer=set()
    for widths in compact:
        for splits in product(*(compositions(w,len(group)) for w,group in zip(widths,parts))):
            result=[0]*k
            for group,vals in zip(parts,splits):
                for j,w in zip(group,vals):result[j]=w
            answer.add(tuple(result))
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exponent", type=int)
    parser.add_argument("--winning-coefficients", action="store_true")
    parser.add_argument("--upper-only", action="store_true")
    args = parser.parse_args()
    p = args.exponent
    mid = (1<<(p-1))*((1<<p)-1)
    hits = [0,0]
    for line in sys.stdin:
        board, *digits = line.split()
        K=next(k for k in range(1,10) if k+(1<<k)==len(digits))
        ports=tuple(map(int,digits[:K]))
        f=list(map(int,digits[K:]))
        printed = False
        for i,target in enumerate([mid-1,mid+1]):
            if args.upper_only and i == 0:
                continue
            widths = solve(f,target)
            if widths:
                hits[i] += len(widths)
                if args.winning_coefficients:
                    if not printed:
                        print(*f,flush=True)
                        printed = True
                else:
                    print(board,*ports,target,sorted(widths),f,flush=True)
    print("Exponent",p,"signed width-solution counts",hits,file=sys.stderr)


if __name__ == "__main__":
    main()
