"""Invert the exact chart along a prescribed sum of its two physical addresses."""
from argparse import ArgumentParser
from itertools import combinations
from math import isqrt
from random import Random
from fermat_receiver import TARGET, address_height
from hyperbola_controls import chart


def possible_x(data, target, total):
    h, q, c, t, e, r, k = data
    linear = q*total+t-c
    lower = linear*linear - 4*q*(target-k-c*total)
    upper = lower + 4*q*(h*total+r+e)
    if upper < 0:
        return
    inner = isqrt(lower-1)+1 if lower>0 else 0
    outer = isqrt(upper)
    step = 2*q
    for start,end in ((-outer,-max(inner,1)),(inner,outer)):
        lo = max(0, -(-(linear+start)//step))
        hi = min(total, (linear+end)//step)
        yield from range(lo,hi+1)


def lower_sum(data,target):
    h,q,c,t,e,r,k=data
    lo,hi=0,1+2*isqrt(target//q+1)
    while lo<hi:
        mid=(lo+hi)//2
        if k+q*(mid*mid//4)+(max(c,t)+h)*mid+r+e>=target:
            hi=mid
        else:
            lo=mid+1
    return lo


def sparse_sums(low,high):
    for width in range(max(1,low.bit_length()),high.bit_length()+1):
        for count in range(4):
            for changes in combinations(range(1,width),count):
                total=(1<<width)-1
                for bit in changes:
                    total^=(1<<bit)-1
                if low<=total<=high:
                    yield total


def search(n, suffix):
    data=chart(n)
    h,q,c,t,e,r,k=data
    low=lower_sum(data,TARGET)
    if suffix is None:
        high=(TARGET-k)//min(c,t)
        sums=sparse_sums(low,high)
    else:
        # For x:y between 1:2 and 2:1 we have xy>=2*(x+y)^2/9.
        high=isqrt(9*(TARGET-k)//(2*q))
        stride=1<<suffix
        first=low+(-1-low)%stride
        sums=range(first,high+1,stride)
    seen=possible=0
    hits=[]
    for total in sums:
        seen+=1
        for x in possible_x(data,TARGET,total):
            y=total-x
            if suffix is not None and max(x,y)>2*min(x,y):
                continue
            residual=TARGET-k-q*x*y-t*x-c*y-e*(y&1)-r*(x&1)
            if residual%(2*h):
                continue
            possible+=1
            difference = address_height(data,x,y)-TARGET
            if suffix is None:
                print("envelope survivor",total,x,y,"height minus target",difference,flush=True)
            if difference==0:
                hits.append((x,y))
    print("spine",n,"terminal ones",suffix,"sums",seen,
          "digital candidates",possible,"hits",hits,flush=True)


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("spine",type=int,nargs="?")
    parser.add_argument("--terminal",type=int)
    args=parser.parse_args()
    if args.spine is None:
        rng=Random(61)
        for n in (4,8):
            data=chart(n)
            for i in range(80):
                total=rng.randrange(1,100)
                x=rng.randrange(total+1)
                target=address_height(data,x,total-x)+i%3-1
                found=[z for z in possible_x(data,target,total)
                       if address_height(data,z,total-z)==target]
                assert found==[z for z in range(total+1) if address_height(data,z,total-z)==target]
        assert sorted(sparse_sums(1,1000))==[v for v in range(1,1001) if len([r for r in bin(v)[2:].split("0") if r])<=2]
        print("Antidiagonal inverse controls passed.")
    else:
        search(args.spine,args.terminal)
