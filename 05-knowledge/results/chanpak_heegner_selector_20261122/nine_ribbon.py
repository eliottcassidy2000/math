"""Choose a common-vertex pair of unit seams with prescribed separation 2*k."""
from argparse import ArgumentParser
from dataclasses import dataclass
from functools import cache
from pathlib import Path
import sys

RESULTS=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(RESULTS/"chanpak_ternary_buckets_20261122"))
from physical_navigation import direct, flipped, star_rows


@cache
def read(x):
    # Iteration keeps very large input addresses independent of Python recursion depth.
    value=0
    while x:
        high,bit=divmod(x,2)
        if bit==0:
            value+=high+(high&1)
        x=high
    return value


def reverse(x,length):
    return int(format(x,f"0{length}b")[::-1],2)


def midpoint(x):
    return 370+388*x+15*read(x)


@dataclass(frozen=True,order=True)
class Template:
    first:int  # least value of (midpoint(x_a)-midpoint(x_b))/2 on this progression
    length:int
    a:int
    b:int
    parity:int
    d:int
    intercept:int


@cache
def bank():
    best={}
    for length in range(1,9):
        modulus=1<<length
        for i in range(modulus):
            a=reverse(i,length)
            for d in (1,3,5,15):
                if i+d>=modulus:
                    continue
                b=reverse(i+d,length)
                if (a^b).bit_count()>4:
                    continue
                for parity in (0,1):
                    top_difference=int(b<modulus//2)-int(a<modulus//2)
                    twice=15*d*parity+388*(a-b)+15*(read(a)-read(b)-parity*top_difference)
                    assert twice%2==0
                    intercept=twice//2
                    step=15*d
                    first=intercept+step*int(parity==0)
                    template=Template(first,length,a,b,parity,d,intercept)
                    key=(step,intercept%step)
                    if key not in best or template<best[key]:
                        best[key]=template
    return best


def progression(r):
    return [bank()[step,r%step] for step in (15,45,75,225) if (step,r%step) in bank()]


def choose(k):
    if k<1:
        raise ValueError("positive half-distance required")
    solutions=[]
    for m in (k,-k):
        for row in progression(m):
            if m<row.first:
                continue
            v,rem=divmod(m-row.intercept,15*row.d)
            assert rem==0
            u=2*v+row.parity
            assert u>=1
            xa=(u<<row.length)+row.a
            xb=(u<<row.length)+row.b
            if m<0:
                xa,xb=xb,xa
            assert midpoint(xa)-midpoint(xb)==2*k
            solutions.append((xa.bit_length(),midpoint(xb),xb,xa))
    return min(solutions)


def shell(x):
    code,cut,apex=0x7103d59dc,134,7
    rows=[0]*9
    from itertools import combinations
    for j,(a,b) in enumerate(combinations(range(9),2)):
        if ((code>>j)^(cut>>a)^(cut>>b))&1==0:
            a,b=b,a
        rows[a]|=1<<b
    ordinary=[i for i in range(9) if i!=apex]
    def project(mask):
        return sum(((mask>>i)&1)<<j for j,i in enumerate(ordinary))
    high=star_rows(tuple(project(rows[i]) for i in ordinary),project(rows[apex]),x,0)
    offset=x.bit_length()
    return flipped(high,offset+1,offset+2),high


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("k",nargs="?",type=int)
    parser.add_argument("--rows",action="store_true")
    args=parser.parse_args()
    # For all m above last, the positive progression suffices in every residue.
    firsts=[min(r.first for r in progression(i)) for i in range(225)]
    last=max(i+225*((firsts[i]-1-i)//225) for i in range(225))
    assert last==3030
    assert all(choose(k) for k in range(1,last+1))
    for x in (1,3,5):
        low,high=shell(x)
        assert (direct(low),direct(high))==(midpoint(x)-1,midpoint(x)+1)
    print("225 residue progressions cover all positive signed half-distances above",last)
    if args.k is not None:
        width,h,x,y=choose(args.k)
        changed=[j for j in range(width) if ((x^y)>>j)&1]
        assert len(changed)<=4
        print("half-distance",args.k,"addresses",x,y,"order",9+width,"bits from right",changed)
        print("heights",h-1,h+1,midpoint(y)-1,midpoint(y)+1)
        if args.rows:
            print("lower and upper outgoing rows at each address")
            for address in (x,y):
                for rows in shell(address):
                    print(*rows)


if __name__=="__main__":
    main()
