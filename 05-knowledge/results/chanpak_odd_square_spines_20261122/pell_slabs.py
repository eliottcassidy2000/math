"""Clear whole E61 rows around the irrational Pell directions using the native inverse."""
from argparse import ArgumentParser
from math import isqrt
from pathlib import Path
from subprocess import run

from rational_carries import chart, TARGET


def interval(n, square_num, square_den, width):
    h,q,c,t,e,r,k=chart(n)
    assert max(c+h,t+h,r+e)<q
    numer=(TARGET-k)**2*square_den
    denom=q*q*square_num
    center=isqrt(isqrt(numer//denom))
    assert denom*center**4 <= numer < denom*(center+1)**4
    low,high=center-2*width-3,center+2*width+3
    # alpha >= 1/2. Throughout this interval y/2 is larger than x/10,
    # even after the largest branch offsets i<=9, j<=1 are subtracted.
    assert 5*(low//2-width-1)>high
    return center,(low-9)//10,high//10


if __name__ == "__main__":
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("binary",type=Path)
    ap.add_argument("spine",type=int,choices=(4,8))
    ap.add_argument("square_numerator",type=int,choices=(1,2))
    ap.add_argument("--bits",type=int,default=36)
    args=ap.parse_args()
    center,low,high=interval(args.spine,args.square_numerator,3-args.square_numerator,1<<args.bits)
    block=40960 if args.spine==4 else 16000
    result=run([args.binary.resolve(),str(args.spine),str(low),str(high),str(block)],
               check=True,capture_output=True,text=True)
    lines=result.stdout.splitlines()
    assert len(lines)==8 and all(line.startswith("branch ") and line.endswith(" hits 0") for line in lines),result.stdout
    assert all(int(line.split()[3])==high-low+1 for line in lines)
    print(args.spine,args.square_numerator,center,low,high,"rows",high-low+1,
          "height tests",sum(int(line.split()[9]) for line in lines),"hits 0")
