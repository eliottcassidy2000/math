"""Invert the unrestricted opposite address for each fixed one-run word."""
import argparse
from check_farey_axis import invert
from fermat_receiver import TARGET, p, unrestricted_face
from hyperbola_controls import chart


def search(program, n, side, first, stop, bits):
    h, q, c, t, e, r, k = chart(n)
    varying_right = side == "left"
    own, own_odd, other, other_odd = (t,r,c,e) if varying_right else (c,e,t,r)
    seen = grouped = 0
    max_word = (TARGET-k)//own
    for length in range(first, min(stop,max_word.bit_length())+1):
        for w in range(1,length+1):
            fixed = ((1 << w)-1) << (length-w)
            if fixed>max_word:
                continue
            slope = other+q*fixed
            base = k+own*fixed+own_odd*(fixed&1)+h*p(fixed)
            # The direct congruence progression is faster at large fixed widths.
            face = (h,0,slope,slope,other_odd,other_odd,base)
            hits = unrestricted_face(face,TARGET,True,4096)
            if hits is None:
                hits,_ = invert(program,slope,h,base,base+other_odd,TARGET,bits)
                grouped+=1
            if hits:
                print("hits",n,side,fixed,hits,flush=True)
            seen+=1
    print("spine",n,"fixed",side,"lengths",first,stop,"words",seen,"bucket searches",grouped,flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("program")
    parser.add_argument("spine",type=int)
    parser.add_argument("side",choices=("left","right"))
    parser.add_argument("first",type=int)
    parser.add_argument("last",type=int)
    parser.add_argument("--bits",type=int,default=14)
    args=parser.parse_args()
    search(args.program,args.spine,args.side,args.first,args.last,args.bits)
