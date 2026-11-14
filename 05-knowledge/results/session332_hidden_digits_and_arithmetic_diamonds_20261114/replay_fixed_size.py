"""Exact size/phase and literal same-cost replacement controls."""
from pathlib import Path
from math import gcd
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from fermat_profiles import glue, ends, strong

F = [(2,3),(3,5),(5,17)]
R = F+[(6,33)]


def require(x, msg):
    if not x:
        raise RuntimeError(msg)


def first_full(modulus, alphabet, goal, limit=700):
    table = [{1}]
    first = None
    for m in range(1,limit):
        table.append(set().union(*(
           {a*v%modulus for v in table[m-k]} for k,a in alphabet if k<=m)))
        if first is None and len(table[-1])==goal:
            first=m
        if first is not None and m==first+3:
            return first,table
    raise RuntimeError("census did not reach the target size")


def physical(word):
    A,v=(0,),0
    for j in word:
        A=glue(A,v,j)
        v=len(A)-1
    require(strong(A),"strong literal terminal output")
    return len(A),sum(ends(A)[-1])


def main():
    for n,expected in enumerate((11,23,47,95,101,119,159,217,323),3):
        f,T=first_full(1<<n,F,1<<(n-1))
        require(f==expected and f%2==1,("first full binary size",n))
        require(all((len(T[m])==1<<(n-1))==(len(T[m-3])==1<<(n-1))
                    for m in range(4,len(T),2)), ("even selector lag",n))
        require(len(T[f+1]) < 1<<(n-1) and
                all(len(T[m]) == 1<<(n-1) for m in (f+2,f+3)),
                ("first isolated gap",n))
        if n==9:
            require(set(range(1,512,2))-T[160]=={87},"unique hole at 160")
    g=27*pow(25,-1,512)%512
    l=17*pow(15,-1,512)%512
    require((g,l,pow(g,16,512))==(83,479,l*l%512),"mod 512 order check")
    for M,miss in [(157,{37}),(159,set())]:
        u=(M-1)//2
        union={x%128 for s in range(16)
               for x in range(-15*s,(u-2*s-2)//3-16*s+1)}
        require(set(range(128))-union==miss,("interval selector",M))

    for bank,goal,expected in [(F,768,319),(R,1536,173)]:
        f,T=first_full(3584,bank,goal)
        require(f==expected,("joint first saturation",bank))
        for m,V in enumerate(T):
            if bank==F:
                require(all(gcd(v,3584)==1 and
                            (pow(v%7,3,7)==1)==((m+(v%4==3))%2==0) for v in V),
                        ("reciprocity constraint",m))
        print("joint phase first saturation",f,goal)
    multi=32*7*13
    multi_full=16*6*12
    for bank,goal in [(F,multi_full//4),(F+[(6,33)],multi_full//2),
                      (F+[(8,129)],multi_full//2),
                      (F+[(6,33),(8,129)],multi_full)]:
        f,T=first_full(multi,bank,goal)
        require(all(len(V)<=goal for V in T),("quadratic rank exact",bank))
        print("joint 7 and 13 sign rank",f,goal)
    require(pow(5,6,49)==43,"primitive seven-power lift")
    require(physical([1,1,1])==(7,27) and physical([5])==(7,33),
            "same-order replacement, opposite reciprocity sign")
    require(physical([1,2,4])==(11,255),"multiplicative ready count")
    print("exact modulo-512 sizes 159 / hole 160 / stable 161; literal same-cost repair")


if __name__=="__main__":
    main()
