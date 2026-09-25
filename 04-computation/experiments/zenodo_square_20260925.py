"""Pythagorean square paths and bounded reflection clocks; standard library."""
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path


def check(ok, label="check failed"):
    if not ok:
        raise RuntimeError(label)


def square(n):
    return n >= 0 and isqrt(n)**2 == n


def reflect(x, s):
    return s*s-x


def ppt(bound):
    ans=[]
    for m in range(2,isqrt(bound)+1):
        for n in range(1,m):
            if gcd(m,n)==1 and (m-n)%2 and m*m+n*n<=bound:
                ans.append((m*m-n*n,2*m*n,m*m+n*n))
    return sorted(ans,key=lambda t:(t[2],t[0]))


def six(t):
    a,b,c=t
    return [a*a,b*b,a*a+2*c+1,b*b-2*c-1,a*a+2,b*b-2]


def path_prefix(t):
    """Independent actual reflection walk, stopped before first invalid vertex."""
    a,b,c=t
    path=[a*a]
    seen=set(path)
    roots=(c,c+1,c,c-1)
    while True:
        y=reflect(path[-1],roots[(len(path)-1)%4])
        if y in seen or not 1<=y<c*c:
            return path,y,"repeat" if y in seen else "outside"
        seen.add(y);path.append(y)


def main():
    lines=[]
    def report(s):lines.append(str(s))
    triples=ppt(10000)
    check(len(triples)==1593)
    direct=[]
    for c in range(2,301):
        for a in range(1,c,2):
            bb=c*c-a*a
            b=isqrt(bb)
            if b>0 and b%2==0 and b*b==bb and gcd(a,b)==1:
                direct.append((a,b,c))
    check(sorted(direct)==sorted(ppt(300)))
    report("ZENODO SQUARE PATH AUDIT 2026-09-25")
    report(f"Universe: all{len(triples)} primitive triples with c<=10000; independent square-root enumeration c<=300; complete reflection prefixes c<=100.")
    extra=[];q_rows=[];rigid=[];inradius=[]
    for a,b,c in triples:
        t=(a,b,c);p=six(t)
        check(a*a+b*b==c*c and gcd(a,b)==1)
        check(b%4==0 and c%4==1)
        check([x%8 for x in p]==[1,0,4,5,3,6])
        check(len(set(p))==6 and min(p)>0 and max(p)<=c*c-5)
        check([p[i]+p[i+1] for i in range(5)]==[c*c,(c+1)**2,c*c,(c-1)**2,c*c])
        actual=[(i,j) for i in range(6) for j in range(i+1,6) if square(p[i]+p[j])]
        expected=[(i,i+1) for i in range(5)]
        if square(2*(a*a+1)):
            expected.append((0,4));extra.append(t)
        check(sorted(actual)==sorted(expected))
        walk=[p[0]]
        for s in (c,c+1,c,c-1,c):walk.append(reflect(walk[-1],s))
        check(walk==p)
        # Exact inverse of the path seed; no general edge-to-PPT claim.
        check((isqrt(p[0]),isqrt(p[1]),isqrt(p[0]+p[1]))==t)
        check(not square(p[4]))
        X,H=Fraction(a*a,c*c),Fraction(a*b,c*c)
        check((X-Fraction(1,2))**2+H*H==Fraction(1,4))
        check((Fraction(1)-X)==Fraction(b*b,c*c))
        r=Fraction(a+b-c,2);s=Fraction(a+b+c,2)
        check(Fraction(a*b,2)==r*s)
        if (c+b)+c==2*(a+b):rigid.append(t)
        if r==1:inradius.append(t)
        M=b*b//2-c
        L=(b*b-a*a-2*c-1)//2
        check(2*L==b*b-a*a-2*c-1 and L%2==0)
        h=L//2
        q=M if h<0 else h
        check(M>=3 and (h<0 or 0<h<M))
        # Endpoint arithmetic tests are separate from complete prefix replay.
        j=q
        A=a*a+2*j;B=b*b-2*j;C=A+2*c+1;D=b*b-2*c-1-2*j
        if h>=0:
            check(B==C==(c+1)**2//2 and A==D)
        else:
            check(B==2*c and C==c*c+1 and D==-1)
        if c<=100:
            prefix,bad,reason=path_prefix(t)
            check(len(prefix)==4*q+2)
            check(reason==("repeat" if h>=0 else "outside"))
            check(bad==(B if h>=0 else c*c+1))
            q_rows.append((t,M,h,q,len(prefix),reason,bad))
    check(rigid==[(3,4,5)] and inradius==[(3,4,5)])
    report(f"Six-vertex path checks:1593; chorded controls={extra}")
    report("Root rigidity and inradius-one controls each exactly[(3,4,5)] in finite universe; all-height proofs in note.")
    report(f"Complete independent prefix controls={q_rows}")
    root=six((3,4,5));prefix,bad,reason=path_prefix((3,4,5))
    check(root==[9,16,20,5,11,14])
    check(prefix==[9,16,20,5,11,14,22,3,13,12,24,1,15,10])
    check(bad==26 and reason=="outside" and reflect(26,5)==-1)
    check(not square(9+5) and square(5+4) and square(11+14))
    report(f"Root path={root}; maximal simple prefix insideQ24={prefix}; next26 then-1.")
    # General second-difference macro, independent direct arithmetic.
    macro_count=0
    for c in range(2,41):
        for d in range(1,c):
            for x in (-5,0,1,c,c*c):
                y=x
                for s in (c,c+d,c,c-d):y=reflect(y,s)
                check(y==x+2*d*d);macro_count+=1
    report(f"General gap-d macro checks={macro_count}; translation2d^2, not multiplication.")
    x=9
    d2walk=[x]
    for s in (5,7,5,3):
        x=reflect(x,s);d2walk.append(x)
    check(d2walk==[9,16,33,-8,17])
    report(f"Gap-d2 positivity hostile={d2walk}")
    check(six((7,24,25))[0]+six((7,24,25))[4]==100)
    report("First induced-path hostile by hypotenuse:(7,24,25), added49--51 edge has sum100.")
    check(extra[0]==(7,24,25))
    report("ALL CHECKS PASSED. A local square path is not Hamiltonicity, and its first macro leaves integer-square PPT seeds.")
    output="\n".join(lines)+"\n"
    out=Path(__file__).resolve().parents[2]/"05-knowledge/results/zenodo_square_20260925.out"
    out.write_text(output,encoding="utf-8")
    print(output,end="")


if __name__=="__main__":main()
