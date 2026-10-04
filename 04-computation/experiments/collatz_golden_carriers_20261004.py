#!/usr/bin/env python3
"""Exact golden-lattice gates and guarded Collatz carrier controls."""
from collections import Counter
from fractions import Fraction
from math import gcd, isqrt


def check(ok, message="check failed"):
    if not ok:
        raise AssertionError(message)


def sign_phi(a,b):
    # Sign of a+b*phi = ((2a+b)+b*sqrt(5))/2, with integer a,b.
    u,v=2*a+b,b
    if v == 0:
        return (u>0)-(u<0)
    if u == 0 or u*v >= 0:
        return 1 if v>0 else -1
    d=u*u-5*v*v
    return (1 if u>0 else -1) * (1 if d>0 else -1)


def beta(state,q):
    a,b=state
    digit=int(sign_phi(b-q,a+b)>0)  # Upper convention: equality -> 0.
    return (b-digit*q,a+b),digit


def trap(q):
    states=set()
    for b in range(-2*q,2*q+1):
        for a in range(-5*q,5*q+1):
            if gcd(gcd(abs(a),abs(b)),q) != 1:
                continue
            if sign_phi(a,b)<0 or sign_phi(a-q,b)>0:
                continue
            # Conjugate x' = (a+b-b*phi)/q, abs(x') <= 3.
            if sign_phi(a+b+3*q,-b)<0 or sign_phi(a+b-3*q,-b)>0:
                continue
            states.add((a,b))
    check(all(beta(s,q)[0] in states for s in states),"trap not invariant")
    return states


def cycles_in(states,q):
    done=set(); cycles=[]
    for start in sorted(states):
        if start in done:
            continue
        path=[]; positions={}; x=start
        while x not in done and x not in positions:
            positions[x]=len(path);path.append(x);x=beta(x,q)[0]
        if x in positions:
            cycles.append(path[positions[x]:])
        done.update(path)
    return cycles


def raw(n):
    return n/2 if n.numerator%2 == 0 else 3*n+1


def decode_cycle(cycle,q):
    A=p=B=0;digits=[]
    for state in cycle:
        d=beta(state,q)[1];digits.append(d)
        if d:
            B=3*B+2**A;p+=1
        else:
            A+=1
    D=2**A-3**p
    check(D != 0)
    root=Fraction(B,D);n=root;orbit=[]
    for d in digits:
        check(n.denominator%2==1 and n.numerator%2==d,"invalid parity")
        orbit.append(n);n=raw(n)
    check(n==root,"cycle failed")
    return dict(A=A,p=p,B=B,D=D,word="".join(map(str,digits)),orbit=orbit)


def odd_step(n):
    m=3*n+1;a=0
    while m%2==0:
        m//=2;a+=1
    return m,a


def summary(word):
    A=p=B=0
    for a in word:
        B=3*B+2**A;A+=a;p+=1
    return A,p,B


def main():
    # A second comparison path uses rational isolating intervals for sqrt(5).
    for a in range(-100,101):
        for b in range(-100,101):
            u,v=2*a+b,b
            lo,hi=Fraction(2236067977499789,10**15),Fraction(2236067977499790,10**15)
            lower=u+v*(lo if v>=0 else hi);upper=u+v*(hi if v>=0 else lo)
            expected=1 if lower>0 else -1 if upper<0 else 0
            check(sign_phi(a,b)==expected,"surds disagree")
    print("Exact surd comparison: 40401 pairs agree with rational sqrt(5) isolation")
    all_cycle_records=[]
    for q in (1,2,11,76):
        states=trap(q);cycles=cycles_in(states,q)
        integral=[]
        for cyc in cycles:
            record=decode_cycle(cyc,q)
            if all(x.denominator==1 for x in record["orbit"]):
                root=min(record["orbit"],key=lambda x:abs(x))
                integral.append(int(root))
            all_cycle_records.append((q,len(cyc),record))
        print("Golden exact denominator",q,"trap states",len(states),"cycles",len(cycles),
              "period counts",sorted(Counter(map(len,cycles)).items()),"integer cycle roots",sorted(integral))
        check(set(integral)==({0,-1} if q==1 else {1} if q==2 else {-5} if q==11 else {-17}))
    # Independent route: enumerate all primitive cyclic legal words through length 18.
    # This uses no lattice points, no beta graph and no quadratic sign tests.
    def pmul(x,y):
        a,b=x;c,d=y
        return (a*c+b*d,a*d+b*c+b*d)
    def pinv(x):
        a,b=x;norm=a*a+a*b-b*b
        return (Fraction(a+b,norm),Fraction(-b,norm))
    def ppow(x,n):
        y=(1,0)
        for _ in range(n):y=pmul(y,x)
        return y
    from math import lcm
    necklace_counts=Counter()
    for L in range(1,19):
        for mask in range(1<<L):
            if mask & (mask<<1) or (mask&1 and mask>>(L-1)&1):
                continue
            word="".join(str(mask>>j&1) for j in range(L))
            if word != min(word[j:]+word[:j] for j in range(L)):
                continue
            if any(L%d==0 and word==word[:d]*(L//d) for d in range(1,L)):
                continue
            num=(0,0)
            for j,d in enumerate(word):
                if d=="1":
                    power=ppow((-1,1),j+1)
                    num=(num[0]+power[0],num[1]+power[1])
            power=ppow((-1,1),L)
            theta=pmul(num,pinv((1-power[0],-power[1])))
            q=lcm(theta[0].denominator,theta[1].denominator)
            if q in (1,2,11,76):necklace_counts[(q,L)]+=1
    lattice_counts=Counter((q,L) for q,L,_ in all_cycle_records)
    check(necklace_counts==lattice_counts,"necklace/lattice census disagreement")
    print("Independent primitive cyclic-word census through length 18 agrees:",sorted(necklace_counts.items()))
    # Integer golden denominators are the Smith factors of G^L-I.
    G=(0,1,1,1)
    def mmul(M,N):
        a,b,c,d=M;e,f,g,h=N
        return (a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h)
    smith=[]
    for L in (2,3,5,18):
        M=(1,0,0,1)
        for _ in range(L):M=mmul(M,G)
        a,b,c,d=M[0]-1,M[1],M[2],M[3]-1
        g=gcd(gcd(abs(a),abs(b)),gcd(abs(c),abs(d)))
        smith.append((L,g,abs(a*d-b*c)//g))
    check(smith==[(2,1,1),(3,2,2),(5,1,11),(18,76,76)])
    print("Smith factors of golden G^L-I:",smith)
    check(pow(4,5,11)==1 and all(pow(4,k,11)!=1 for k in range(1,5)))
    check(pmul((-1,7),(-4,1))==(11,-22))
    check(pmul((1,4),(-4,1))==(0,-11))
    print("Same norm-11 denominator ideal: theta(-5)=(-1+7phi)/11 and theta(1/13)=(1+4phi)/11")
    # Exact theta roots, then propagate backwards through ordinary integer paths.
    roots={1:(0,1,2),-1:(1,0,1),-5:(-1,7,11),-17:(9,41,76),0:(0,0,1)}
    memo={}
    cycle_table=[]
    for root,t in roots.items():
        n=root;state=t;orbit=[]
        while n not in orbit:
            orbit.append(n);memo[n]=(state,root)
            a,b,q=state
            nxt,d=beta((a,b),q)
            check(d==n%2,"golden root digit mismatch")
            state=(*nxt,q);n=n//2 if n%2==0 else 3*n+1
        check(n==root and state==t)
        odd=sum(x%2 for x in orbit);halves=len(orbit)-odd
        cycle_table.append((root,len(orbit),odd,halves,t))
    counts={r:0 for r in roots};max_prefix=0
    for n0 in range(-10000,10001):
        n=n0;path=[]
        while n not in memo and len(path)<10000:
            path.append(n);n=n//2 if n%2==0 else 3*n+1
        check(n in memo,"unresolved orbit in bounded census")
        max_prefix=max(max_prefix,len(path))
        state,root=memo[n]
        for source in reversed(path):
            a,b,q=state;d=source%2
            state=(b-a-d*q,a+d*q,q)
            check(gcd(gcd(abs(state[0]),abs(state[1])),q)==1,"denominator changed")
            check(beta(state[:2],q)[0]==(a,b),"forward/backward golden mismatch")
            memo[source]=(state,root)
        counts[root]+=1
    print("Root/period/odd/halvings/theta(a,b,q):",cycle_table)
    print("Integer census -10000..10000 root counts:",sorted(counts.items()))
    # Check the parallel carry polynomial and signed reverse decoding.
    from itertools import product
    word_tests=0
    for length in range(1,5):
        for word in product(range(1,5),repeat=length):
            A,p,B=summary(word)
            partial=0;terms=[]
            for i,a in enumerate(word):
                terms.append((p-1-i,partial));partial+=a
            check(B==sum(3**i*2**j for i,j in terms))
            residue=((2**A-B)*pow(3**p,-1,2**(A+1)))%2**(A+1)
            for n0 in (residue,residue+2**(A+1),residue-2**(A+1)):
                n=n0
                for a in word:
                    n,actual=odd_step(n);check(actual==a)
                check(2**A*n==3**p*n0+B)
            word_tests+=1
    print("Ordered carry polynomial and exact signed cylinders:",word_tests,"words, 3 representatives each")
    shadow_tests=0
    for root,word in [(-1,(1,)),(-5,(1,2)),(-17,(1,1,1,2,1,1,4))]:
        A,p,B=summary(word);check((2**A-3**p)*root==B)
        for k in range(1,5):
            for b in (1,3,-1,-3):
                n0=root+2**(A*k+1)*b;n=n0
                for _ in range(k):
                    for a in word:
                        n,actual=odd_step(n);check(actual==a)
                check(n==root+2*3**(p*k)*b)
                shadow_tests+=1
    print("Three rational-anchor register families:",shadow_tests,"exact repeated-block controls")
    # All-equal interchange ray; equality with the same-atom pair occurs only at t=1.
    matches=[]
    for t in range(1,501):
        H,L=8*t*t+4*t,4*t*t+4*t
        if H%6==0:
            N=H//6
            if (odd_step(4*N-1)[0],odd_step(4*N+1)[0])==(H-1,L-1):
                matches.append(t)
    check(matches==[1]);print("Uniform interchange-ray / same-atom match through t=500:",matches)
    # The inherited square-window prime characterization.
    exceptional=[]
    for n in range(2,10001):
        if all(n%d for d in range(2,isqrt(n)+1)):
            j=1
            while (2*j+1)**2<=n:j+=1
            if 2*n<(2*j+1)**2:exceptional.append(n)
    check(exceptional==[2,3,11]);print("Odd-square-window prime exceptions through 10000:",exceptional)
    # Clocks alone forget the ordered carry. The exact golden coordinate with
    # its boundary convention determines the word; its denominator does not.
    check(summary((1,2))==(3,2,5) and summary((2,1))==(3,2,7))
    print("Carry-order hostile: clocks (A,p)=(3,2), carries 5 and 7")
    # Build concrete guarded ribbons; all fields have finite independent checks.
    from json import dumps
    examples=[]
    for n0 in (7,9,27,-3,-9,-27,-7,-25,151,64,54,-18):
        odd_source=n0;sleeve=0
        while odd_source%2==0:
            odd_source//=2;sleeve+=1
        n=odd_source;word=[];seen=[]
        while n not in seen:
            seen.append(n);n,a=odd_step(n);word.append(a)
        split=seen.index(n);prefix=word[:split];cycleword=word[split:];endpoint=n
        A,p,B=summary(prefix);CA,Cp,CB=summary(cycleword)
        check((2**CA-3**Cp)*endpoint==CB)
        check(2**A*endpoint- B == 3**p*odd_source)
        backward=endpoint
        for a in reversed(prefix):
            num=2**a*backward-1
            check(num%3==0 and (num//3)%2==1)
            backward=num//3
        check(2**sleeve*backward==n0)
        # Raw parity polynomial, and its exact eventual-period seal.
        prefixbits="0"*sleeve+"".join("1"+"0"*a for a in prefix)
        cyclebits="".join("1"+"0"*a for a in cycleword)
        polynomial=Counter()
        for j,d in enumerate(prefixbits):
            if d=="1":polynomial[j]+=1;polynomial[j+len(cyclebits)]-=1
        for j,d in enumerate(cyclebits):
            if d=="1":polynomial[j+len(prefixbits)]+=1
        binary_coordinate=Fraction(sum(c*2**j for j,c in polynomial.items()),1-2**len(cyclebits))
        actual_bits=0;current=n0;precision=64
        for j in range(precision):
            actual_bits+=(current%2)*2**j
            current=current//2 if current%2==0 else 3*current+1
        modulus=2**precision
        check(binary_coordinate.numerator*pow(binary_coordinate.denominator,-1,modulus)%modulus==actual_bits,
              "sealed polynomial/binary parity mismatch")
        numerator=(Fraction(0),Fraction(0))
        for j,coefficient in polynomial.items():
            power=ppow((-1,1),j)
            numerator=(numerator[0]+coefficient*power[0],numerator[1]+coefficient*power[1])
        power=ppow((-1,1),len(cyclebits))
        theta=pmul(pmul((-1,1),numerator),pinv((1-power[0],-power[1])))
        a,b,q=memo[n0][0]
        check(theta==(Fraction(a,q),Fraction(b,q)),"sealed polynomial/golden mismatch")
        # The incoming two-register Fibonacci reader is a second, finite path
        # to the golden prefix. Its input here is the TIME word, not digits(n0).
        X=Y=0
        for digit in prefixbits:
            d=int(digit);X,Y=Y+d,X+Y+2*d
        golden_prefix=(2*Y-3*X,2*X-Y)
        direct_prefix=(0,0)
        for digit in prefixbits:
            u,v=pmul((0,1),direct_prefix)
            direct_prefix=(u+int(digit),v)
        check(golden_prefix==direct_prefix,"Fibonacci reader/golden prefix mismatch")
        ta,tb,tq=memo[endpoint][0]
        advanced=pmul(ppow((0,1),len(prefixbits)),theta)
        check(advanced==(golden_prefix[0]+Fraction(ta,tq),golden_prefix[1]+Fraction(tb,tq)),
              "golden prefix/tail bridge mismatch")
        examples.append(dict(integer=n0,binary_sleeve=sleeve,prefix=prefix,endpoint=endpoint,cycle_word=cycleword,
                             affine=dict(A=A,p=p,B=B),
                             exact_residue=((2**A-B)*pow(3**p,-1,2**(A+1)))%2**(A+1),
                             modulus=2**(A+1),golden=dict(a=a,b=b,q=q),
                             parity_coordinate=dict(numerator=binary_coordinate.numerator,
                                                    denominator=binary_coordinate.denominator),
                             time_word_fibonacci_reader=dict(X=X,Y=Y,golden_a=golden_prefix[0],
                                                             golden_b=golden_prefix[1]),
                             seal_polynomial={str(j):v for j,v in sorted(polynomial.items()) if v}))
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument("--write-certificates")
    args=parser.parse_args()
    if args.write_certificates:
        from pathlib import Path
        Path(args.write_certificates).write_text(dumps(examples,indent=2)+"\n")
    print("Verified ribbon examples:",[(e["integer"],e["endpoint"],e["golden"]["q"]) for e in examples])
    print("Fibonacci digit reader / golden prefix-tail bridge:",len(examples),"exact controls")
    print("ALL CHECKS PASSED")


if __name__=="__main__":main()
