"""Exact selected-descent certificates for recreated dyadic precision."""
from fractions import Fraction


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def v2(n):
    return (n & -n).bit_length()-1


def U(n):
    z = 3*n+1
    return z >> v2(z)


def iterate(n, k):
    for _ in range(k):
        n = U(n)
    return n


def core_certificate(q):
    n = 3*q
    A = J = 0
    earlier_slopes = []
    while n >= 2*q:
        need(J < 1000, ('core cap', q))
        P = A
        z = 3*n+1
        A += v2(z)
        n = z >> v2(z)
        J += 1
        earlier_slopes.append(3**(J+1) < 2**(A+1))
    b = n
    clipped = max(P+1, A-(((2*q-1)//b).bit_length()-1))
    for R in range(P+1, A+1):
        D = 2**(R+1)-3**(J+1)
        B = 3*2**(A-R)*b-6*q+1
        if D > 0 and D > B:
            break
    else:
        raise RuntimeError(('missing threshold',q))
    need(R <= clipped <= A, ('threshold ordering',q))
    # This is a finite-bank observation, not used by the universal proof.
    need(not any(earlier_slopes[:-1]), ('earlier slope in bank',q))
    K = R+1
    residue = ((6*q-1)*pow(3,-1,1 << K)) % (1 << K)
    return dict(q=q,J=J,A=A,P=P,b=b,R=R,clip=clipped,K=K,residue=residue)


def main():
    rows = [core_certificate(q) for q in range(1,342,2)]
    need(max(r['J'] for r in rows)==40 and max(r['A'] for r in rows)==66, 'bank maxima')
    need(max(r['R'] for r in rows)==64, 'optimized maximum')
    full = []
    for q in range(1,342,2):
        c = 3*q
        A = J = 0
        while c != 1:
            need(J < 1000, ('full core cap',q))
            z = 3*c+1
            A += v2(z)
            c = z >> v2(z)
            J += 1
        full.append((q,J,A))
    need(max(r[1] for r in full)==56 and max(r[2] for r in full)==99, 'full route maxima')

    hard_cases = selected = collision_controls = clipped_controls = 0
    for row in rows:
        q,J,A,P,b = (row[k] for k in ('q','J','A','P','b'))
        for R in range(1,13):
            for v in range(1,64,2):
                if (2**(R+1)*v-1) % 3:
                    continue
                n = (2**(R+1)*v+6*q-1)//3
                need(U(n)==v*2**R+3*q, ('swap',q,R,v))
                hard_cases += 1
        for R in sorted({row['R'],row['R']+1,A,A+1,A+8,A+40}):
            for v in range(1,80,2):
                if (2**(R+1)*v-1) % 3:
                    continue
                n = (2**(R+1)*v+6*q-1)//3
                endpoint = iterate(n,J+1)
                if R > A:
                    expected = 3**J*v*2**(R-A)+b
                elif R == A:
                    raw = 3**J*v+b
                    expected = raw >> v2(raw)
                    collision_controls += 1
                else:
                    need(R>P, ('prior exactness',q,R))
                    expected = 3**J*v+2**(A-R)*b
                    clipped_controls += 1
                need(endpoint == expected and endpoint < n, ('selected descent',q,R,v))
                selected += 1

    exceptions = 0
    for row in rows:
        q,J,K,r = (row[k] for k in ('q','J','K','residue'))
        need(r%4==3, ('residue parity',q))
        for n in range(r,2*q+1,1 << K):
            need(n>1 and iterate(n,J+1)<n, ('small residue exception',q,n))
            exceptions += 1
    need(exceptions==1431,'exception count')

    minimal = []
    for row in sorted(rows,key=lambda r:(r['K'],r['q'])):
        if not any(row['residue']%(1 << old['K'])==old['residue'] for old in minimal):
            minimal.append(row)
    for i,a in enumerate(minimal):
        for b in minimal[i+1:]:
            K = min(a['K'],b['K'])
            need(a['residue']%(1 << K)!=b['residue']%(1 << K), 'cylinders disjoint')
    density = sum((Fraction(1,1 << r['K']) for r in minimal),Fraction())
    need(len(minimal)==65 and density==Fraction(6985206796614369409,36893488147419103232),'density certificate')

    # Independent direct-prefix coverage check, using only selected class data.
    coverage = 0
    for n in range(3,100000,2):
        matches = [r for r in minimal if n%(1 << r['K'])==r['residue']]
        need(len(matches)<=1, ('unique cylinder',n))
        if matches:
            need(iterate(n,matches[0]['J']+1)<n, ('direct cylinder',n))
            coverage += 1

    for j in range(1,65):
        q = (64**j-1)//9
        need(3*(3*q)+1==2**(6*j), ('unbounded q core',j))
        for R in (1,2,3,4,6*j,6*j+1):
            for v in range(1,40,2):
                if (2**(R+1)*v-1)%3:
                    continue
                n=(2**(R+1)*v+6*q-1)//3
                endpoint=iterate(n,2)
                need((endpoint<n)==(R>=3), ('q family threshold',j,R,v))

    first_descent = {}
    for c in (2,4):
        for k in range(1,129):
            start = c*8**k-5
            need(v2(3*(start-2)+1)==2, ('fixed low precision',c,k))
            need(c*8**(k+1)-5==8*start+35, ('affine source recursion',c,k))
            n = start
            for j in range(k):
                need(n==c*9**j*8**(k-j)-5, ('pair exact',c,k,j))
                middle = U(n)
                need(middle==(3*c//2)*9**j*8**(k-j)-7 and middle>start, ('middle exact',c,k,j))
                n=U(middle)
                need(n>start, ('no first descent',c,k,j))
            need(n==c*9**k-5, ('pair final',c,k))
        table = []
        for k in range(1,9):
            start=c*8**k-5
            n=start
            for t in range(1,10001):
                n=U(n)
                if n<start:
                    break
            else:
                raise RuntimeError(('first descent cap',c,k))
            need(t>2*k, ('first descent lower bound',c,k,t))
            table.append((k,start,t,n))
        first_descent[c]=table
    for k in range(1,129):
        boundary=4*8**k-5
        for X,expected in ((boundary-1,k-1),(boundary,k),(4*8**(k+1)-6,k)):
            count=0
            while 4*8**(count+1)-5<=X:
                count+=1
            need(count==expected, ('27 family counting',k,X))
    need(iterate(27,3)==47>27, 'smallR short-core hostile')

    print('reset_20260926_swaplift: exact certificates')
    print('171 core coefficients q odd1..341; maxJ40,maxA66 atq9; max optimizedR64')
    print('full-to1 comparison: maxJ56,maxA99 atq339')
    print('hard-swap identities:',hard_cases,'; selected source controls:',selected)
    print('final-collision controls:',collision_controls,'; final-swap clipping controls:',clipped_controls)
    print('small progression exceptions checked:',exceptions)
    print('minimal disjoint cylinders:',len(minimal),'; natural density:',density)
    print('density amongoddintegers:',2*density,'; direct covered sources below100000:',coverage)
    print('unbounded q=(64^j-1)/9 controls j1..64: precision3 suffices,1/2 fail two-step descent')
    print('fixed q1,R1 positive shadow controls c2/4,k1..128: all first2k iterates exceed the source')
    print('first descent table columns k,source,odd_steps,endpoint; cap10000:')
    print('c2:',first_descent[2])
    print('c4, family containing27:',first_descent[4])
    print('27-family counting boundary controls k1..128 pass; no monotonic first-descent inference')
    print('certificate table: q J A Aprev b clippedR optimizedR residue modulus_exponent')
    for r in rows:
        print(*(r[k] for k in ('q','J','A','P','b','clip','R','residue','K')))
    print('minimal-cylinder coefficient labels:',','.join(str(r['q']) for r in minimal))
    print('PASS: all checks active under -O')


if __name__=='__main__':
    main()
