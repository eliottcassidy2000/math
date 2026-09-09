"""Exact checksum/projective-tournament transport audit, standard library only.

No trigonometric approximation: determinant signs use exact cyclic residues.
No repository producers are imported. Exceptions remain active with python -O.
"""
from collections import Counter
from itertools import product
from math import gcd
from pathlib import Path

GATES = 0
LINES = []


def check(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def emit(line):
    LINES.append(line)


def checksum(word):
    return sum((i + 1) * b for i, b in enumerate(word)) % len(word)


def rotate(word, r=1):
    r %= len(word)
    return word[-r:] + word[:-r] if r else word


def determinant_sign(a, b, q):
    delta = (b-a) % (2*q)
    return 0 if delta in (0, q) else (1 if delta < q else -1)


def section_matrix(sigma):
    q = len(sigma)
    return tuple(tuple(0 if x == y else sigma[x]*sigma[y]*(1 if y>x else -1)
                       for y in range(q)) for x in range(q))


def signed_step(sigma, k=1):
    q = len(sigma)
    result = [0]*q
    for x in range(q):
        wrap, y = divmod(x+k, q)
        result[y] = sigma[x] * (-1 if wrap%2 else 1)
    return tuple(result)


def canonical(sigma):
    return tuple(sigma[0]*v for v in sigma)


def iterate(sigma, r, k=1):
    for _ in range(r):
        sigma = signed_step(sigma, k)
    return sigma


def strongly_connected(matrix):
    q=len(matrix)
    for start in range(q):
        seen={start}
        stack=[start]
        while stack:
            x=stack.pop()
            for y in range(q):
                if matrix[x][y] > 0 and y not in seen:
                    seen.add(y)
                    stack.append(y)
        if len(seen) != q:
            return False
    return True


def det(a,b):
    return a[0]*b[1]-a[1]*b[0]


def pfaffian4(h):
    return h[0][1]*h[2][3]-h[0][2]*h[1][3]+h[0][3]*h[1][2]


def main():
    words=0
    valid=0
    for m in range(2,13,2):
        local=0
        for word in product((0,1),repeat=m):
            j=sum(word)
            if j in (0,m):
                continue
            words+=1
            s=checksum(word)
            g=gcd(j,m)
            check(checksum(rotate(word)) == (s+j)%m, "literal cyclic response")
            orbit={checksum(rotate(word,r)) for r in range(m)}
            check(len(orbit) == m//g, "response orbit size")
            if (m//g)%2 == 0:
                q=m//(2*g)
                k=j//g
                check(gcd(k,2*q)==1 and k%2==1,"response normalization")
                s0=s%g
                r=(s-s0)//g
                check((r<q)==(s<m//2),"half-open decoder")
                moved=checksum(rotate(word,q))
                check(moved==(s+m//2)%m,"q rotations hit antipode")
                check((moved<m//2)!=(s<m//2),"actual output reverses")
                check(((moved-s0)//g)%q==r%q,"quotient address returns")
                actual_period=next(d for d in range(1,m+1) if rotate(word,d)==word)
                check(actual_period%(2*q)==0,"physical orbit retains response order")
                valid+=1
                local+=1
            else:
                check((s+m//2)%m not in orbit,"hostile missing antipode")
        emit(f"TAILS m={m}: antipode-compatible words={local}")
    emit(f"LITERAL UNIVERSE: {words} nonconstant words at every even m=2..12; {valid} admit the response antipode")

    total_sections=0
    for q in list(range(1,13))+[16]:
        states=[(1,)+tail for tail in product((-1,1),repeat=q-1)]
        state_set=set(states)
        total_sections+=len(states)
        step={s:canonical(signed_step(s)) for s in states}
        check(set(step.values())==state_set,"section transition permutation")
        histogram=Counter()
        unseen=set(states)
        while unseen:
            first=next(iter(unseen))
            cycle=[]
            x=first
            while x not in cycle:
                cycle.append(x)
                x=step[x]
            check(x==first,"honest section orbit")
            unseen.difference_update(cycle)
            histogram[len(cycle)]+=1
        for r in range(1,q+1):
            if q%r:
                continue
            fixed=sum(canonical(iterate(s,r))==s for s in states)
            expected=2**(r-1) if (q//r)%2 else 0
            check(fixed==expected,"fixed-section formula")
        two_part=q & -q
        check(min(histogram)==two_part,"sharp minimum section period")
        if q&(q-1)==0:
            check(dict(histogram)=={q:2**(q-1)//q},"every dyadic section has full period")
        if q%2:
            alternating=tuple((-1)**x for x in range(q))
            check(canonical(signed_step(alternating))==alternating,"odd carousel fixed section")
            S=section_matrix(alternating)
            check(all(sum(v>0 for v in row)==(q-1)//2 for row in S),"odd regular tournament")
        # Sample every lawful response step for small q, with every section.
        if q<=8:
            for k in range(1,2*q):
                if gcd(k,2*q)>1:
                    continue
                for sigma in states:
                    nxt=signed_step(sigma,k)
                    S=section_matrix(sigma)
                    Sn=section_matrix(nxt)
                    check(iterate(sigma,q,k)==tuple(-v for v in sigma),"negative monodromy")
                    for x in range(q):
                        for y in range(q):
                            check(Sn[(x+k)%q][(y+k)%q]==S[x][y],"determinant tournament transport")
                            rx=x+(q if sigma[x]<0 else 0)
                            ry=y+(q if sigma[y]<0 else 0)
                            check(S[x][y]==determinant_sign(rx,ry,q),"exact angular sign")
        emit(f"SECTIONS q={q}: count={len(states)}, period:orbit_count={dict(sorted(histogram.items()))}")
    emit(f"SECTION UNIVERSE: {total_sections} labelled sections modulo global sign")

    # Strong noncosmetic hostile: the full determinant section returns while
    # the actual weight-one checksum verdict flips under four rotations.
    sigma=(1,-1,1,-1)
    S=section_matrix(sigma)
    check(strongly_connected(S),"strong four-vertex control")
    check(iterate(sigma,4)==tuple(-v for v in sigma),"explicit negative return")
    check(section_matrix(iterate(sigma,4))==S,"complete edge field sign returns")
    word=(0,0,0,0,0,0,0,1)
    other=rotate(word,4)
    check(checksum(word)==0 and checksum(other)==4,"literal H/T return witness")
    emit(f"STRONG HOSTILE m=8,j=1,q=4: {''.join(map(str,word))} -> {''.join(map(str,other))}; checksum 0 heads -> 4 tails; same quotient address and section tournament")
    # A common non-diagonal frame change conjugates the rotation. At q=2,
    # A=[[1,1],[0,1]] and A*Rot(pi/2)*A^-1=[[1,-2],[1,-1]].
    frame=((1,0),(1,1))
    conjugated=lambda v:(v[0]-2*v[1],v[0]-v[1])
    check(conjugated(frame[0])==frame[1],"sheared frame first transport")
    check(conjugated(frame[1])==tuple(-v for v in frame[0]),"sheared frame wrap transport")
    for v in frame:
        check(conjugated(conjugated(v))==tuple(-z for z in v),"nondiagonal negative return")
    check(det(*frame)==det(*(conjugated(v) for v in frame)),"nondiagonal determinant covariance")
    emit("GAUGE q=2: sheared quarter rotation [[1,-2],[1,-1]] has square -I and preserves the full determinant")
    # One root sign reconstructs every orientation section from its tournament.
    for q in range(1,9):
        for sigma in product((-1,1),repeat=q):
            S=section_matrix(sigma)
            rebuilt=(sigma[0],)+tuple(sigma[0]*S[0][x] for x in range(1,q))
            check(rebuilt==sigma,"one-bit sufficient decoder")

    # Sharp four-vertex amplitude obstruction, all 64 tournament sign fields.
    values=Counter()
    pairs=[(i,j) for i in range(4) for j in range(i+1,4)]
    for signs in product((-1,1),repeat=6):
        H=[[0]*4 for _ in range(4)]
        for (i,j),s in zip(pairs,signs):
            H[i][j]=s
            H[j][i]=-s
        p=pfaffian4(H)
        check(p%2==1,"odd sign Pfaffian")
        values[p]+=1
    vectors=((1,0),(1,1),(0,1),(-1,1))
    H=[[det(a,b) for b in vectors] for a in vectors]
    check(pfaffian4(H)==0,"rank-two amplitude Plucker identity")
    check(all(H[i][j]>0 for i,j in pairs),"tie-free integer amplitude witness")
    emit(f"AMPLITUDE HOSTILE: six determinants={[H[i][j] for i,j in pairs]}, exact Pfaffian=0; tournament sign Pfaffians={dict(sorted(values.items()))}")
    emit("SCOPE: exact checksum response/section transport; one output sheet and edge amplitudes remain necessary; no improved extractor deadline or LRC closure")
    emit(f"PASS {GATES} always-active exact gates")


if __name__=='__main__':
    main()
    transcript='\n'.join(LINES)+'\n'
    print(transcript,end='')
    out=Path(__file__).resolve().parents[1]/'05-knowledge/results/tournament_orientation_transport_20260908.out'
    out.write_text(transcript,encoding='utf-8',newline='\n')
