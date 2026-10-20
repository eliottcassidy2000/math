#!/usr/bin/env python3
from itertools import combinations, permutations, product
from verify import tournament, h_dp, h_perm, h_word, strong, need


def invs(word,bit):
    return [(1<<i)|(1<<j) for i in range(len(word)) for j in range(i+1,len(word))
            if word[i]>>bit&1 and not(word[j]>>bit&1)]


def cross(word):
    a,b=invs(word,0),invs(word,1)
    return bool(a and b) and all(x&y for x in a for y in b)


def templates(a,b):
    return [(1,)*a+(2,0),(2,)*a+(1,0),
            (1,)*a+(2,)+(1,)*b,(2,)*a+(1,)+(2,)*b,
            (1,)*a+(3,)*b+(0,),(2,)*a+(3,)*b+(0,),
            (3,)+(0,)*a+(1,)*b,(3,)+(0,)*a+(2,)*b,
            (3,1)+(2,)*a,(3,2)+(1,)*a,(3,)+(0,)*a,(3,)*a+(0,)]


def cycle_counts(out,m):
    c=[0,0,0]
    for length in range(3,len(out)+1,2):
        for sub in combinations(range(len(out)),length):
            for rest in permutations(sub[1:]):
                p=(sub[0],)+rest
                if all(out[p[i]]>>p[(i+1)%length]&1 for i in range(length)):
                    x,y=m in sub,m+1 in sub
                    need(x or y,('cycle outside feedback',p))
                    c[2 if x and y else 0 if x else 1]+=1
    return tuple(c)


def main():
    bank={(3,1,0),(3,2,0)}
    for a in range(1,9):
        for b in range(1,9):bank.update(templates(a,b))
    counts=[]
    for m in range(2,9):
        good=set()
        for w in product(range(4),repeat=m):
            if w[0]==0 or w[-1]==3:continue
            if cross(w):good.add(w)
        claimed={w for w in bank if len(w)==m}
        need(good==claimed,('reduced grammar',m,good^claimed))
        counts.append((m,len(good)))
    print('exact reduced grammar controls by chain length:',counts)
    cases=0
    for a in range(1,5):
        families=[((1,)*a+(2,0),(3*2**(a-1)-1,1,3*2**(a-1)),6*2**a+1),
                  ((2,)*a+(1,0),(1,3*2**(a-1)-1,3*2**a-2),9*2**a-3)]
        for b in range(1,5):
            if a+b>5:continue
            families += [((2,)*a+(1,)+(2,)*b,
                          (2**(b-1),2**(a-1),3*2**(a+b-1)-2**(a-1)-2**(b-1)),3*2**(a+b)+1),
                         ((2,)*a+(3,)*b+(0,),
                          (2**(b-1),2**(a+b-1),2**(a+b)-2**(b-1)),3*2**(a+b)+1)]
        for w,expected,h in families:
            out=tournament(w,2,1)
            need(cross(w) and cycle_counts(out,len(w))==expected,('direct odd cycles',w))
            need(h_dp(out)==h==1+2*sum(expected),('path count',w))
            need(strong(out),('strong ray',w))
            if len(out)<=7:need(h_perm(out)==h,('permutation',w))
            cases+=1
    print('direct cycle-table/path checks:',cases)
    for w,h in [((3,1,0),13),((3,2,0),15)]:
        out=tournament(w,2,1);c=cycle_counts(out,len(w))
        need(h_perm(out)==h_dp(out)==h==1+2*sum(c),('exception',w))
        print('short exception:',w,'cycles',c,'H',h)
    out=(50,16,3,7,12,30)
    triangles=[set(s) for s in combinations(range(6),3)
               if all(sum(bool(out[v]>>u&1) for u in s if u!=v)==1 for v in s)]
    pair=next((sorted(a),sorted(b)) for a,b in combinations(triangles,2) if a.isdisjoint(b))
    need(h_perm(out)==h_dp(out)==29 and strong(out),'H29 witness')
    print('H29 witness:',out,'disjoint triangles:',pair,'all720 permutations: PASS')
    w=tuple(map(int,'1100100'));out=tournament(w)
    need(h_perm(out)==h_dp(out)==h_word(w)==87 and strong(out),'H87 hostile')
    print('H87 one-apex hostile:',out,'all40320 permutations: PASS')
    # Direct finite validation of root's reflection-filter counting proof.
    for L in range(5,13):
        X=1<<L;values={1};removed=0;fixed=0;words=0
        for m in range(2,L+2):
            for mid in product(range(2),repeat=m-2):
                w=(1,)+mid+(0,);rev=tuple(1-x for x in reversed(w));words+=1
                need(h_word(w)==h_word(rev),('involution',w))
                fixed+=w==rev
                if m==L+1 and w[:2]==(1,1) and w[-2:]==(0,0):
                    need(h_word(w)>X,('outer four inversions',w));removed+=1
                h=h_word(w)
                if h<=X:values.add(h)
        bound=7*X//16+2**((L+1)//2-1)
        need(words==X-1 and removed==X//8 and fixed==2**((L+1)//2)-1,'bank counts')
        need(len(values)<=bound,('density bound',L))
    print('dyadic reflection-filter controls L=5..12: PASS')
    print('ALL CHECKS PASSED')

if __name__=='__main__':main()
