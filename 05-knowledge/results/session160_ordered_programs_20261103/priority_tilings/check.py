from fractions import Fraction
from math import lcm
from itertools import product, permutations
from functools import reduce
P=[Fraction(17,143),Fraction(133,85),Fraction(17,19),Fraction(23,17),Fraction(1015,69),Fraction(23,29),Fraction(31,23),Fraction(38073,62),Fraction(31,37),Fraction(41,31),Fraction(86,287),Fraction(41,43),Fraction(13,41)]
D=[q.denominator for q in P]; m=len(P)
def feasible(bits):
    L=lcm(*(D[i] for i,x in enumerate(bits) if x))
    return all(x or L % D[i] for i,x in enumerate(bits))
A=list(product((0,1),repeat=m)); F=[x for x in A if feasible(x)]; both=sum(feasible(tuple(1-b for b in x)) for x in F)
if len(F)!=2592 or both!=512: raise AssertionError((len(F),both))
print('denominators',D,'patterns feasible',len(F),'of',len(A),'native-complement-preserved ordered patterns',both)
for x in F:
    if not feasible(tuple(1-b for b in x)):
        print('native complement exits image:', ''.join(map(str,x)), 'witness L=',lcm(*(D[i] for i,b in enumerate(x) if b)));break
else: raise AssertionError('expected a native complement exit')
for n in range(1,10001):
    e=tuple(int(n%d==0) for d in D)
    # independently execute two-state transfer
    state=0; marker=None
    for i,b in enumerate(e):
        if state==0 and b: state=1; marker=i
    expected=next((i for i,q in enumerate(P) if (n*q).denominator==1),None)
    if marker != expected: raise AssertionError((marker, expected))
print('checked 10000 integer scans')
# n=4, three chords (0,2),(0,3),(1,3)
ch=[(0,2),(0,3),(1,3)]
def adjacency(bits):
    B=dict(zip(ch,bits))
    return {(i,j) if (j==i+1 or B[i,j]) else (j,i) for i in range(4) for j in range(i+1,4)}
def iso(a,b):
    A=adjacency(a); B=adjacency(b)
    return next((p for p in permutations(range(4)) if {(p[i],p[j]) for i,j in A}==B),None)
for a in product((0,1),repeat=3):
 for b in product((0,1),repeat=3):
    fa=next((i for i,x in enumerate(a) if x),None); fb=next((i for i,x in enumerate(b) if x),None)
    p=iso(a,b)
    if p and fa!=fb:
        print('unmarked quotient collapses different scan outputs:',a,b,'map',p,'first',fa,fb);quit()

raise AssertionError('expected an unmarked collision')
