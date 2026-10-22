from itertools import product
from math import factorial
# Independent coefficient of q0*q1*q2*q3*x01^2*x02*x23 in the
# connected signed forest gas, using multivariate clone multiplicities.
# Edges: double 01 loop plus 02,23 tail. Odd runs cancel.
edges={(0,1):0,(0,2):1,(2,3):2}
target=(2,1,1)
def term(word,cuts):
    sig=[0]*3;sgn=1
    for start,end in cuts:
        for p in range(start,end-1):
            a,b=word[p:p+2]
            if a==b:return None
            edge=tuple(sorted((a,b)))
            if edge not in edges:return None
            sig[edges[edge]]+=1
            sgn*=1 if a<b else -1
    if tuple(sig)!=target:return None
    count=[word.count(i) for i in range(4)]
    if not all(count):return None
    mob=1
    for c in count:mob*=(-1)**(c-1)*factorial(c-1)
    return sgn*mob
p4=sum(term(w,[(0,5)]) or 0 for w in product(range(4),repeat=5))
p22=sum(term(w,[(0,3),(3,6)]) or 0 for w in product(range(4),repeat=6))
assert p22%2==0
print({'single_4_run':p4,'two_2_runs':p22//2,'connected_lollipop_coefficient':p4+p22//2,'closed_formula_expected':-4})
assert p4+p22//2==-4
