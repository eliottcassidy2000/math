"""Minimal ternary carry modules and the previous radical's literal 9-clock."""
from math import comb
import json
from support import HERE,PREV,need,save
def mul(a,b):
    r=len(a);out=[0]*r
    for i,x in enumerate(a):
        for j,y in enumerate(b[:r-i]):out[i+j]=(out[i+j]+x*y)%3
    return tuple(out)
def row(v,M):return tuple(sum(x*M[i][j] for i,x in enumerate(v))%3 for j in range(len(M[0])))
def substitution(r):
    t=(0,2,1)+((0,)*(r-3)) if r>=3 else (0,2)
    result=[];power=(1,)+(0,)*(r-1)
    for i in range(r):result.append(power);power=mul(power,t)
    return result
def main():
    checks=[]
    for k in range(1,5):
        N=3**k;r=3**(k-1)+1;D=substitution(r)
        orbit=[tuple(comb(j,i)%3 if i<=j else 0 for i in range(r)) for j in range(N)]
        u4=tuple(comb(4,i)%3 if i<=4 else 0 for i in range(r))
        need(len(set(orbit))==N,'faithful orbit of exact ternary clock length')
        need(all(mul(orbit[j],u4)==orbit[(j+4)%N] and row(orbit[j],D)==orbit[2*j%N] for j in range(N)),'literal translation/doubling intertwiner')
        need(pow(2,2*3**(k-1),N)==1 and (k==1 or pow(2,2*3**(k-2),N)!=1),'full doubling order')
        checks.append({'k':k,'cycle':N,'linear_dimension':r,'affine_dimension':r-1,'doubling_order':2*3**(k-1),'group_order':2*3**(2*k-1)})
    a=json.loads((PREV/'radical_clocks.json').read_text())
    from fractions import Fraction
    M=[[int((q:=Fraction(x)).numerator)*pow(q.denominator,-1,3)%3 for x in rr] for rr in a['nilpotent_matrix']]
    E=None
    for i in range(27):
        v=tuple(int(j==i) for j in range(27));chain=[v]
        for _ in range(4):chain.append(row(chain[-1],M))
        if any(chain[3]):E=chain[:4];need(not any(chain[4]),'four-step inherited radical chain');break
    need(E is not None,'primitive actual matrix chain modulo3')
    states=[]
    for label in range(9):
        c=[comb(label,j)%3 if j<=label else 0 for j in range(4)]
        states.append(tuple(sum(c[i]*E[i][j] for i in range(4))%3 for j in range(27)))
    need(len(set(states))==9,'nine distinct full observer vectors')
    U=[[int(i==j)+M[i][j] for j in range(27)] for i in range(27)]
    need(all(row(states[j],U)==states[(j+1)%9] for j in range(9)),'literal inherited unipotent clock')
    labels=[4,8,3,7,2,6,10,5,9,4]
    need(all((b-a)%9==4 for a,b in zip(labels,labels[1:])),'original displayed cycle')
    save('clock.json',{'status':'FINITE_EXACT_MINIMAL_CLOCK_MODELS','levels':checks,'four_state_substitution':substitution(4),
       'inherited_row_index':i,'chain_basis':E,'states_by_residue':states,'displayed_labels':labels,
       'scope':'formal radical/module intertwiner; no assertion these vectors or invertible clocks are legal physical source words'})
    print('PASS minimal linear dimensions2,4,10,28; inherited radical realizes the9-cycle',flush=True)
if __name__=='__main__':main()
