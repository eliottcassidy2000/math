"""Independent 26 exact replacement identities and terminal 2-saturation certificates."""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib

def need(b,m):
    if not b:raise ArithmeticError(m)
def rank(rows):
    basis={}
    for x in rows:
        while x:
            p=x.bit_length()-1
            if p in basis:x^=basis[p]
            else:basis[p]=x;break
    return len(basis)
def plus(P,Q,A,B):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and y+v==0:return None
    m=(3*x*x+2*A*x+B)/(2*y) if x==u else (v-y)/(u-x)
    z=m*m-A-x-u
    return z,-y+m*(x-z)
def determinant(M):
    a=[row[:] for row in M];value=F(1)
    for j in range(len(a)):
        k=next((k for k in range(j,len(a)) if a[k][j]),None)
        need(k is not None,'invertible transport')
        if k!=j:a[j],a[k]=a[k],a[j];value=-value
        pivot=a[j][j];value*=pivot
        for i in range(j+1,len(a)):
            factor=a[i][j]/pivot
            for k in range(j,len(a)):a[i][k]-=factor*a[j][k]
    return value

def main():
    root=Path(__file__).resolve().parents[1]/'halving';source=root/'halving_certificate.json';data=json.loads(source.read_text())
    result=[]
    for item in data['curves']:
        curve=item['curve'];ident=curve['id'];z=json.loads((root/f'phase2/extension_{ident}.json').read_text())
        need(hashlib.sha256(source.read_bytes()).hexdigest()==z['source_sha256'],'source hash')
        a1,a2,a3,a4,a6=map(F,curve['ainvs']);A=a2+a1*a1/4;B=a4+a1*a3/2;C=a6+a3*a3/4
        disc=16*(A*A*B*B-4*B**3-4*A**3*C-27*C*C+18*A*B*C)
        points=[(F(x),F(y)+(a1*F(x)+a3)/2) for x,y in curve['points']]
        def on(P):return P is None or P[1]**2==P[0]**3+A*P[0]**2+B*P[0]+C
        need(all(on(P) for P in points),'original points')
        specials={'O':0,'T':0}
        def rowscheck(rows):
            masks=[]
            for row in rows:
                p,r=row['p'],row['root'];need(p>2 and all(p%d for d in range(2,int(p**.5)+1)),'prime')
                red=lambda q:q.numerator*pow(q.denominator,-1,p)%p
                need(red(disc)!=0,'good prime');need((r**3+red(A)*r*r+red(B)*r+red(C))%p==0,'cubic root')
                mask=0
                for k,(x,y) in enumerate(points):
                    if x.denominator%p==0:
                        bit=0;specials['O']+=1
                    else:
                        v=(red(x)-r)%p
                        if v==0:
                            need(red(y)==0,'reduction to named2torsion')
                            v=(3*r*r+2*red(A)*r+red(B))%p;specials['T']+=1
                        need(v!=0,'nonzero character argument')
                        legendre=pow(v,(p-1)//2,p);need(legendre in [1,p-1],'Euler sign');bit=int(legendre==p-1)
                    mask|=bit<<k
                need(mask==row['bits'],'recomputed character row');masks.append(mask)
            return rank(masks)
        M=[[F(i==j) for j in range(30)] for i in range(30)]
        for number,step in enumerate(z['steps'],1):
            need(rowscheck(step['rows_before'])==step['rank_before']==16+number,'successive binary rank')
            S=step['support'];i=step['replace_index'];need(i in S and len(set(S))==len(S),'legal support')
            terms=[points[j] for j in S]
            while len(terms)>1:terms=[plus(terms[k],terms[k+1],A,B) if k+1<len(terms) else terms[k] for k in range(0,len(terms),2)]
            total=terms[0];half=tuple(map(F,step['half_completed_square']))
            need(on(half) and plus(half,half,A,B)==total==tuple(map(F,step['sum_completed_square'])),'exact half identity')
            new=[sum(M[j][k] for j in S)/2 for k in range(30)]
            need(new==list(map(F,step['new_basis_row'])),'exact source coefficient transport')
            M[i]=new;points[i]=half
        need(points==[tuple(map(F,P)) for P in z['final_points_completed_square']],'named final points')
        need(M==[list(map(F,row)) for row in z['final_basis_matrix']],'named final matrix')
        need(rowscheck(z['final_character_rows'])==30,'full final rank')
        p=z['two_torsion_exclusion_prime'];red=lambda q:q.numerator*pow(q.denominator,-1,p)%p
        need(p>2 and all(p%d for d in range(2,int(p**.5)+1)) and red(disc)!=0,'good torsion exclusion prime')
        need(all((x**3+red(A)*x*x+red(B)*x+red(C))%p for x in range(p)),'no rational2torsion')
        need(len(z['steps'])==13 and determinant(M)==F(1,8192),'index8192')
        need(all(q.denominator in [1,2] for row in M for q in row),'elementary quotient exponent2')
        transport=json.loads((root/f'phase2/inverse_transport_{ident}.json').read_text())
        N=[list(map(F,row)) for row in transport['original_basis_in_final_basis']]
        need(all(q.denominator==1 for row in N for q in row),'integral inverse transport')
        need(all(sum(N[i][k]*M[k][j] for k in range(30))==F(i==j) for i in range(30) for j in range(30)),'inverse exact')
        result.append({'curve':ident,'verified_replacements':13,'index':8192,'final_binary_rank':30,
                       'quotient':'C2^13','two_saturated':True,'special_character_evaluations':specials})
    print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
