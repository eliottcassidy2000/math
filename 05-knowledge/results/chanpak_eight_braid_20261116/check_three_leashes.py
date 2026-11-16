"""Exact braid matrices and literal span-3 modular false negative."""
from pathlib import Path
import sys
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

HERE=Path(__file__).resolve().parent
for folder in (HERE,HERE.parent/'target_ribbon_root',
               HERE.parent/'chanpak_target_ribbon_audit_20261116'):
    if (folder/'verify_rank_alias.py').exists():
        sys.path.insert(0,str(folder))
        break
from verify_rank_alias import O,B,h,v0,product,rank,ow,rw
def need(cond):
    if not cond:
        raise RuntimeError('certificate failure')
def main():
    U=Matrix([[1,1],[0,1]])
    V=Matrix([[1,0],[-1,1]])
    need(U*V*U==V*U*V)
    C=U**8*V**8*U**(-8)*V**(-8)
    need(C==Matrix([[4033,512],[512,65]]) and C!=Matrix.eye(2)
         and all(x%8==0 for x in C-Matrix.eye(2)))
    u=(0,2,7,0,2)
    v=(0,1,7,0,2)
    row_and_cols=[([0,2,3],[0,2,5],36),([0,2,5],[0,2,5],32),
                  ([0,2,5],[0,2,5],-108),([0,2],[0,2],72)]
    for w,hh,expected,minor in zip((u,v,u+(1,),v+(1,)),(27,17,51,33),
                                    ((1,2,2),(1,4,4),(1,2,6),(1,4)),
                                    row_and_cols):
        P=product(w)
        K=O*P*B
        need(int((Matrix([list(h)])*P*Matrix(v0))[0])==hh)
        i,j,d=minor
        need(K.extract(i,j).det()==d)
        rr=rank(K)
        dd=smith_normal_form(K,domain=ZZ)
        need(tuple(abs(dd[i,i]) for i in range(rr))==expected)
        if w==v:
            need(K[0,0]==17)
            for k in range(1,K.rows):
                m=int(K[k,0])%4
                need(all((K[k,l]-m*K[0,l])%4==0 for l in range(K.cols)))
        # rw excludes the mandatory first vertex, which is included in p.
        ss=sorted(map(len,ow),reverse=True)[:rr]
        pp=sorted((len(ww)+1 for ww in rw),reverse=True)[:rr]
        print(w,hh,expected,'worst minor exponent',sum(ss)+sum(pp))
    print('Longest suffix and prefix length (including initial vertex):',
        max(map(len,ow)),1+max(map(len,rw)))
    print('PASS: braid commutator and literal mod-8 false negative.')
if __name__=='__main__':
    main()
