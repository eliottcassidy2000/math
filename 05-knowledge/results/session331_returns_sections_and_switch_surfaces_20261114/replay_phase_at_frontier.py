from replay_matching_sections import ctx_coeff, replace, coefficients, labelled
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from profile_tools import strong, hp_all

bank={
    'a':((14,28,8,16,5),(9,10,2)),
    'b':((14,12,24,16,3),(11,17,4)),
    'c':((6,28,8,17,5),(13,10,2)),
    'd':((14,4,24,18,3),(13,12,2)),
    'e':((22,12,24,1,10),(15,11,2)),
    'g':((6,12,24,17,3),(15,13,2)),
    '+':((22,12,24,17,2),(13,6,0)),
    '-':((6,28,24,17,1),(11,6,0)),
}
for R, row in bank.values():
    assert strong(R) and tuple(ctx_coeff(R,k) for k in (1,2,3))==row
    assert ctx_coeff(R,4)==0
assert {tuple(ctx_coeff(R,k) for k in (1,2,3))
        for R in labelled(5) if strong(R) and ctx_coeff(R,3)}=={
        row for R,row in bank.values() if row[2]}
def det(M):
    (a,b,c),(d,e,f),(g,h,i)=M
    return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
for e,names,dd in [(0,'bcd',60),(1,'ade',-16)]:
    rows=[bank[k][1] for k in names]
    assert det(rows)==dd and dd%11
    assert all((x-e*y)%3 for x,y,z in rows)
    for p in range(11):
        for q in range(11):
            assert any((x+y*p+z*q)%11 for x,y,z in rows)
for U,ex in [((2,12,9,17,7),'bcd'),((2,16,3,7,45,15),'ade')]:
    P=coefficients(U,3); h,p,t=P; q=p+2*h
    Y={}
    for k,(R,row) in bank.items():
        G=replace(R,4,U)
        assert strong(G)
        Y[k]=hp_all(G)[-1]
        assert Y[k]==sum(a*b for a,b in zip(row,P))
    assert (Y['-'],Y['+'])==(6*q-h,6*q+h) and 2*Y['g']-Y['b']==9*q+h
    assert any(Y[k]%3 and Y[k]%11 for k in ex)
    assert any(Y[k]%3 and Y[k]%11 for k in '+-')
    print(len(U),P,[(k,Y[k]) for k in ex+'+-'])
print('same-order 99-unit frontier: PASS')
