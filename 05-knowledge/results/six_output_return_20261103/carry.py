"""Arithmetic branch words and the two distinct integer lattices behind six."""
import json,itertools
from pathlib import Path
from sympy import Matrix,ZZ
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).resolve().parent
def need(ok,message):
    if not ok:raise ValueError(message)
def run():
    for h in range(1,16,2):
      for c2 in range(1,9):
       for c in (2,4):
        h2=h+2*c2;v1=(9*h-3*h2)//2+c
        need((v1%2==1)==(c2%2==0),'clone-splice intermediate parity')
        if c2%2==0:
            v2=3*c2//2
            for m in range(1,8):need(v1+(2**(m+1)-2)*v2==3*(h+(2**m-2)*c2)+c,'independent corrected splice filter')
    levels=[]
    for k in range(10):
        values={};mod={}
        for word in itertools.product((2,4),repeat=k):
            c=0
            for digit in word:c=3*c+digit
            D=sum(((digit-2)//2)*3**(k-1-i) for i,digit in enumerate(word))
            need(c==3**k-1+2*D,'ordered ternary carry')
            heap=1;ternary='1'
            for digit in word:heap=2*heap+(digit-2)//2;ternary+=str(digit-2)
            need(int(ternary,3)==c+1 and heap.bit_length()-1==k,'carry length and binary-prefix address')
            need(c not in values,'distinct words have distinct integer carries');values[c]=word
            mod.setdefault(c%7,[]).append(word)
        levels.append(dict(length=k,integer_carries=len(values),mod7_carries=len(mod)))
    c242=3*(3*2+4)+2;c424=3*(3*4+2)+4
    need((c242,c424)==(32,46) and c242%7==c424%7,'first fixed-length phase collision')
    need(levels[0]['integer_carries']==1 and all(x['integer_carries']==x['mod7_carries'] for x in levels[:3]),'no earlier fixed-length collision')
    # Formal branch-character and carry group; physical full-word realization remains open.
    seen={(1,0)};todo=[(1,0)]
    for s,c in todo:
        for d in (2,4):
            z=(3*s%7,(3*c+d)%7)
            if z not in seen:seen.add(z);todo.append(z)
    need(len(seen)==42,'formal multiplier/carry group')
    squares={1,2,4}
    need({3*x%7 for x in squares}=={3,5,6},'ternary branch reverses the Paley phase relation')
    # Q_r = ordered nonempty r-path covers; F_r allows empty slots.
    G=Matrix([[1,0,0],[2,2,0],[0,3,3]])
    P=Matrix([[1,0,0],[2,1,0],[3,3,1]])
    need(P.det()==1 and P*G*P.inv()==Matrix.diag(1,2,3),'physical ordinal-padding diagonalization over Z')
    powers=[]
    for t in range(1,9):
        D=smith_normal_form(G**t,domain=ZZ);diagonal=[abs(int(D[i,i])) for i in range(3)]
        need(diagonal==[1,1,6**t],'one cyclic padding inverse defect')
        powers.append(dict(t=t,smith=diagonal))
    E=Matrix([[1,0,0],[1,2,0],[1,6,6]])
    D=smith_normal_form(E,domain=ZZ);evaluation_smith=[abs(int(D[i,i])) for i in range(3)]
    need(evaluation_smith==[1,2,6] and E.det()==12,'three-query coefficient lattice is a different object')
    out=dict(scope='arithmetic words and formal integer response modules; not universal physical branch realization',levels=levels,
             carry_code='C+1 has ternary digits 1 followed by 0/2; k is its digit length minus one. Binary heap address replaces those suffix digits by 0/1.',
             splice_parity_checks=128,splice_identity_checks=448,
             first_phase_collision=dict(words=[[2,4,2],[4,2,4]],carries=[32,46]),formal_group_size=42,
             padding_matrix=[list(map(int,G.row(i))) for i in range(3)],empty_slot_basis=[list(map(int,P.row(i))) for i in range(3)],
             padding_smith=powers,evaluation_smith=evaluation_smith)
    (HERE/'carry.json').write_text(json.dumps(out,indent=2)+'\n');print('carry and lattice checks passed',flush=True);return out
if __name__=='__main__':run()
