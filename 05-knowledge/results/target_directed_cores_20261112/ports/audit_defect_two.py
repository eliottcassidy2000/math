"""Independent small-core, literal-row, and zero-prefix audit of the defect-two source."""
from pathlib import Path
import json
import module_lift as m
import intertwiner_probe as i

HERE=Path(__file__).resolve().parent


def word(x):
    return tuple(map(int,bin(x)[2:][::-1])) if x else ()


def construct(core,J,x):
    w=word(x);n=len(core);p=n+len(w);r=p+1
    T=[0]*(r+1)
    for u,row in enumerate(core):
        T[u]=row|sum(1<<v for v in range(n,p))
        if u in J:T[p]|=1<<u
        else:T[u]|=1<<p
    for j,bit in enumerate(w):
        v=n+j
        T[v]|=sum(1<<u for u in range(v+1,p))
        if bit:T[v]|=1<<p
        else:T[p]|=1<<v
    T[p]|=1<<r;T[r]=((1<<r)-1)^(1<<p)
    assert all(((T[a]>>b)&1)+((T[b]>>a)&1)==1 for a in range(len(T)) for b in range(a+1,len(T)))
    return tuple(T),p,r


def main():
    source=json.loads((HERE.parent/'constructors/defect_two_inverse_receipt.json').read_text())
    core=(6,4,8,3);J=(1,2);x=1133276
    H,ends=m.tables(core);_,starts=m.tables(m.opposite(core));full=(1<<len(core))-1
    E=[sum(row[j] for j in J) for row in ends]
    S=[sum(row[j] for j in J) for row in starts]
    moments=dict(h=H[-1],p2=sum(H[a]*H[full^a] for a in range(full+1)),e=E[-1],
                 F=sum(E[a]*H[full^a] for a in range(full+1)),
                 G=sum(S[a]*H[full^a] for a in range(full+1)),budget=S[-1],
                 C=sum(E[a]*S[full^a] for a in range(full+1)))
    assert moments==dict(h=5,p2=32,e=3,F=16,G=11,budget=2,C=4)
    assert moments==source['parameters']
    T,p,r=construct(core,J,x)
    assert T==tuple(source['rows']) and (p,r)==(source['root_parent'],source['root'])
    g=i.coeff(word(x),1)[2]
    z=(5*g+43*x+3*(x&1)+18,5*g+11*x+3*(x&1)+2,
       5*g+11*x+3*(x&1)+4,5*g+43*x+3*(x&1)+20)
    assert z==(53084231,16819383,16819385,53084233)
    assert z==tuple(source['root_state_H_s_q_F'])
    small=[]
    for y in range(16):
        U,a,b=construct(core,J,y);gy=i.coeff(word(y),1)[2]
        hy=5*gy+43*y+3*(y&1)+18;qy=5*gy+11*y+3*(y&1)+4
        actual=m.profile(U,b)
        assert actual==(hy,qy-2,qy,hy+2)
        V=m.birth(U,b);active=len(U);V=m.birth(V,active)
        actual_nr=m.profile(V,active)
        assert actual_nr==(3*hy+4,hy,hy,3*hy+4)
        small.append(dict(outer_word=y,state=actual,NR_state=actual_nr))
    extension=[]
    for ell in range(6):
        U,a,b=construct(core,J,x<<ell)
        old_names=list(range(4))+list(range(4+ell,len(U)))
        assert m.induced(U,old_names)==T
        gy=i.coeff(word(x<<ell),1)[2]
        assert gy==g+((1<<ell)-1)*x
        hy=5*gy+43*(x<<ell)+18;qy=5*gy+11*(x<<ell)+4
        assert hy==54397248*(1<<ell)-1313017
        assert qy==18132416*(1<<ell)-1313031
        extension.append(dict(zero_prefix_length=ell,H=hy,q=qy,root=b,predecessor=a,
                              retained_old_vertex_names=old_names))
    adjacent_x=45157534;gx=i.coeff(word(adjacent_x),1)[2]
    A,ap,ar=construct((0,),(0,),adjacent_x)
    h=gx+3*adjacent_x+(adjacent_x&1);q=adjacent_x+gx+(adjacent_x&1)
    assert h==159252695
    result=dict(status='PASS',scope='Independent subset DP on the four-vertex core; all16 short-word full root and NR profiles; literal27-vertex row comparison; six retained-vertex prefix controls; direct word-automaton positive neighbor159252695. All-length formulas use the proved funnel law.',
                moments=moments,target_state=z,binary_word=word(x),g=g,short_controls=small,prefix_controls=extension,
                adjacent_ready_control=dict(target=h,word_integer=adjacent_x,g=gx,rows=A,root=ar,predecessor=ap,
                                            state=(h,q-1,q,h+1)))
    (HERE/'defect_two_audit_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],moments=moments,target_state=z,adjacent_ready_state=result['adjacent_ready_control']['state'])))


if __name__=='__main__':main()
