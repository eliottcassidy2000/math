"""Exact retained-predecessor reroot and all-clock primitive-exit certificate."""
from pathlib import Path
from math import gcd, isqrt
import json
import module_lift as m


def run():
    controls=[]; coeff=[]
    for i,recipe in enumerate(m.RECIPES):
        values=[]
        for name,M in [('K1',(0,)),('TT2',m.tt(2)),('TT3',m.tt(3)),('C3',(2,4,1))]:
            T,block,x,y=m.build(recipe,M)
            z=m.profile(T,x)
            HM,_=m.tables(M); full=(1<<len(M))-1
            h=HM[-1];c=sum(HM[a]*HM[full^a] for a in range(1,full))
            values.append((name,h,c,z))
            control=dict(chart=i+1,module=name,state_at_x=z)
            if i==1 and name=='K1':
                U=m.birth(T,x)
                U=m.birth(U,x);root=len(U)-1
                U=m.birth(U,root)
                actual=m.profile(U,root)
                expected=(z[0]+5*z[3]-3*z[2],z[0]+z[3]-z[2],z[0]+z[3]-z[2],z[0]+5*z[3]-3*z[2])
                assert actual==expected
                control.update(RNR_rows=U,RNR_root=root,RNR_state=actual)
            controls.append(control)
        a=values[0][3];b=tuple((v-u)//2 for u,v in zip(a,values[1][3]))
        assert all(tuple(aa*h+bb*c for aa,bb in zip(a,b))==z for _,h,c,z in values)
        coeff.append(dict(chart=i+1,h_coeff=a,C2_coeff=b))
    assert coeff[1]['h_coeff']==(2759,15,15,5055)
    assert coeff[1]['C2_coeff']==(1386,0,0,2592)
    assert 7799*14346-3978*27989==108*5039
    A,B,C,D=2689128,-5495,9697896,-24605
    determinant=A*D-B*C
    factors={2:4,3:3,5:1,7:1,13:2,5039:1}
    product=1
    for p,e in factors.items():
        assert all(p%d for d in range(2,isqrt(p)+1))
        product*=p**e
    assert abs(determinant)==product
    residue=(-B*pow(A,-1,5039))%5039
    assert residue==3365
    assert pow(2,2519,5039)==1
    assert pow(residue,2519,5039)==5038
    # Every excluded small prime either never divides u or requires X=0.
    small=[]
    for p in (2,3,5,7,13):
        powers=set();x=1
        while x not in powers:
            powers.add(x);x=x*2%p
        assert all(gcd(A*x+B,C*x+D)%p for x in powers)
        small.append(dict(prime=p,powers_of_two=sorted(powers)))
    samples=[]
    for k in range(31):
        X=1<<k
        H=936936*X-455;F=1752192*X-4515;q=525
        u=A*X+B;b=C*X+D
        assert u==H+F-q and b==H+5*F-3*q
        assert gcd(u,b)==1 and (u%3,b%3)==(1,1)
        assert gcd(u,q,2*F-q)==1
        samples.append(dict(clock=k,input_state_at_x=(H,q,q,F),balanced=(b,u,u,b),existing_pair_content=gcd(u,q,2*F-q)))
    tt_samples=[]
    assert 3978*(-703)-(-157)*14346==-4*27*5039
    assert 157*pow(3978,-1,5039)%5039==2365
    assert pow(2365,2519,5039)==5038
    for k in range(1,32):
        X=1<<k;H=1386*X-13;F=2592*X-129;q=15
        u=3978*X-157;b=14346*X-703
        assert (H,q,q,F)==tuple(a+z*(X-2) for a,z in zip(coeff[1]['h_coeff'],coeff[1]['C2_coeff']))
        assert (u,b)==(H+F-q,H+5*F-3*q)
        assert gcd(u,b)==1 and (u%3,b%3)==(2,2)
        assert gcd(u,q,2*F-q)==1
        tt_samples.append(dict(clock=k,input_state_at_x=(H,q,q,F),balanced=(b,u,u,b),existing_pair_content=gcd(u,q,2*F-q)))
    growing=[]
    residue_powers={str(v):pow(v,2519,5039) for v in (2,35,676,2365)}
    assert residue_powers=={'2':1,'35':1,'676':1,'2365':5038}
    assert gcd(35,2*676)==1
    for a in range(9):
        for bcount in range(7):
            if a+bcount==0:continue
            h=35**bcount;P=676**bcount*2**a;c=P-2*h
            H=1386*P-13*h;F=2592*P-129*h;q=15*h
            u=3978*P-157*h;exit_b=14346*P-703*h
            assert c>=0 and gcd(h,c)==1 and h%3
            assert (u,exit_b)==(H+F-q,H+5*F-3*q)
            assert gcd(u,exit_b)==1 and (exit_b-2*u)%3==h%3
            assert pow(P*pow(h,-1,5039)%5039,2519,5039)==1
            growing.append(dict(a=a,b=bcount,source_order=7*bcount+a+11,
                                exit_order=7*bcount+a+14,
                                module_H=str(h),module_P2=str(P),module_C2=str(c),
                                source_state=list(map(str,(H,q,q,F))),balanced=list(map(str,(exit_b,u,u,exit_b)))))
    result=dict(status='PASS',scope='Three inherited quotients; complete rooted subset-DP at predecessor x for four small modules; one literal RNR replay; exact finite factor and subgroup certificates proving all clocks; initial clock values are controls only',
                coefficients=coeff,controls=controls,
                primitive_certificate=dict(u=(A,B),b=(C,D),resultant=determinant,absolute_resultant_factorization=factors,
                    small_prime_controls=small,remaining_prime=5039,required_residue=residue,exponent=2519,
                    power_of_two=pow(2,2519,5039),power_of_required_residue=pow(residue,2519,5039)),samples=samples,
                transitive_clock_certificate=dict(u=(3978,-157),b=(14346,-703),resultant=-544212,
                    absolute_resultant_factorization={2:2,3:3,5039:1},remaining_prime=5039,required_residue=2365,
                    exponent=2519,power_of_required_residue=pow(2365,2519,5039)),transitive_clock_samples=tt_samples,
                growing_module_certificate=dict(domain='a,b>=0; a+b>0',module='D^b ordinal TT_a',
                    module_base_rows=(8,5,9,18,71,31,47),base_H=35,base_P2=676,
                    prime=5039,exponent=2519,residue_powers=residue_powers,
                    coefficient_coordinates='h=35^b, P=676^b*2^a',
                    source_H='1386P-13h',source_q_s='15h',source_F='2592P-129h',
                    exit_u='3978P-157h',exit_B='14346P-703h',
                    direct_C3_guard='B-2u = h mod3, nonzero',
                    scalar_control_encoding='All coordinates are decimal strings'),growing_module_controls=growing)
    assert result['primitive_certificate']['u']==(2689128,-5495)
    assert result['primitive_certificate']['b']==(9697896,-24605)
    Path(__file__).with_name('port_escape_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(coefficients=coeff,primitive_certificate=result['primitive_certificate'],first_sample=samples[0]),indent=2))


if __name__=='__main__':run()
