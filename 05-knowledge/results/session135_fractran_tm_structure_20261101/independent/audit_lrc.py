#!/usr/bin/env python3
"""Independent literal machines and Farey-gap equality-sidecar checks."""
from fractions import Fraction as F
from pathlib import Path
import sys
sys.dont_write_bytecode=True
import hashlib,importlib.util,itertools,json,math
ROOT=Path(__file__).resolve().parents[1]/'lrc_machine'
def check(ok,msg):
    if not ok:raise RuntimeError(msg)
def run(n,rows):
    steps=0
    while True:
        for a,b in rows:
            if n%b==0:n=n//b*a;steps+=1;break
        else:return n,steps
        check(steps<100000,'literal budget')
def parse(xs):return [(F(x).numerator,F(x).denominator) for x in xs]
def main():
    m=json.loads((ROOT/'machines.json').read_text());rows=parse(m['five_fraction_primitive']);primitive=0
    for a,b,c,d in itertools.product(range(7),repeat=4):
        end,cost=run(11*2**a*3**b*5**c*7**d,rows)
        check((end%17==0)==(a>=b and c>=d),'primitive verdict')
        check(cost==min(a,b)+min(c,d)+1,'primitive cost');primitive+=1
    table={(r['state'],r['read']):(r['next'],r['write'],r['move']) for r in m['tm_transitions']}
    check(len(table)==12,'TM transitions');tmcases=0
    for a,b in itertools.product(range(41),repeat=2):
        tape=dict(enumerate('1'*a+'#'+'1'*b));q='H';head=steps=0
        while q not in ('A','R'):
            q,w,move=table[q,tape.get(head,'_')];tape[head]=w;head+=move;steps+=1
            check(steps<=4*(a+b+2)*(min(a,b)+2),'TM bound')
        check((q=='A')==(a>=b),'TM verdict');tmcases+=1
    whole=[]
    for case in m['whole_row_programs']:
        layout=case['layout'];state=layout['control'];pairs=case['inequalities']
        for (a,b),(p,q) in zip(pairs,layout['register_primes']):state*=p**a*q**b
        check(state.bit_length()==case['input_bit_length'],'input size')
        check(hashlib.sha256(state.to_bytes((state.bit_length()+7)//8,'big')).hexdigest()==case['input_sha256_big_endian'],'input pin')
        end,cost=run(state,parse(case['fractions']))
        t=F(case['t']);delta=F(case['delta']);den=t.denominator
        verdict=all(F(min((t.numerator*v)%den,(-t.numerator*v)%den),den)>=delta for v in case['speeds'])
        check((end%layout['accept']==0)==verdict,'whole verdict');check(cost==sum(min(a,b) for a,b in pairs)+1,'whole cost')
        whole.append({'name':case['name'],'accepted':verdict,'steps':cost,'input_bits':state.bit_length()})
    spec=importlib.util.spec_from_file_location('producer_lrc',ROOT/'check.py');p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
    harmonic=[]
    for k in range(1,121):
        # Independent Farey interval construction: nearest reduced fractions
        # a/b,c/d leave [a/b+lambda/b,c/d-lambda/d], including equality.
        fs=sorted({F(a,b) for b in range(1,k+1) for a in range(b+1)})
        result=[];widths=[]
        for left,right in zip(fs,fs[1:]):
            b,d=left.denominator,right.denominator
            check(right.numerator*b-left.numerator*d==1,'Farey determinant')
            check(b+d>k,'Farey neighbor maximality')
            l=left+F(1,(k+2)*b);r=right-F(1,(k+2)*d)
            if l<=r:result.append((l,r));widths.append(b+d)
        check(result==p.unit_harmonic_windows(k),('harmonic',k))
        check(all(s in (k+1,k+2) for s in widths),'missing denominator stratum')
        # All retained endpoints must satisfy every original integer-distance inequality.
        for l,r in result:
            for t in (l,r,(l+r)/2):
                check(all(F(min((t.numerator*j)%t.denominator,(-t.numerator*j)%t.denominator),t.denominator)>=F(1,k+2) for j in range(1,k+1)),('endpoint',k,t))
        harmonic.append({'k':k,'intervals':sum(l<r for l,r in result),'singletons':sum(l==r for l,r in result)})
    print(json.dumps({'status':'PASS','primitive_literal_cases':primitive,'literal_TM_cases':tmcases,'whole_row_runs':whole,'harmonic_Farey_range':[1,120],'harmonic_selected':[harmonic[i-1] for i in (1,3,12,80,120)]},indent=2))
if __name__=='__main__':main()
