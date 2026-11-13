"""Which small prime-lane surplus failures are repaired by actual frame moves?"""
from pathlib import Path
from math import isqrt
import gzip,json

HERE=Path(__file__).resolve().parent

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def atom(n):return prime(n) and n!=7 or n%7==0 and prime(n//7) and n//7>=5

def run():
    head=json.loads(gzip.decompress((HERE/'surplus.json.gz').read_bytes()))['head']
    frames=json.loads(gzip.decompress((HERE/'frames.json.gz').read_bytes()))
    rows=[]
    for rec in head:
        h=rec['h']
        if h<3:continue
        for digit in (2,4):
            N=3*h+digit
            if not atom(N) or str(digit//2) in rec['witnesses']:continue
            hits=[]
            for t in frames['transitions']:
                if t['digit']!=digit:continue
                a,b=t['odd_amplitude'],t['offset'];e=t['min_exponent']
                while a*2**e+b<=h:
                    if a*2**e+b==h:hits.append(dict(amplitude=a,offset=b,exponent=e,source=t['source'],target=t['target']))
                    e+=1
            rows.append(dict(h=h,target=N,digit=digit,frame_families=hits))
    result=dict(status='FINITE-EXACT',scope='Parent counts3..73 present in the complete order8 bank; prime and7p target lanes; absence of direct surplus roots requires the separate all-order padding argument. Frame repairs are certified infinite identities, but this is not complete across all frames.',rows=rows)
    suffix='' if __debug__ else '.optimized'
    (HERE/('triage'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print([(r['h'],r['target'],[(t['amplitude'],t['offset'],t['exponent']) for t in r['frame_families']]) for r in rows],flush=True)

if __name__=='__main__':run()
