"""Exact membership and actual-source inverse in the two quadratic 7p lanes."""
import argparse,json
from inverse import source_certificate

SEEDS={'r29':(-1,13,529),'r43':(1,13,501)}

def invert(p,seed):
    r,a,b=SEEDS[seed];k=0
    if p<13:return None
    while a<p:
        a,b=b,40*b-13*a-10*r;k+=1
    if a!=p:return None
    return (seed,)+('r21',)*k+('b',)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('p',type=int)
    parser.add_argument('--seed',choices=tuple(SEEDS));args=parser.parse_args()
    answers={}
    for seed in (args.seed,) if args.seed else SEEDS:
        word=invert(args.p,seed)
        if word is not None:answers[seed]=source_certificate(word)
    print(json.dumps({'p':args.p,'target_H':7*args.p,'status':'CONSTRUCTED' if answers else 'ABSENT_IN_TWO_LANES',
                      'witnesses':answers,'primality':'not inferred from lane membership'}))
if __name__=='__main__':main()
