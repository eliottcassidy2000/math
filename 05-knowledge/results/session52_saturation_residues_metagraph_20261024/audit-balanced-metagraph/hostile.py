#!/usr/bin/env python3
"""Reject plausible gauge, endpoint, parity, and quotient mistakes."""
from audit import matrix,scores,C_matrix,R_matrix,extend,cycle_flip,require,literal
from collections import Counter
from pathlib import Path
import argparse,copy,json

def must_reject(label,test):
    try:test()
    except RuntimeError as e:return {'test':label,'status':'REJECTED','reason':str(e)}
    raise RuntimeError('hostile accepted: '+label)

def eq(a,b,msg):require(a==b,msg)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='hostile.json');args=ap.parse_args()
    proof=literal(4);D=proof['internal_ports'];B=proof['blue_endpoint']['HL'];tests=[]
    candidates=[matrix(4,w) for w in range(8)]
    good=[A for A in candidates if set(scores(A))<={1,2} and scores(A)[0]==2 and scores(A)[-1]==1][0]
    bad=[A for A in candidates if set(scores(A))<={1,2} and scores(A)[0]==1 and scores(A)[-1]==2][0]
    tests.append(must_reject('full converse substituted for C',lambda:eq(list(map(list,zip(*good))),C_matrix(good),'path arcs must remain directed forward')))
    tests.append(must_reject('endpoint criterion reversed to LH',lambda:require(set(scores(C_matrix(bad)))<={1,2},'LH partner leaves balanced layer')))
    tests.append(must_reject('start/end transport correction signs swapped',lambda:eq(scores(C_matrix(good)),tuple(3-d-(i==0)+(i==3) for i,d in enumerate(scores(good))),'actual endpoint correction has plus at start, minus at end')))
    wrong=extend(good);n=4
    for i in range(n):wrong[i][n],wrong[n][i]=wrong[n][i],wrong[i][n]
    tests.append(must_reject('regular extension arc directions inverted',lambda:eq(scores(wrong),(2,)*5,'only z to high and low to z gives regular extension')))
    # Retaining all added-vertex arcs is wrong for internal vertices even though cycle closes.
    correct=cycle_flip(extend(good));wrong=cycle_flip(extend(good))
    for i in range(1,3):wrong[i][4]=extend(good)[i][4];wrong[4][i]=extend(good)[4][i]
    tests.append(must_reject('lift reverses old nonpath arcs but not z interior arcs',lambda:eq(wrong,correct,'all off-cycle arcs, including z arcs, must reverse')))
    A=matrix(2,0)
    tests.append(must_reject('free C pairing incorrectly includes n2',lambda:require(C_matrix(A)!=A,'n2 has zero free bits and C is identity')))
    tests.append(must_reject('odd internal-port count accepted under free C',lambda:require((D+1)%2==0,'n>=4 internal carrier is C-free')))
    tests.append(must_reject('odd blue-internal-port count accepted',lambda:require((B+1)%2==0,'blue carrier also has free C pairing')))
    proof6=literal(6)
    tests.append(must_reject('R edge quotient substituted for native C-lines',lambda:eq(proof6['optional_R_edge_orbits'],proof6['native_C_lines'],'n6 has 18 optional R edge-orbits but 30 native C-lines')))
    result={'status':'PASS','hostile_tests':tests};Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
