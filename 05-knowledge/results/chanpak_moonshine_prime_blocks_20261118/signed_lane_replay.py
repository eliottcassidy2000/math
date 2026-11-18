"""Literal checks for signed defective-parent repairs."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'moonshine_prime_blocks'))
from replay import profile
P=(
 (17,19,(254,252,120,240,224,192,128,4),(2,5)),
 (73,75,(254,92,248,112,96,194,128,26),(4,6)),
 (297,299,(238,220,184,112,225,66,132,40),(1,2)),
 (133,147,(190,252,248,112,224,192,1,72),(0,4)),
 (287,301,(190,124,120,240,160,192,17,70),(1,3)),
 (581,595,(62,76,56,240,98,194,133,23),(4,7)),
)
def run():
    for h,k,A,(i,j) in P:
        assert profile(A)[0]==h
        B=list(A);B[i]^=1<<j;B[j]^=1<<i
        assert profile(tuple(B))[0]==k
        print(h,k,'reversed pair',(i,j))
    assert [7*(4**u-1)//3 for u in (3,5,7,9)]==[147,77*31,301*127,147*3*19*73]
if __name__=='__main__':run()
