"""Four-black R4 dissolution before any context ending at a singleton."""
from pathlib import Path
import json
from moves import need, join, is_path, blue, complement
from moves_five import tilings

B = (6,4,8,3)
D = (10,0,3,6)
WORDS = [(0,1,2,3,'W'),(3,'W',2,0,1),
         (3,1,'W',2,0),(2,0,3,1,'W')]


def expand(word,k):
    return tuple(v for t in word for v in (range(4,4+k) if t == 'W' else [t]))


def verify(w):
    state = join([B,w])
    n = len(state)
    differences = []
    for word in WORDS:
        p = expand(word,len(w))
        need(is_path(state,p), 'end-dissolution word is not a path')
        need(not blue(state,p), 'end-dissolution word is blue')
        differences.append(state[p[0]].bit_count()+state[p[-1]].bit_count()-(n-1))
        state = complement(state,p)
    need(differences == [-1,4-n,n-3,1], 'endpoint formula mismatch')
    need(state == join([D,w]), 'literal labelled target mismatch')


def main():
    # W=(L,R) is the formal two-port boundary identity; the one-port
    # collapse W=(R) is checked separately, without invoking that lift.
    verify((2,0))
    verify((0,))
    rows = []
    for k in range(1,8):
        count = 0
        for prefix in tilings(k-1):
            w = join([prefix,(0,)])
            verify(w)
            count += 1
        rows.append(dict(context_order=k,all_fixed_path_sink_contexts=count,
                         black_steps=4*count))
    report = dict(B=B,D=D,words=WORDS,
                  endpoint_differences=['-1','4-n','n-3','1'],controls=rows)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
