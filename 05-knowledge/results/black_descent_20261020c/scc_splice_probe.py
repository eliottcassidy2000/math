"""Search even black macros with a shared path edge for context extension."""
from itertools import permutations
import json
from pathlib import Path
from scc_probe import word_tournament, paths, blue, reverse_path


def relabel(out, p):
    result = [0]*len(out)
    for u, row in enumerate(out):
        for v in range(len(out)):
            if (row >> v) & 1:
                result[p[u]] |= 1 << p[v]
    return tuple(result)


def steps(out, edge):
    for p in paths(out):
        if edge in {tuple(sorted(e)) for e in zip(p,p[1:])} and not blue(out,p):
            yield reverse_path(out,p), p


def ball_two(out, edge):
    result = {}
    for middle,p in steps(out,edge):
        for end,q in steps(middle,edge):
            result.setdefault(end,(p,q))
    return result


def main():
    records = []
    for n in range(6,11):
        source = word_tournament((1,)*(n-3)+(3,))
        target = word_tournament((1,)*(n-4)+(3,1))
        edge=(0,1)
        forward=ball_two(source,edge)
        answer=None
        # The prefix is fixed pointwise; only the four changed vertices move.
        for perm in permutations(range(n-4,n)):
            p=tuple(range(n-4))+perm
            t=relabel(target,p)
            backward=ball_two(t,edge)
            common=forward.keys() & backward.keys()
            if common:
                mid=min(common)
                cert=forward[mid]+tuple(tuple(reversed(q))
                                       for q in reversed(backward[mid]))
                actual=source
                for q in cert:
                    if blue(actual,q):
                        raise RuntimeError("blue witness")
                    actual=reverse_path(actual,q)
                if actual != t:
                    raise RuntimeError("wrong endpoint")
                answer={"target_relabelling":p,"paths":cert,"final_out":t}
                break
        record={"n":n,"two_step_states":len(forward),"four_step_macro":answer}
        records.append(record)
        print(json.dumps(record),flush=True)
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(records,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
