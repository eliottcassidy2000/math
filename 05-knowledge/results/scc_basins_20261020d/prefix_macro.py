"""Universal prefix interchange: exact positional tests and class census."""
import csv
import json
from itertools import combinations
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent/"black_descent_20261020c"))
from scc_probe import paths, blue, reverse_path, need, decode


def induced(out, vertices):
    return tuple(sum(1 << j for j,v in enumerate(vertices) if out[u] >> v & 1)
                 for u in vertices)


def macro_source_target(A):
    m = len(A)
    x,a,b,c = range(m,m+4)
    tail = sum(1 << v for v in (x,a,b,c))
    source = tuple(r | tail for r in A) + (
        (1 << a) | (1 << b) | (1 << c), 1 << b, 1 << c, 1 << a)
    target = tuple(r | tail for r in A) + (
        (1 << a) | (1 << c), (1 << b) | (1 << c),
        (1 << x) | (1 << c), 0)
    return source,target


def macro(A,p):
    m = len(A)
    x,a,b,c = range(m,m+4)
    source,target = macro_source_target(A)
    words = (p+(x,c,a,b),p[::-1]+(c,x,b,a),
             (x,a,b,c)+p,(c,b,a,x)+p[::-1])
    actual = source
    colours = []
    for word in words:
        need(sorted(word) == list(range(m+4)), "not spanning")
        colours.append(blue(actual,word))
        actual = reverse_path(actual,word)
    need(actual == target, "incorrect literal endpoint")
    return tuple(colours),words


def has_join(out, left, right):
    return all(out[u] >> v & 1 for u in left for v in right)


def cycle3(out, vertices):
    return all(sum(out[u] >> v & 1 for v in vertices if v != u) == 1
               for u in vertices)


def predicted_boundary_blues(A,p):
    """SCC form of the first and last possible blue step."""
    m = len(A)
    first = ((m == 3 and cycle3(A,p)) or
             (m >= 4 and cycle3(A,p[:3]) and
              has_join(A,p[:3],p[3:]) and
              has_join(A,p[3:4],p[4:]) and
              blue(induced(A,p[4:]),tuple(range(m-4)))))
    last = ((m == 1) or
            (m >= 4 and has_join(A,p[:1],p[1:]) and
             cycle3(A,p[1:4]) and has_join(A,p[1:4],p[4:]) and
             blue(induced(A,p[4:]),tuple(range(m-4)))))
    return first,last


def standard_path_tilings(m):
    pairs = [(i,j) for i,j in combinations(range(m),2) if j>i+1]
    for code in range(1 << len(pairs)):
        out = [sum(1 << j for j in range(i+1,m)) for i in range(m)]
        for k,(i,j) in enumerate(pairs):
            if code >> k & 1:
                out[i] ^= 1 << j
                out[j] ^= 1 << i
        yield tuple(out)


def local_audit():
    receipt = []
    for m in range(1,8):
        counts = [0]*4
        total = 0
        for A in standard_path_tilings(m):
            p = tuple(range(m))
            colours,_ = macro(A,p)
            first,last = predicted_boundary_blues(A,p)
            need((colours[0],colours[3]) == (first,last), "SCC blue description")
            need(not colours[1], "second step blue")
            if m>=2 and A[p[-1]].bit_count() != m-3:
                need(not colours[2], "endpoint-score criterion failed")
            for j,c in enumerate(colours):
                counts[j] += c
            total += 1
        receipt.append({"prefix_order":m,"fixed_path_tilings":total,
                        "blue_counts_by_step":counts})
    return receipt


def class_audit():
    receipt = []
    for m in range(3,9):
        table = ROOT/f"04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n{m}.tsv"
        with table.open() as f:
            rows = list(csv.DictReader(f,delimiter="\t"))
        excluded = []
        success = 0
        max_trials = 0
        hardest = None
        for row in rows:
            A = decode(m,int(row["canon"]))
            first_all = last_all = True
            witness = None
            trials = 0
            for p in paths(A):
                colours,words = macro(A,p)
                first_all &= colours[0]
                last_all &= colours[3]
                trials += 1
                if not any(colours):
                    witness = words
                    break
            if witness is not None:
                success += 1
                if trials > max_trials:
                    max_trials = trials
                    hardest = {"prefix_code":int(row["canon"]),"paths":witness}
            else:
                # Every source/target HP is obtained from an A path and one of
                # the three cyclic rotations of its C3. These rotations are
                # automorphisms, so first_all/last_all certify actual purity.
                need(first_all or last_all, "active endpoints without a black macro")
                excluded.append({"prefix_code":int(row["canon"]),
                                 "pure_source":first_all,"pure_target":last_all,
                                 "all_prefix_paths":trials})
        receipt.append({"prefix_order":m,"ordinary_classes":len(rows),
                        "successful_classes":success,"excluded_pure_endpoint":excluded,
                        "maximum_paths_tried":max_trials,"hardest_witness":hardest})
    return receipt


def main():
    result = {"fixed_path_audit":local_audit(),"all_class_existence":class_audit()}
    HERE.joinpath("prefix_macro.json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8")
    summary = {**result,"all_class_existence":[
        {k:v for k,v in r.items() if k not in ("hardest_witness","excluded_pure_endpoint")}
        | {"excluded_classes":len(r["excluded_pure_endpoint"])}
        for r in result["all_class_existence"]]}
    text = json.dumps(summary,indent=2)+"\n"
    HERE.joinpath("prefix_macro.out").write_text(text,encoding="utf-8")
    print(text,end="")


if __name__ == "__main__":
    main()
