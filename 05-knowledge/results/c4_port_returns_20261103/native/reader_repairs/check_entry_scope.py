"""Ten cheap literal controls for the exact end/start scope boundary."""
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def tables(rows):
    n = len(rows)
    size = 1 << n
    ends = [[0]*n for _ in range(size)]
    h = [1]+[0]*(size-1)
    for mask in range(1, size):
        for last in range(n):
            if mask >> last & 1:
                rest = mask ^ (1 << last)
                ends[mask][last] = (sum(ends[rest][u] for u in range(n)
                                        if rows[u] >> last & 1) if rest else 1)
        h[mask] = sum(ends[mask])
    return h, ends


def profile(rows):
    n = len(rows)
    full = (1 << n)-1
    h, ends = tables(rows)
    opposite = [full ^ row ^ (1 << v) for v, row in enumerate(rows)]
    hc, starts = tables(opposite)
    need(h == hc, "converse preserves subset path counts")
    f = [sum(ends[m][v]*h[full ^ m] for m in range(full+1)) for v in range(n)]
    s = [sum(starts[m][v]*h[full ^ m] for m in range(full+1)) for v in range(n)]
    return {"H": h[-1], "E2": f, "S2": s,
            "end_defects": [x-h[-1] for x in f],
            "start_defects": [x-h[-1] for x in s]}


def reader(rows, mark, word):
    rows = list(rows)
    for bit in word:
        n = len(rows)
        rows = [row | (0 if bit == "0" and v == mark else 1 << n)
                for v, row in enumerate(rows)]+[(1 << mark) if bit == "0" else 0]
    return rows


def attach(rows, mark, clock):
    n = len(rows)
    result = list(rows)+[0]*(clock+1)
    root = n+clock
    result[mark] |= 1 << root
    for i in range(clock):
        result[n+i] = ((1 << n)-1) | sum(1 << (n+j) for j in range(i+1, clock))
    result[root] = ((1 << root)-1) ^ (1 << mark)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE/"entry_scope.json")
    args = parser.parse_args()
    seed = [8,1,3,6]
    first = attach(seed,0,1)
    need(first == [40,1,3,6,15,30], "exact H19 literal matrix")
    a = profile(first)
    need(a["H"] == 19 and a["S2"][1] == 20 and a["start_defects"][1] == 1,
         "retained vertex1 is start-ready1")
    need(not any(b in (1,2) for b in a["end_defects"]), "end-only theorem at seed order4")
    triangle = attach([2,4,1],0,1)
    b = profile(triangle)
    need(b["H"] == 9 and b["end_defects"].count(2) == 2,
         "C3 clock1 has exactly two end-ready2 vertices")
    component = reader(seed,0,"00")
    third = attach(component,0,2)
    c = profile(third)
    need(c["start_defects"][5] == 1, "last vertex of nonempty00 reader is start-ready1")
    need(not any(d in (1,2) for d in c["end_defects"]), "end-only theorem for nonempty reader")
    records = [
        {"case":"strong4_empty_reader_clock1", "rows":first, "names":[0,1,2,3,"a1","r"], **a},
        {"case":"C3_empty_reader_clock1_exception", "rows":triangle, "names":[0,1,2,"a1","r"], **b},
        {"case":"strong4_reader00_clock2", "rows":third,
         "names":[0,1,2,3,"z1","z2","a1","a2","r"], **c},
    ]
    for mark, word, q in ((0,"10100",1),(1,"100",3)):
        component = reader(seed,mark,word)
        rows = attach(component,mark,1)
        observed = profile(rows)
        last = len(component)-1
        need(observed["start_defects"][last] == q, "arbitrary reader ending00 gives sigma=q")
        records.append({"case":f"strong4_mark{mark}_reader{word}_clock1", "rows":rows,
                        "seed_q":q,"tested_vertex_index":last, **observed})
    for word, expected_h, penultimate in (("11100",469,168),("01100",475,182)):
        rows = attach(reader([4,1,2],0,word),0,2)
        observed = profile(rows)
        need(observed["H"] == expected_h and observed["start_defects"][7] == 1 and
             observed["start_defects"][6] == penultimate, "actual alternate reader entry")
        records.append({"case":f"bank0_reader{word}_clock2", "rows":rows,
                        "last_reader_index":7,"penultimate_reader_index":6,**observed})
    seed_ends = tables(seed)[1][-1]
    for mark, word, expected in ((0,"0",1),(1,"0",2),(0,"10",0)):
        component = reader(seed,mark,word)
        rows = attach(component,mark,1)
        observed = profile(rows)
        last = len(component)-1
        if word == "0":
            need(seed_ends[mark] == expected, "seed endpoint control")
        need(observed["start_defects"][last] == expected, "complete terminal-reader law")
        records.append({"case":f"terminal_law_mark{mark}_word{word}", "rows":rows,
                        "seed_endpoint_count":seed_ends[mark],"tested_vertex_index":last,**observed})
    report = {"status":"FINITE-EXACT ten explicit scope controls; no search",
              "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "records":records}
    args.output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(records,indent=2))


if __name__ == "__main__":
    main()
