"""Literal same-count alternative source141 ->215 ->425 and a scoped ready family."""
import json
from pathlib import Path

from audit_compiler_selector import construct, g
from check_five_cores import CORES, check, converse, paths_dp, rooted_profiles, strong


def attach_out(rows, out):
    n = len(rows)
    return tuple(row | (0 if out >> v & 1 else 1 << n)
                 for v, row in enumerate(rows)) + (out,)


def endpoint_table(rows):
    n = len(rows)
    full = (1 << n)-1
    table = [[0]*n for _ in range(full+1)]
    for v in range(n):
        table[1 << v][v] = 1
    for mask in range(1,full+1):
        for v, count in enumerate(table[mask]):
            if not count:
                continue
            available = rows[v] & ~mask
            while available:
                bit = available & -available
                available ^= bit
                table[mask | bit][bit.bit_length()-1] += count
    return [1]+[sum(row) for row in table[1:]], table


def ready_word_controls():
    cases = queries = 0
    for m in range(2,9):
        transitive = tuple(sum(1 << w for w in range(v+1,m)) for v in range(m))
        for interior in range(1 << (m-2)):
            integer = (1 << (m-1)) | (interior << 1)
            bits = [int(bit) for bit in bin(integer)[2:]]
            rows = attach_out(transitive,sum(bit << v for v,bit in enumerate(bits)))
            hp,_ = endpoint_table(rows)
            _,starts = endpoint_table(converse(rows))
            full = (1 << len(rows))-1
            gaps = [sum(starts[mask][v]*hp[full ^ mask]
                        for mask in range(1,full+1))-hp[full] for v in range(len(rows))]
            expected = set()
            if bits[-2] == 0:
                expected.add(m-1)
            if m>=3 and all(bits[:-1]):
                expected.add(m-2)
            check({v for v,gap in enumerate(gaps) if gap == 1} == expected,
                  "complete one-apex converse1-ready root classification")
            check(hp[full] == g(integer), "one-apex binary count mismatch")
            cases += 1
            queries += len(rows)
    return {"all_strong_normalized_cut_words_base_orders_2_through_8":cases,
            "all_root_queries":queries,"scope":"Finite controls of the written all-order proof"}


def main():
    parent = construct(CORES[0], 0, 1, 61)
    middle = attach_out(parent, 128)
    child = attach_out(middle, 512)
    check([paths_dp(t) for t in (parent,middle,child)] == [141,215,425], "count mismatch")
    check(all(strong(t) for t in (parent,middle,child)), "strongness mismatch")
    profile = rooted_profiles(converse(parent))[7]
    check(profile["ready_gap"] == 1 and profile["s_start"] == 68
          and profile["q_delete"] == 69, "converse readiness mismatch")
    family = []
    for m in range(4,9):
        b = (1 << m)-3
        source = construct(CORES[0],0,1,b)
        root = len(source)-2
        first = attach_out(source,1 << root)
        second = attach_out(first,1 << len(source))
        expected = 9*(1 << (m-2))-3
        check(paths_dp(source) == expected, "source family count")
        check(paths_dp(second) == 3*expected+2, "source family branch")
        family.append({"m":m,"B":b,"source_H":expected,"output_H":3*expected+2})
    intersections = []
    for k in range(1,9):
        target = 34*(1 << k)+5
        us = [u for u in range(1,(target-1)//3+1) if g(u)+3*u+(u&1) == target]
        intersections.append({"k":k,"H":target,"one_apex_ready_parameters_u":us,
                              "cut_words":[bin(4*u)[2:] for u in us]})
    report = {"status":"PROVED literal alternative source and scoped family; no general source-return theorem",
              "parent":{"rows":list(parent),"H":141,"order":len(parent),
                        "five_core_index":0,"out_mask":0,"A":1,"B":61},
              "middle":{"rows":list(middle),"H":215,"new_outmask":128},
              "child":{"rows":list(child),"H":425,"new_outmask":512},
              "converse_parent_root":7,"converse_ready_profile":profile,
              "family_m_ge_4":family,"neutral_C4_value_intersections":intersections,
              "complete_one_apex_ready_root_controls":ready_word_controls(),
              "scope":"Alternative same-H realizer; no H-neutral native transformation from either hostile C4 source is supplied. Two births, not a one-birth certificate."}
    Path(__file__).with_name("alternate141.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"counts":[141,215,425],"converse_profile":profile,
                      "intersections":intersections},indent=2))


if __name__ == "__main__":
    main()
