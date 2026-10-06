"""Unbounded initial-one runs with a source-bounded induction expense.

The production recognizer uses actual block guards, no stored ROOT routes.
The family is rooted at 5 with the single terminal word (4).
"""
from dataclasses import dataclass
from fractions import Fraction as F
import json
import collatz_floor_transport_deadlines_20261005 as floors
import collatz_inductive_floor_receipts_20261005 as prior

CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def natural(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError("exact integer in range required")


def unit(n):
    natural(n, 3)
    if n % 2 != 1 or n % 3 == 0:
        raise ValueError("nonroot odd unit required")


def valuation(n):
    natural(n, 1)
    return (n & -n).bit_length()-1


@dataclass(frozen=True)
class Block:
    source: int
    parent: int
    run: int
    terminal: int

    @property
    def word(self):
        return (1,)*self.run+(self.terminal,)

    @property
    def expense(self):
        return self.terminal-1


def exponent_phase(parent, run):
    """Unique exponent class modulo 2*3**run; digit lifting, not a long search."""
    unit(parent)
    natural(run)
    a = 2 if parent % 3 == 1 else 1
    period, modulus = 2, 3
    for _ in range(run):
        modulus *= 3
        a = next(a+d*period for d in range(3)
                 if (pow(2,a+d*period-1,modulus)*parent+1) % modulus == 0)
        period *= 3
    minimum = run+2
    if a < minimum:
        a += ((minimum-a+period-1)//period)*period
    return a, period


def inverse_block(parent, run, terminal):
    """All phases are accepted, with the strong terminal>=run+2 guard."""
    unit(parent)
    natural(run)
    natural(terminal, run+2)
    denominator = 3**(run+1)
    numerator = (1 << (run+terminal))*parent+(1 << (run+1))
    if numerator % denominator:
        raise ValueError("ternary source guard failed")
    source = numerator//denominator-1
    unit(source)
    if source <= parent or valuation(source+1) != run+1:
        raise ValueError("source rank or maximal-run guard failed")
    return Block(source,parent,run,terminal)


def palette(parent, run, group=0):
    """Two unit children in any three successive admissible exponent phases."""
    natural(group)
    a0, period = exponent_phase(parent,run)
    result = []
    for digit in range(3):
        a = a0+period*(3*group+digit)
        denominator = 3**(run+1)
        numerator = (1 << (run+a))*parent+(1 << (run+1))
        if numerator % denominator:
            raise ValueError("phase construction failed")
        if (numerator//denominator-1) % 3:
            result.append(inverse_block(parent,run,a))
    if len(result) != 2:
        raise ValueError("unit phase split failed")
    return tuple(result)


def parse(source):
    """An O(log source)-length arithmetic block; no missing suffix exploration."""
    unit(source)
    run = valuation(source+1)-1
    checkpoint = 3**run*((source+1) >> run)-1
    terminal = valuation(3*checkpoint+1)
    parent = (3*checkpoint+1) >> terminal
    if terminal < run+2 or parent <= 1 or parent % 3 == 0:
        return None
    result = inverse_block(parent,run,terminal)
    if result.source != source:
        raise ValueError("inverse/forward source mismatch")
    return result


def recognize(source):
    """Total canonical membership test for the proper tree rooted at 5."""
    natural(source,1)
    if source % 2 != 1:
        raise ValueError("odd source required")
    current, blocks = source, []
    while current != 5:
        if current <= 1 or current % 3 == 0:
            return {"member":False,"frontier":current,"blocks":tuple(blocks)}
        block = parse(current)
        if block is None:
            return {"member":False,"frontier":current,"blocks":tuple(blocks)}
        blocks.append(block)
        current = block.parent
    return {"member":True,"frontier":5,"blocks":tuple(blocks)}


def expense_ceiling(source):
    """Largest E>=0 with 4*(4/3)**E <= source-1, using integers only."""
    unit(source)
    if source < 5:
        raise ValueError("family source at least five")
    left, right, count = 16, 3*(source-1), 0
    while left <= right:
        count += 1
        left *= 4
        right *= 3
    return count


def receipt(source):
    """Analytic source and designated-leaf floors, then a measurement degree."""
    report = recognize(source)
    if not report["member"]:
        return None
    ceiling = expense_ceiling(source)
    kmax = ceiling//2+1
    leaf,a = prior.least_leaf(source)
    epsilon = floors.weight(ceiling+1,kmax+2)
    return {"source":source,"expense_ceiling":ceiling,
            "odd_step_deadline":ceiling+1,
            "source_floor":floors.weight(ceiling,kmax),
            "leaf":leaf,"leaf_valuation":a,"index":(leaf-3)//6,
            "leaf_floor":epsilon,"degree":prior.localization_degree(epsilon),
            "readout_floor":epsilon/2,
            "moment_values_consumed":0,"target_ROOT_words_consumed":0}


def export_word(source):
    """Export only after structural membership is established."""
    report = recognize(source)
    if not report["member"]:
        raise ValueError("source outside the proved family")
    return sum((block.word for block in report["blocks"]), ())+(4,)


def rejected(call):
    try:
        call()
    except (ValueError,TypeError):
        return True
    return False


def main():
    report = {}
    phase_controls = 0
    for parent in (5,13,35,739):
        for run in range(7):
            a0,period = exponent_phase(parent,run)
            need(a0 >= run+2 and a0-period < run+2, "least strong exponent in phase")
            need((pow(2,a0-1,3**(run+1))*parent+1) % 3**(run+1) == 0,
                 "full ternary phase")
            for group in range(3):
                children = palette(parent,run,group)
                for block in children:
                    need(parse(block.source) == block, "canonical block parser")
                    current = block.source
                    for i,a in enumerate(block.word):
                        current,actual = floors.step(current)
                        need(actual == a, "independent literal valuation")
                        need(current > block.source if i+1 < len(block.word)
                             else current == parent < block.source,
                             "exact first descent at terminal")
                    rho = F(2**(run+block.terminal),3**(run+1))
                    need(block.source-1 >= rho*(parent-1), "multiplicative height")
                    need(rho >= F(4,3)**block.expense, "expense pays for growth")
                    phase_controls += 1
            if run < 2:
                old = prior.one_step_palette(parent) if run == 0 else prior.palette(parent)
                need(tuple(x.source for x in palette(parent,run)) ==
                     tuple(x.source for x in old), "old palettes recovered exactly")
    report["phase_controls"] = phase_controls

    levels, sources, counts = [5], {5}, []
    for depth in range(4):
        counts.append(len(levels))
        for source in levels:
            recognized = recognize(source)
            need(recognized["member"], "generated source recognized")
            word = export_word(source)
            need(word == prior.literal_root(source), "independent strict ROOT replay")
            expense = sum(block.expense for block in recognized["blocks"])
            cap = expense_ceiling(source)
            need(expense <= cap, "source-only expense budget")
            need(4*4**cap <= (source-1)*3**cap, "ceiling lower endpoint")
            need(4*4**(cap+1) > (source-1)*3**(cap+1), "ceiling upper endpoint")
            l,k = len(word)-1,sum((a-1)//2 for a in word)
            need(l <= expense and k <= expense//2+1, "counter bounds")
            need(floors.weight(cap,cap//2+1) <= floors.word_weight(word),
                 "analytic source floor")
        if depth < 3:
            nxt=[child.source for parent in levels for run in range(4)
                 for child in palette(parent,run)]
            need(len(set(nxt)) == len(nxt) and not sources.intersection(nxt),
                 "unique parents and disjoint generations")
            sources.update(nxt)
            levels=nxt
    report["bounded_generation"] = {"runs":[0,1,2,3],"group":0,
                                    "levels":counts,"sources":len(sources),
                                    "max_bits":max(x.bit_length() for x in sources)}

    old_levels, old_count = [5],0
    for depth in range(5):
        for n in old_levels:
            need(recognize(n)["member"], "all old hybrid controls retained")
            old_count += 1
        if depth < 4:
            old_levels=[r.source for n in old_levels for r in prior.hybrid_palette(n)]
    report["old_hybrid_inclusion_controls"] = old_count
    need(export_word(23) == (1,1,5,4), "strict new member")
    need(not prior.recognize_hybrid(23)["member"], "23 outside prior mixed tree")
    need(all(not recognize(n)["member"] or prior.recognize_hybrid(n)["member"]
             for n in range(1,23,2)), "23 least new source below itself")

    samples=[]
    for n in (5,23,35,739,palette(5,5)[0].source):
        token=receipt(n)
        word=export_word(n)
        leafword=(token["leaf_valuation"],)+word
        need(floors.replay(token["leaf"],leafword) == 1, "actual leaf identity")
        need(floors.word_weight(leafword) >= token["leaf_floor"], "analytic leaf floor")
        d=token["degree"]
        need(8*F(16,25)**d <= token["readout_floor"], "positive localized readout")
        need(d == 0 or 8*F(16,25)**(d-1) > token["readout_floor"], "least sufficient degree")
        samples.append({"source":n,"expense_ceiling":token["expense_ceiling"],
                        "odd_step_deadline":token["odd_step_deadline"],
                        "leaf":token["leaf"],"degree":d,
                        "leaf_floor_denominator_bits":token["leaf_floor"].denominator.bit_length()})
    report["measurement_controls"] = samples

    # The strong guard is sufficient, not necessary for a paid first reset.
    need(floors.replay(55,(1,1,3)) == 47, "weaker paid block exists")
    need(47 < 55 and F(54,46) < F(4,3)**2, "expense-growth bound fails without strong guard")
    need(parse(55) is None, "strong grammar does not silently use the weaker row")
    need(not recognize(7)["member"] and not recognize(27)["member"], "proper-domain hostiles")
    need(not recognize(21)["member"], "other direct ROOT ray is outside chosen root-five grammar")
    need(not recognize(85)["member"], "unit direct ROOT ray still needs a separate root boundary rule")
    report["weaker_guard_hostile"] = {"source":55,"parent":47,"word":[1,1,3],
                                    "block_length":3,"expense":2,
                                    "height_ratio":"27/23",
                                    "scope":"A paid block outside the strong grammar, not an unrooted source."}
    hostiles=(lambda:recognize(True),lambda:recognize(5.0),lambda:parse(3),
              lambda:inverse_block(5,True,5),lambda:inverse_block(5,2,3),
              lambda:inverse_block(5,2,4),lambda:palette(5,2,False),
              lambda:exponent_phase(3,2),lambda:expense_ceiling(1),
              lambda:export_word(7))
    for call in hostiles:
        need(rejected(call), "typed or guarded hostile rejected")
    report["hostiles"]=len(hostiles)
    report["checks"]=CHECKS
    report["scope"]="Unbounded run and exponent grammar with explicit floors; universal entry remains OPEN."
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
