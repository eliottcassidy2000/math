"""Source-owned delayed paid returns; exact guards, no universal acceptance claim."""
from dataclasses import dataclass, replace
from fractions import Fraction
from pathlib import Path
import ast
import hashlib
import collatz_charged_terminal_family_20261012a as old


def need(ok, message):
    if not ok:
        raise ValueError(message)


def E(n):
    old.odd(n)
    return (n+1).bit_length()-1


def rate_value(rate):
    need(type(rate) is int and rate >= 1, 'positive exact rate')
    return rate


def deadline(source, budget):
    """Least inherited ROOT-or-budget deadline, computed by exact binary search."""
    old.odd(source), old.natural(budget)
    if source == 1:
        return 0
    def enough(t):
        return ((source-1)*pow(3, t) << budget) < (1 << (2*t+1))
    lo, hi = -1, 4*(budget+source.bit_length())
    need(enough(hi), 'proved finite ceiling')
    while hi-lo > 1:
        mid = (hi+lo)//2
        if enough(mid):
            hi = mid
        else:
            lo = mid
    return budget+hi


def paid(source, endpoint, charge, rate):
    old.odd(source), old.odd(endpoint), old.natural(charge), rate_value(rate)
    return 2*(endpoint+1) <= source+1 and charge <= rate*(E(source)**2-E(endpoint)**2)


@dataclass(frozen=True)
class Return:
    source: int
    rate: int
    word: tuple
    status: str


def audit_return(result):
    need(type(result) is Return, 'exact Return')
    old.odd(result.source), rate_value(result.rate), old.word(result.word)
    need(type(result.status) is str and result.status in ('ROOT', 'PAID', 'PENDING', 'EXHAUSTED'),
         'known exact status')
    if result.source == 1:
        need(result.word == () and result.status == 'ROOT', 'unique ROOT receipt')
        return 1, 0, None
    need(result.status != 'ROOT', 'ROOT status only for initial ROOT')
    n, charge, first_half = result.source, 0, None
    budget = result.rate*(E(result.source)**2-1)
    for i, a in enumerate(result.word):
        n = old.replay(n, (a,))
        charge += a % 2
        if first_half is None and 2*(n+1) <= result.source+1:
            first_half = (i+1, n, charge)
        if i+1 < len(result.word):
            need(not paid(result.source, n, charge, result.rate) and charge <= budget,
                 'first decisive event was not discarded')
    is_paid = paid(result.source, n, charge, result.rate)
    expected = 'PAID' if is_paid else ('EXHAUSTED' if charge > budget else 'PENDING')
    need(result.status == expected, 'status bound to actual prefix')
    need(len(result.word) <= deadline(result.source, budget), 'finite source deadline')
    return n, charge, first_half


def search(source, rate=3, step_cap=None):
    old.odd(source), rate_value(rate)
    if step_cap is not None:
        old.natural(step_cap)
    initial = Return(source, rate, (), 'ROOT' if source == 1 else 'PENDING')
    return resume(initial, step_cap)


def resume(result, step_cap=None):
    """Resume authenticated debt at the SAME source; cap counts new odd steps."""
    n, charge, _ = audit_return(result)
    if step_cap is not None:
        old.natural(step_cap)
    if result.status != 'PENDING':
        return result
    budget = result.rate*(E(result.source)**2-1)
    end = deadline(result.source, budget)
    if step_cap is not None:
        end = min(end, len(result.word)+step_cap)
    word = list(result.word)
    status = 'PENDING'
    while len(word) < end:
        a = old.valuation(3*n+1)
        n = (3*n+1) >> a
        word.append(a)
        charge += a % 2
        if paid(result.source, n, charge, result.rate):
            status = 'PAID'
            break
        if charge > budget:
            status = 'EXHAUSTED'
            break
    need(status != 'PENDING' or len(word) < deadline(result.source, budget), 'deadline forced a decisive event')
    answer = Return(result.source, result.rate, tuple(word), status)
    audit_return(answer)
    return answer


def ceil_div(a, b):
    return -((-a)//b)


@dataclass(frozen=True)
class Rule:
    word: tuple
    rate: int
    drop: int


def rule_data(rule):
    need(type(rule) is Rule, 'exact Rule')
    old.word(rule.word), rate_value(rule.rate)
    need(bool(rule.word) and type(rule.drop) is int and rule.drop >= 1, 'nonempty word and positive exact drop')
    p, q, b = old.carrier(rule.word)
    scale = 1 << rule.drop
    gap = q-scale*p
    need(gap > 0, 'strict leading contraction')
    charge = sum(a % 2 for a in rule.word)
    minimum_E = max(rule.drop+1, ceil_div(charge+rule.rate*rule.drop**2, 2*rule.rate*rule.drop))
    cut = max(3, (1 << minimum_E)-1, ceil_div(scale*b+(scale-1)*q, gap))
    for i in range(len(rule.word)):
        pp, qq, bb = old.carrier(rule.word[:i])
        cut = max(cut, ceil_div(3*qq-bb, pp))
    modulus = 2*q
    residue = ((q-b)*pow(p, -1, modulus)) % modulus
    first = residue+max(0, ceil_div(cut-residue, modulus))*modulus
    return dict(P=p, Q=q, B=b, charge=charge, modulus=modulus, residue=residue,
                first=first, cut=cut, minimum_E=minimum_E)


def compile_return(receipt):
    """Compile the exact word into a tail AP; may exclude the discovery point."""
    audit_return(receipt)
    need(receipt.status == 'PAID', 'supplied paid receipt')
    p, q, _ = old.carrier(receipt.word)
    candidates = []
    drop = 1
    while (p << drop) < q:
        rule = Rule(receipt.word, receipt.rate, drop)
        data = rule_data(rule)
        candidates.append((data['first'], drop, rule))
        drop += 1
    need(bool(candidates), 'actual half exit forces strict leading half contraction')
    return min(candidates, key=lambda row: row[:2])[2]


def apply_rule(rule, source):
    old.odd(source)
    data = rule_data(rule)
    if source < data['first'] or source % data['modulus'] != data['residue']:
        return None
    endpoint = old.replay(source, rule.word)
    need((endpoint+1) << rule.drop <= source+1, 'full carry payment')
    need(data['charge'] <= rule.rate*(E(source)**2-E(endpoint)**2), 'source-owned quadratic invoice')
    return endpoint, data['charge']


def lift_data(receipt):
    """One authenticated point plus the continuous tail covers its entire native AP."""
    endpoint, charge, _ = audit_return(receipt)
    need(receipt.status == 'PAID', 'paid discovery point')
    rule = compile_return(receipt)
    data = rule_data(rule)
    need(receipt.source % data['modulus'] == data['residue'], 'same native source cylinder')
    need(data['first'] <= receipt.source+data['modulus'], 'no unverified native source between seed and tail')
    return rule, data, endpoint, charge


def apply_lift(receipt, source):
    old.odd(source)
    rule, data, endpoint, charge = lift_data(receipt)
    if source < receipt.source or (source-receipt.source) % data['modulus']:
        return None
    if source == receipt.source:
        return endpoint, charge
    return apply_rule(rule, source)


def even_absorption_cut(old_drop, extra_drop, even_steps):
    """Minimum integer n+1 from the exact carry-aware even-block estimate."""
    need(all(type(v) is int and v >= 1 for v in (old_drop, extra_drop, even_steps)), 'positive exact parameters')
    q, p = 4**even_steps, 3**even_steps
    gap = q-(p << extra_drop)
    need(gap > 0, 'even block contracts enough')
    return ceil_div((q-p) << (old_drop+extra_drop+1), gap)


def record_sources():
    # Read only the inherited literal candidate table, never execute that producer.
    path = Path(__file__).with_name('collatz_band_budget_hostile_20261012.py')
    tree = ast.parse(path.read_text(encoding='utf-8'))
    rows = next(ast.literal_eval(stmt.value) for stmt in tree.body
                if isinstance(stmt, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'GLIDE_RECORDS' for t in stmt.targets))
    need(type(rows) is list and len(rows) == 34, 'frozen candidate-table cardinality')
    for source, _ in rows:
        old.odd(source)
    return tuple(source for source, _ in rows)


# Finite discovery only: on 2026-10-13, each of these nine endpoints was
# iterated to first ROOT with a declared 10,000-odd-step cap. Production never
# runs that discovery: it decodes and strictly authenticates the frozen word.
# A hexadecimal digit is one valuation (not a byte pair).
ROOT_SUFFIXES = (
    (14157451, '12314112222321141122143311111212332223511114121121112311212111113111422431154', 77, 146, 50, '9e494369e0b2be6638538ade5c9501c2ff67b0c77569870e7a063713eaad7a39'),
    (10445855, '11112211112211112112131c211112112112114411211122114711112115131111111112311111211322212321215122312242111111243331234', 117, 209, 81, '48eb55188dffb34568511758600f5a3febd8f7967f5b0b020301fb0c9224ac52'),
    (249622248673, '221111161213212312221212161311252441311111151343113123111261112221642131234', 75, 157, 47, '09f77208f3d6a9a3c95d570d99cda14e5cb5ef4ab58e00221120dbe15caa551d'),
    (522082008095, '11112211311111211311212711122117112212121512421222211521236222211152212312114611121112131171114121121112311212111113111422431154', 128, 242, 84, '86d2dd97548a81c2e25ccb402d9f8e8b1bd4530594a756310810fd723c2e8f16'),
    (4994812878239, '1111212123121611132522111372213111151211331121211432121312231122211311112121131133113121212116213123341221224321231214131234', 124, 239, 83, '876e31fe91b56c14887b0e6807c275111acfe1a910288e9963b6ec10139058f2'),
    (141295968112997, '41412322111311113121311221172113322121131125113311215361311211215131211132631211112122121211311222211113212111121223312153325154', 128, 250, 90, 'c53df30ccf624946eafa9e2fbc3656290cd6d48e35a97aa1cab1622179eeda48'),
    (22391847502773035, '122321214131111232111221113312132341134222123122212361151211241223441312432122112112121211351131112121251421212113222223321121112311212111113111422431154', 153, 297, 93, '1b28654e64068d361b3d1f4da49b05471ba98a985f9adb641ee96a2a96ac6c45'),
    (29997399390080243, '1321135132122113253213114521111943152114211251212112111222312115242211421121122311213121213112122126328', 103, 218, 64, 'c350048c5b24bef0d12c7a9b576ef99c2b77d967b42ebdff8833813da38eacca'),
    (889585424873412547, '14111231322421111151325441116221231213326113131232112332225d121112131171114121121112311212111113111422431154', 108, 231, 73, '64a766d512dd228c11da3e351674dff70b684624c227099ce6a052661e79ab32'),
)


def root_suffix(endpoint):
    """Decode a supplied finite bank entry; no orbit discovery or near-source reuse."""
    old.odd(endpoint)
    row = next((row for row in ROOT_SUFFIXES if row[0] == endpoint), None)
    need(row is not None, 'endpoint absent from finite ROOT bank')
    _, digits, rank, cost, charge, digest = row
    need(type(digits) is str and all(c in '123456789abcdef' for c in digits), 'exact positive valuation digits')
    w = tuple(int(c, 16) for c in digits)
    need((len(w), sum(w), sum(a % 2 for a in w)) == (rank, cost, charge), 'ROOT word metadata')
    need(hashlib.sha256(bytes(w)).hexdigest() == digest, 'ROOT word digest')
    need(old.replay(endpoint, w) == 1, 'strict supplied endpoint ROOT word')
    return w


def discharge(receipt, supplied_child_word):
    """Substitute an actual child ROOT word at exactly its retained source port."""
    endpoint, _, _ = audit_return(receipt)
    old.word(supplied_child_word)
    need(receipt.status in ('PAID', 'ROOT'), 'ground only a completed dependency')
    need(old.replay(endpoint, supplied_child_word) == 1, 'supplied exact-child first-hit ROOT')
    full = receipt.word+supplied_child_word
    need(old.replay(receipt.source, full) == 1, 'strict original-source ROOT splice')
    return full


def main():
    checks = 0
    def check(ok, why):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(why)
    def reject(fn):
        try:
            fn()
        except ValueError:
            check(True, 'hostile rejected')
        else:
            check(False, 'hostile accepted')

    for source in range(3, 128, 2):
        for budget in range(9):
            check(deadline(source, budget) == old.charge_deadline(source, budget), 'independent deadline path')
    rows, rooted = [], []
    for source in record_sources():
        receipt = search(source)
        endpoint, charge, first = audit_return(receipt)
        check(receipt.status == 'PAID', 'declared 34-source exact outcome')
        check(old.replay(source, receipt.word) == endpoint, 'independent strict receipt replay')
        first_length, first_endpoint, first_charge = first
        if not paid(source, first_endpoint, first_charge, 3):
            prefix = search(source, step_cap=first_length)
            check(prefix.status == 'PENDING' and resume(prefix) == receipt, 'unpaid first exit retained and resumed')
            rule, data, _, _ = lift_data(receipt)
            check(data['first'] <= source+data['modulus'], 'no omitted native source')
            for t in (0, 1, 2, 17, 2**30+1):
                n = source+data['modulus']*t
                target, cost = apply_lift(receipt, n)
                check(cost == charge and 2*(target+1) <= n+1, 'all-height return receipt control')
                check(data['P']*n+data['B'] == data['Q']*target, 'independent affine endpoint')
                back = target
                for a in reversed(receipt.word):
                    numerator = (back << a)-1
                    check(numerator % 3 == 0 and numerator//3 > 1, 'strict inverse replay')
                    back = numerator//3
                check(back == n, 'backward source identity')
            rows.append((source, first_length, len(receipt.word), first_charge, charge,
                         E(source)-E(endpoint), rule.drop, data['Q'].bit_length(),
                         data['first'] > source,
                         hashlib.sha256(bytes(receipt.word)).hexdigest()))
            suffix = root_suffix(endpoint)
            full = discharge(receipt, suffix)
            back = 1
            for a in reversed(full):
                numerator = (back << a)-1
                check(numerator % 3 == 0 and numerator//3 > 1, 'ROOT independent backward replay')
                back = numerator//3
            check(back == source, 'grounded original-source identity')
            rooted.append((source, endpoint, len(suffix), sum(suffix), sum(a % 2 for a in suffix),
                           len(full), sum(full), sum(a % 2 for a in full),
                           hashlib.sha256(bytes(full)).hexdigest()))
    check(len(rows) == 9, 'exact hostile subuniverse')
    check(sum(r[-2] for r in rows) == 3, 'exact native-jump count')
    check([r[2]-r[1] for r in rows] == [2,1,2,24,8,1,3,1,1], 'extra return clocks')
    least = search(56924955)
    check(least.word[-2:] == (2,2) and audit_return(least)[:2] == (14157451,151), 'charge-free least repair')
    check(search(1008932249296231).word[-1] == 3, 'C4 hostile paid without raising coefficient')
    # Slope alone loses the ordered carry at a legal actual source.
    hostile = Rule((2,1,1,2,3), 3, 1)
    hd = rule_data(hostile)
    check((hd['P'],hd['Q'],hd['B']) == (243,512,373) and hd['cut'] == 49, 'leading projection carries a finite debt')
    check(old.replay(9,hostile.word) == 5 and 2*(5+1)>9+1 and apply_rule(hostile,9) is None, 'source9 false linearized exit')
    check(even_absorption_cut(1,1,3) == 30, 'three-free-letter repair cut')
    for source in range(29, 512, 2):
        for x in range(3, (source-1)//2+1, 2):
            y, used, all_even = x, 0, True
            while y != 1 and used < 3:
                a = old.valuation(3*y+1)
                if a % 2:
                    all_even = False
                    break
                y = (3*y+1) >> a
                used += 1
            if all_even:
                check(4*(y+1) <= source+1, 'zero-charge absorption with first-hit boundary')
    for source in range(3,1024,2):
        result = search(source, step_cap=2)
        completed = resume(result)
        check(completed == search(source), 'cap and resumption preserve exact first paid event')
    check(search(27,1).status == 'EXHAUSTED', 'failed invoice is a retained finite outcome')
    check(search(1) == Return(1,3,(),'ROOT'), 'unique initial ROOT interface')
    check(search(9,step_cap=0).status == 'PENDING', 'zero-work cap retains source')
    good = search(56924955, step_cap=198)
    for fn in (lambda: search(True),lambda: search(9,3.0),lambda: resume(replace(good,source=56924955.0)),
               lambda: resume(replace(good,rate=True)),lambda: resume(replace(good,status='PAID')),
               lambda: resume(Return(1,3,(2,),'ROOT')),lambda: rule_data(replace(hostile,drop=True)),
               lambda: apply_rule(hostile,9.0),lambda: even_absorption_cut(1,1,2),
               lambda: compile_return(good),lambda: search(9,step_cap=-1),
               lambda: paid(9,5,True,3),lambda: root_suffix(14157451.0),
               lambda: root_suffix(14157453),lambda: discharge(least, root_suffix(14157451)+(2,)),
               lambda: discharge(least, root_suffix(10445855)),lambda: discharge(good, ())):
        reject(fn)
    check(discharge(search(1), ()) == (), 'unique grounded ROOT empty word')
    print('Adaptive original-source return; global acceptance remains OPEN')
    print('Candidate universe: 34 inherited glide sources; coefficient3 first-exit failures:9')
    for row in rows:
        print('source,first_length,paid_length,first_charge,paid_charge,bit_drop,rule_drop,guard_bits,lattice_jump,word_sha256:',row)
    print('All9 exact native APs are paid: one checked seed plus generic tail;3 require native-lattice jump')
    print('Leading-only hostile:9 --(2,1,1,2,3)-->5;243/512<1/2 but half-band fails; cut49')
    print('Three subsequent even letters (or earlier ROOT) absorb another bit at source>=29')
    print('Finite ROOT discovery: exactly9 named endpoints, cap10000 each; frozen suffixes only in production')
    for row in rooted:
        print('source,child,suffix_rank,suffix_cost,suffix_charge,ROOT_rank,ROOT_cost,ROOT_charge,ROOT_word_sha256:', row)
    print('All9 original points are grounded; no infinite-family child ROOT or universal charge bound is inferred')
    print('checks=',checks,sep='')


if __name__ == '__main__':
    main()
