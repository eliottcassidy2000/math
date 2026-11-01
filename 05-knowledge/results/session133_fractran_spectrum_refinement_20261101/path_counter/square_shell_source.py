def square_shell(n):
    # q=n>=1 parametrizes one pointed self-converse one-apex class.
    if n == 0:
        return 0
    q = n
    value = 1
    while q > 1:
        parent = q // 2
        if q % 2 == 0:
            value += parent + parent % 2
        q = parent
    return n * n + n % 2 + 2 * value - 1
