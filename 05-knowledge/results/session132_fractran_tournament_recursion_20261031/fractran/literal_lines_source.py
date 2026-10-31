# Input t=n-3, output L_n for n>=3.
def count(n):
    output = 1
    growth = 2
    t = 0
    while t < n:
        j = 0
        while j < growth:
            output = output * 2
            j += 1
        growth += 1
        t += 1
    return output
