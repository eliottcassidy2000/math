"""Local no-op stand-in for numba (mac-mini has no numba; no third-party install). njit returns the function unchanged."""
def njit(*args, **kwargs):
    if args and callable(args[0]) and len(args) == 1 and not kwargs:
        return args[0]
    def deco(f):
        return f
    return deco
prange = range
