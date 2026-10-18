"""Helper: disable numba on-disk caching for dispatchers of imported session modules (so no session file is touched)."""
from numba.core.caching import NullCache
from numba.core.dispatcher import Dispatcher


def nocache(module):
    n = 0
    for name in dir(module):
        obj = getattr(module, name)
        if isinstance(obj, Dispatcher):
            obj._cache = NullCache()
            n += 1
    return n
