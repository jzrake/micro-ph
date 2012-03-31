
import quantities as pq
from functools import wraps
from collections import OrderedDict


class LimitedSizeDict(OrderedDict):
    """
    Used by the cacheing decorator to limit the size of the cache.
    """
    def __init__(self, *args, **kwargs):
        self.size_limit = kwargs.pop("size_limit", None)
        OrderedDict.__init__(self, *args, **kwargs)
        self._check_size_limit()

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        self._check_size_limit()

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)


class memoized(object):
   """
   Decorator that caches a function's return value each time it is called.  If
   called later with the same arguments, the cached value is returned, and not
   re-evaluated.
   """
   def __init__(self, size_limit=10000):
       self.cache = LimitedSizeDict(size_limit=size_limit)
       
   def __call__(self, func):
       @wraps(func)
       def new_func(*args):
           hashed_args = tuple(self._hash_quantity(a) for a in args)
           try:
               res = self.cache[hashed_args]
               print "using hashed value", hashed_args
               return res
           except KeyError:
               print "creating new cached value for", hashed_args
               value = func(*args)
               self.cache[hashed_args] = value
               return value
       return new_func

   def _hash_quantity(self, q):
       if isinstance(q, pq.Quantity):
           return hash((float(q.magnitude), q.dimensionality.string))
       elif isinstance(q, list):
           return self._hash_quantity(tuple(q))
       elif isinstance(q, tuple):
           return tuple(self._hash_quantity(a) for a in q)
       else:
           return hash(q)
