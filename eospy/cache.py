
import quantities as pq

class memoized(object):
   """
   Decorator that caches a function's return value each time it is called.  If
   called later with the same arguments, the cached value is returned, and not
   re-evaluated.
   """
   def __init__(self):
       self.cache = {}
       print "creating a new instance of memoize"
       
   def __call__(self, func):
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
