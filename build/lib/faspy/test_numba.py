from numba import njit, jit
import time

abc = 'abc'
defa = 'defa'

@njit
def test(a,b):
    return "_".join([a, b])

start0  = time.perf_counter()
result = test(abc,defa)
end0 = time.perf_counter()
print(f"Time required for all process to complete: {end0 - start0:0.8f} seconds")
print(result)