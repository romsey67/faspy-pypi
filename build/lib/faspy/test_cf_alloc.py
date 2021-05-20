from risk import var_utils2 as vu2, var_utils as vu
import time

start0 = time.perf_counter()
result1 = vu2.get_cfalloc2(0.01,0.01,0.1, 0.1)
print(result1)
end0 = time.perf_counter()
print(f"Time required for all process to complete: {end0 - start0:0.6f} seconds")

start1 = time.perf_counter()
result2 = vu.get_cfalloc(0.01,0.01,0.1, 0.1)
print(result2)
end1 = time.perf_counter()
print(f"Time required for all process to complete: {end1 - start1:0.6f} seconds")

start2 = time.perf_counter()
result1 = vu.get_cfweight(1,3,2 )
print(result1)
end2 = time.perf_counter()
print(f"Time required for all process to complete: {end2 - start2:0.8f} seconds")


start3 = time.perf_counter()
result1 = vu2.get_cfweight2(1,3,2 )
print(result1)
end3 = time.perf_counter()
print(f"Time required for all process to complete: {end3 - start3:0.8f} seconds")