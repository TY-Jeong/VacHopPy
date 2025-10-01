import os
import time
import inspect
import functools
import tracemalloc
import numpy as np


def monitor_performance(func):
    """
    A decorator that measures and prints the execution time 
    and peak memory usage of a function.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        is_verbose = False
        try:
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()
            is_verbose = bound_args.arguments.get('verbose', False)
        except TypeError:
            pass

        if not is_verbose:
            return func(*args, **kwargs)

        tracemalloc.start()
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        _, peak_mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        print(f"Execution Time: {elapsed_time:.3f} seconds")
        print(f"Peak RAM Usage: {peak_mem / 1024**3:.3f} GB")
        
        return result
    return wrapper