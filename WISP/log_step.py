
import contextlib
import sys
import time
import os

@contextlib.contextmanager
def suppress_output():
    """
    Context manager that temporarily redirects stdout and stderr to os.devnull,
    silencing any print or logging output within its block.

    Usage:
        with suppress_output():
            # code whose stdout/stderr you want to suppress
            ...
    """
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr

class LogStep:
    def __init__(self, name):
        self.name = name
        self.start_time = None

    def start(self):
        self.start_time = time.time()
        print(f"Starting step: {self.name}")
        sys.stdout.flush()

    def end(self):
        end_time = time.time()
        duration = end_time - self.start_time
        print(f"Ending step: {self.name} [{duration:.2f} seconds]")
        sys.stdout.flush()

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.end()