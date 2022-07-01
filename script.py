# import os.path
# from multiprocessing import Pool
# import sys
# import time
# import re


# import glob
# import os
# import multiprocessing
# def process(file, string):
#     with open(file) as f:   # open file
#         for line in f:       # process line by line
#             if string in line:    # search for string
#                 print(f'found string in file {file}')
#                 break
# def main():
#     folder = r"C:\Users\Julien\Documents\University\Ete2022\Stage\Code\PhysiCell_V.1.10.1"
#
#
#     p = multiprocessing.Pool()
#     for f in glob.glob(folder+"*."):
#         # launch a process for each file (ish).
#         # The result will be approximately one process per CPU core available.
#         if os.path.isfile(f):
#             p.apply_async(lambda file: process(file,'Cell_Definition'), [f])
#         else:
#             print(f)
#
#     p.close()
#     p.join() # Wait for all child processes to close.
#
# from multiprocessing import Process, freeze_support
# if __name__ == '__main__':
#     freeze_support()
#     main()
# import os.path
# from multiprocessing import Pool
# import sys
# import time
#
# def process_file(name):
#     ''' Process one file: count number of lines and words '''
#     linecount=0
#     wordcount=0
#     with open(name, 'r') as inp:
#         for line in inp:
#             linecount+=1
#             wordcount+=len(line.split(' '))
#
#     return name, linecount, wordcount
#
# def process_files_parallel(arg, dirname, names):
#     ''' Process each file in parallel via Poll.map() '''
#     pool=Pool()
#     results=pool.map(process_file, [os.path.join(dirname, name) for name in names])
#
# def process_files(arg, dirname, names):
#     ''' Process each file in via map() '''
#     results=map(process_file, [os.path.join(dirname, name) for name in names])
#
# from multiprocessing import Process, freeze_support
# if __name__ == '__main__':
#     freeze_support()
#     start=time.time()
#     os.path.walk('C:\Users\Julien\Documents\Python\MultiTreadStringSearch\', process_files, None)
#     print("process_files()", time.time()-start)
#     start=time.time()
#     os.path.walk('C:\Users\Julien\Documents\Python\MultiTreadStringSearch\', process_files_parallel, None)
#     print("process_files_parallel()", time.time()-start)
#
# from pathlib import Path
# from concurrent.futures import ProcessPoolExecutor, as_completed
# from typing import Iterable
#
# def get_files(directory: Path) -> Iterable[Path]:
#     # Using glob simplifies the code in these cases.
#     return (file for file in directory.glob("**/*") if file.is_file())
#
# def foo(file: Path):
#     ...
#
# # This can be easily called from a test
# def process_files(directory: Path):
#     # I am using sum, so that the result is computed lazily and we
#     # do not need to build a list of all files. If the directory is
#     # very large, this could save a lot of memory.
#     # Since get_files returns a one-shot generator, we cannot
#     # save it to a variable and reuse it here and below.
#     file_count = sum(1 for _ in get_files(directory))
#
#     with concurrent.futures.ProcessPoolExecutor() as executor:
#         futures = executor.map(foo, get_files(directory))
#         # Reading the value returned by the executor is very
#         # important, because exceptions happening in the `foo`
#         # function will be captured and not be raised until the
#         # value is read, thus obfuscating the errors.
#         #
#         # Another nice solution for the progress would be
#         # using the library tqdm.
#         for i, _ in enumerate(as_completed(futures)):
#             print(f"Processed file {i+1} / {file_count}")
#
# if __name__ == "__main__":
#     process_files(Path(""))
#







# def process_file(name, string):
#     ''' Process one file: count number of lines and words '''
#     linecount=0
#     wordcount=0
#
#     if re.search('mandy', 'Mandy Pande', re.IGNORECASE):
#         return name, linecount, wordcount
#
# def process_files_parallel(arg, dirname, names):
#     ''' Process each file in parallel via Poll.map() '''
#     pool=Pool()
#     results=pool.map(process_file, [os.path.join(dirname, name) for name in names])
#
# def process_files(arg, dirname, names):
#     ''' Process each file in via map() '''
#     results=map(process_file, [os.path.join(dirname, name) for name in names])
#
# if __name__ == '__main__':
#     start=time.time()
#     os.path.walk('input/', process_files, None)
#     print "process_files()", time.time()-start
#
#     start=time.time()
#     os.path.walk('input/', process_files_parallel, None)
#     print "process_files_parallel()", time.time()-start







from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Iterable
from multiprocessing import Process, freeze_support
def get_files(directory: Path) -> Iterable[Path]:
    # Using glob simplifies the code in these cases.
    return (file for file in directory.glob("**/*") if file.is_file())

def foo(file: Path):
    ...

# This can be easily called from a test
def process_files(directory: Path):
    freeze_support()

    # I am using sum, so that the result is computed lazily and we
    # do not need to build a list of all files. If the directory is
    # very large, this could save a lot of memory.
    # Since get_files returns a one-shot generator, we cannot
    # save it to a variable and reuse it here and below.
    file_count = sum(1 for _ in get_files(directory))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = executor.map(foo, get_files(directory))
        # Reading the value returned by the executor is very
        # important, because exceptions happening in the `foo`
        # function will be captured and not be raised until the
        # value is read, thus obfuscating the errors.
        #
        # Another nice solution for the progress would be
        # using the library tqdm.
        for i, _ in enumerate(as_completed(futures)):
            print(f"Processed file {i+1} / {file_count}")

if __name__ == "__main__":
    freeze_support()
    process_files(Path(""))
