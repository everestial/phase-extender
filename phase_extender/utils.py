import resource
import collections


def print_authorship():

    """printing authorship"""

    print("#######################################################################")
    print("        Welcome to phase-extender version %i       " % 1.1)
    print("  Author: kiranNbishwa (bkgiri@uncg.edu, kirannbishwa01@gmail.com) ")
    print("#######################################################################")
    print()


""" Merge ordered dictionaries, mapping each key to lists of all values mapped to in order of occurence."""


def accumulate(data):
    acc = collections.OrderedDict()
    for dataum in data:
        for kg, vg in dataum.items():
            acc.setdefault(kg, []).append(vg)
    return acc


""" to monitor memory """


def current_mem_usage():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
