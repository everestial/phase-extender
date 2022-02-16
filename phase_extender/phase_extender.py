import os
import sys
import time
import argparse

from phase_extender.phaser import phase_converter
from phase_extender.arg_builders import get_args
from phase_extender.val_extractor import args_to_val
from phase_extender.utils import print_authorship


def main():
    parser = argparse.ArgumentParser()
    args_namespace = get_args(parser)
    print_authorship()
    parsed_args = args_to_val(args_namespace)
    phase_converter(*parsed_args)


if __name__ == "__main__":
    main()
