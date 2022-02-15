import os
import sys
import time
import argparse

from phaser import phase_converter
from arg_builders import get_args
from val_extractor import args_to_val
from utils import print_authorship


def main():
    parser = argparse.ArgumentParser()
    args_namespace = get_args(parser)
    print_authorship()
    parsed_args = args_to_val(args_namespace)
    phase_converter(*parsed_args)


if __name__ == "__main__":
    main()
