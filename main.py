#!/usr/bin/env python3
# -*- coding:utf-8 -*-
u"""
Created at 2020.09.04 by Zhang Yiming
"""
import warnings
warnings.filterwarnings("ignore")

import click

from cli.climb import climb
from cli.diff import diff

@click.group()
def main():
    pass

main.add_command(climb)
main.add_command(diff)


if __name__ == '__main__':
    main()
