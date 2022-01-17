import sys

import numpy as np

from System.MSLogging import log_get_error
from System.MSSystem import LOG_PRIME_NUMBER, TOLERANCE_TYPE


def toolGetWord(input_string, index, d):
    if d is not None and input_string[0] != d:
        input_string = d + input_string
    if d is not None and input_string[-1] != d:
        input_string = input_string + d
    return input_string.split(d)[index + 1]


def toolStr2List(input_string, input_separator):
    if input_string[-1] != input_separator:
        input_string = input_string + input_separator
    return [float(item.strip()) for item in input_string.split(input_separator) if item.strip() != '']

