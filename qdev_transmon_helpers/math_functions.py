import numpy as np
from math import sqrt

# TODO: docstrings


def qubit_from_push(g, push, bare_res):
    """
    Get estimated qubit location given coupling, push on resonator
    and bare resonator position.

    Args:
        g: coupling in Hz
        push: push in Hz
        bare_res: resonator position in Hz

    Returns:
        estimated qubit frequency
    """
    delta = g**2 / push
    return bare_res - delta


def g_from_qubit(qubit, bare_res, push):
    """
    Get estimated coupling strength given qubit frequency, push
    on resonator and bare resonator position.

    Args:
        g: coupling in Hz
        push: push in Hz
        bare_res: resonator position in Hz

    Returns:
        estimated qubit frequency
    """
    delta = bare_res - qubit
    return sqrt(delta * push)


def exp_decay_sin(x, a, b, c, d, e):
    return a * np.exp(-x / b) * np.sin(c * x + d) + e


def exp_decay(x, a, b, c):
    return a * np.exp(-x / b) + c
