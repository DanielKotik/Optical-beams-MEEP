import math


def critical(n1, n2):
    """Calculate critical angle in degrees."""
    assert n1 > n2, "\nWarning: Critical angle is not defined, since n1 <= n2!"
    return math.degrees(math.asin(n2/n1))


def brewster(n1, n2):
    """Calculate Brewster angle in degrees."""
    return math.degrees(math.atan(n2/n1))
