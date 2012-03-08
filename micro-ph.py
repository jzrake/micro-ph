#!/usr/bin/env python

from eospy import evaluate_term

if __name__ == "__main__":
    print evaluate_term("number_density", +1, 1, 1)
    print evaluate_term("number_density", -1, 1, 1, massless=False)
