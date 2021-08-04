"""
Simple experiment with Scipy optimisation.
We want to see if constraints that use absolute values are solved correctly.
"""
from typing import Sequence

import numpy as np
from scipy.optimize import NonlinearConstraint, minimize


def loss(x: np.ndarray) -> np.ndarray:
    """ Reasonably representative 1D problem with minimum at 0.5. """
    return (x[0] - 1.)**2 + np.abs(x[0])


def constraint(x: np.ndarray) -> Sequence[np.ndarray]:
    return [np.abs(x[0])]


if __name__ == '__main__':
    nonlinear_constraint_1 = NonlinearConstraint(constraint, -np.inf, 1.)
    x0 = np.array([-0.2])  # Start on the wrong side of zero
    res = minimize(loss, x0, method='trust-constr', constraints=[nonlinear_constraint_1], options={'verbose': 1})
    print(res)

    nonlinear_constraint_2 = NonlinearConstraint(constraint, -np.inf, 0.25)
    x0 = np.array([-0.2])  # Start on the wrong side of zero
    res = minimize(loss, x0, method='trust-constr', constraints=[nonlinear_constraint_2], options={'verbose': 1})
    print(res)
