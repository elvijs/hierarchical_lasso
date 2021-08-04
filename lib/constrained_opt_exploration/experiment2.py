"""
Example 5 from Hugo's paper.

Just checking the interfaces and making sure we can optimise a nontrivial example.
"""
from typing import Sequence

import numpy as np
from scipy.optimize import NonlinearConstraint, minimize


A = np.array(
    [
        [1., 0., 0., -1., 0.],
        [1., 0., 0., 0., -1.],
        [0., 1., 0., -1., 0.],
        [0., 0., 1., 0., -1.],
    ]
)
Y = np.array([1., 2., 3., 4., 5., 6., 7., 8.]).reshape(8, 1)
X = np.arange(1., 5 * 8 + 1, 1.).reshape((8, 5))
lambda_ = 1
theta0 = np.zeros(shape=(5,))


def loss(theta: np.ndarray) -> np.ndarray:
    # print(f"Y={Y}\nX={X}\ntheta={theta}")
    return 0.5*np.linalg.norm(Y - X@theta, ord=2) + lambda_*np.linalg.norm(theta, ord=1)


def constraint(theta: np.ndarray) -> Sequence[float]:
    return [A[i, :] @ np.abs(theta) for i in range(4)]


if __name__ == '__main__':
    nonlinear_constraint = NonlinearConstraint(constraint, 0., np.inf)
    res = minimize(loss, theta0, method='trust-constr', constraints=[nonlinear_constraint], options={'verbose': 1})
    print(res)

    res = minimize(loss, theta0, method='SLSQP', constraints=[nonlinear_constraint], options={'verbose': 1})
    print(res)
