""" This module implements the `HierarchicalLasso` class. """
from typing import Optional, Dict

import numpy as np
import scipy.optimize

Array = np.ndarray


class HierarchicalLasso:
    """
    This class implements the constrained lasso described in section 3.1 of
    "Lasso for hierarchical polynomial models", by Hugo Maruri-Aguilar,
    Simon Lunagomez, see https://arxiv.org/abs/2001.07778.

    This is a linear model::
        y = Xw + e

    TODO: describe how we obtain e

    We fit the parameter w by minimising the objective function::

        0.5 * ||y - Xw||^2_2 + lambda * ||w||_1

    """
    # TODO: classic Lasso in Scikit-learn weighs the L2 term by 1/(2 * n_samples).
    #  Perhaps we should be doing the same.

    def __init__(
            self,
            *,
            _lambda=1.0,
            max_iter=1000,
            optimisation_method: str = "trust-constr",
            optimisation_kwargs: Optional[Dict] = None,
    ):
        # TODO document the arguments
        self._lambda = _lambda
        self._max_iter = max_iter
        self._w = None
        self._e = None

        self._optimisation_method = optimisation_method
        self._optimisation_kwargs = optimisation_kwargs or dict()
        # TODO: implement other scikit-learn convenience/speedup methods

    def fit(self, X: Array, y: Array) -> "HierarchicalLasso":
        """Fit Hierarchical Lasso.

        That is, optimise the w and e parameters.

        Parameters
        ----------
        X : {ndarray} of (n_samples, n_features)
            Data.

        y : {ndarray} of shape (n_samples, n_targets)
            Target. Will be cast to X's dtype if necessary.

        TODO: do scikit-learn style input checking
        TODO: admit (n_samples,) in y
        TODO: admit sparse arrays
        """
        n_samples, n_features = X.shape
        _lambda = self._lambda

        n_targets = y.shape[1]

        def objective(w: np.ndarray) -> float:
            # TODO there might be ways of doing the matrix multiplication without the reshape;
            #  this function will be called repeatedly, so perhaps better to create a flatter X outside
            #  One for later once the correctness has been proved.
            w_ = w.reshape((n_features, n_targets))
            return 0.5 * np.linalg.norm(y - X @ w_, ord=2) + _lambda * np.linalg.norm(w_, ord=1)

        # TODO: rescale X and y.
        #  Use the offsets to compute the intercept and fit using the scaled versions
        self._e = 42.

        # TODO add the Hessian and the Jacobian

        w0 = np.zeros(shape=(n_features * n_targets,))  # The objective function
        result = scipy.optimize.minimize(
            objective,
            w0,
            method=self._optimisation_method,
            **self._optimisation_kwargs,
        )
        # TODO: check convergence
        self._w = result.x.reshape((n_features, n_targets))

        # return self for chaining fit and predict calls - this is consistent with scikit-learn
        return self

    def predict(self, X):
        """
        Predict using the linear model Y = Xw + e.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Samples.

        Returns
        -------
        C : array, shape (n_samples, n_targets)
            Returns predicted values.
        """
        self._check_have_been_fit()
        return X @ self._w + self._e

    def _check_have_been_fit(self) -> None:
        assert self._w is not None and self._e is not None
