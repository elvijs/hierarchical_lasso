""" This module implements the `HierarchicalLasso` class. """
import numpy as np

Array = np.ndarray


class HierarchicalLasso:
    """
    This class implements the constrained lasso described in section 3.1 of
    "Lasso for hierarchical polynomial models", by Hugo Maruri-Aguilar,
    Simon Lunagomez, see https://arxiv.org/abs/2001.07778.

    This is a linear model::
        y = Xw + e

    We fit the parameters w and e by minimising the objective function::

        0.5 * ||y - Xw||^2_2 + lambda * ||w||_1

    """
    # TODO: classic Lasso in Scikit-learn weighs the L2 term by 1/(2 * n_samples).
    #  Perhaps we should be doing the same.

    def __init__(
            self,
            *,
            _lambda=1.0,
            max_iter=1000,
    ):
        self._lambda = _lambda
        self._max_iter = max_iter
        self._w = None
        self._e = None
        # TODO: implement other scikit-learn convenience/speedup methods

    def fit(self, X: Array, y: Array) -> None:
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

        # TODO: implement

        # return self for chaining fit and predict calls - this is consistent with scikit-learn
        return self
