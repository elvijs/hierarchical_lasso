""" This module implements the `HierarchicalLasso` class. """
from typing import Optional, Dict, Tuple

import numpy as np
import scipy.optimize
from sklearn.preprocessing import normalize

Array = np.ndarray


class HierarchicalLasso:
    """
    This class implements the constrained lasso described in section 3.1 of
    "Lasso for hierarchical polynomial models", by Hugo Maruri-Aguilar,
    Simon Lunagomez, see https://arxiv.org/abs/2001.07778.

    This is a linear model::
        y = X @ w + e, where @ denotes the matrix multiplication

    TODO: describe how we obtain e

    We fit the parameter w by minimising the objective function::

        0.5 * ||y - Xw||^2_2 + lambda * ||w||_1

    """
    # TODO: classic Lasso in Scikit-learn weighs the L2 term by 1/(2 * n_samples).
    #  Perhaps we should be doing the same.
    _LAMBDA_VALUE_MSG = (
        "We have seen some numerical instability with smaller lambda, please avoid for now. "
        "See the skipped model tests for details."
    )

    def __init__(
            self,
            *,
            lambda_=1.0,
            max_iter=1000,
            optimisation_method: str = "trust-constr",
            optimisation_kwargs: Optional[Dict] = None,
            normalise: bool = False,
    ):
        # TODO document the arguments
        # assert lambda_ >= 0.4, self._LAMBDA_VALUE_MSG

        self._lambda = lambda_
        self._max_iter = max_iter
        self._normalise = normalise
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

        try:
            n_targets = y.shape[1]
            if len(y.shape) > 2:
                raise IndexError
        except IndexError:
            raise ValueError(f"Expected y to have shape (n_samples, n_targets), instead got {y.shape}")

        # Description of the algorithm:
        # 1. We'll normalise the data, so that y_n = y - mean(y) and X_n = (X - mean(X)) / SD(X)
        # TODO: variance or SD?
        # 2. solve y_n = X_n @ w_n and then de-normalise to get
        #    y - mean(y) = (X - mean(X))/SD(X) @ w_n, or
        #    y = X @ (w_n/SD(X)) + mean(y) - mean(X).w_n/SD(X) giving us
        #    w = w_n/SD(X) and e = mean(y) - mean(X).w_n/SD(X)

        X_normalised, y_normalised, X_offset, y_offset, X_scale = self._normalise_data(X, y, scale=self._normalise)

        def objective(w: np.ndarray) -> float:
            # TODO there might be ways of doing the matrix multiplication without the reshape;
            #  this function will be called repeatedly, so perhaps better to create a flatter X outside.
            w_ = w.reshape((n_features, n_targets))
            return 0.5 * np.linalg.norm(y_normalised - X_normalised @ w_, ord=2) + _lambda * np.linalg.norm(w_, ord=1)

        w0 = np.zeros(shape=(n_features * n_targets,))
        # Origin seems like the canonical starting point.
        # Scipy expects a vector instead of a matrix, so that's what we provide.

        # TODO add the Hessian and the Jacobian
        result = scipy.optimize.minimize(
            objective,
            w0,
            method=self._optimisation_method,
            **self._optimisation_kwargs,
        )
        # TODO: check convergence

        normalised_w = result.x.reshape((n_features, n_targets))
        self._w = normalised_w / X_scale.reshape((n_features, 1))
        self._e = y_offset - np.dot(X_offset, self._w)  # Note that self._w is already divided by X_scale

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
        assert (self._w is not None and self._e is not None), "Please call .fit() before attempting to predict."

    @staticmethod
    def _normalise_data(
            X: Array,
            y: Array,
            scale: bool,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Center and optionally scale datd.

        normalised_X = (X - X_offset) / X_scale
        normalised_y = y - y_offset

        X - X_offset and y - y_offset have 0 mean along axis 0.
        X_scale is the L2 norm of X - X_offset.

        This function also systematically makes y consistent with X.dtype.

        :return: normalised_X, normalised_y, X_offset, y_offset, X_scale
        """
        y_copy = y.copy()
        y_copy = np.asarray(y_copy, dtype=X.dtype)
        X_copy = X.copy()

        X_offset = np.average(X_copy, axis=0)
        X_copy -= X_offset
        if scale:
            X_copy, X_scale = normalize(X_copy, axis=0, copy=False, return_norm=True)
        else:
            X_scale = np.ones_like(X_offset)

        y_offset = np.average(y_copy, axis=0)
        y_copy -= y_offset

        return X_copy, y_copy, X_offset, y_offset, X_scale
