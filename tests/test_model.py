""" Tests for the hierarchical lasso model. """
import numpy as np
from sklearn.linear_model import Lasso

from hierarchical_lasso import HierarchicalLasso


def test_can_construct() -> None:
    HierarchicalLasso()


def test_can_fit() -> None:
    n_features, n_samples, n_outputs = 4, 3, 2
    X = np.ones(shape=(n_samples, n_features))
    y = np.ones(shape=(n_samples, n_outputs))

    HierarchicalLasso().fit(X, y)


def test_can_predict() -> None:
    n_features, n_samples, n_outputs = 4, 3, 2
    X = np.ones(shape=(n_samples, n_features))
    y = np.ones(shape=(n_samples, n_outputs))

    model = HierarchicalLasso().fit(X, y)
    model.predict(X)


def test_equivalence_with__scikit_learn_lasso() -> None:
    n_features, n_samples, n_outputs = 4, 3, 2
    X = np.ones(shape=(n_samples, n_features))
    y = np.ones(shape=(n_samples, n_outputs))

    scikit_model = Lasso().fit(X, y)
    model = HierarchicalLasso().fit(X, y)

    # TODO compare the intercept and w directly
    scikit_prediction = scikit_model.predict(X)
    prediction = model.predict(X)

    np.testing.assert_array_equal(scikit_prediction, prediction)
