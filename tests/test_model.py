""" Tests for the hierarchical lasso model. """
import numpy as np
import pytest
from sklearn.linear_model import Lasso

from hierarchical_lasso import HierarchicalLasso
from hierarchical_lasso.model import Array


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


def test_y_vector__throws_value_error() -> None:
    X = np.ones(shape=(2, 1))
    y = np.ones(shape=(2,))  # this is not allowed - we always want y to be a matrix

    with pytest.raises(ValueError):
        HierarchicalLasso().fit(X, y)


@pytest.mark.parametrize("X, y, lambda_, true_w, true_e", [
    # TODO: add a single sample case?
    (
        np.array([0., 1.]).reshape(2, 1),
        np.array([0., 1.]).reshape(2, 1),
        0.,
        np.array([[1.]]),
        np.array([0.]),
    ),
    (
        np.array([1., 2.]).reshape(2, 1),
        np.array([0., 1.]).reshape(2, 1),
        0.,
        np.array([[1.]]),
        np.array([-1.]),
    ),
    (
        np.array([0., 1.]).reshape(2, 1),
        np.array([0., 1.]).reshape(2, 1),
        5.,
        np.array([[3.]]),
        np.array([-1.]),
    ),
])
def test_analytically_solvable_cases(X: Array, y: Array, lambda_: float, true_w: Array, true_e: Array) -> None:
    # The true values above can be easily verified by pen and paper
    model = HierarchicalLasso(lambda_=lambda_).fit(X, y)
    np.testing.assert_array_almost_equal(model._w, true_w)
    np.testing.assert_array_almost_equal(model._e, true_e)


@pytest.mark.parametrize("lambda_", np.linspace(0., 5., 10))
@pytest.mark.parametrize("X, y", [
    # Degenerate data
    (np.ones(shape=(4, 3)), np.ones(shape=(4, 2))),
    (np.zeros(shape=(10, 3)), np.ones(shape=(10, 2))),
    # Random data with a few different shapes
    (np.random.rand(4, 3), np.random.rand(4, 2)),
    (np.random.rand(10, 5), np.random.rand(10, 1)),
    (np.random.rand(4, 3) * 12., np.random.rand(4, 2) * 34.),  # check we can do with bigger numbers
    (np.random.rand(10, 5) + 56., np.random.rand(10, 1) + 78.),
])
def test_equivalence_with__scikit_learn_lasso(lambda_: float, X: Array, y: Array) -> None:
    scikit_model = Lasso(alpha=lambda_).fit(X, y)
    model = HierarchicalLasso(lambda_=lambda_).fit(X, y)

    new_X = np.random.rand(*X.shape)
    scikit_prediction = scikit_model.predict(new_X)
    prediction = model.predict(new_X)

    # We use slightly different shape conventions to scikit-learn,
    # so flattening arrays for comparison
    np.testing.assert_array_equal(scikit_model.coef_.T.flatten(), model._w.flatten())
    np.testing.assert_array_equal(scikit_model.intercept_.flatten(), model._e.flatten())
    np.testing.assert_array_equal(scikit_prediction.flatten(), prediction.flatten())


# Fixtures

@pytest.fixture(scope="module", autouse=True)
def fix_numpy_random_seed() -> None:
    np.random.seed(42)
