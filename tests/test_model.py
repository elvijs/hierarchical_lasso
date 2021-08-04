""" Tests for the hierarchical lasso model. """
from typing import Tuple, Sequence

import numpy as np
import pytest
from sklearn.linear_model import Lasso

from hierarchical_lasso import HierarchicalLasso
from hierarchical_lasso.model import Array, unconstrained_lasso_objective

DECIMALS_TOLERANCE = 3
RELATIVE_TOLERANCE = 1e-6  # TODO tighten the bounds
ABSOLUTE_TOLERANCE = 1e-8

X_and_y = Tuple[np.ndarray, np.ndarray]


@pytest.fixture(scope="function", autouse=True)
def fix_numpy_random_seed() -> None:
    np.random.seed(42)


@pytest.fixture
def degenerate_data_1() -> X_and_y:
    return np.ones(shape=(4, 3)), np.ones(shape=(4, 2))


@pytest.fixture
def degenerate_data_2() -> X_and_y:
    return np.zeros(shape=(10, 3)), np.ones(shape=(10, 2))


@pytest.fixture
def big_degenerate_data() -> X_and_y:
    return np.array([0., 1.]).reshape(-1, 1), np.array([0., 1.]).reshape(-1, 1) * 3


@pytest.fixture
def non_degenerate_data_1() -> X_and_y:
    return np.random.rand(6, 3), np.random.rand(6, 2)


@pytest.fixture
def non_degenerate_data_2() -> X_and_y:
    return np.random.rand(10, 5), np.random.rand(10, 1)


@pytest.fixture
def big_random_data_1() -> X_and_y:
    return np.random.rand(2, 1), np.random.rand(2, 1) * 1_000_000.  # TODO make the shapes a bit bigger


@pytest.fixture
def big_random_data_2() -> X_and_y:
    return np.random.rand(10, 5) + 56., np.random.rand(10, 1) + 78.


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
        np.array([[0.]]),
        np.array([.5]),
    ),
])
def test_analytically_solvable_cases(X: Array, y: Array, lambda_: float, true_w: Array, true_e: Array) -> None:
    # The true values above can be easily verified by pen and paper
    model = HierarchicalLasso(lambda_=lambda_).fit(X, y)
    np.testing.assert_array_almost_equal(model._w, true_w)
    np.testing.assert_array_almost_equal(model._e, true_e)

    n_samples, n_features = X.shape
    n_outputs = y.shape[1]

    # Let's also check the scikit learn implementation agrees
    scikit_model = Lasso(alpha=lambda_).fit(X, y)
    np.testing.assert_allclose(
        scikit_model.coef_.reshape((n_features, n_outputs)),
        model._w,
        atol=ABSOLUTE_TOLERANCE,
    )
    np.testing.assert_allclose(
        scikit_model.intercept_.reshape((n_outputs,)),
        model._e,
        atol=ABSOLUTE_TOLERANCE,
    )


# Note that we cannot generate random data in the pytest decorator
# as we want to always generate the same random numbers.
# TODO uncomment the right things here
# @pytest.mark.parametrize("lambda_", np.linspace(1., 10., 10))  # np.linspace(0., 10., 11))
# @pytest.mark.parametrize("lambda_", np.linspace(0., 10., 11))  # np.linspace(0., 10., 11))
@pytest.mark.parametrize("data_fixture_id", [
    # "degenerate_data_1",
    # "degenerate_data_2",
    # "non_degenerate_data_1",
    # "non_degenerate_data_2",
    "big_degenerate_data",
    # "big_random_data_1",
    # "big_random_data_2",
])
def test_equivalence_with__scikit_learn_lasso(data_fixture_id: str, request) -> None:
    lambda_ = 1.
    X, y = request.getfixturevalue(data_fixture_id)
    n_samples, n_features = X.shape
    n_outputs = y.shape[1]
    another_X, another_y = X.copy(), y.copy()

    scikit_model = Lasso(alpha=lambda_, normalize=True).fit(X, y)
    model = HierarchicalLasso(lambda_=lambda_, normalise=True).fit(another_X, another_y)

    new_X = np.random.rand(*X.shape)
    new_X_copy = new_X.copy()
    scikit_prediction = scikit_model.predict(new_X)
    prediction = model.predict(new_X_copy)

    def loss(coeff: np.ndarray) -> float:
        return np.linalg.norm(y - X @ coeff, ord=2) / (2 * n_samples) + lambda_ * np.linalg.norm(coeff, ord=1)

    # We use slightly different shape conventions to scikit-learn, hence the reshapes
    assert loss(scikit_model.coef_.reshape((n_features, n_outputs))) == loss(model._w)
    np.testing.assert_allclose(
        scikit_model.coef_.reshape((n_features, n_outputs)),
        model._w,
        rtol=RELATIVE_TOLERANCE,
    )
    np.testing.assert_allclose(
        scikit_model.intercept_.reshape((n_outputs,)),
        model._e,
        rtol=RELATIVE_TOLERANCE,
    )
    np.testing.assert_allclose(
        scikit_prediction.reshape(prediction.shape),
        prediction,
        rtol=RELATIVE_TOLERANCE,
    )


@pytest.mark.parametrize(
    "X, y, lambda_, ws, expected_losses",
    [
        # (
        #     np.array([[0.]]).reshape(-1, 1),
        #     np.array([[0.]]).reshape(-1, 1),
        #     1.,
        #     # objective = |w|_1
        #     list(np.array([float(i)]) for i in range(-10, 10, 1)),
        #     [abs(i) for i in range(-10, 10, 1)],
        # ),
        (
            np.array([[1.]]).reshape(-1, 1),
            np.array([[0.]]).reshape(-1, 1),
            1.,
            # objective = |w|^2_2 + |w|_1
            list(np.array([float(i)]) for i in range(-10, 10, 1)),
            [i**2 + abs(i) for i in range(-10, 10, 1)],
        ),
        # TODO: multi output
    ]
)
def test_unconstrained_lasso_objective(
        X: np.ndarray,
        y: np.ndarray,
        lambda_: float,
        ws: Sequence[np.ndarray],
        expected_losses: Sequence[float],
) -> None:
    loss = unconstrained_lasso_objective(X, y, lambda_)
    assert [loss(w) for w in ws] == expected_losses
