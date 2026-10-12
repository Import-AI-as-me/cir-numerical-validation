"""Independent numerical and output-protection checks (standard unittest).

The Riccati ODE and Gaussian quadratic-form benchmarks do not reuse the
closed-form CIR bond formula or the discrete affine recursion under test.
"""
import contextlib
import io
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import study


def riccati_rk4(parameters, load, steps=4000):
    """Integrate B'=load-kB-sigma^2 B^2/2, a'=k theta B from zero."""
    kappa = parameters['kappa']
    theta = parameters['theta']
    sigma = parameters['sigma']
    dt = parameters['T'] / steps
    state = np.zeros(2)

    def rhs(y):
        b = y[0]
        return np.array([load-kappa*b-0.5*sigma*sigma*b*b,
                         kappa*theta*b])

    for _ in range(steps):
        k1 = rhs(state)
        k2 = rhs(state + dt*k1/2)
        k3 = rhs(state + dt*k2/2)
        k4 = rhs(state + dt*k3)
        state += dt*(k1+2*k2+2*k3+k4)/6
    b, a = state
    return math.exp(-a-b*parameters['r0']), float(b), float(-a)


def gaussian_square_grid(parameters, steps, dimension, load):
    """Independent CIR benchmark when 4*kappa*theta/sigma^2 is an integer.

    CIR is the sum of squares of independent OU coordinates. The discount
    expectation is then the transform of a Gaussian quadratic form. Matrix
    covariance and a determinant replace the scalar affine recurrence.
    """
    p = parameters
    dt = p['T']/steps
    times = np.arange(1, steps+1, dtype=float)*dt
    weights = np.full(steps, load*dt)
    weights[-1] /= 2
    covariance = p['sigma']**2/(4*p['kappa'])*(
        np.exp(-p['kappa']*np.abs(times[:, None]-times[None, :])/2)
        - np.exp(-p['kappa']*(times[:, None]+times[None, :])/2))
    sqrt_weights = np.sqrt(weights)
    matrix = np.eye(steps)+2*sqrt_weights[:, None]*covariance*sqrt_weights[None, :]
    weighted_mean = sqrt_weights*math.sqrt(p['r0'])*np.exp(-p['kappa']*times/2)
    sign, logdet = np.linalg.slogdet(matrix)
    if sign <= 0:
        raise ArithmeticError('OU covariance transform must be positive definite')
    logprice = (-load*dt*p['r0']/2-dimension*logdet/2
                - weighted_mean @ np.linalg.solve(matrix, weighted_mean))
    return math.exp(float(logprice))


class IndependentNumericalValidation(unittest.TestCase):
    def test_continuous_price_and_coefficients_against_riccati_ode(self):
        for name, p in study.SCENARIOS.items():
            for load in (0.35, 1.0, 2.0):
                with self.subTest(scenario=name, load=load):
                    expected = riccati_rk4(p, load)
                    actual = study.continuous(p, load)
                    np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-11)

    def test_grid_transform_against_gaussian_quadratic_form(self):
        for dimension in (1, 3, 6):
            p = dict(kappa=0.6, sigma=0.3, r0=0.03, T=2.0,
                     theta=dimension*0.3**2/(4*0.6))
            for steps in (1, 5, 17):
                for load in (0.4, 2.0):
                    with self.subTest(dimension=dimension, steps=steps, load=load):
                        expected = gaussian_square_grid(p, steps, dimension, load)
                        actual = study.grid_affine(p, steps, load)[0]
                        self.assertAlmostEqual(actual, expected, delta=2e-13)

    def test_zero_load_is_one_for_any_initial_rate(self):
        for p in study.SCENARIOS.values():
            for function, arguments in ((study.continuous, (p,)),
                                        (study.grid_affine, (p, 13))):
                price, b, loga = function(*arguments, load=0.0)
                self.assertAlmostEqual(price, 1.0, delta=2e-14)
                self.assertAlmostEqual(b, 0.0, delta=2e-14)
                self.assertAlmostEqual(loga, 0.0, delta=2e-14)

    def test_initial_rate_derivatives_by_finite_difference(self):
        for p in study.SCENARIOS.values():
            for function, suffix in ((study.continuous, ()),
                                     (study.grid_affine, (7,))):
                price, b, _ = function(p, *suffix)
                epsilon = 1e-5
                high = function(dict(p, r0=p['r0']+epsilon), *suffix)[0]
                low = function(dict(p, r0=p['r0']-epsilon), *suffix)[0]
                finite_difference = (high-low)/(2*epsilon)
                self.assertAlmostEqual(finite_difference, -b*price, delta=2e-9)

    def test_price_cancellation_does_not_cancel_sensitivity(self):
        # Frozen, independently calculated example from the research report.
        p = dict(study.SCENARIOS['regular'], r0=0.04008534689970658)
        continuous_price, continuous_b, _ = study.continuous(p)
        grid_price, grid_b, _ = study.grid_affine(p, 5)
        self.assertAlmostEqual(continuous_price, grid_price, delta=1e-13)
        self.assertGreater(grid_b/continuous_b-1, 0.03)
        self.assertLess(grid_b/continuous_b-1, 0.031)

    def test_unknown_simulation_scheme_is_rejected(self):
        with self.assertRaises(ValueError):
            study.experiment(study.SCENARIOS['regular'], 'typo', 2, 2, 4, 1)


class OutputProtectionValidation(unittest.TestCase):
    def test_existing_directory_is_rejected_and_untouched(self):
        with tempfile.TemporaryDirectory(prefix='.validation-test-', dir=study.ROOT) as folder:
            path = Path(folder)
            sentinel = path/'results.json'
            sentinel.write_text('preserve archived results', encoding='utf-8')
            with self.assertRaises(ValueError):
                study.prepare_output_dir(path)
            self.assertEqual(sentinel.read_text(encoding='utf-8'), 'preserve archived results')

    def test_repository_root_is_never_an_output_destination(self):
        with self.assertRaises(ValueError):
            study.prepare_output_dir(study.ROOT)

    def test_verify_only_writes_no_result_files(self):
        with tempfile.TemporaryDirectory(prefix='.validation-test-', dir=study.ROOT) as folder:
            destination = Path(folder)/'must-not-be-created'
            # This test checks CLI side effects; mathematical checks are tested above.
            with patch.object(study, 'verify', return_value={'independent_cli_test': True}):
                with contextlib.redirect_stdout(io.StringIO()) as captured:
                    study.main(['--verify-only', '--output-dir', str(destination)])
            self.assertFalse(destination.exists())
            self.assertIn('passed', captured.getvalue())


if __name__ == '__main__':
    unittest.main(verbosity=2)
