from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import TestCase, run_module_suite, \
                          assert_array_almost_equal, assert_almost_equal, \
                          assert_allclose

from scipy.signal import disc2continuous as d2c
from scipy.signal import dlsim, ss2tf, ss2zpk, lsim2

# Author: Jens Wurm <wurm.jens@gmail.com>
# October 01, 2015


class TestD2C(TestCase):
    def test_zoh(self):
        ac = np.eye(2)
        bc = 0.5 * np.ones((2, 1))
        cc = np.array([[0.75, 1.0], [1.0, 1.0], [1.0, 0.25]])
        dc = np.array([[0.0], [0.0], [-0.33]])

        ad_truth = 1.648721270700128 * np.eye(2)
        bd_truth = 0.324360635350064 * np.ones((2, 1))
        # c and d in discrete should be equal to their continuous counterparts
        dt_requested = 0.5

        ad, bd, cd, dd, dt = d2c((ac, bc, cc, dc), dt_requested, method='zoh')

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cc, cd)
        assert_array_almost_equal(dc, dd)
        assert_almost_equal(dt_requested, dt)

    def test_gbt(self):
        ac = np.eye(2)
        bc = 0.5 * np.ones((2, 1))
        cc = np.array([[0.75, 1.0], [1.0, 1.0], [1.0, 0.25]])
        dc = np.array([[0.0], [0.0], [-0.33]])

        dt_requested = 0.5
        alpha = 1.0 / 3.0

        ad_truth = 1.6 * np.eye(2)
        bd_truth = 0.3 * np.ones((2, 1))
        cd_truth = np.array([[0.9, 1.2],
                             [1.2, 1.2],
                             [1.2, 0.3]])
        dd_truth = np.array([[0.175],
                             [0.2],
                             [-0.205]])

        ad, bd, cd, dd, dt = d2c((ac, bc, cc, dc), dt_requested,
                                 method='gbt', alpha=alpha)

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cd_truth, cd)
        assert_array_almost_equal(dd_truth, dd)

    def test_euler(self):
        ac = np.eye(2)
        bc = 0.5 * np.ones((2, 1))
        cc = np.array([[0.75, 1.0], [1.0, 1.0], [1.0, 0.25]])
        dc = np.array([[0.0], [0.0], [-0.33]])

        dt_requested = 0.5

        ad_truth = 1.5 * np.eye(2)
        bd_truth = 0.25 * np.ones((2, 1))
        cd_truth = np.array([[0.75, 1.0],
                             [1.0, 1.0],
                             [1.0, 0.25]])
        dd_truth = dc

        ad, bd, cd, dd, dt = d2c((ac, bc, cc, dc), dt_requested,
                                 method='euler')

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cd_truth, cd)
        assert_array_almost_equal(dd_truth, dd)
        assert_almost_equal(dt_requested, dt)

    def test_backward_diff(self):
        ac = np.eye(2)
        bc = 0.5 * np.ones((2, 1))
        cc = np.array([[0.75, 1.0], [1.0, 1.0], [1.0, 0.25]])
        dc = np.array([[0.0], [0.0], [-0.33]])

        dt_requested = 0.5

        ad_truth = 2.0 * np.eye(2)
        bd_truth = 0.5 * np.ones((2, 1))
        cd_truth = np.array([[1.5, 2.0],
                             [2.0, 2.0],
                             [2.0, 0.5]])
        dd_truth = np.array([[0.875],
                             [1.0],
                             [0.295]])

        ad, bd, cd, dd, dt = d2c((ac, bc, cc, dc), dt_requested,
                                 method='backward_diff')

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cd_truth, cd)
        assert_array_almost_equal(dd_truth, dd)

    def test_bilinear(self):
        ac = np.eye(2)
        bc = 0.5 * np.ones((2, 1))
        cc = np.array([[0.75, 1.0], [1.0, 1.0], [1.0, 0.25]])
        dc = np.array([[0.0], [0.0], [-0.33]])

        dt_requested = 0.5

        ad_truth = (5.0 / 3.0) * np.eye(2)
        bd_truth = (1.0 / 3.0) * np.ones((2, 1))
        cd_truth = np.array([[1.0, 4.0 / 3.0],
                             [4.0 / 3.0, 4.0 / 3.0],
                             [4.0 / 3.0, 1.0 / 3.0]])
        dd_truth = np.array([[0.291666666666667],
                             [1.0 / 3.0],
                             [-0.121666666666667]])

        ad, bd, cd, dd, dt = d2c((ac, bc, cc, dc), dt_requested,
                                 method='bilinear')

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cd_truth, cd)
        assert_array_almost_equal(dd_truth, dd)
        assert_almost_equal(dt_requested, dt)

        # Same continuous system again, but change sampling rate

        ad_truth = 1.4 * np.eye(2)
        bd_truth = 0.2 * np.ones((2, 1))
        cd_truth = np.array([[0.9, 1.2], [1.2, 1.2], [1.2, 0.3]])
        dd_truth = np.array([[0.175], [0.2], [-0.205]])

        dt_requested = 1.0 / 3.0

        ad, bd, cd, dd, dt = c2d((ac, bc, cc, dc), dt_requested,
                                 method='bilinear')

        assert_array_almost_equal(ad_truth, ad)
        assert_array_almost_equal(bd_truth, bd)
        assert_array_almost_equal(cd_truth, cd)
        assert_array_almost_equal(dd_truth, dd)
        assert_almost_equal(dt_requested, dt)

    def test_transferfunction(self):
        numc = np.array([0.25, 0.25, 0.5])
        denc = np.array([0.75, 0.75, 1.0])

        numd = np.array([[1.0 / 3.0, -0.427419169438754, 0.221654141101125]])
        dend = np.array([1.0, -1.351394049721225, 0.606530659712634])

        dt_requested = 0.5

        num, den, dt = d2c((numc, denc), dt_requested, method='zoh')

        assert_array_almost_equal(numd, num)
        assert_array_almost_equal(dend, den)
        assert_almost_equal(dt_requested, dt)

    def test_zerospolesgain(self):
        zeros_c = np.array([0.5, -0.5])
        poles_c = np.array([1.j / np.sqrt(2), -1.j / np.sqrt(2)])
        k_c = 1.0

        zeros_d = [1.23371727305860, 0.735356894461267]
        polls_d = [0.938148335039729 + 0.346233593780536j,
                   0.938148335039729 - 0.346233593780536j]
        k_d = 1.0

        dt_requested = 0.5

        zeros, poles, k, dt = d2c((zeros_c, poles_c, k_c), dt_requested,
                                  method='zoh')

        assert_array_almost_equal(zeros_d, zeros)
        assert_array_almost_equal(polls_d, poles)
        assert_almost_equal(k_d, k)
        assert_almost_equal(dt_requested, dt)

    def test_gbt_with_sio_tf_and_zpk(self):
        """Test method='gbt' with alpha=0.25 for tf and zpk cases."""
        # State space coefficients for the continuous SIO system.
        A = -1.0
        B = 1.0
        C = 1.0
        D = 0.5

        # The continuous transfer function coefficients.
        cnum, cden = ss2tf(A, B, C, D)

        # Continuous zpk representation
        cz, cp, ck = ss2zpk(A, B, C, D)

        h = 1.0
        alpha = 0.25

        # Explicit formulas, in the scalar case.
        Ad = (1 + (1 - alpha) * h * A) / (1 - alpha * h * A)
        Bd = h * B / (1 - alpha * h * A)
        Cd = C / (1 - alpha * h * A)
        Dd = D + alpha * C * Bd

        # Convert the explicit solution to tf
        dnum, dden = ss2tf(Ad, Bd, Cd, Dd)

        # Compute the discrete tf using cont2discrete.
        c2dnum, c2dden, dt = d2c((cnum, cden), h, method='gbt', alpha=alpha)

        assert_allclose(dnum, c2dnum)
        assert_allclose(dden, c2dden)

        # Convert explicit solution to zpk.
        dz, dp, dk = ss2zpk(Ad, Bd, Cd, Dd)

        # Compute the discrete zpk using cont2discrete.
        c2dz, c2dp, c2dk, dt = d2c((cz, cp, ck), h, method='gbt', alpha=alpha)

        assert_allclose(dz, c2dz)
        assert_allclose(dp, c2dp)
        assert_allclose(dk, c2dk)

    def test_discrete_approx(self):
        """
        Test that the solution to the discrete approximation of a continuous
        system actually approximates the solution to the continuous system.
        This is an indirect test of the correctness of the implementation
        of cont2discrete.
        """

        def u(t):
            return np.sin(2.5 * t)

        a = np.array([[-0.01]])
        b = np.array([[1.0]])
        c = np.array([[1.0]])
        d = np.array([[0.2]])
        x0 = 1.0

        t = np.linspace(0, 10.0, 101)
        dt = t[1] - t[0]
        u1 = u(t)

        # Use lsim2 to compute the solution to the continuous system.
        t, yout, xout = lsim2((a, b, c, d), T=t, U=u1, X0=x0,
                              rtol=1e-9, atol=1e-11)

        # Convert the continuous system to a discrete approximation.
        dsys = d2c((a, b, c, d), dt, method='bilinear')

        # Use dlsim with the pairwise averaged input to compute the output
        # of the discrete system.
        u2 = 0.5 * (u1[:-1] + u1[1:])
        t2 = t[:-1]
        td2, yd2, xd2 = dlsim(dsys, u=u2.reshape(-1, 1), t=t2, x0=x0)

        # ymid is the average of consecutive terms of the "exact" output
        # computed by lsim2.  This is what the discrete approximation
        # actually approximates.
        ymid = 0.5 * (yout[:-1] + yout[1:])

        assert_allclose(yd2.ravel(), ymid, rtol=1e-4)


if __name__ == "__main__":
    run_module_suite()