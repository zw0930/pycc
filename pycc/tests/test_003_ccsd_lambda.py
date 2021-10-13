"""
Test CCSD Lambda equation solution using various molecule test cases.
"""

# Import package, test suite, and other packages as needed
import psi4
import pycc
import pytest
from ..data.molecules import *
from .setup_psi4_calculations import setup_h2o_sto3g, setup_h2o_cc_pvdz


def test_lambda_ccsd_h2o(setup_h2o_sto3g, setup_h2o_cc_pvdz):

    rhf_e, rhf_wfn = setup_h2o_sto3g

    maxiter = 75
    e_conv = 1e-12
    r_conv = 1e-12

    ccsd = pycc.ccenergy(rhf_wfn)
    eccsd = ccsd.solve_ccsd(e_conv, r_conv)
    hbar = pycc.cchbar(ccsd)
    cclambda = pycc.cclambda(ccsd, hbar)
    lccsd = cclambda.solve_lambda(e_conv, r_conv)
    epsi4 = -0.070616830152761
    lpsi4 = -0.068826452648939
    assert (abs(epsi4 - eccsd) < 1e-11)
    assert (abs(lpsi4 - lccsd) < 1e-11)

    # cc-pVDZ basis set
    #psi4.set_options({'basis': 'cc-pVDZ'})
    #rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)
    rhf_e, rhf_wfn = setup_h2o_cc_pvdz
    ccsd = pycc.ccenergy(rhf_wfn)
    eccsd = ccsd.solve_ccsd(e_conv, r_conv, maxiter)
    hbar = pycc.cchbar(ccsd)
    cclambda = pycc.cclambda(ccsd, hbar)
    lccsd = cclambda.solve_lambda(e_conv, r_conv)
    epsi4 = -0.222029814166783
    lpsi4 = -0.217838951550509
    assert (abs(epsi4 - eccsd) < 1e-11)
    assert (abs(lpsi4 - lccsd) < 1e-11)
