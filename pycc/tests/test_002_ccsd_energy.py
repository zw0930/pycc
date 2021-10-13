"""
Test CCSD equation solution using various molecule test cases.
"""

# Import package, test suite, and other packages as needed
import psi4
import pycc
from .setup_psi4_calculations import setup_h2o_sto3g

def test_ccsd_h2o(setup_h2o_sto3g):

    rhf_e, rhf_wfn = setup_h2o_sto3g

    maxiter = 75
    e_conv = 1e-12
    r_conv = 1e-12
    ccsd = pycc.ccenergy(rhf_wfn)
    eccsd = ccsd.solve_ccsd(e_conv, r_conv, maxiter)
    epsi4 = -0.070616830152761
    assert (abs(epsi4 - eccsd) < 1e-11)

    # cc-pVDZ basis set
    psi4.set_options({'basis': 'cc-pVDZ'})
    rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)
    ccsd = pycc.ccenergy(rhf_wfn)
    eccsd = ccsd.solve_ccsd(e_conv, r_conv, maxiter)
    epsi4 = -0.222029814166783
    assert (abs(epsi4 - eccsd) < 1e-11)
