"""
Psi4 Setup
"""

# Import package, test suite, and other packages as needed
import psi4
import pytest
from ..data.molecules import *


@pytest.fixture()
def setup_h2o_sto3g():
    # Psi4 Setup for h2o in the STO-3G basis set
    psi4.set_memory('2 GB')
    psi4.core.set_output_file('output.dat', False)
    psi4.set_options({'basis': 'STO-3G',
                      'scf_type': 'pk',
                      'mp2_type': 'conv',
                      'freeze_core': 'true',
                      'e_convergence': 1e-12,
                      'd_convergence': 1e-12,
                      'r_convergence': 1e-12,
                      'diis': 1})
    mol = psi4.geometry(moldict["H2O"])
    rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)
    return rhf_e, rhf_wfn

@pytest.fixture()
def setup_h2o_cc_pvdz():
    # Psi4 Setup for h2o in the cc-pVDZ basis set
    psi4.set_memory('2 GB')
    psi4.core.set_output_file('output.dat', False)
    #psi4.set_options({'basis': 'cc-pVDZ'})
    psi4.set_options({'basis': 'cc-pVDZ',
                      'scf_type': 'pk',
                      'mp2_type': 'conv',
                      'freeze_core': 'true',
                      'e_convergence': 1e-12,
                      'd_convergence': 1e-12,
                      'r_convergence': 1e-12,
                      'diis': 1})


    mol = psi4.geometry(moldict["H2O"])
    rhf_e, rhf_wfn = psi4.energy('SCF', return_wfn=True)
    return rhf_e, rhf_wfn
