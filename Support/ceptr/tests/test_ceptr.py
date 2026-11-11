"""Tests for cepter."""

import math
import os
import pathlib
import random
import shutil
import subprocess
import sys

import cantera as ct

import ceptr.converter as converter
import ceptr.qssa as qssa
import ceptr.utilities as cu
from ceptr import __version__


def mechanism_path(mname):
    """Determine mechanism path."""
    this_file_dir = pathlib.Path(__file__).parent.resolve()
    model_path = "Mechanisms"
    return this_file_dir.parents[2] / model_path / mname


def test_version():
    """Test the version."""
    assert __version__ == "0.1.0"


def test_air():
    """Test mechanism generation of air."""
    mech_path = mechanism_path("air")
    fname = mech_path / "mechanism.yaml"
    mechanism = ct.Solution(fname)
    interface = None
    chemistry = "homogeneous"
    conv = converter.Converter(mechanism, interface, chemistry)
    conv.writer()
    conv.formatter()


def test_lidryer():
    """Test mechanism generation of LiDryer."""
    mech_path = mechanism_path("LiDryer")
    fname = mech_path / "mechanism.yaml"
    mechanism = ct.Solution(fname)
    interface = None
    chemistry = "homogeneous"
    conv = converter.Converter(mechanism, interface, chemistry)
    conv.writer()
    conv.formatter()


def test_dodecane_lu():
    """Test mechanism generation of dodecane_lu."""
    mech_path = mechanism_path("dodecane_lu")
    fname = mech_path / "mechanism.yaml"
    mechanism = ct.Solution(fname)
    interface = None
    chemistry = "homogeneous"
    conv = converter.Converter(mechanism, interface, chemistry)
    conv.writer()
    conv.formatter()


def test_ch4_lean_qss(tmp_path):
    """Test QSS reduction and mechanism generation of CH4_lean_qss.

    Runs the same two stages as ``qssa -lq`` and ``convert -lq`` (QSS
    reduction of the skeletal mechanism, then conversion with the symbolic
    QSS Jacobian) on a copy of the mechanism inputs, so the committed files
    are left untouched.
    """
    mech_path = mechanism_path("CH4_lean_qss")
    for name in ["skeletal.yaml", "non_qssa_list.yaml", "qssa_input.toml"]:
        shutil.copy(mech_path / name, tmp_path)
    qssa.process_qss(
        tmp_path / "skeletal.yaml", tmp_path / "non_qssa_list.yaml", False, 2
    )
    mechanism = ct.Solution(tmp_path / "qssa.yaml")
    conv = converter.Converter(
        mechanism,
        None,
        "homogeneous",
        jacobian=True,
        qss_format_input=tmp_path / "qssa_input.toml",
        qss_symbolic_jacobian=True,
    )
    conv.writer()
    conv.formatter()
    header = (tmp_path / "mechanism.H").read_text()
    assert "aJacobian" in header
    assert (tmp_path / "mechanism.cpp").exists()


def convert_with_hash_seed(fname, out_dir, seed):
    """Convert a copy of a mechanism in a fresh interpreter with a hash seed."""
    out_dir.mkdir()
    shutil.copy(fname, out_dir / "mechanism.yaml")
    script = (
        "import sys, cantera as ct, ceptr.converter as converter;"
        "conv = converter.Converter(ct.Solution(sys.argv[1]), None, 'homogeneous');"
        "conv.writer();"
        "conv.formatter()"
    )
    env = dict(os.environ, PYTHONHASHSEED=str(seed))
    subprocess.run(
        [sys.executable, "-c", script, str(out_dir / "mechanism.yaml")],
        env=env,
        check=True,
        capture_output=True,
    )
    return [(out_dir / f).read_text() for f in ["mechanism.H", "mechanism.cpp"]]


def test_hash_seed_independence(tmp_path):
    """Check that the generated code does not depend on PYTHONHASHSEED.

    posf10325 has reactions with a species on both sides with different
    coefficients (e.g. CH3 + POSF10325 => ... + 0.7 CH3 + ...), which used to
    give the production-rate writer tied sort keys, so the order of the wdot
    updates followed the hash-randomised set order.
    """
    fname = mechanism_path("posf10325") / "mechanism.yaml"
    reference = convert_with_hash_seed(fname, tmp_path / "seed0", 0)
    for seed in [1, 2, 3]:
        assert (
            convert_with_hash_seed(fname, tmp_path / f"seed{seed}", seed) == reference
        ), f"output with PYTHONHASHSEED={seed} differs from PYTHONHASHSEED=0"


def test_kc_exp():
    """Check the product form of the equilibrium constants against exp."""
    mech_path = mechanism_path("dodecane_lu")
    mechanism = ct.Solution(mech_path / "mechanism.yaml")
    conv = converter.Converter(mechanism, None, "homogeneous")
    species_info = conv.species_info
    rng = random.Random(0)
    g_rt = [rng.uniform(-50.0, 50.0) for _ in range(species_info.n_species)]
    env = {
        "g_RT": g_rt,
        "ge": [math.exp(-g) for g in g_rt],
        "inv_ge": [math.exp(g) for g in g_rt],
        "exp": math.exp,
    }
    for reaction in mechanism.reactions():
        arg = eval(cu.sorted_kc_exp_arg(mechanism, species_info, reaction), env)
        for sign in (-1, 1):
            product, exp_form = cu.kc_exp_parts(mechanism, species_info, reaction, sign)
            assert product is not None and "exp" not in product
            assert math.isclose(eval(exp_form, env), math.exp(sign * arg))
            assert math.isclose(eval(product, env), math.exp(sign * arg), rel_tol=1e-12)


def test_kc_product_no_overflow():
    """Check that Kc products cannot overflow where they are used.

    Regression test: in LiDryer at 79 K, H + H2O2 <=> H2O + OH has a finite
    Kc (about 1e190), but a naive factor order overflows on its first two
    factors. The product form must be off there (exp fallback), and on the
    range where it is on, every left-to-right partial product must be finite
    and the result must match exp of the summed exponent.
    """
    mech_path = mechanism_path("LiDryer")
    mechanism = ct.Solution(mech_path / "mechanism.yaml")
    conv = converter.Converter(mechanism, None, "homogeneous")
    species_info = conv.species_info
    t_low, t_high, _ = cu.kc_product_plan(mechanism, species_info)
    assert t_low > 79.0
    for temp in [t_low, 250.0, 300.0, 1000.0, 3000.0, t_high]:
        g_rt = {
            s: float(cu.gibbs_rt(mechanism, s, [temp])[0])
            for s in species_info.nonqssa_species_list
        }
        g_list = [g_rt[s] for s in species_info.nonqssa_species_list]
        env = {
            "g_RT": g_list,
            "ge": [math.exp(-g) for g in g_list],
            "inv_ge": [math.exp(g) for g in g_list],
            "exp": math.exp,
        }
        for reaction in mechanism.reactions():
            arg = eval(cu.sorted_kc_exp_arg(mechanism, species_info, reaction), env)
            for sign in (-1, 1):
                product, _ = cu.kc_exp_parts(mechanism, species_info, reaction, sign)
                assert product is not None
                value = 1.0
                for factor in product.split(" * "):
                    value *= eval(factor, env)
                    assert math.isfinite(value) and value > 0.0
                assert math.isclose(value, math.exp(sign * arg), rel_tol=1e-10)
