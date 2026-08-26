"""Solve a 1D premixed flame in Cantera and save in Pele-readable format."""

###############################################################
# ADIABATIC_FLAME - A freely-propagating, premixed flat flame
#                   of fuel in AIR (O2 + N2)
###############################################################

# import :

import argparse
import csv
import os
import re

import numpy as np
import yaml
from cantera import Solution, FreeFlame

#################################################################
# Parse arguments
#################################################################
parser = argparse.ArgumentParser(
    prog="Cantera PMF Generator",
    description="Use Cantera to solve a 1D premixed flame and save in a Pele-readable format",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    fromfile_prefix_chars='@',
)

parser.add_argument(
    "-m",
    "--mechanism",
    default="drm19",
    help=(
        "Name of PelePhysics mechanism from Mechanisms, i.e. the directory "
        "holding it relative to Mechanisms. Sub-mechanisms of a collection are "
        "given as a path (C3MechLite_v401/H2-NH3_25sp_noHeAr), or by their "
        "innermost directory alone when that name is unambiguous"
    ),
)
parser.add_argument(
    "-pp", "--pp_home", default="../../../", help="Path to PelePhysics directory"
)
parser.add_argument(
    "-f", "--fuel", default="CH4:1", help="Fuel stream mole-basis Cantera composition"
)
parser.add_argument(
    "-ox",
    "--oxidizer",
    default="O2:1, N2:3.76",
    help="Oxidizer stream mole-basis Cantera composition",
)
parser.add_argument(
    "-Y",
    "--massfrac",
    default=None,
    help="String of mass fractions (overrides -f and -ox options)",
)
parser.add_argument(
    "-T",
    "--temperature",
    default=300,
    type=float,
    help="Unburned mixture temperature [K]",
)
parser.add_argument(
    "-p", "--pressure", default=101325, type=float, help="Pressure [Pa]"
)
parser.add_argument("-phi", "--phi", default=1.0, type=float, help="Equivalence Ratio")
parser.add_argument(
    "-tr",
    "--transport",
    default="mixture-averaged",
    choices=[
        "mixture-averaged",
        "multicomponent",
        "unity-Lewis-number",
        "mixture-averaged-CK",
        "multicomponent-CK",
        # legacy Cantera (<3.0) spellings, kept for convenience
        "Mix",
        "Multi",
        "UnityLewis",
        "CK_Mix",
        "CK_Multi",
    ],
    help="Cantera transport model used for the flame solution",
)
parser.add_argument(
    "-s",
    "--soret",
    action="store_true",
    help=(
        "Include thermal diffusion (Soret effect). Support depends on the "
        "transport model and the Cantera version in use (older Cantera releases "
        "only evaluate thermal diffusion coefficients for the multicomponent "
        "models); Cantera raises an error if the combination is unsupported"
    ),
)
parser.add_argument("-d", "--domain", default=0.04, type=float, help="Domain width [m]")
parser.add_argument("-v", "--verbose", default=1, type=int, help="Verbosity level")
parser.add_argument(
    "-o",
    "--output",
    default=None,
    help=(
        "Path to directory where flames will be saved, if not specified files are saved to PelePhysics/Mechanisms/<mechanism>/PMFs/"
    ),
)
args = parser.parse_args()

#################################################################
# Prepare your run
#################################################################
# Parameter values :

# Mixture
mechanism = args.mechanism
fuel_species = args.fuel
ox_species = args.oxidizer
Y_species = args.massfrac

# General
p = args.pressure  # pressure [Pa]
tin = args.temperature  # unburned gas temperature [K]
phi = args.phi  # Eq. ratio [-]

# Cantera transport model. The legacy spellings are deprecated in Cantera >= 3.0,
# so map them onto the current names before handing them to Cantera.
transport_aliases = {
    "Mix": "mixture-averaged",
    "Multi": "multicomponent",
    "UnityLewis": "unity-Lewis-number",
    "CK_Mix": "mixture-averaged-CK",
    "CK_Multi": "multicomponent-CK",
}
transport = transport_aliases.get(args.transport, args.transport)
multicomponent_models = ("multicomponent", "multicomponent-CK")
soret = args.soret

# Recent Cantera versions also evaluate thermal diffusion coefficients outside of
# the multicomponent models, so no combination is screened out here: if the
# installed Cantera cannot enable Soret for the requested transport model, it
# raises the error itself when f.soret_enabled is set below.

# Multicomponent transport is expensive and less robust from a cold start, so the
# refinement ladder below is always run mixture-averaged and the final solution is
# continued onto the requested model.
ladder_transport = (
    "mixture-averaged" if transport in multicomponent_models else transport
)

# Refined grid at inlet and outlet, 6 points in x-direction :
domain_size = args.domain  # Domain size [m]
initial_grid = (
    domain_size * np.array([0.0, 0.001, 0.01, 0.02, 0.029, 0.03], "d") / 0.03
)  # m

# Set tolerance properties
tol_ss = [1.0e-8, 1.0e-9]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-8, 1.0e-9]  # [rtol atol] for time stepping
loglevel = args.verbose  # amount of diagnostic output
refine_grid = True  # True to enable refinement

#################
# Find mechanism in PelePhysics
pp_path = os.path.join(args.pp_home, "Mechanisms")
if not (os.path.exists(pp_path)):
    raise RuntimeError("Invalid path to PelePhysics: " + args.pp_home)

mech_paths = [
    re.sub(r'\s--plog=\S+', '', name).strip()
    for name in open(os.path.join(pp_path, "list_mech")).readlines()
    if not name.startswith("#")
]
# A mechanism is named by the directory holding it, relative to Mechanisms.
# Collections such as C3MechLite keep their sub-mechanisms one level deeper, so
# that name may itself be a path, e.g. C3MechLite_v401/H2-NH3_25sp_noHeAr.
mech_names = [os.path.dirname(name) for name in mech_paths]
mech_paths = dict(zip(mech_names, mech_paths))

# QSS: we will solve with skeletal mechanism, then eliminate QSS species
qss_data = [
    name.split()
    for name in open(os.path.join(pp_path, "list_qss_mech")).readlines()
    if not name.startswith("#")
]
qss_names = [os.path.dirname(name[0]) for name in qss_data]
qss_paths = dict(zip(qss_names, [name[2] for name in qss_data]))
qss_nonqss = dict(zip(qss_names, [name[3] for name in qss_data]))

# The innermost directory alone is accepted as a shorthand for a nested
# mechanism whenever it is unambiguous, so that -m H2-NH3_25sp_noHeAr works as
# well as the full -m C3MechLite_v401/H2-NH3_25sp_noHeAr.
all_names = mech_names + qss_names
basenames = [os.path.basename(name) for name in all_names]
shorthands = {
    base: name
    for base, name in zip(basenames, all_names)
    if basenames.count(base) == 1 and base not in all_names
}
mechanism = shorthands.get(mechanism, mechanism)

if mechanism in mech_names:
    mech_has_qssa = False
    mech_path = mech_paths[mechanism]
elif mechanism in qss_names:
    mech_has_qssa = True
    mech_path = qss_paths[mechanism]
    with open(os.path.join(pp_path, qss_nonqss[mechanism])) as f:
        nonqss_spec_list = yaml.safe_load(f)["species"]
else:
    raise RuntimeError(
        "Requested mechanism ("
        + mechanism
        + ") found in neither list_mech or list_qssa_mech"
    )
mech_path = os.path.join(pp_path, mech_path)

#################
# Print information
# A nested mechanism name carries a directory prefix that must not leak into the
# file names, so the innermost directory is used to tag the solutions.
label_pre = "pmf-" + os.path.basename(mechanism)
if Y_species is not None:
    num_input_species = len([item for item in Y_species.split(',') if ':' in item])
    label = "Y" + str(num_input_species) + "_T" + str(tin) + "_P" + str(p)
else:
    label = fuel_species.split(":")[0] + "_PHI" + str(phi) + "_T" + str(tin) + "_P" + str(p)

# Tag the file name with the transport model so that solutions obtained with
# different models do not overwrite each other. The default (mixture-averaged,
# no Soret) keeps the historical file names untouched.
transport_tags = {
    "mixture-averaged": "Mix",
    "multicomponent": "Multi",
    "unity-Lewis-number": "UnityLe",
    "mixture-averaged-CK": "MixCK",
    "multicomponent-CK": "MultiCK",
}
if transport != "mixture-averaged" or soret:
    label += "_" + transport_tags[transport] + ("Soret" if soret else "")

#################
# Directory the solutions are written to
if args.output is None:
    outdir = os.path.join(pp_path, mechanism, "PMFs")
else:
    outdir = args.output
if not os.path.exists(outdir):
    os.makedirs(outdir)

#################
# Create and Run Flame:

# Set gas state to that of the unburned gas
gas = Solution(mech_path, "gas")
if Y_species is None:
    gas.TP = tin, p
    gas.set_equivalence_ratio(phi, fuel_species, ox_species, basis="mole")
else:
    gas.TPY = tin, p, Y_species
    print("\nMass fractions read into Cantera:")
    for spec, massfrac in zip(gas.species_names, gas.Y):
        print(f"  {spec}: {massfrac}")

# Create the free laminar premixed flame
f = FreeFlame(gas, initial_grid)

f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

f.transport_model = ladder_transport

# No energy for starters
f.energy_enabled = False

# Refinement criteria
f.set_refine_criteria(ratio=7.0, slope=1, curve=1)

# Max number of times the Jacobian will be used before it must be re-evaluated
f.set_max_jac_age(10, 10)

# Set time steps whenever Newton convergence fails
f.set_time_step(1.0e-06, [1, 2, 5, 10])  # s
f.max_time_step_count = 3000
f.max_grid_points = 1000

# Calculation
f.solve(loglevel, refine_grid)

#################
# Second flame:

# Energy equation enabled
f.energy_enabled = True

# Refinement criteria when energy equation is enabled
f.set_refine_criteria(ratio=5.0, slope=0.5, curve=0.5)

# Calculation
f.solve(loglevel, refine_grid)

#################
# Third flame and so on ...:
f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.1)

f.solve(loglevel, refine_grid)

##################
# Fourth flame and so on ...:
f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05, prune=0.01)

f.solve(loglevel, refine_grid)

##################
# Fifth flame and so on ...
f.set_refine_criteria(ratio=2.0, slope=0.02, curve=0.02, prune=0.01)

f.solve(loglevel, refine_grid)

##################
# Continuation onto the requested transport model, then onto Soret:

if transport != ladder_transport:
    f.transport_model = transport
    f.solve(loglevel, refine_grid)

if soret:
    f.soret_enabled = True
    f.solve(loglevel, refine_grid)

print(
    transport
    + (" + Soret" if soret else "")
    + " flamespeed = ",
    f.velocity[0],
)

#################################################################
# Save your results
#################################################################

# Always save with Mole Fractions of all species and CGS units
nz = f.flame.n_points
csv_file = os.path.join(outdir, str(label_pre + "-" + label + "-X.dat"))
with open(csv_file, "w") as outfile:
    writer = csv.writer(
        outfile, delimiter=" ", quotechar=" ", quoting=csv.QUOTE_MINIMAL
    )
    writer.writerow(
        ['VARIABLES = "X" "temp" "u" "rho"']
        + ['"X_' + x + '"' for x in gas.species_names]
    )
    writer.writerow(["ZONE I=" + str(nz) + " FORMAT=POINT" + " SPECFORMAT=MOLE"])
    for kk in range(len(f.grid)):
        f.set_gas_state(kk)
        writer.writerow(
            [f.grid[kk] * 1e2, gas.T, f.velocity[kk] * 1e2, gas.density * 1e-3]
            + list(gas.X)
        )

# QSS Mechanisms: save again with QSS species removed and others renormalized to unity sum
if mech_has_qssa:
    csv_file = os.path.join(
        outdir, str(label_pre + "-" + label + "-X-qssa-removed.dat")
    )
    with open(csv_file, "w") as outfile:
        writer = csv.writer(
            outfile, delimiter=" ", quotechar=" ", quoting=csv.QUOTE_MINIMAL
        )
        writer.writerow(
            ['VARIABLES = "X" "temp" "u" "rho"']
            + ['"X_' + x + '"' for x in gas.species_names if x in nonqss_spec_list]
        )
        writer.writerow(["ZONE I=" + str(nz) + " FORMAT=POINT" + " SPECFORMAT=MOLE"])
        max_total_qss = 0.0
        for kk in range(len(f.grid)):
            f.set_gas_state(kk)
            zeroqss_mole_fracs = np.array(
                [
                    X if (gas.species_name(ii) in nonqss_spec_list) else 0.0
                    for ii, X in enumerate(gas.X)
                ]
            )
            total_qss = 1.0 - np.sum(zeroqss_mole_fracs)
            max_total_qss = max(max_total_qss, total_qss)
            # this will modify rho so that it is consistent with the right T/p with the modified mole fractions
            gas.TPX = gas.T, gas.P, zeroqss_mole_fracs
            writer.writerow(
                [f.grid[kk] * 1e2, gas.T, f.velocity[kk] * 1e2, gas.density * 1e-3]
                + [
                    X
                    for ii, X in enumerate(gas.X)
                    if (gas.species_name(ii) in nonqss_spec_list)
                ]
            )
    print("Maximum local total mole fraction of QSS species: ", max_total_qss)
