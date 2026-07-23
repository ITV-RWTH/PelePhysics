#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include <turbinflow.H>
#include <mechanism.H>

// Standalone regression test for the species-projection capability of the
// TurbInflow utility. It reads a turbulence file that carries species mass
// fractions (produced by the TurbInflowGenerator with the `species` option),
// projects it onto a low-Z inflow boundary through TurbInflow::add_turb, and
// verifies that:
//   * the requested species are recognized and reported,
//   * their mass fractions are written into the correct state components,
//   * the projected values match the file (bounded, spatially varying, and
//     O2 + N2 == 1 as encoded by the generator input), and
//   * velocity fluctuations are still applied.
//
// It does not depend on a flow solver, so it exercises the utility directly.

using namespace amrex;

int
main(int argc, char* argv[])
{
  Initialize(argc, argv);
  int status = 0;
  {
    // Domain matching the turbulence file scaling used in the input
    // (0.01 m cube, 32^3, periodic in x,y and inflow on low-Z).
    const int ncell = 32;
    Box domain(
      IntVect(AMREX_D_DECL(0, 0, 0)),
      IntVect(AMREX_D_DECL(ncell - 1, ncell - 1, ncell - 1)));
    RealBox rb(
      {AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(0.01, 0.01, 0.01)});
    Array<int, AMREX_SPACEDIM> is_per = {AMREX_D_DECL(1, 1, 0)};
    Geometry geom(domain, rb, 0, is_per);

    // Initialize the TurbInflow utility from the input file.
    pele::physics::turbinflow::TurbInflow turb_inflow;
    turb_inflow.init(geom);
    AMREX_ALWAYS_ASSERT(turb_inflow.is_initialized());

    const int nspec_file = turb_inflow.max_turb_species();
    Print() << "TurbInflow reports " << nspec_file
            << " projected species in the file\n";
    if (nspec_file <= 0) {
      Print() << "FAIL: expected the turb file to carry species\n";
      Finalize();
      return 1;
    }

    // Which mechanism species are provided by the file?
    auto flags = turb_inflow.turb_species_flags(NUM_SPECIES);
    AMREX_ALWAYS_ASSERT(flags[O2_ID] == 1 && flags[N2_ID] == 1);

    // Build a state FAB covering the low-Z ghost plane, laid out as
    // [u, v, w, <species...>]; here spec_comp = AMREX_SPACEDIM.
    const int spec_comp = AMREX_SPACEDIM;
    const int ncomp = AMREX_SPACEDIM + NUM_SPECIES;
    const int dir = 2;
    const auto side = Orientation::low;

    Box gbx = domain;
    gbx.growLo(dir, 2); // include the low-Z ghost cells
    FArrayBox data(gbx, ncomp, The_Async_Arena());

    // Single ghost plane just outside the low-Z face.
    Box bndryBox = adjCellLo(domain, dir, 1);

    // Mimic PeleLM::fillTurbInflow: zero velocity, sentinel species.
    data.setVal<RunOn::Device>(0.0, gbx, 0, AMREX_SPACEDIM);
    data.setVal<RunOn::Device>(-1.0, gbx, spec_comp, NUM_SPECIES);

    turb_inflow.add_turb(bndryBox, data, 0, geom, 0.0, dir, side, spec_comp);
    Gpu::streamSynchronize();

    // Inspect the projected plane.
    const auto& a = data.const_array();
    Real o2min = 1.0e30, o2max = -1.0e30;
    Real max_sum_err = 0.0, max_abs_vel = 0.0;
    int ncovered = 0, nsentinel = 0;
    const int k = bndryBox.smallEnd(dir);
    for (int j = bndryBox.smallEnd(1); j <= bndryBox.bigEnd(1); ++j) {
      for (int i = bndryBox.smallEnd(0); i <= bndryBox.bigEnd(0); ++i) {
        const Real yO2 = a(i, j, k, spec_comp + O2_ID);
        const Real yN2 = a(i, j, k, spec_comp + N2_ID);
        if (yO2 < 0.0) {
          nsentinel++;
          continue; // cell not covered by the patch
        }
        ncovered++;
        o2min = std::min(o2min, yO2);
        o2max = std::max(o2max, yO2);
        max_sum_err = std::max(max_sum_err, std::abs(yO2 + yN2 - 1.0));
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          max_abs_vel = std::max(max_abs_vel, std::abs(a(i, j, k, d)));
        }
      }
    }

    Print() << "Covered cells: " << ncovered << ", sentinel cells: "
            << nsentinel << "\n";
    Print() << "Projected Y(O2): min = " << o2min << ", max = " << o2max
            << "\n";
    Print() << "max |Y(O2)+Y(N2) - 1| = " << max_sum_err << "\n";
    Print() << "max |velocity fluctuation| = " << max_abs_vel << "\n";

    auto check = [&](bool ok, const char* msg) {
      Print() << (ok ? "  PASS: " : "  FAIL: ") << msg << "\n";
      if (!ok) {
        status = 1;
      }
    };
    check(ncovered > 0, "at least one inflow cell received projected species");
    check(
      o2min >= 0.05 - 1.0e-6 && o2max <= 0.40 + 1.0e-6,
      "projected Y(O2) stays within the file's [0.05, 0.40] range");
    check(
      o2max - o2min > 0.05,
      "projected Y(O2) varies in space (matches the file pattern, not a "
      "uniform default)");
    check(max_sum_err < 1.0e-6, "projected Y(O2) + Y(N2) == 1");
    check(max_abs_vel > 0.0, "velocity fluctuations were also projected");

    if (status == 0) {
      Print() << "TurbInflow species projection test PASSED\n";
    } else {
      Print() << "TurbInflow species projection test FAILED\n";
    }
  }
  Finalize();
  return status;
}
