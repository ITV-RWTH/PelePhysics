#ifndef TURBINFLOW_H
#define TURBINFLOW_H

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

namespace pele::physics::turbinflow {

AMREX_ENUM(TurbInterpType, linear, quadratic);

struct TurbParm
{
  // Turbulence data file
  std::string m_turb_file;

  // Verbosity
  int verbose = 1;

  // Geometry of the turbulent data
  amrex::GpuArray<int, 3> npboxcells = {{0}};
  amrex::GpuArray<amrex::Real, 3> pboxlo = {{0.0}};
  amrex::GpuArray<amrex::Real, 3> dx = {{0.0}};
  amrex::GpuArray<amrex::Real, 3> dxinv = {{0.0}};
  amrex::GpuArray<amrex::Real, 3> pboxsize = {{0.0}};
  int dir;                       // Direction where to use this TurbParm
  amrex::Orientation::Side side; // Side on which to use this TurbParm
  amrex::Real time_shift =
    0.0; // Enable to offset the location in the turb file

  int nplane = 32; // Number of turb planes stored at once in memory
  amrex::FArrayBox* sdata = nullptr; // Currently loaded chunk of turb data
  int izlo = -1;                     // Currently loaded low plane
  int izhi = -1;                     // Currently loaded high plane
  amrex::Real szlo = -1.0e12;        // Position of the current chunk low plane
  amrex::Real szhi = -1.0e11;        // Position of the current chunk high plane

  amrex::IntVect periodicity{AMREX_D_DECL(1, 1, 1)};
  bool istimeplanes = false; // false indicates turb data is periodic normal to
                             // flow dir, true requires time array in turb data

  bool extrap_nonperiodic =
    false; // if true, allow interpolation where stencil touches ghost cells for
           // nonperiodic directions
  bool tile_periodic = false; // if true, turbulence is patched over entire face
                              // with a periodic tiled pattern
  bool time_periodic = false; // if true, turbulence data is looped in time
  TurbInterpType interp_type =
    TurbInterpType::quadratic; // linear required for non-periodic directions

  amrex::Real turb_scale_loc =
    1.; // Spatial scale relating turbfile geom -> actual case scale
  amrex::Real turb_scale_vel = 1.; // Velocity scaling
  amrex::Real turb_conv_vel =
    1.; // Mean inlet velocity, used to move through data planes

  amrex::Real convected_distance =
    -1.; // If >= 0, used directly as the through-plane sampling position
         // (physical convected distance, e.g. the time integral of the inlet
         // velocity), bypassing turb_conv_vel * time. Negative => disabled.

  amrex::Vector<amrex::Real> planeTimes;

  amrex::Vector<long> offset;
  int kmax; // Number of plane in Turbfile

  // Optional species mass-fraction planes carried by the turb file, in
  // addition to the 3 velocity components. These are stored in the file
  // after the velocity components (component-major order). Species are
  // projected onto the inflow just like velocity, but without the velocity
  // scaling / high-side sign flip (they are passive scalars).
  int n_species = 0; // Number of species planes stored in the turb file
  int ncomp =
    AMREX_SPACEDIM; // Total number of stored components (velocity + species)
  amrex::Vector<std::string> spec_names; // Species names, in file order
  amrex::Vector<int>
    spec_idx; // Mechanism index for each file species (host, file order)
  amrex::Gpu::DeviceVector<int> spec_idx_d; // ... same mapping on device
};

struct TurbInflow
{
public:
  TurbInflow() = default;

  ~TurbInflow() = default;

  void init(amrex::Geometry const& geom);

  void add_turb(
    amrex::Box const& bx,
    amrex::FArrayBox& data,
    const int dcomp,
    amrex::Geometry const& geom,
    const amrex::Real time,
    const int dir,
    const amrex::Orientation::Side& side,
    const int spec_comp = -1);

  bool is_initialized() const { return turbinflow_initialized; }

  // Total number of species-mass-fraction planes carried by any turb patch.
  // Zero for legacy (velocity-only) turb files.
  int max_turb_species() const
  {
    int n = 0;
    for (const auto& t : tp) {
      n = amrex::max<int>(n, t.n_species);
    }
    return n;
  }

  // For each mechanism species index in [0, num_species), returns 1 if any
  // turb patch provides that species (and hence seeds it into the state at
  // the inflow), 0 otherwise. Used by the flow solver to know which species
  // to pick up from the seeded inflow data.
  amrex::Vector<int> turb_species_flags(int num_species) const
  {
    amrex::Vector<int> flags(num_species, 0);
    for (const auto& t : tp) {
      for (int idx : t.spec_idx) {
        if (idx >= 0 && idx < num_species) {
          flags[idx] = 1;
        }
      }
    }
    return flags;
  }

  // Override the through-plane sampling position with a supplied physical
  // convected distance (e.g. the time-integrated inlet velocity). Applies to
  // all loaded turb planes. Pass a negative value to revert to the default
  // turb_conv_vel * time behavior.
  void set_convected_distance(amrex::Real dist)
  {
    for (auto& t : tp) {
      t.convected_distance = dist;
    }
  }

  static void read_turb_planes(TurbParm& a_tp, amrex::Real z);

  static void read_one_turb_plane(TurbParm& a_tp, int iplane, int k);

  static void fill_turb_plane(
    TurbParm& a_tp,
    const amrex::Vector<amrex::Real>& x,
    const amrex::Vector<amrex::Real>& y,
    amrex::Real z,
    amrex::FArrayBox& v);

  static void set_turb(
    int normDir,
    int transDir1,
    int transDir2,
    amrex::FArrayBox& v,
    amrex::FArrayBox& data,
    const int dcomp);

  // Copy the species mass-fraction planes accumulated in v (stored at
  // components [AMREX_SPACEDIM, AMREX_SPACEDIM + a_tp.n_species)) into the
  // state FAB, at component spec_comp + mechanism_index for each file species.
  static void set_turb_species(
    int normDir,
    int transDir1,
    int transDir2,
    amrex::FArrayBox& v,
    amrex::FArrayBox& data,
    const int spec_comp,
    const TurbParm& a_tp);

private:
  amrex::Vector<TurbParm> tp;
  bool turbinflow_initialized = false;
};
} // namespace pele::physics::turbinflow

#endif
