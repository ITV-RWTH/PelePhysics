#include <turbinflow.H>
#include <mechanism.H>

namespace pele::physics::turbinflow {
void
TurbInflow::init(amrex::Geometry const& /*geom*/)
{
  amrex::ParmParse ppr;

  int n_tp = 0;
  n_tp = ppr.countval("turbinflows");
  amrex::Vector<std::string> tp_list;
  if (n_tp > 0) {
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      AMREX_SPACEDIM == 3, "TurbInflow::init(): TurbInflows are only supported "
                           "for 3 dimensional simulations for now");
    tp.resize(n_tp);
    tp_list.resize(n_tp);
    for (int n = 0; n < n_tp; n++) {
      ppr.get("turbinflows", tp_list[n], n);
    }
  }

  for (int n = 0; n < n_tp; n++) {

    amrex::ParmParse pp("turbinflow." + tp_list[n]);
    if (pp.countval("turb_file") > 0) {

      // Query data
      pp.query("turb_file", tp[n].m_turb_file);
      tp[n].dir = -1;
      pp.query("dir", tp[n].dir);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        tp[n].dir >= 0 && tp[n].dir < AMREX_SPACEDIM,
        "Injection direction is needed: 0, 1 or 2");
      std::string side;
      pp.query("side", side);
      if (side == "low") {
        tp[n].side = amrex::Orientation::low;
      } else if (side == "high") {
        tp[n].side = amrex::Orientation::high;
      } else {
        amrex::Abort("turbinflow.side can only be low or high");
      }
      pp.query("time_offset", tp[n].time_shift);
      pp.query("turb_scale_loc", tp[n].turb_scale_loc);
      pp.query("turb_scale_vel", tp[n].turb_scale_vel);
      pp.query("verbose", tp[n].verbose);
      pp.query("extrap_nonperiodic", tp[n].extrap_nonperiodic);
      pp.query("tile_periodic", tp[n].tile_periodic);
      pp.query("time_periodic", tp[n].time_periodic);
      pp.query("interp_type", tp[n].interp_type);
      if (tp[n].verbose > 0) {
        amrex::Print() << "Initializing turbInflow " << tp_list[n]
                       << " with file " << tp[n].m_turb_file
                       << " (location coordinates in will be scaled by "
                       << tp[n].turb_scale_loc
                       << " and velocity out to be scaled by "
                       << tp[n].turb_scale_vel << ") \n";
      }

      // Get the turbcenter on the injection face
      amrex::Vector<amrex::Real> turb_center(AMREX_SPACEDIM - 1, 0);
      pp.getarr("turb_center", turb_center);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        turb_center.size() == AMREX_SPACEDIM - 1,
        "turb_center must have AMREX_SPACEDIM-1 elements");
      for (amrex::Real& tc : turb_center) {
        tc *= tp[n].turb_scale_loc;
      }

      pp.query("turb_nplane", tp[n].nplane);
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        tp[n].nplane > 3, "need at least 4 turb planes for 3 point "
                          "interpolation stencil + 1 extra");
      pp.query("turb_conv_vel", tp[n].turb_conv_vel);

      // Set other stuff
      std::string turb_header = tp[n].m_turb_file + "/HDR";
      std::ifstream is(turb_header.c_str());
      if (!is.is_open()) {
        amrex::Abort("Unable to open input file " + turb_header);
      }
      amrex::Array<int, AMREX_SPACEDIM> npts = {{0}};
      amrex::Array<amrex::Real, AMREX_SPACEDIM> probsize = {{0}};

      AMREX_D_TERM(is >> npts[0], >> npts[1], >> npts[2]);
      AMREX_D_TERM(is >> probsize[0], >> probsize[1], >> probsize[2]);
      AMREX_D_TERM(
        is >> tp[n].periodicity[0], >> tp[n].periodicity[1],
        >> tp[n].periodicity[2]); // Will use zperiodicity to single whether
                                  // using periodic or time per plane mode

      tp[n].istimeplanes = AMREX_D_PICK(, false, tp[n].periodicity[2] == 0);
      if (tp[n].periodicity[0] == 0 || tp[n].periodicity[1] == 0) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          tp[n].interp_type == TurbInterpType::linear,
          "Linear interpolation required for turbinflow with nonperiodic "
          "directions");
      }

      for (int idim = 0; idim < 2; ++idim) {
        tp[n].dx[idim] = probsize[idim] / amrex::Real(npts[idim] - 1);
        tp[n].dxinv[idim] = 1.0 / tp[n].dx[idim];
      }
      AMREX_D_TERM(, , tp[n].dx[2] = probsize[2] / amrex::Real(npts[2]);)
      AMREX_D_TERM(, , tp[n].dxinv[2] = 1.0 / tp[n].dx[2];)

      // The following is relative to the injection face:
      // 0 and 1 are transverse directions, 2 is normal
      // one ghost point on each side, tangential to inflow face
      AMREX_D_TERM(
        tp[n].pboxsize[0] = probsize[0] - 2.0 * tp[n].dx[0];
        , tp[n].pboxsize[1] = probsize[1] - 2.0 * tp[n].dx[1];
        , tp[n].pboxsize[2] = probsize[2];)

      AMREX_D_TERM(
        tp[n].npboxcells[0] = npts[0] - 3;, tp[n].npboxcells[1] = npts[1] - 3;
        , tp[n].npboxcells[2] = npts[2];)

      // Center the turbulence
      AMREX_D_TERM(
        tp[n].pboxlo[0] = turb_center[0] - 0.5 * tp[n].pboxsize[0];
        , tp[n].pboxlo[1] = turb_center[1] - 0.5 * tp[n].pboxsize[1];
        , tp[n].pboxlo[2] = 0.0;)

      // Swirl type: we can't load more planes than are available
      if (tp[n].istimeplanes) {
        tp[n].nplane = AMREX_D_PICK(
          tp[n].nplane, tp[n].nplane, amrex::min<int>(tp[n].nplane, npts[2]));
      }

      AMREX_D_TERM(, , tp[n].kmax = npts[2];)

      // Optional species header (backward compatible). Immediately after the
      // periodicity line, the turb file may declare species mass-fraction
      // planes as: "SPECIES <n> <name1> ... <nameN>". Legacy (velocity-only)
      // files jump straight to the integer plane offsets, so we peek the next
      // whitespace-delimited token to tell the two apart.
      std::string first_tok;
      is >> first_tok;
      const bool has_species_hdr = (first_tok == "SPECIES");
      if (has_species_hdr) {
        is >> tp[n].n_species;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          tp[n].n_species >= 0, "Negative species count in turb file HDR");
        tp[n].spec_names.resize(tp[n].n_species);
        for (int s = 0; s < tp[n].n_species; ++s) {
          is >> tp[n].spec_names[s];
        }
      }
      tp[n].ncomp = AMREX_SPACEDIM + tp[n].n_species;

      // Resolve file species names to mechanism species indices
      if (tp[n].n_species > 0) {
        amrex::Vector<std::string> mech_names(NUM_SPECIES);
        CKSYMS_STR(mech_names);
        tp[n].spec_idx.resize(tp[n].n_species);
        for (int s = 0; s < tp[n].n_species; ++s) {
          int found = -1;
          for (int m = 0; m < NUM_SPECIES; ++m) {
            if (mech_names[m] == tp[n].spec_names[s]) {
              found = m;
              break;
            }
          }
          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            found >= 0,
            "TurbInflow: a species listed in the turb file was not found in "
            "the active mechanism");
          tp[n].spec_idx[s] = found;
        }
        tp[n].spec_idx_d.resize(tp[n].n_species);
        amrex::Gpu::copy(
          amrex::Gpu::hostToDevice, tp[n].spec_idx.begin(),
          tp[n].spec_idx.end(), tp[n].spec_idx_d.begin());
        if (tp[n].verbose > 0) {
          amrex::Print() << "  turbInflow " << tp_list[n] << " projects "
                         << tp[n].n_species << " species:";
          for (int s = 0; s < tp[n].n_species; ++s) {
            amrex::Print() << " " << tp[n].spec_names[s];
          }
          amrex::Print() << "\n";
        }
      }

      amrex::Box sbx(
        amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
        amrex::IntVect(
          AMREX_D_DECL(npts[0] - 1, npts[1] - 1, tp[n].nplane - 1)));

      tp[n].sdata =
        new amrex::FArrayBox(sbx, tp[n].ncomp, amrex::The_Async_Arena());

      // Offset for each plane in Binary TurbFile. There are ncomp * kmax
      // planes: the AMREX_SPACEDIM velocity components first, then any species.
      tp[n].offset.resize(tp[n].kmax * tp[n].ncomp);
      if (has_species_hdr) {
        for (auto& off : tp[n].offset) {
          is >> off;
        }
      } else {
        // first_tok was actually the first plane offset for a legacy file
        tp[n].offset[0] = std::stol(first_tok);
        for (long i = 1; i < static_cast<long>(tp[n].offset.size()); ++i) {
          is >> tp[n].offset[i];
        }
      }

      if (tp[n].istimeplanes) {
        tp[n].planeTimes.resize(tp[n].kmax);
        for (int i = 0; i < tp[n].kmax; i++) {
          is >> tp[n].planeTimes[i]; // Time for each plane
        }
      }
      is.close();
    }
    turbinflow_initialized = true;
  }
}

void
TurbInflow::add_turb(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const int dir,
  const amrex::Orientation::Side& side,
  const int spec_comp)
{
  AMREX_ALWAYS_ASSERT(turbinflow_initialized);

  // Sentinel marking an inflow cell not covered by a turb patch: species are
  // only projected where the patch provides data, so uncovered cells keep this
  // value and are skipped by set_turb_species (mass fractions are always >= 0).
  constexpr amrex::Real spec_sentinel = -1.0;

  // Box on which we will access data
  amrex::Box bvalsBox = bx;
  int planeLoc =
    (side == amrex::Orientation::low ? geom.Domain().smallEnd()[dir] - 1
                                     : geom.Domain().bigEnd()[dir] + 1);
  bvalsBox.setSmall(dir, planeLoc);
  bvalsBox.setBig(dir, planeLoc);

  // Define box that we will fill with turb: need to be z-normal
  // Get transverse directions
  int tdir1 = (dir != 0) ? 0 : 1;
  int tdir2 = (dir != 0) ? ((dir == 2) ? 1 : 2) : 2;
  int tr1Lo = bvalsBox.smallEnd()[tdir1];
  int tr1Hi = bvalsBox.bigEnd()[tdir1];
  int tr2Lo = bvalsBox.smallEnd()[tdir2];
  int tr2Hi = bvalsBox.bigEnd()[tdir2];
  const amrex::IntVect lo(AMREX_D_DECL(tr1Lo, tr2Lo, planeLoc));
  const amrex::IntVect hi(AMREX_D_DECL(tr1Hi, tr2Hi, planeLoc));
  amrex::Box turbBox(lo, hi);

  // Number of species-mass-fraction planes to project on this face. Species
  // are only handled when the caller supplies a target component (spec_comp)
  // in the state FAB.
  int nSpecFace = 0;
  if (spec_comp >= 0) {
    for (auto& tpn : tp) {
      if (tpn.dir == dir && tpn.side == side) {
        nSpecFace = amrex::max<int>(nSpecFace, tpn.n_species);
      }
    }
  }

  amrex::FArrayBox v(turbBox, 3 + nSpecFace, amrex::The_Async_Arena());
  v.setVal<amrex::RunOn::Device>(0.0, turbBox, 0, 3); // velocity accumulator
  if (nSpecFace > 0) {
    v.setVal<amrex::RunOn::Device>(spec_sentinel, turbBox, 3, nSpecFace);
  }

  // Add turbulence from all the tp acting on this face
  for (auto& tpn : tp) {

    if (tpn.dir == dir && tpn.side == side) {

      // 0 and 1 are the two transverse directions
      amrex::Vector<amrex::Real> x(turbBox.size()[0]), y(turbBox.size()[1]);
      for (int i = turbBox.smallEnd()[0]; i <= turbBox.bigEnd()[0]; ++i) {
        x[i - turbBox.smallEnd()[0]] =
          (geom.ProbLo()[tdir1] + (i + 0.5) * geom.CellSize(tdir1)) *
          tpn.turb_scale_loc;
      }
      for (int j = turbBox.smallEnd()[1]; j <= turbBox.bigEnd()[1]; ++j) {
        y[j - turbBox.smallEnd()[1]] =
          (geom.ProbLo()[tdir2] + (j + 0.5) * geom.CellSize(tdir2)) *
          tpn.turb_scale_loc;
      }

      // Get the turbulence
      amrex::Real z;
      if (tpn.istimeplanes) {
        z = time + tpn.time_shift;
      } else if (tpn.convected_distance >= 0.0) {
        // Through-plane position supplied as a physical convected distance
        // (e.g. the time-integrated inlet velocity); turb_conv_vel is bypassed.
        z = tpn.convected_distance * tpn.turb_scale_loc;
      } else {
        z = (time + tpn.time_shift) * tpn.turb_conv_vel * tpn.turb_scale_loc;
      }
      fill_turb_plane(tpn, x, y, z, v);

      // Project this patch's species onto the state right away (species use
      // SET, not superposition), then reset the species slots so a subsequent
      // overlapping patch starts from the sentinel again.
      if (spec_comp >= 0 && tpn.n_species > 0) {
        set_turb_species(dir, tdir1, tdir2, v, data, spec_comp, tpn);
        v.setVal<amrex::RunOn::Device>(spec_sentinel, turbBox, 3, tpn.n_species);
      }
    }
  }

  // Moving the velocity fluctuations into data
  set_turb(dir, tdir1, tdir2, v, data, dcomp);
}

void
TurbInflow::set_turb_species(
  int normDir,
  int transDir1,
  int transDir2,
  amrex::FArrayBox& v,
  amrex::FArrayBox& data,
  const int spec_comp,
  const TurbParm& a_tp)
{
  const auto& box = v.box(); // z-normal plane
  const auto& v_in = v.array();
  const auto& v_out = data.array();
  const int nspec = a_tp.n_species;
  const int* sidx = a_tp.spec_idx_d.data();

  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // From z-normal box index to data box index
    int idx[3] = {0};
    idx[transDir1] = i;
    idx[transDir2] = j;
    idx[normDir] = k;
    for (int c = 0; c < nspec; ++c) {
      const amrex::Real val = v_in(i, j, k, AMREX_SPACEDIM + c);
      // Skip cells not covered by this patch (sentinel < 0)
      if (val >= 0.0) {
        v_out(idx[0], idx[1], idx[2], spec_comp + sidx[c]) = val;
      }
    }
  });
}

void
TurbInflow::set_turb(
  int normDir,
  int transDir1,
  int transDir2,
  amrex::FArrayBox& v,
  amrex::FArrayBox& data,
  const int dcomp)
{
  // copy velocity fluctuations from plane into data
  const auto& box = v.box(); // z-normal plane
  const auto& v_in = v.array();
  const auto& v_out = data.array(dcomp);

  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // From z-normal box index to data box index
    int idx[3] = {0};
    idx[transDir1] = i;
    idx[transDir2] = j;
    idx[normDir] = k;
    v_out(idx[0], idx[1], idx[2], transDir1) =
      v_in(i, j, k, 0); // transverse velocity 1
    v_out(idx[0], idx[1], idx[2], transDir2) =
      v_in(i, j, k, 1); // transverse velocity 2
    v_out(idx[0], idx[1], idx[2], normDir) =
      v_in(i, j, k, 2); // normal velocity
  });
}

void
TurbInflow::read_one_turb_plane(TurbParm& a_tp, int iplane, int k)
{
  // There are AMREX_SPACEDIM * kmax planes of FABs.
  // The first component are in the first kmax planes,
  // the second component in the next kmax planes, ....
  // Note also that both (*plane) and (*ncomp) start from
  // 1 not 0 since they're passed from Fortran.

  std::string turb_data = a_tp.m_turb_file + "/DAT";
  std::ifstream ifs(turb_data.c_str());
  if (!ifs.is_open()) {
    amrex::Abort("Unable to open input file " + turb_data);
  }

  amrex::Box dstBox = a_tp.sdata->box();
  dstBox.setSmall(AMREX_SPACEDIM - 1, iplane);
  dstBox.setBig(AMREX_SPACEDIM - 1, iplane);

  for (int n = 0; n < a_tp.ncomp; ++n) {

    const long offset_idx = k + (n * a_tp.kmax);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      offset_idx < a_tp.offset.size(), "Bad turb fab offset idx");

    const long start = a_tp.offset[offset_idx];
    ifs.seekg(start, std::ios::beg);

    if (!ifs.good()) {
      amrex::Abort("getplane(): seekg() failed");
    }

    amrex::FArrayBox tmp;
    tmp.readFrom(ifs);
    if (a_tp.verbose > 2) {
      amrex::Print() << "   for d = " << n << " and k = " << k
                     << ": minval = " << tmp.min<amrex::RunOn::Device>(0)
                     << ", maxval = " << tmp.max<amrex::RunOn::Device>(0)
                     << std::endl;
    }
    amrex::Box srcBox = tmp.box();
    a_tp.sdata->copy<amrex::RunOn::Device>(tmp, srcBox, 0, dstBox, n, 1);
  }
  ifs.close();
}

void
TurbInflow::read_turb_planes(TurbParm& a_tp, amrex::Real z)
{
  if (a_tp.istimeplanes) {
    // If time_periodic is enabled, wrap the time value to be within bounds
    if (a_tp.time_periodic) {
      const amrex::Real t_start = a_tp.planeTimes[0];
      const amrex::Real t_end = a_tp.planeTimes[a_tp.kmax - 2];
      const amrex::Real period = t_end - t_start;
      if (period > 0.0) {
        // Wrap z to be within [t_start, t_end)
        amrex::Real z_wrapped = z - t_start;
        z_wrapped = z_wrapped - std::floor(z_wrapped / period) * period;
        z = z_wrapped + t_start;
      }
    }

    if (z < a_tp.planeTimes[0] || z >= a_tp.planeTimes[a_tp.kmax - 2]) {
      amrex::Error(
        "TurbInflow::read_turb_planes(): Requested time (" + std::to_string(z) +
        ") is outside bounds of turbulence data [" +
        std::to_string(a_tp.planeTimes[0]) + ", " +
        std::to_string(a_tp.planeTimes[a_tp.kmax - 2]) +
        ")"); // Need one turbplane forward for interpolation
    }
    for (a_tp.izlo = 0; (a_tp.izlo <= (a_tp.kmax - a_tp.nplane + 1)) &&
                        (a_tp.planeTimes[a_tp.izlo] <= z);
         ++a_tp.izlo) {
    } // Stop when first plane later than time=z
    a_tp.izlo -= 2; // read one extra earlier tplane to prevent rereading
    a_tp.izlo = amrex::max<int>(a_tp.izlo, 0);
    a_tp.izhi = a_tp.izlo + a_tp.nplane - 1;
    a_tp.szlo = a_tp.planeTimes[a_tp.izlo];
    a_tp.szhi = a_tp.planeTimes[a_tp.izhi - 1]; // need one plane forward in
                                                // time for interpolating
  } else {
    a_tp.izlo = (int)(floor(z * a_tp.dxinv[2] - 0.5)) -
                1; // read one extra earlier tplane to prevent rereading
    a_tp.izhi = a_tp.izlo + a_tp.nplane - 1;
    a_tp.szlo = (static_cast<amrex::Real>(a_tp.izlo) + 0.5) * a_tp.dx[2];
    a_tp.szhi = (static_cast<amrex::Real>(a_tp.izhi) - 0.5) *
                a_tp.dx[2]; // need one plane forward in time for interpolating
  }
  if (a_tp.verbose > 1) {
    std::string varname = a_tp.istimeplanes ? "t" : "z";
    amrex::Print() << "read_turb_planes filling " << a_tp.izlo << " to "
                   << a_tp.izhi << std::endl
                   << " --> now have interp data for range [" << a_tp.szlo
                   << ", " << a_tp.szhi << ") with current " << varname << " = "
                   << z << std::endl;
  }

  for (int iplane = 0; iplane < a_tp.nplane; ++iplane) {
    int k = a_tp.izlo + iplane;
    if (!a_tp.istimeplanes) {
      k = (a_tp.npboxcells[2] + k) %
          a_tp.npboxcells[2]; // "wrap" planes if data is periodic
    }
    read_one_turb_plane(a_tp, iplane, k);
  }
}

void
TurbInflow::fill_turb_plane(
  TurbParm& a_tp,
  const amrex::Vector<amrex::Real>& x,
  const amrex::Vector<amrex::Real>& y,
  amrex::Real z,
  amrex::FArrayBox& v)
{
  // If time_periodic is enabled and istimeplanes, wrap the time value
  if (a_tp.istimeplanes && a_tp.time_periodic && a_tp.kmax > 0) {
    const amrex::Real t_start = a_tp.planeTimes[0];
    const amrex::Real t_end = a_tp.planeTimes[a_tp.kmax - 2];
    const amrex::Real period = t_end - t_start;
    if (period > 0.0) {
      // Wrap z to be within [t_start, t_end)
      amrex::Real z_wrapped = z - t_start;
      z_wrapped = z_wrapped - std::floor(z_wrapped / period) * period;
      z = z_wrapped + t_start;
    }
  }

  const amrex::Real tplanes_lo = a_tp.szlo;
  const amrex::Real tplanes_hi = a_tp.szhi;

  if ((z < tplanes_lo) || (z >= tplanes_hi)) {
    if (a_tp.verbose > 1) {
      std::string varname = a_tp.istimeplanes ? "t" : "z";
      amrex::Print() << "Reading new data because " << varname << " = " << z
                     << " is outside the range [" << tplanes_lo << ", "
                     << tplanes_hi << ")" << std::endl;
    }
    read_turb_planes(a_tp, z);
  }

  const auto& bx = v.box();
  const auto& vd = v.array();

  amrex::Gpu::DeviceVector<amrex::Real> x_dev(x.size());
  amrex::Gpu::DeviceVector<amrex::Real> y_dev(y.size());
  amrex::Gpu::copyAsync(
    amrex::Gpu::hostToDevice, x.begin(), x.end(), x_dev.begin());
  amrex::Gpu::copyAsync(
    amrex::Gpu::hostToDevice, y.begin(), y.end(), y_dev.begin());
  amrex::Real* xd = x_dev.data();
  amrex::Real* yd = y_dev.data();

  amrex::Real velScale = (a_tp.side == amrex::Orientation::high)
                           ? -a_tp.turb_scale_vel
                           : a_tp.turb_scale_vel;
  const auto& npboxcells = a_tp.npboxcells;
  const auto& pboxlo = a_tp.pboxlo;
  const auto& pboxsize = a_tp.pboxsize;
  const auto& szlo = a_tp.szlo;
  const auto& dxinv = a_tp.dxinv;
  const auto& dx = a_tp.dx;
  const auto& sd = a_tp.sdata->array();
  const auto& ext_nonper = a_tp.extrap_nonperiodic;
  const auto& tile_per = a_tp.tile_periodic;
  const auto& periodicity = a_tp.periodicity;
  const bool lininterp = a_tp.interp_type == TurbInterpType::linear;
  // Number of components to project into v. sdata always holds a_tp.ncomp
  // (velocity + species), but v may be velocity-only (e.g. a velocity-only
  // add_turb call, or a caller that did not request species), so never write
  // past v's component count.
  const int lncomp = amrex::min<int>(a_tp.ncomp, v.nComp());
  amrex::Real cz[3];
  int k0 = -1;
  if (a_tp.istimeplanes) {
    AMREX_ALWAYS_ASSERT(
      z >= a_tp.planeTimes[a_tp.izlo] && z <= a_tp.planeTimes[a_tp.izhi]);
    for (k0 = 1; k0 < a_tp.nplane - 2 && a_tp.planeTimes[a_tp.izlo + k0] <= z;
         ++k0) {
    } // Stop when first plane later than time=z, then go back one
    k0 -= 1;
    const auto& t0 = a_tp.planeTimes[a_tp.izlo + k0];
    const auto& t1 = a_tp.planeTimes[a_tp.izlo + k0 + 1];
    const auto& t2 = a_tp.planeTimes[a_tp.izlo + k0 + 2];
    AMREX_ALWAYS_ASSERT(z >= t0 && z <= t2);
    cz[0] = (z - t1) * (z - t2) / ((t0 - t1) * (t0 - t2));
    cz[1] = (z - t0) * (z - t2) / ((t1 - t0) * (t1 - t2));
    cz[2] = (z - t0) * (z - t1) / ((t2 - t0) * (t2 - t1));
  } else {
    amrex::Real zz =
      (z - szlo) * dxinv[2];    // How many dz away from the left side ?
    k0 = (int)(std::floor(zz)); // What's the closest point ?
    zz -= amrex::Real(k0);
    cz[0] =
      lininterp ? 1.0 - zz : 0.5 * (zz - 1.0) * (zz - 2.0); // Weight of k0
    cz[1] = lininterp ? zz : zz * (2.0 - zz);               // Weight of k0 + 1
    cz[2] = lininterp ? 0.0 : 0.5 * zz * (zz - 1.0);        // Weight of k0 + 2
  }

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real cx[3], cy[3], ydata[3];
    amrex::Real zdata[3][3];

    const amrex::Real x_from_pboxlo = xd[i - bx.smallEnd(0)] - pboxlo[0];
    amrex::Real xx =
      (x_from_pboxlo)*dxinv[0] +
      0.5; // from 1st cell center (ghost cell so 1/2 dx outside of pbox)
    const amrex::Real y_from_pboxlo = yd[j - bx.smallEnd(1)] - pboxlo[1];
    amrex::Real yy =
      (y_from_pboxlo)*dxinv[1] +
      0.5; // from 1st cell center (ghost cell so 1/2 dx outside of pbox)
    int i0 = (int)(std::floor(xx));
    int j0 = (int)(std::floor(yy));
    xx -= amrex::Real(i0);
    yy -= amrex::Real(j0);
    // Wrap if needed
    if (tile_per) {
      i0 = (i0 % npboxcells[0] + npboxcells[0]) % npboxcells[0];
      j0 = (j0 % npboxcells[1] + npboxcells[1]) % npboxcells[1];
    }
    cx[0] = lininterp ? 1.0 - xx : 0.5 * (xx - 1.0) * (xx - 2.0);
    cy[0] = lininterp ? 1.0 - yy : 0.5 * (yy - 1.0) * (yy - 2.0);
    cx[1] = lininterp ? xx : xx * (2.0 - xx);
    cy[1] = lininterp ? yy : yy * (2.0 - yy);
    cx[2] = lininterp ? 0.0 : 0.5 * xx * (xx - 1.0);
    cy[2] = lininterp ? 0.0 : 0.5 * yy * (yy - 1.0);

    if (
      (x_from_pboxlo >= 0.0 && x_from_pboxlo < pboxsize[0]) ||
      (tile_per && periodicity[0] != 0)) {
      if (
        (x_from_pboxlo < 0.5 * dx[0] || x_from_pboxlo > pboxsize[0] - dx[0]) &&
        (periodicity[0] == 0) && !ext_nonper) {
        amrex::Error(
          "TurbInflow interp stencil touches ghost cell for nonperiodic "
          "direction and extrap_nonperiodic option not used");
      }
      if (
        (y_from_pboxlo >= 0.0 && y_from_pboxlo < pboxsize[1]) ||
        (tile_per && periodicity[1] != 0)) {
        if (
          (y_from_pboxlo < 0.5 * dx[1] ||
           y_from_pboxlo > pboxsize[1] - dx[1]) &&
          (periodicity[1] == 0) && !ext_nonper) {
          amrex::Error(
            "TurbInflow interp stencil touches ghost cell for nonperiodic "
            "direction and extrap_nonperiodic option not used");
        }

        for (int n = 0; n < lncomp; ++n) {
          for (int ii = 0; ii <= 2; ++ii) {
            for (int jj = 0; jj <= 2; ++jj) {
              zdata[ii][jj] = cz[0] * sd(i0 + ii, j0 + jj, k0, n) +
                              cz[1] * sd(i0 + ii, j0 + jj, k0 + 1, n) +
                              cz[2] * sd(i0 + ii, j0 + jj, k0 + 2, n);
            }
          }
          for (int ii = 0; ii <= 2; ++ii) {
            ydata[ii] = cy[0] * zdata[ii][0] + cy[1] * zdata[ii][1] +
                        cy[2] * zdata[ii][2];
          }
          const amrex::Real interp =
            cx[0] * ydata[0] + cx[1] * ydata[1] + cx[2] * ydata[2];
          if (n < AMREX_SPACEDIM) {
            // Velocity fluctuation: scaled, sign-flipped on the high side,
            // and superimposed across overlapping patches.
            vd(i, j, k, n) += velScale * interp;
          } else {
            // Species mass fraction: passive scalar, projected as-is (SET).
            vd(i, j, k, n) = interp;
          }
        }
      }
    }
  });
  amrex::Gpu::synchronize(); // Ensure that DeviceVector's don't leave scope
                             // early
}
} // namespace pele::physics::turbinflow
