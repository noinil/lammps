/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "bond_harmonic.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
/* ------------------------------------------------ //
//      _                         _            _    //
//   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
//  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
// | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
//  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
//                       |_|                        //
// ------------------ tc_test --------------------- */
#include <fstream>
#include <iomanip>
// ===================================================

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondHarmonic::BondHarmonic(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondHarmonic::~BondHarmonic()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
  }
}

/* ---------------------------------------------------------------------- */

void BondHarmonic::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr,rk;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  /* ------------------------------------------------ //
  //      _                         _            _    //
  //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
  //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
  // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
  //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
  //                       |_|                        //
  // ------------------ tc_test --------------------- */
  double ftan[nlocal][3];
  double etan = 0;
  double etan_total = 0;
  std::ofstream forces_file("n_force_bd_harmonic.dat");
  std::ofstream energy_file("p_energy_bd_harmonic.dat");
  forces_file << " stacking forces: " << std::endl;
  energy_file << " stacking energy: " << std::endl;
  energy_file << std::setw(6) << "bd2_i"
              << std::setw(6) << "i1"
              << std::setw(6) << "i2"
              << std::setw(11) << "E_bd2"
              << std::endl;
  energy_file << " ---------------------------------------------"
              << std::endl;
  // ===================================================

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    dr = r - r0[type];
    rk = k[type] * dr;

    // force & energy

    if (r > 0.0) fbond = -2.0*rk/r;
    else fbond = 0.0;

    if (eflag) ebond = rk*dr;

    /* ------------------------------------------------ //
    //      _                         _            _    //
    //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
    //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
    // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
    //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
    //                       |_|                        //
    // ------------------ tc_test --------------------- */
    etan = rk*dr;
    etan_total += etan;
    energy_file << std::setw(6) << n + 1
                << std::setw(6) << atom->tag[i1]
                << std::setw(6) << atom->tag[i2] << " "
                << std::setw(10) << etan
                << std::endl;
    // ===================================================

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
        /* ------------------------------------------------ //
        //      _                         _            _    //
        //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
        //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
        // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
        //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
        //                       |_|                        //
        // ------------------ tc_test --------------------- */
        ftan[i1][0] += delx*fbond;
        ftan[i1][1] += dely*fbond;
        ftan[i1][2] += delz*fbond;

        f[i1][0] += delx*fbond;
        f[i1][1] += dely*fbond;
        f[i1][2] += delz*fbond;
    }

    if (newton_bond || i2 < nlocal) {
        /* ------------------------------------------------ //
        //      _                         _            _    //
        //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
        //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
        // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
        //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
        //                       |_|                        //
        // ------------------ tc_test --------------------- */
        ftan[i2][0] -= delx*fbond;
        ftan[i2][1] -= dely*fbond;
        ftan[i2][2] -= delz*fbond;

        f[i2][0] -= delx*fbond;
        f[i2][1] -= dely*fbond;
        f[i2][2] -= delz*fbond;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
  /* ------------------------------------------------ //
  //      _                         _            _    //
  //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
  //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
  // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
  //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
  //                       |_|                        //
  // ------------------ tc_test --------------------- */
  energy_file << "Total bond_harmonic energy: " << etan_total << std::endl;
  energy_file << " ================================================== "
              << std::endl;
  forces_file << std::setw(6) << "i"
              << std::setw(10) << "fx"
              << std::setw(12) << "fy"
              << std::setw(12) << "fz"
              << std::setw(12) << "|f|"
              << std::endl;
  forces_file << " -----------------------------------------------------"
              << std::endl;
  for (int mm = 1; mm < nlocal + 1; mm++) {
      int nn = atom->map(mm);
      double ftc0 = ftan[nn][0] * ftan[nn][0] > 1e-8 ? ftan[nn][0] : 0;
      double ftc1 = ftan[nn][1] * ftan[nn][1] > 1e-8 ? ftan[nn][1] : 0;
      double ftc2 = ftan[nn][2] * ftan[nn][2] > 1e-8 ? ftan[nn][2] : 0;
      double ffftc = ftan[nn][0] * ftan[nn][0] + ftan[nn][1] * ftan[nn][1] + ftan[nn][2] * ftan[nn][2];
      ffftc = sqrt(ffftc);
      forces_file << std::setw(6) << mm
                  << std::setprecision(2) << std::setw(10) << ftc0 << "  "
                  << std::setw(10) << ftc1 << "  "
                  << std::setw(10) << ftc2 << "  "
                  << std::setw(10) << ffftc
                  << std::endl;
  }
  forces_file << " ================================================== "
              << std::endl;
  // ===================================================

}

/* ---------------------------------------------------------------------- */

void BondHarmonic::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(r0,n+1,"bond:r0");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondHarmonic::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double r0_one = force->numeric(FLERR,arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondHarmonic::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondHarmonic::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondHarmonic::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],r0[i]);
}

/* ---------------------------------------------------------------------- */

double BondHarmonic::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double rk = k[type] * dr;
  fforce = 0;
  if (r > 0.0) fforce = -2.0*rk/r;
  return rk*dr;
}
