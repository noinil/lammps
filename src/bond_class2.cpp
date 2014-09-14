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

/* ----------------------------------------------------------------------
   Contributing author: Eric Simon (Cray)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "bond_class2.h"
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

BondClass2::BondClass2(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondClass2::~BondClass2()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(r0);
    memory->destroy(k2);
    memory->destroy(k3);
    memory->destroy(k4);
  }
}

/* ---------------------------------------------------------------------- */

void BondClass2::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr,dr2,dr3,dr4,de_bond;

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
  std::ofstream forces_file("n_force_bd_class2.dat");
  std::ofstream energy_file("p_energy_bd_class2.dat");
  forces_file << " stacking forces: " << std::endl;
  energy_file << " stacking energy: " << std::endl;
  energy_file << std::setw(6) << "bd2_i"
              << std::setw(6) << "i1"
              << std::setw(6) << "i2"
              << std::setw(11) << "E_stk"
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
    dr2 = dr*dr;
    dr3 = dr2*dr;
    dr4 = dr3*dr;

    // force & energy

    de_bond = 2.0*k2[type]*dr + 3.0*k3[type]*dr2 + 4.0*k4[type]*dr3;
    if (r > 0.0) fbond = -de_bond/r;
    else fbond = 0.0;

    if (eflag) ebond = k2[type]*dr2 + k3[type]*dr3 + k4[type]*dr4;

    /* ------------------------------------------------ //
    //      _                         _            _    //
    //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
    //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
    // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
    //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
    //                       |_|                        //
    // ------------------ tc_test --------------------- */
    etan = k2[type]*dr2 + k3[type]*dr3 + k4[type]*dr4;
    etan_total += etan;
    energy_file << std::setw(6) << n + 1
                << std::setw(6) << atom->tag[i1]
                << std::setw(6) << atom->tag[i2] << " "
                << std::setprecision(3) << std::setw(10) << etan
                << std::endl;
    // ===================================================

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fbond;
      f[i1][1] += dely*fbond;
      f[i1][2] += delz*fbond;
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
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fbond;
      f[i2][1] -= dely*fbond;
      f[i2][2] -= delz*fbond;
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
  energy_file << "Total bond class2 energy: " << etan_total << std::endl;
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

void BondClass2::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(r0,n+1,"bond:r0");
  memory->create(k2,n+1,"bond:k2");
  memory->create(k3,n+1,"bond:k3");
  memory->create(k4,n+1,"bond:k4");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs from one line in input script or data file
------------------------------------------------------------------------- */

void BondClass2::coeff(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double r0_one = force->numeric(FLERR,arg[1]);
  double k2_one = force->numeric(FLERR,arg[2]);
  double k3_one = force->numeric(FLERR,arg[3]);
  double k4_one = force->numeric(FLERR,arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    r0[i] = r0_one;
    k2[i] = k2_one;
    k3[i] = k3_one;
    k4[i] = k4_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondClass2::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondClass2::write_restart(FILE *fp)
{
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&k2[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&k3[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&k4[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondClass2::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&k2[1],sizeof(double),atom->nbondtypes,fp);
    fread(&k3[1],sizeof(double),atom->nbondtypes,fp);
    fread(&k4[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k2[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k3[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k4[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondClass2::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,r0[i],k2[i],k3[i],k4[i]);
}

/* ---------------------------------------------------------------------- */

double BondClass2::single(int type, double rsq, int i, int j, double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double dr2 = dr*dr;
  double dr3 = dr2*dr;
  double dr4 = dr3*dr;
  double de_bond = 2.0*k2[type]*dr + 3.0*k3[type]*dr2 + 4.0*k4[type]*dr3;
  if (r > 0.0) fforce = -de_bond/r;
  else fforce = 0.0;
  return (k2[type]*dr2 + k3[type]*dr3 + k4[type]*dr4);
}
