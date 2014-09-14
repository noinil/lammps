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
#include "angle_harmonic.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
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
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleHarmonic::AngleHarmonic(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleHarmonic::~AngleHarmonic()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(theta0);
  }
}

/* ---------------------------------------------------------------------- */

void AngleHarmonic::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double dtheta,tk;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;

  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
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
  std::ofstream forces_file("n_force_angle.dat");
  std::ofstream energy_file("p_energy_angle.dat");
  forces_file << " angle forces: " << std::endl;
  energy_file << " angle energy: " << std::endl;
  energy_file << std::setw(6) << "ang_i"
              << std::setw(6) << "i1"
              << std::setw(6) << "i2"
              << std::setw(6) << "i3"
              << std::setprecision(3) << std::setw(11) << "E_ang"
              << std::endl;
  energy_file << " ---------------------------------------------"
              << std::endl;
  // ===================================================

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // angle (cos and sin)

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    // force & energy

    dtheta = acos(c) - theta0[type];
    tk = k[type] * dtheta;

    if (eflag) eangle = tk*dtheta;
    /* ------------------------------------------------ //
    //      _                         _            _    //
    //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
    //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
    // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
    //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
    //                       |_|                        //
    // ------------------ tc_test --------------------- */
    etan = tk*dtheta;
    etan_total += etan;
    energy_file << std::setw(6) << n + 1
                << std::setw(6) << atom->tag[i1]
                << std::setw(6) << atom->tag[i2]
                << std::setw(6) << atom->tag[i3] << " "
                << std::setw(10) << etan
                << std::endl;
    // ===================================================

    a = -2.0 * tk * s;
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
      /* ------------------------------------------------ //
      //      _                         _            _    //
      //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
      //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
      // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
      //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
      //                       |_|                        //
      // ------------------ tc_test --------------------- */
      ftan[i1][0] -= f1[0];
      ftan[i1][1] -= f1[1];
      ftan[i1][2] -= f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
      /* ------------------------------------------------ //
      //      _                         _            _    //
      //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
      //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
      // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
      //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
      //                       |_|                        //
      // ------------------ tc_test --------------------- */
      ftan[i2][0] -= f1[0] + f3[0];
      ftan[i2][1] -= f1[1] + f3[1];
      ftan[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
      /* ------------------------------------------------ //
      //      _                         _            _    //
      //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
      //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
      // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
      //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
      //                       |_|                        //
      // ------------------ tc_test --------------------- */
      ftan[i3][0] += f3[0];
      ftan[i3][1] += f3[1];
      ftan[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
  /* ------------------------------------------------ //
  //      _                         _            _    //
  //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
  //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
  // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
  //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
  //                       |_|                        //
  // ------------------ tc_test --------------------- */
  energy_file << "Total angle energy: " << etan_total << std::endl;
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

void AngleHarmonic::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k,n+1,"angle:k");
  memory->create(theta0,n+1,"angle:theta0");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleHarmonic::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nangletypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double theta0_one = force->numeric(FLERR,arg[2]);

  // convert theta0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    theta0[i] = theta0_one/180.0 * MY_PI;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleHarmonic::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleHarmonic::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&theta0[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nangletypes,fp);
    fread(&theta0[1],sizeof(double),atom->nangletypes,fp);
  }
  MPI_Bcast(&k[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta0[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleHarmonic::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],theta0[i]/MY_PI*180.0);
}

/* ---------------------------------------------------------------------- */

double AngleHarmonic::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double dtheta = acos(c) - theta0[type];
  double tk = k[type] * dtheta;
  return tk*dtheta;
}
