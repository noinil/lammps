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
   Contributing author: Daniel Hinckley (Wisconsin/UChicago) dhinckley@wisc.edu
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "dihedral_list.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "math_const.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

#define TOLERANCE 0.05
#define SMALL     0.001
/* ---------------------------------------------------------------------- */

DihedralList::DihedralList(LAMMPS *lmp) : Dihedral(lmp) {}

/* ---------------------------------------------------------------------- */

DihedralList::~DihedralList()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(kperiodic);
    memory->destroy(kgauss);
    memory->destroy(phi);
    memory->destroy(sigm);
    memory->destroy(style);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralList::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,i,m,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double edihedral,f1[3],f2[3],f3[3],f4[3];
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
  double df,df1,ddf1,fg,hg,fga,hgb,gaa,gbb;
  double dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;
  double c,s,p,sx2,sy2,sz2;
  double cosphi, aphi, dtor, expt;

  edihedral = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  int *tag = atom->tag;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double etorsi, virixx, virixy, virixz, viriyx, viriyy, viriyz,
      virizx, virizy, virizz, xlen, ylen, zlen, xhlf, yhlf, zhlf,
      xlni, ylni, zlni, dabx, daby, dabz, dbcx, dbcy, dbcz, dcdx,
      dcdy, dcdz, drab, irab, eabx, eaby, eabz, drbc, irbc, ebcx,
      ebcy, ebcz, drcd, ircd, ecdx, ecdy, ecdz, pacx, pacy, pacz,
      cosb, isb2, pbdx, pbdy, pbdz, cosc, isc2, isnc, padx, pady,
      padz, cnum, isnb, ctau, fra1, frb1, frb2, frc1,
      frc2, frd1, fabx, faby, fabz, fbcx, fbcy, fbcz, fcdx, fcdy,
      fcdz;
  double vb1,vb2, vb3, eax,eay,eaz,ebx,eby,ebz;
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
  for (int nn = 0; nn < nlocal; nn++) {
      ftan[nn][0] = 0;
      ftan[nn][1] = 0;
      ftan[nn][2] = 0;
  }
  std::ofstream forces_file("n_force_dihedral.dat");
  std::ofstream energy_file("p_energy_dihedral.dat");
  forces_file << " dihedral forces: " << std::endl;
  energy_file << " dihedral energy: " << std::endl;
  energy_file << std::setw(6) << "dih_i"
              << std::setw(6) << "i1"
              << std::setw(6) << "i2"
              << std::setw(6) << "i3"
              << std::setw(6) << "i4"
              << std::setw(11) << "E_dih"
              << std::endl;
  energy_file << " ---------------------------------------------"
              << std::endl;
  // ===================================================

  for (n = 0; n < ndihedrallist; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = tag[i1] + tag[i2] + tag[i3] + tag[i4];

    // 1st bond
    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];
    vb1 = sqrt(vb1x * vb1x + vb1y * vb1y  + vb1z * vb1z );

    // 2nd bond
    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];
    vb2 = sqrt(vb2x * vb2x + vb2y * vb2y  + vb2z * vb2z );

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;

    // 3rd bond
    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];
    vb3 = sqrt(vb3x * vb3x + vb3y * vb3y  + vb3z * vb3z );

    // c,s calculation

    // a and b are cross-products
    ax = vb1y*vb2zm - vb1z*vb2ym;
    ay = vb1z*vb2xm - vb1x*vb2zm;
    az = vb1x*vb2ym - vb1y*vb2xm;
    eax = ax / (vb1 * vb2);
    eay = ay / (vb1 * vb2);
    eaz = az / (vb1 * vb2);
    bx = vb3y*vb2zm - vb3z*vb2ym;
    by = vb3z*vb2xm - vb3x*vb2zm;
    bz = vb3x*vb2ym - vb3y*vb2xm;
    ebx = bx / (vb2 * vb3);
    eby = by / (vb2 * vb3);
    ebz = bz / (vb2 * vb3);


    // Geting the squared length of a and b
    rasq = ax*ax + ay*ay + az*az;
    rbsq = bx*bx + by*by + bz*bz;
    rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm; //
    rg = sqrt(rgsq);

    rginv = ra2inv = rb2inv = 0.0;
    if (rg > 0) rginv = 1.0/rg;
    if (rasq > 0) ra2inv = 1.0/rasq;
    if (rbsq > 0) rb2inv = 1.0/rbsq;
    rabinv = sqrt(ra2inv*rb2inv);


    double cx = eay * ebz - eaz * eby;
    double cy = eaz * ebx - eax * ebz;
    double cz = eax * eby - eay * ebx;
    c = (ax*bx + ay*by + az*bz)*rabinv;

    double cnum = (cx * vb2x + cy * vb2y + cz * vb2z)/vb2;

    s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      int me;
      MPI_Comm_rank(world,&me);
      if (screen) {
        char str[128];
        sprintf(str,"Dihedral problem: %d " BIGINT_FORMAT " %d %d %d %d",
                me,update->ntimestep,
                atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
        error->warning(FLERR,str,0);
        fprintf(screen,"  1st atom: %d %g %g %g\n",
                me,x[i1][0],x[i1][1],x[i1][2]);
        fprintf(screen,"  2nd atom: %d %g %g %g\n",
                me,x[i2][0],x[i2][1],x[i2][2]);
        fprintf(screen,"  3rd atom: %d %g %g %g\n",
                me,x[i3][0],x[i3][1],x[i3][2]);
        fprintf(screen,"  4th atom: %d %g %g %g\n",
                me,x[i4][0],x[i4][1],x[i4][2]);
      }
    }

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    // Added specifically for torsion dihedral
    cosphi = c;
    aphi = acos(c);
    //printf("cnum=%lf\n",cnum);
    if (cnum < 0.0)
    {
        aphi = -aphi;
    }

    dtor = aphi - phi[type]; // Change made to get things to "work"
    //dtor = aphi + phi[type];

    if (dtor < -MY_PI)
    {
        dtor += 2.0 * MY_PI;
    }
    else if (dtor > MY_PI)
    {
        dtor -= 2.0 * MY_PI;
    }
    df = 0.0;
    edihedral = 0.0;
    if (style[type] == 1) {
        // Periodic torsion
        df = -kperiodic[type] * sin(dtor);
        if (eflag) edihedral = kperiodic[type] * (1-cos(dtor));
        // ==================== tc_test ====================--------------------
        etan = kperiodic[type] * (1-cos(dtor));
    } else if (style[type] == 2) {
        // Gaussian well + Periodic torsion
        expt = exp(-dtor * dtor/(2.0*sigm[type]*sigm[type]));
        df = -kgauss[type] * dtor * expt / (sigm[type] * sigm[type]);
        df += -kperiodic[type] * sin(dtor);
        if (eflag) edihedral = -kgauss[type] * expt;
        if (eflag) edihedral += kperiodic[type] * (1-cos(dtor));
        // ==================== tc_test ====================--------------------
        etan = -kgauss[type] * expt + kperiodic[type] * (1-cos(dtor));
    } else if (style[type] == 3) {
        df = -kperiodic[type] * sin(dtor) - sigm[type] * 3.0 * sin(3.0*dtor);
        if (eflag) edihedral = kperiodic[type] * (1-cos(dtor)) + sigm[type] * (1-cos(3.0 *dtor));
        // ==================== tc_test ====================--------------------
        etan = kperiodic[type] * (1-cos(dtor)) + sigm[type] * (1-cos(3.0 *dtor));
    }
    /* ------------------------------------------------ //
    //      _                         _            _    //
    //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
    //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
    // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
    //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
    //                       |_|                        //
    // ------------------ tc_test --------------------- */
    etan_total += etan;
    energy_file << std::setw(6) << n + 1
                << std::setw(6) << atom->tag[i1]
                << std::setw(6) << atom->tag[i2]
                << std::setw(6) << atom->tag[i3]
                << std::setw(6) << atom->tag[i4] << "   "
                << std::setw(10) << etan
                << std::endl;
    // ===================================================


    fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
    hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
    fga = fg*ra2inv*rginv;
    hgb = hg*rb2inv*rginv;
    gaa = -ra2inv*rg;
    gbb = rb2inv*rg;

    dtfx = gaa*ax;
    dtfy = gaa*ay;
    dtfz = gaa*az;
    dtgx = fga*ax - hgb*bx;
    dtgy = fga*ay - hgb*by;
    dtgz = fga*az - hgb*bz;
    dthx = gbb*bx;
    dthy = gbb*by;
    dthz = gbb*bz;

    sx2 = df*dtgx;
    sy2 = df*dtgy;
    sz2 = df*dtgz;

    f1[0] = df*dtfx;
    f1[1] = df*dtfy;
    f1[2] = df*dtfz;

    f2[0] = sx2 - f1[0];
    f2[1] = sy2 - f1[1];
    f2[2] = sz2 - f1[2];

    f4[0] = df*dthx;
    f4[1] = df*dthy;
    f4[2] = df*dthz;

    f3[0] = -sx2 - f4[0];
    f3[1] = -sy2 - f4[1];
    f3[2] = -sz2 - f4[2];

    // apply force to each of 4 atoms

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
      ftan[i1][0] += f1[0];
      ftan[i1][1] += f1[1];
      ftan[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
      /* ------------------------------------------------ //
      //      _                         _            _    //
      //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
      //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
      // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
      //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
      //                       |_|                        //
      // ------------------ tc_test --------------------- */
      ftan[i2][0] += f2[0];
      ftan[i2][1] += f2[1];
      ftan[i2][2] += f2[2];
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

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
      /* ------------------------------------------------ //
      //      _                         _            _    //
      //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
      //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
      // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
      //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
      //                       |_|                        //
      // ------------------ tc_test --------------------- */
      ftan[i4][0] += f4[0];
      ftan[i4][1] += f4[1];
      ftan[i4][2] += f4[2];
    }

    if (evflag)
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,edihedral,f1,f3,f4,
               vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
  }

  /* ------------------------------------------------ //
  //      _                         _            _    //
  //   __| |_   _ _ __ ___  _ __   | |_ ___  ___| |_  //
  //  / _` | | | | '_ ` _ \| '_ \  | __/ _ \/ __| __| //
  // | (_| | |_| | | | | | | |_) | | ||  __/\__ \ |_  //
  //  \__,_|\__,_|_| |_| |_| .__/   \__\___||___/\__| //
  //                       |_|                        //
  // ------------------ tc_test --------------------- */
  energy_file << "Total dihedral energy: " << etan_total << std::endl;
  energy_file << " ================================================== "
              << std::endl;
  forces_file << std::setw(6) << "dih_i"
              << std::setw(12) << "fx"
              << std::setw(12) << "fy"
              << std::setw(12) << "fz"
              << std::setw(12) << "|f|"
              << std::endl
              << " -------------------------------------------------------------"
              << std::endl;
  for (int mm = 1; mm < nlocal + 1; mm++) {
      int nn = atom->map(mm);
      // double ftc0 = ftan[nn][0] * ftan[nn][0] > 1e-8 ? ftan[nn][0] : 0;
      // double ftc1 = ftan[nn][1] * ftan[nn][1] > 1e-8 ? ftan[nn][1] : 0;
      // double ftc2 = ftan[nn][2] * ftan[nn][2] > 1e-8 ? ftan[nn][2] : 0;
      double ftc0 = ftan[nn][0];
      double ftc1 = ftan[nn][1];
      double ftc2 = ftan[nn][2];
      double ffftc = ftan[nn][0] * ftan[nn][0] + ftan[nn][1] * ftan[nn][1] + ftan[nn][2] * ftan[nn][2];
      ffftc = sqrt(ffftc);
      forces_file << std::setw(6) << mm << "  "
                  << std::setprecision(5) << std::setw(10) << ftc0 << "  "
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

void DihedralList::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  //memory->create(k,n+1,"dihedral:k");
  //memory->create(phi,n+1,"dihedral:phi");
  //memory->create(sigm,n+1,"dihedral:sigm");
  memory->create(setflag,n+1,"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void DihedralList::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal dihedral_style command");

  FILE *fp = fopen(arg[0],"r");
  char line[1024];
  if (fp == NULL)
    error->all(FLERR,"Cannot open dihedral list file");

  int num = 1;
  while(fgets(line,1024,fp)) ++num;
  rewind(fp);
  int array_size = 20*(num+2); // Need to optimize this value...
  memory->create(kperiodic,array_size,"dihedral:kperiodic");
  memory->create(kgauss,array_size,"dihedral:kgauss");
  memory->create(phi,array_size,"dihedral:phi");
  memory->create(sigm,array_size,"dihedral:sigm");
  memory->create(style,array_size,"dihedral:style");

  char *ptr;
  int idx, id1, id2, id3, id4;

  // Allocate arrays that are 4*ndihedrals long
  while(fgets(line,1024,fp)) {
    ptr = strtok(line," \t\n\r\f");

    // skip empty lines
    if (!ptr) continue;

    // skip comment lines starting with #
    if (*ptr == '#') continue;

    id1 = atoi(ptr);

    // The second site
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted dihedral list file");
    id2 = atoi(ptr);

    // The third site
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted dihedral list file");
    id3 = atoi(ptr);

    // The fourth site
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted dihedral list file");
    id4 = atoi(ptr);

    // Setting the idx in the base array
    idx = id1 + id2 + id3 + id4;
    if ((idx-1) > array_size)
        error->all(FLERR,"Parameter array in dihedral_list.cpp is too short!");

    // phi
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted dihedral list file");
    phi[idx] = force->numeric(FLERR,ptr);
    phi[idx] = phi[idx] * MY_PI/180.0 + MY_PI;

    // kperiodic
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted dihedral list file");
    kperiodic[idx] = force->numeric(FLERR,ptr);

    // kgauss
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted dihedral list file");
    kgauss[idx] = force->numeric(FLERR,ptr);

    // sigm or K3
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted dihedral list file");
    sigm[idx] = force->numeric(FLERR,ptr);

    // Types [0] - 1-  Cosine, 2 - Cosine+Gaussian, 3 - Cafemol
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted dihedral list file");
    style[idx] = atoi(ptr);

    //printf("%d %f %f %f\n",idx,kperiodic[idx],phi[idx],sigm[idx]);
  }

  fclose(fp);

}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void DihedralList::coeff(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Incorrect args for dihedral coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->ndihedraltypes,ilo,ihi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    //printf("phi_one=%lf (%lf)\n",phi[i],phi[i] * 180.0 / MY_PI);
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void DihedralList::write_restart(FILE *fp)
{
  fwrite(&kperiodic[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&kgauss[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&phi[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&sigm[1],sizeof(double),atom->ndihedraltypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void DihedralList::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&kperiodic[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&kgauss[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&phi[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&sigm[1],sizeof(double),atom->ndihedraltypes,fp);
  }
  MPI_Bcast(&kperiodic[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kgauss[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&phi[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigm[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    setflag[i] = 1;
  }
}
