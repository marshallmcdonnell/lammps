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
   Contributing author: Paul Crozier, Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "compute_widom.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_hybrid.h"
#include "molecule.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "comm.h"
#include "compute.h"
#include "group.h"
#include "domain.h"
#include "region.h"
#include "random_park.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "thermo.h"
#include "output.h"
#include "neighbor.h"
#include <iostream>

using namespace std;
using namespace LAMMPS_NS;
using namespace MathConst;

// large energy value used to signal overlap

#define MAXENERGYSIGNAL 1.0e100

// this must be lower than MAXENERGYSIGNAL
// by a large amount, so that it is still
// less than total energy when negative
// energy contributions are added to MAXENERGYSIGNAL

#define MAXENERGYTEST 1.0e50

// ensemble options
enum{NVT,NPT,NVE};

/* ---------------------------------------------------------------------- */

ComputeWidom::ComputeWidom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  idregion(NULL), full_flag(0), ngroups(0), groupstrings(NULL), ngrouptypes(0), grouptypestrings(NULL),
  grouptypebits(NULL), grouptypes(NULL), random_equal(NULL), random_unequal(NULL),
  ensemble(0), tflag(0), pflag(0),
  wtest(0), wtest_sq(0), beta_mu_ex(0), beta_mu_ex_sigma(0),
  chemical_potential(0), chemical_potential_sigma(0)

{
  if (narg < 9) error->all(FLERR,"Illegal compute widom command");

  if (atom->molecular == 2)
    error->all(FLERR,"Compute widom does not (yet) work with atom_style template");

  vector_flag = 1;
  size_vector = 2;
  extvector = 0;

  dynamic_group_allow = 1;

  // required args

  nevery = force->inumeric(FLERR,arg[3]);
  ninsertions = force->inumeric(FLERR,arg[4]);
  nwidom_type = force->inumeric(FLERR,arg[5]);
  seed = force->inumeric(FLERR,arg[6]);
  reservoir_temperature = force->numeric(FLERR,arg[7]);

  if (strcmp(arg[8],"nvt") == 0) {
    ensemble = NVT;
  } else if (strcmp(arg[8],"npt") == 0) {
    ensemble = NPT;
  } else {
    error->all(FLERR,"Unsupporte ensemble in widom command");
  }


  if (nevery <= 0) error->all(FLERR,"Illegal compute widom command");
  if (ninsertions < 0) error->all(FLERR,"Illegal compute widom command");
  if (seed <= 0) error->all(FLERR,"Illegal compute widom command");
  if (reservoir_temperature < 0.0)
    error->all(FLERR,"Illegal compute widom command");

  // read options from end of input line

  options(narg-9,&arg[9]);

  // check that pressure is set if NPT
    
  if (ensemble == NPT && !pflag) 
    error->all(FLERR,"Must specifiy pressure if using NPT in compute widom command");

  // random number generator, same for all procs

  random_equal = new RanPark(lmp,seed);

  // random number generator, not the same for all procs

  random_unequal = new RanPark(lmp,seed);

  // error checks on region and its extent being inside simulation box

  region_xlo = region_xhi = region_ylo = region_yhi =
    region_zlo = region_zhi = 0.0;
  if (regionflag) {
    if (domain->regions[iregion]->bboxflag == 0)
      error->all(FLERR,"Compute widom region does not support a bounding box");
    if (domain->regions[iregion]->dynamic_check())
      error->all(FLERR,"Compute widom region cannot be dynamic");

    region_xlo = domain->regions[iregion]->extent_xlo;
    region_xhi = domain->regions[iregion]->extent_xhi;
    region_ylo = domain->regions[iregion]->extent_ylo;
    region_yhi = domain->regions[iregion]->extent_yhi;
    region_zlo = domain->regions[iregion]->extent_zlo;
    region_zhi = domain->regions[iregion]->extent_zhi;

    if (region_xlo < domain->boxlo[0] || region_xhi > domain->boxhi[0] ||
        region_ylo < domain->boxlo[1] || region_yhi > domain->boxhi[1] ||
        region_zlo < domain->boxlo[2] || region_zhi > domain->boxhi[2])
      error->all(FLERR,"Compute widom region extends outside simulation box");

    // estimate region volume using MC trials

    double coord[3];
    int inside = 0;
    int attempts = 10000000;
    for (int i = 0; i < attempts; i++) {
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      if (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) != 0)
        inside++;
    }

    double max_region_volume = (region_xhi - region_xlo)*
     (region_yhi - region_ylo)*(region_zhi - region_zlo);

    region_volume = max_region_volume*static_cast<double> (inside)/
     static_cast<double> (attempts);
  }

  if (charge_flag && atom->q == NULL)
    error->all(FLERR,"Compute widom atom has charge, but atom style does not");

  // Setup temperature if using NVE
  // create a new compute temp style
  // id = fix-ID + temp and
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all
  if ( ensemble == NVE ) {

    // temperature
    int n = strlen(id) + 6;
    id_temp = new char[n];
    strcpy(id_temp,id);
    strcat(id_temp,"_temp");

    char **newarg = new char*[3];
    newarg[0] = id_temp;
    newarg[1] = (char *) "all";
    newarg[2] = (char *) "temp";
    modify->add_compute(3,newarg);
    delete [] newarg;
    tflag = 1;
  }


  // compute the number of MC cycles that occur nevery timesteps

  ncycles = ninsertions;

  //ncycles = ninsertions + nrotations;

  widom_nmax = 0;

  vector = new double[size_vector];
    
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void ComputeWidom::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal compute widom command");

  // defaults

  //max_rotation_angle = 10*MY_PI/180;
  regionflag = 0;
  iregion = -1;
  region_volume = 0;
  max_region_attempts = 1000;
  exclusion_group = 0;
  exclusion_group_bit = 0;
  charge = 0.0;
  charge_flag = false;
  full_flag = false;
  ngroups = 0;
  int ngroupsmax = 0;
  groupstrings = NULL;
  ngrouptypes = 0;
  int ngrouptypesmax = 0;
  grouptypestrings = NULL;
  grouptypes = NULL;
  grouptypebits = NULL;
  energy_intra = 0.0;
  tfac_insert = 1.0;
  overlap_cutoffsq = 0.0;
  overlap_flag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute widom command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for compute widom does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      regionflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"maxangle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute widom command");
      //max_rotation_angle = force->numeric(FLERR,arg[iarg+1]);
      //max_rotation_angle *= MY_PI/180;
      iarg += 2;
    } else if (strcmp(arg[iarg],"pressure") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute widom command");
      pressure = force->numeric(FLERR,arg[iarg+1]);
      pflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"charge") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute widom command");
      charge = force->numeric(FLERR,arg[iarg+1]);
      charge_flag = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"full_energy") == 0) {
      full_flag = true;
      iarg += 1;
    } else if (strcmp(arg[iarg],"group") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute widom command");
      if (ngroups >= ngroupsmax) {
        ngroupsmax = ngroups+1;
        groupstrings = (char **)
          memory->srealloc(groupstrings,
                           ngroupsmax*sizeof(char *),
                           "compute_widom:groupstrings");
      }
      int n = strlen(arg[iarg+1]) + 1;
      groupstrings[ngroups] = new char[n];
      strcpy(groupstrings[ngroups],arg[iarg+1]);
      ngroups++;
      iarg += 2;
    } else if (strcmp(arg[iarg],"grouptype") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal compute widom command");
      if (ngrouptypes >= ngrouptypesmax) {
        ngrouptypesmax = ngrouptypes+1;
        grouptypes = (int*) memory->srealloc(grouptypes,ngrouptypesmax*sizeof(int),
                         "compute_widom:grouptypes");
        grouptypestrings = (char**)
          memory->srealloc(grouptypestrings,
                           ngrouptypesmax*sizeof(char *),
                           "compute_widom:grouptypestrings");
      }
      grouptypes[ngrouptypes] = atoi(arg[iarg+1]);
      int n = strlen(arg[iarg+2]) + 1;
      grouptypestrings[ngrouptypes] = new char[n];
      strcpy(grouptypestrings[ngrouptypes],arg[iarg+2]);
      ngrouptypes++;
      iarg += 3;
    } else if (strcmp(arg[iarg],"intra_energy") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute widom command");
      energy_intra = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"tfac_insert") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute widom command");
      tfac_insert = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"overlap_cutoff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute widom command");
      double rtmp = force->numeric(FLERR,arg[iarg+1]);
      overlap_cutoffsq = rtmp*rtmp;
      overlap_flag = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal compute widom command");
  }
}

/* ---------------------------------------------------------------------- */

ComputeWidom::~ComputeWidom()
{
  if (regionflag) delete [] idregion;
  delete random_equal;
  delete random_unequal;
  delete [] vector;

  if (ngroups > 0) {
    for (int igroup = 0; igroup < ngroups; igroup++)
      delete [] groupstrings[igroup];
    memory->sfree(groupstrings);
  }

  if (ngrouptypes > 0) {
    memory->destroy(grouptypes);
    memory->destroy(grouptypebits);
    for (int igroup = 0; igroup < ngrouptypes; igroup++)
      delete [] grouptypestrings[igroup];
    memory->sfree(grouptypestrings);
  }
  if (full_flag && group) {
    int igroupall = group->find("all");
    neighbor->exclusion_group_group_delete(exclusion_group,igroupall);
  }

  if (tflag) modify->delete_compute(id_temp);

}


/* ---------------------------------------------------------------------- */

void ComputeWidom::init()
{

  triclinic = domain->triclinic;

  // decide whether to switch to the full_energy option

  if (!full_flag) {
    if ((force->kspace) ||
        (force->pair == NULL) ||
        (force->pair->single_enable == 0) ||
        (force->pair_match("hybrid",0)) ||
        (force->pair_match("eam",0)) ||
        (force->pair->tail_flag)
        ) {
      full_flag = true;
      if (comm->me == 0)
        error->warning(FLERR,"Compute widom using full_energy option");
    }
  }

  if (full_flag) {
    char *id_pe = (char *) "thermo_pe";
    int ipe = modify->find_compute(id_pe);
    c_pe = modify->compute[ipe];
  }

  int *type = atom->type;

  if (nwidom_type <= 0 || nwidom_type > atom->ntypes)
    error->all(FLERR,"Invalid atom type in compute widom command");

  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use compute widom in a 2d simulation");

  // create a new group for interaction exclusions
  // used for attempted atom or molecule deletions
  // skip if already exists from previous init()

  if (full_flag && !exclusion_group_bit) {
    char **group_arg = new char*[4];

    // create unique group name for atoms to be excluded

    int len = strlen(id) + 30;
    group_arg[0] = new char[len];
    sprintf(group_arg[0],"ComputeWidom:gcmc_exclusion_group:%s",id);
    group_arg[1] = (char *) "subtract";
    group_arg[2] = (char *) "all";
    group_arg[3] = (char *) "all";
    group->assign(4,group_arg);
    exclusion_group = group->find(group_arg[0]);
    if (exclusion_group == -1)
      error->all(FLERR,"Could not find compute widom exclusion group ID");
    exclusion_group_bit = group->bitmask[exclusion_group];

    // neighbor list exclusion setup
    // turn off interactions between group all and the exclusion group

    int narg = 4;
    char **arg = new char*[narg];;
    arg[0] = (char *) "exclude";
    arg[1] = (char *) "group";
    arg[2] = group_arg[0];
    arg[3] = (char *) "all";
    neighbor->modify_params(narg,arg);
    delete [] group_arg[0];
    delete [] group_arg;
    delete [] arg;
  }

  // get mass
  gas_mass = atom->mass[nwidom_type];

  if (gas_mass <= 0.0)
    error->all(FLERR,"Illegal compute widom gas mass <= 0");

  // compute beta, lambda, sigma, and the zz factor
  // For LJ units, lambda=1
  beta = 1.0/(force->boltz*reservoir_temperature);

  sigma = sqrt(force->boltz*reservoir_temperature*tfac_insert/gas_mass/force->mvv2e);

  imagezero = ((imageint) IMGMAX << IMG2BITS) |
             ((imageint) IMGMAX << IMGBITS) | IMGMAX;

  // construct group bitmask for all new atoms
  // aggregated over all group keywords

  groupbitall = 1 | groupbit;
  for (int igroup = 0; igroup < ngroups; igroup++) {
    int jgroup = group->find(groupstrings[igroup]);
    if (jgroup == -1)
      error->all(FLERR,"Could not find specified compute widom group ID");
    groupbitall |= group->bitmask[jgroup];
  }

  // construct group type bitmasks
  // not aggregated over all group keywords

  if (ngrouptypes > 0) {
    memory->create(grouptypebits,ngrouptypes,"compute_widom:grouptypebits");
    for (int igroup = 0; igroup < ngrouptypes; igroup++) {
      int jgroup = group->find(grouptypestrings[igroup]);
      if (jgroup == -1)
        error->all(FLERR,"Could not find specified compute widom group ID");
      grouptypebits[igroup] = group->bitmask[jgroup];
    }
  }

}

/* ----------------------------------------------------------------------
   attempt Monte Carlo insertions to calculate chemical potential
   before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void ComputeWidom::compute_vector()
{
  invoked_vector = update->ntimestep;

  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];
  if (triclinic) {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  } else {
    sublo = domain->sublo;
    subhi = domain->subhi;
  }

  if (regionflag) volume = region_volume;
  else volume = domain->xprd * domain->yprd * domain->zprd;

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  // Zero out the Boltzmann factor variables
  wtest = 0.0; 
  wtest_sq = 0.0; 

  if (full_flag) {
    energy_stored = energy_full();
    if (overlap_flag && energy_stored > MAXENERGYTEST)
        error->warning(FLERR,"Energy of old configuration in "
                       "compute widom is > MAXENERGYTEST.");

    for (int i = 0; i < ncycles; i++) {
      atomic_insertion_full();
    }

    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    comm->exchange();
    atom->nghost = 0;
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);

  } else {

    for (int i = 0; i < ncycles; i++) {
      atomic_insertion();
    }

    // Do ensemble-dependent calculations
    if ( ensemble == NVT ) {
        //get_nvt_mu_id();
    } else if ( ensemble == NPT ) {
      //get_npt_mu_id();
      wtest *= beta * ( pressure * volume) / (atom->natoms+1);
    } else if ( ensemble == NVE ) {
      // get temperature
      int icompute = modify->find_compute(id_temp);
      if (icompute < 0)
        error->all(FLERR,"Temperature ID for compute widom does not exist");
      temperature = modify->compute[icompute]->compute_scalar();

      //wtest *= avg_t ^-3/2 * wtest*T_i^3/2
    }

  }

  // Reduce the accumulated boltzmann factors sampled from test insertions
  double wtest_all = 0.0;
  double wtest_sq_all = 0.0;
  MPI_Allreduce(&wtest,&wtest_all,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&wtest_sq,&wtest_sq_all,1,MPI_DOUBLE,MPI_SUM,world);

  // Get the statistics of wtest
  double wtest_sum = wtest_all / (double) ncycles;
  double wtest_sq_sum = wtest_sq_all / (double) ncycles;
  double var_wtest = wtest_sq_sum - wtest_sum * wtest_sum;
  
  // Collect chemical potential & stats
  beta_mu_ex = -1.0 * log( wtest_all / (double) ncycles );
  beta_mu_ex_sigma = sqrt( var_wtest ) / wtest_all;

  chemical_potential = beta_mu_ex / beta;
  chemical_potential_sigma = beta_mu_ex_sigma / beta;

  // Output vector
  vector[0] = chemical_potential;
  vector[1] = chemical_potential_sigma;

}



/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void ComputeWidom::atomic_insertion()
{
  double lamda[3];

  // pick coordinates for insertion point

  double coord[3];
  if (regionflag) {
    int region_attempt = 0;
    coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
    coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
    coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
    while (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) == 0) {
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      region_attempt++;
      if (region_attempt >= max_region_attempts) return;
    }
    if (triclinic) domain->x2lamda(coord,lamda);
  } else {
    if (triclinic == 0) {
      coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
      coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
      coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
    } else {
      lamda[0] = random_equal->uniform();
      lamda[1] = random_equal->uniform();
      lamda[2] = random_equal->uniform();

      // wasteful, but necessary

      if (lamda[0] == 1.0) lamda[0] = 0.0;
      if (lamda[1] == 1.0) lamda[1] = 0.0;
      if (lamda[2] == 1.0) lamda[2] = 0.0;

      domain->lamda2x(lamda,coord);
    }
  }

  // Flag if the coordinate is on this proc
  int proc_flag = 0;
  if (triclinic == 0) {
    domain->remap(coord);
    if (!domain->inside(coord))
      error->one(FLERR,"Compute widom put atom outside box");
    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) proc_flag = 1;
  } else {
    if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
        lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
        lamda[2] >= sublo[2] && lamda[2] < subhi[2]) proc_flag = 1;
  }

  // Test insertions if on this proc
  if (proc_flag) {
    int ii = -1;
    if (charge_flag) {
      ii = atom->nlocal + atom->nghost;
      if (ii >= atom->nmax) atom->avec->grow(0);
      atom->q[ii] = charge;
    }
    double insertion_energy = energy(ii,nwidom_type,-1,coord);

    // collect energy sum for insertion
    double expE_kbT  = exp( -beta * insertion_energy);
    wtest    += expE_kbT;
    wtest_sq += expE_kbT * expE_kbT;

  }

}

/* ----------------------------------------------------------------------
-------------------------------------------------------------------------*/ 

void ComputeWidom::atomic_insertion_full()
{
  double lamda[3];

  double energy_before = energy_stored;

  double coord[3];
  if (regionflag) {
    int region_attempt = 0;
    coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
    coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
    coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
    while (domain->regions[iregion]->match(coord[0],coord[1],coord[2]) == 0) {
      coord[0] = region_xlo + random_equal->uniform() * (region_xhi-region_xlo);
      coord[1] = region_ylo + random_equal->uniform() * (region_yhi-region_ylo);
      coord[2] = region_zlo + random_equal->uniform() * (region_zhi-region_zlo);
      region_attempt++;
      if (region_attempt >= max_region_attempts) return;
    }
    if (triclinic) domain->x2lamda(coord,lamda);
  } else {
    if (triclinic == 0) {
      coord[0] = xlo + random_equal->uniform() * (xhi-xlo);
      coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
      coord[2] = zlo + random_equal->uniform() * (zhi-zlo);
    } else {
      lamda[0] = random_equal->uniform();
      lamda[1] = random_equal->uniform();
      lamda[2] = random_equal->uniform();

      // wasteful, but necessary

      if (lamda[0] == 1.0) lamda[0] = 0.0;
      if (lamda[1] == 1.0) lamda[1] = 0.0;
      if (lamda[2] == 1.0) lamda[2] = 0.0;

      domain->lamda2x(lamda,coord);
    }
  }

  int proc_flag = 0;
  if (triclinic == 0) {
    domain->remap(coord);
    if (!domain->inside(coord))
      error->one(FLERR,"Compute widom put atom outside box");
    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
        coord[2] >= sublo[2] && coord[2] < subhi[2]) proc_flag = 1;
  } else {
    if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
        lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
        lamda[2] >= sublo[2] && lamda[2] < subhi[2]) proc_flag = 1;
  }

  if (proc_flag) {
    atom->avec->create_atom(nwidom_type,coord);
    int m = atom->nlocal - 1;

    // add to groups
    // optionally add to type-based groups

    atom->mask[m] = groupbitall;
    for (int igroup = 0; igroup < ngrouptypes; igroup++) {
      if (nwidom_type == grouptypes[igroup])
        atom->mask[m] |= grouptypebits[igroup];
    }

    if (charge_flag) atom->q[m] = charge;
    modify->create_attribute(m);
  }

  atom->natoms++;
  if (atom->tag_enable) {
    atom->tag_extend();
    if (atom->map_style) atom->map_init();
  }
  atom->nghost = 0;
  if (triclinic) domain->x2lamda(atom->nlocal);
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  if (force->kspace) force->kspace->qsum_qsq();
  if (force->pair->tail_flag) force->pair->reinit();
  double energy_after = energy_full();

  double insertion_energy = energy_after - energy_before;

  // collect energy sum for insertion
  double expE_kbT  = exp( -beta * insertion_energy);
  wtest    += expE_kbT;
  wtest_sq += expE_kbT * expE_kbT;

  // Reset the insertion
  atom->natoms--;
  if (proc_flag) atom->nlocal--;
  if (force->kspace) force->kspace->qsum_qsq();
  if (force->pair->tail_flag) force->pair->reinit();
}

/* ----------------------------------------------------------------------
   compute particle's interaction energy with the rest of the system
------------------------------------------------------------------------- */

double ComputeWidom::energy(int i, int itype, tagint imolecule, double *coord)
{
  double delx,dely,delz,rsq;

  double **x = atom->x;
  int *type = atom->type;
  tagint *molecule = atom->molecule;
  int nall = atom->nlocal + atom->nghost;
  pair = force->pair;
  cutsq = force->pair->cutsq;

  double fpair = 0.0;
  double factor_coul = 1.0;
  double factor_lj = 1.0;

  double total_energy = 0.0;

  for (int j = 0; j < nall; j++) {

    if (i == j) continue;

    delx = coord[0] - x[j][0];
    dely = coord[1] - x[j][1];
    delz = coord[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    int jtype = type[j];

    // if overlap check requested, if overlap,
    // return signal value for energy

    if (overlap_flag && rsq < overlap_cutoffsq)
      return MAXENERGYSIGNAL;

    if (rsq < cutsq[itype][jtype])
      total_energy +=
        pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
  }

  return total_energy;
}

/* ----------------------------------------------------------------------
   compute system potential energy
-------------------------------------------------------------------------*/ 

double ComputeWidom::energy_full()
{

  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build(1);
  int eflag = 1;
  int vflag = 0;

  // if overlap check requested, if overlap,
  // return signal value for energy

  if (overlap_flag) {
    int overlaptestall;
    int overlaptest = 0;
    double delx,dely,delz,rsq;
    double **x = atom->x;
    tagint *molecule = atom->molecule;
    int nall = atom->nlocal + atom->nghost;
    for (int i = 0; i < atom->nlocal; i++) {
      for (int j = i+1; j < nall; j++) {

        delx = x[i][0] - x[j][0];
        dely = x[i][1] - x[j][1];
        delz = x[i][2] - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq < overlap_cutoffsq) {
          overlaptest = 1;
          break;
        }
      }
      if (overlaptest) break;
    }
    MPI_Allreduce(&overlaptest, &overlaptestall, 1,
                  MPI_INT, MPI_MAX, world);
    if (overlaptestall) return MAXENERGYSIGNAL;
  }

  // clear forces so they don't accumulate over multiple
  // calls within compute widom timestep, e.g. for fix shake

  size_t nbytes = sizeof(double) * (atom->nlocal + atom->nghost);
  if (nbytes) memset(&atom->f[0][0],0,3*nbytes);

  if (modify->n_pre_force) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) force->kspace->compute(eflag,vflag);

  // unlike Verlet, not performing a reverse_comm() or forces here
  // b/c Widom does not care about forces
  // don't think it will mess up energy due to any post_force() fixes

  if (modify->n_post_force) modify->post_force(vflag);
  if (modify->n_end_of_step) modify->end_of_step();

  // NOTE: all fixes with THERMO_ENERGY mask set and which
  //   operate at pre_force() or post_force() or end_of_step()
  //   and which user has enable via fix_modify thermo yes,
  //   will contribute to total MC energy via pe->compute_scalar()

  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void ComputeWidom::toggle_intramolecular(int i)
{
  if (atom->avec->bonds_allow)
    for (int m = 0; m < atom->num_bond[i]; m++)
      atom->bond_type[i][m] = -atom->bond_type[i][m];

  if (atom->avec->angles_allow)
    for (int m = 0; m < atom->num_angle[i]; m++)
      atom->angle_type[i][m] = -atom->angle_type[i][m];

  if (atom->avec->dihedrals_allow)
    for (int m = 0; m < atom->num_dihedral[i]; m++)
      atom->dihedral_type[i][m] = -atom->dihedral_type[i][m];

  if (atom->avec->impropers_allow)
    for (int m = 0; m < atom->num_improper[i]; m++)
      atom->improper_type[i][m] = -atom->improper_type[i][m];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double ComputeWidom::memory_usage()
{
  double bytes = widom_nmax * sizeof(int);
  return bytes;
}
