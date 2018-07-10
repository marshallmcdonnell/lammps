/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(widom,ComputeWidom)

#else

#ifndef LMP_FIX_WIDOM_H
#define LMP_FIX_WIDOM_H

#include <cstdio>
#include "compute.h"

namespace LAMMPS_NS {

class ComputeWidom : public Compute {
 public:
  ComputeWidom(class LAMMPS *, int, char **);
  ~ComputeWidom();
  void init();
  virtual void compute_vector();
  void atomic_insertion();
  void atomic_insertion_full();
  double energy(int, int, tagint, double *);
  double energy_full();
  void toggle_intramolecular(int);
  double memory_usage();

 private:
  int exclusion_group,exclusion_group_bit;
  int nwidom_type,nevery,seed;
  int ncycles,ninsertions;
  int regionflag;           // 0 = anywhere in box, 1 = specific region
  int iregion;              // widom region
  char *idregion;           // widom region id
  bool pressure_flag;       // true if user specified reservoir pressure
  bool charge_flag;         // true if user specified atomic charge
  bool full_flag;           // true if doing full system energy calculations

  int ensemble;             // ensemble type to determine mu_ex equation
  int tflag, pflag;         // flags for temperature and pressure
  char *id_temp, *id_press; // IDs for temperature and pressure computes

  int groupbitall;          // group bitmask for inserted atoms
  int ngroups;              // number of group-ids for inserted atoms
  char** groupstrings;      // list of group-ids for inserted atoms
  int ngrouptypes;          // number of type-based group-ids for inserted atoms
  char** grouptypestrings;  // list of type-based group-ids for inserted atoms
  int* grouptypebits;       // list of type-based group bitmasks
  int* grouptypes;          // list of type-based group types

  int widom_nmax;
  int max_region_attempts;
  double gas_mass;
  double reservoir_temperature;
  double tfac_insert;
  double wtest, wtest_sq;
  double beta_mu_ex, beta_mu_ex_sigma;
  double chemical_potential, chemical_potential_sigma;
  double beta,sigma,volume;
  double pressure,temperature,charge;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double region_xlo,region_xhi,region_ylo,region_yhi,region_zlo,region_zhi;
  double region_volume;
  double energy_stored;  // full energy of old/current configuration
  double *sublo,*subhi;
  int *local_gas_list;
  double **cutsq;
  imageint imagezero;
  double overlap_cutoffsq; // square distance cutoff for overlap
  int overlap_flag;

  double energy_intra;

  class Pair *pair;

  class RanPark *random_equal;
  class RanPark *random_unequal;

  class Atom *model_atom;

  class Compute *c_pe;

  class Fix *fixrigid, *fixshake;
  int rigidflag, shakeflag;
  char *idrigid, *idshake;
  int triclinic;                         // 0 = orthog box, 1 = triclinic

  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute widom does not (yet) work with atom_style template

Self-explanatory.

E: Compute widom region does not support a bounding box

Not all regions represent bounded volumes.  You cannot use
such a region with the fix widom command.

E: Compute widom region cannot be dynamic

Only static regions can be used with fix widom.

E: Compute widom region extends outside simulation box

Self-explanatory.

E: Compute widom molecule must have coordinates

The defined molecule does not specify coordinates.

E: Compute widom molecule must have atom types

The defined molecule does not specify atom types.

E: Atom type must be zero in fix widom mol command

Self-explanatory.

E: Compute widom molecule has charges, but atom style does not

Self-explanatory.

E: Compute widom molecule template ID must be same as atom_style template ID

When using atom_style template, you cannot insert molecules that are
not in that template.

E: Compute widom atom has charge, but atom style does not

Self-explanatory.

E: Cannot use fix widom rigid and not molecule

UNDOCUMENTED

E: Cannot use fix widom shake and not molecule

Self-explanatory.

E: Cannot use fix widom rigid and shake

UNDOCUMENTED

E: Cannot use fix widom rigid with MC moves

UNDOCUMENTED

E: Cannot use fix widom shake with MC moves

UNDOCUMENTED

E: Molecule template ID for fix widom does not exist

Self-explanatory.

W: Molecule template for fix widom has multiple molecules

The fix widom command will only create molecules of a single type,
i.e. the first molecule in the template.

E: Region ID for fix widom does not exist

Self-explanatory.

W: Compute widom using full_energy option

Compute widom has automatically turned on the full_energy option since it
is required for systems like the one specified by the user. User input
included one or more of the following: kspace, a hybrid
pair style, an eam pair style, tail correction,
or no "single" function for the pair style.

E: Invalid atom type in fix widom command

The atom type specified in the widom command does not exist.

E: Compute widom cannot exchange individual atoms belonging to a molecule

This is an error since you should not delete only one atom of a
molecule.  The user has specified atomic (non-molecular) gas
exchanges, but an atom belonging to a molecule could be deleted.

E: All mol IDs should be set for fix widom group atoms

The molecule flag is on, yet not all molecule ids in the fix group
have been set to non-zero positive values by the user. This is an
error since all atoms in the fix widom group are eligible for deletion,
rotation, and translation and therefore must have valid molecule ids.

E: Compute widom molecule command requires that atoms have molecule attributes

Should not choose the widom molecule feature if no molecules are being
simulated. The general molecule flag is off, but widom's molecule flag
is on.

E: Compute widom rigid fix does not exist

UNDOCUMENTED

E: Compute widom and fix rigid/small not using same molecule template ID

UNDOCUMENTED

E: Compute widom shake fix does not exist

Self-explanatory.

E: Compute widom and fix shake not using same molecule template ID

Self-explanatory.

E: Cannot use fix widom in a 2d simulation

Compute widom is set up to run in 3d only. No 2d simulations with fix widom
are allowed.

E: Could not find fix widom exclusion group ID

Self-explanatory.

E: Could not find fix widom rotation group ID

Self-explanatory.

E: Illegal fix widom gas mass <= 0

The computed mass of the designated gas molecule or atom type was less
than or equal to zero.

E: Cannot do Widom on atoms in atom_modify first group

This is a restriction due to the way atoms are organized in a list to
enable the atom_modify first command.

W: Compute widom is being applied to the default group all

This is allowed, but it will result in Monte Carlo moves
being performed on all the atoms in the system, which is
often not what is intended.

E: Could not find specified fix widom group ID

Self-explanatory.

E: fix widom does currently not support full_energy option with molecules on more than 1 MPI process.

UNDOCUMENTED

W: Energy of old configuration in fix widom is > MAXENERGYTEST.

This probably means that a pair of atoms are closer than the
overlap cutoff distance for keyword overlap_cutoff.

E: Compute widom put atom outside box

This should not normally happen.  Contact the developers.

E: Compute widom ran out of available molecule IDs

See the setting for tagint in the src/lmptype.h file.

E: Compute widom ran out of available atom IDs

See the setting for tagint in the src/lmptype.h file.

E: Too many total atoms

See the setting for bigint in the src/lmptype.h file.

U: Compute widom can not currently be used with fix rigid or fix rigid/small

Self-explanatory.

*/
