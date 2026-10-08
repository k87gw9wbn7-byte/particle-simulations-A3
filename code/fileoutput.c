#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "memory.h"
#include "structs.h"

static FILE *checked_open(const char *name, const char *mode)
{
  FILE *f=fopen(name,mode);
  if (!f) { perror(name); exit(EXIT_FAILURE); }
  return f;
}

// Write the particle positions to a pdb file, a text format that visualization
// programs such as VMD and OVITO read directly.
// The filename (without extension) is given by p_parameters->filename_pdb.
// If reset = 1 the data is written to the file deleting data it possibly contained.
// If reset = 0 the data is appended.
void record_trajectories_pdb(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
{
  FILE *fp_traj;
  char filename[1029];
  double rs = p_parameters->rescale_output;

  snprintf(filename, sizeof(filename), "%s%s", p_parameters->filename_pdb, ".pdb");
  if (reset == 1)
  {
    fp_traj = checked_open(filename, "w");
  }
  else
  {
    fp_traj = checked_open(filename, "a");
  }

  fprintf(fp_traj, "MODEL\n");
  fprintf(fp_traj, "REMARK TIME = %f\n", time);
  fprintf(fp_traj, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-10s%-3s\n", rs*p_parameters->L.x, rs*p_parameters->L.y, rs*p_parameters->L.z, 90.0, 90.0, 90.0, "P 1", "1");
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    /* C and O are visualization labels for A and B, not chemical atoms. */
    char element = p_vectors->type[i] == 0 ? 'C' : 'O';
    fprintf(fp_traj, "HETATM%5u  %c   UNK A   1    %8.3f%8.3f%8.3f  1.00  0.00           %c\n", (unsigned int)(i+1) % 100000, element, rs*p_vectors->r[i].x, rs*p_vectors->r[i].y, rs*p_vectors->r[i].z, element);
  }
  fprintf(fp_traj, "ENDMDL\n");

  fclose(fp_traj);
}

// Write the particle positions to a xyz file
// The filename (without extension) is given by p_parameters->filename_xyz.
// If reset = 1 the data is written to the file deleting data it possibly contained.
// If reset = 0 the data is appended.
void record_trajectories_xyz(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
{
  FILE *fp_traj;
  char filename[1029];
  double rs = p_parameters->rescale_output;

  snprintf(filename, sizeof(filename), "%s%s", p_parameters->filename_xyz, ".xyz");
  if (reset == 1)
  {
    fp_traj = checked_open(filename, "w");
  }
  else
  {
    fp_traj = checked_open(filename, "a");
  }

  fprintf(fp_traj, "%lu\n", p_parameters->num_part);
  fprintf(fp_traj, "time = %f\n", time);
  struct Vec3D *r = p_vectors->r;
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    fprintf(fp_traj, "  %c        %10.5f %10.5f %10.5f\n", p_vectors->type[i]==0 ? 'C' : 'O', rs*r[i].x, rs*r[i].y, rs*r[i].z);
  }

  fclose(fp_traj);
}

// Save the state of the simulation to a binary restart file: the number of
// particles followed by the position, velocity and force arrays. Parameters,
// types and topology are NOT stored; they are set up again when restarting
// (see main.c), so the restart file stays valid when settings change.
void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
  FILE* p_file = checked_open( p_parameters->restart_out_filename, "wb");
  size_t num_part = p_parameters->num_part;
  size_t sz = num_part*sizeof(struct Vec3D);

  fwrite(&num_part, sizeof(size_t), 1, p_file);
  fwrite(p_vectors->r, sz, 1, p_file);
  fwrite(p_vectors->v ,sz, 1, p_file);
  fwrite(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}

// Keep the A2 binary layout (size_t, r, v, f); require the configured count
// to match. Types, chains, box and RNG state are not in this legacy format.
void load_restart(struct Parameters *p, struct Vectors *v)
{
  FILE *f=checked_open(p->restart_in_filename,"rb");
  size_t n;
  if (fread(&n,sizeof(n),1,f)!=1 || n!=p->num_part || n>v->size) {
    fprintf(stderr,"Restart bead count does not match settings.\n"); exit(EXIT_FAILURE);
  }
  if (fread(v->r,sizeof(Vec3D),n,f)!=n || fread(v->v,sizeof(Vec3D),n,f)!=n ||
      fread(v->f,sizeof(Vec3D),n,f)!=n) {
    fprintf(stderr,"Truncated restart file.\n"); exit(EXIT_FAILURE);
  }
  fclose(f);
}

void record_to_csv(int reset, struct Parameters *p_parameters, size_t step, double time, double Epot, double Ekin, double T, const double *T_type)
{
  char filename[1029];
  
  snprintf(filename, sizeof(filename), "%s%s", p_parameters->filename_csv, ".csv");

  FILE *fp = checked_open(filename, reset == 1 ? "w" : "a");
  if (reset == 1)
  {
    fprintf(fp, "step,time,Epot,Ekin,Etot,T");
    for (int t = 0; t < TYPES; t++)
      fprintf(fp, ",T_type%d", t);
    fprintf(fp, "\n");
    fclose(fp);
    return;
  }
  fprintf(fp, "%lu,%.10g,%.10g,%.10g,%.10g,%.10g",
          (long unsigned)step, time, Epot, Ekin, Epot + Ekin, T);
  for (int t = 0; t < TYPES; t++)
    fprintf(fp, ",%.10g", T_type[t]);
  fprintf(fp, "\n");
  fclose(fp);
}

