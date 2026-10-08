#ifndef FILEOUTPUT_H_
#define FILEOUTPUT_H_

/**
 * @brief Output particle positions to a pdb file for visualization in tools such as OVITO.
 * @param reset 1: write new file or overwrite existing file, 0: append data
 * @param p_parameters used members: filename_pdb
 * @param p_vectors used members: r
 * @param time time stamp
 */
void record_trajectories_pdb(int reset, struct Parameters * p_parameters, struct Vectors * p_vectors, double time);

/**
 * @brief Output particle positions to xyz file
 * @param reset 1: write new file or overwrite existing file, 0: append data
 * @param p_parameters used members: filename_xyz
 * @param p_vectors used members: r
 * @param time time stamp
 */
void record_trajectories_xyz(int reset, struct Parameters * p_parameters, struct Vectors * p_vectors, double time);

/**
 * @brief Save a restart file
 * 
 * @param p_parameters 
 * @param p_vectors 
 */
void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Load a restart file
 * 
 * @param p_parameters 
 * @param p_vectors 
 */
void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Output one line of thermodynamic data (step, time, energies, temperature) to a csv file.
 * @param reset 1: write new file or overwrite existing file (and write the header line), 0: append data
 * @param p_parameters used members: filename_csv
 * @param step time step number
 * @param time time stamp
 * @param Epot potential energy
 * @param Ekin kinetic energy
 * @param T Kinetic Temp
 */
void record_to_csv(int reset, struct Parameters * p_parameters, size_t step, double time, double Epot, double Ekin, double T, const double *T_type);

#endif /* FILEOUTPUT_H_ */
