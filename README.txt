chem_rate_fit
Fortran code for least-squares fits of rate coefficients, initial concentrations, 
 and (less so) initial times to chemical kinetic data.
The user supplies the mechanism (a set of elementary reactions) and initial
 values for all concentrations and rate coefficients.  Flagged parameters
 are then fit.
Propagated uncertainties and correlation coefficients are reported.
The least squares routine is limdif from Minpack.  

The input file is assumed to be named "rate.dat", final output is sent to "rate.log".
