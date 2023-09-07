# GeoCentroidAnalysis
Compares structures and excludes dublicates by a geometric centroid and E_tot analysis

In my computational chemistry projects, I performed the analysis of hundreds of conformers, each with an explicit water molecule positioned around a central fragment using Davor Sakic's "kick" procedure, as described at http://sw.pharma.hr/kick/.

The script I developed compares all optimized structures (in solution, SMD(H2O)) and identifies potential duplicates and mirror structures. This is achieved by comparing the geometric centroids of the conformers and applying an energy cutoff value.


Copies the 10 best by E_tot to a subfolder.

Feel free to reuse some code snippets or contact me.
