# GeoCentroidAnalysis
Compares structures and sorts out dublicates by geeomtric centroid and E_tot analysis

During my computational projects I had to analyse hundrets of conformers with one explicit water molecule around a central fragment created by Davor Sakic's "kick" procedure: http://sw.pharma.hr/kick/

This skript compares all (solution, SMD(H2O)) optimized structures and sorts out possible doublets as well as mirror structures by comparing the geometric centroids and by using a energy cutoff value.

Copies the 50 best by E_tot to a subfolder.

Feels free to reuse some code snippets.
