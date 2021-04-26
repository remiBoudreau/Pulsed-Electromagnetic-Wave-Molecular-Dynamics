# pulsed-electromagnetic-wave-molecular-dynamics
Addon module for performing born-oppenheimer molecular dynamics in the presence of a time-varying electric field  at the tddft level of theory via the GAMESS-US package.

GAMESS-US is required to be installed and available using the command 'rungms'.

To execute:  
1. Replace the GAMESS-US input file in the 'PEW_MD_INPUT' directory with your own molecule. Note it must be described at the TDDFT level of theory. See bR.inp as an example.
2. Edit the parameters in the 'PEW_MD.sh' file to those for your simulation. 
3. Navigate to the 'PEW-MD' directory and execute with './PEW_MD.sh'

For details of this code, refer to the following chapter of my [thesis](https://tspace.library.utoronto.ca/bitstream/1807/103164/1/Boudreau_Jean-Michel_Remi_202011_MSc_thesis.pdf): 
1. Necessary theory, Chapter 2.
2. Code Structure, Chapter 3.

Documentation in progress.
