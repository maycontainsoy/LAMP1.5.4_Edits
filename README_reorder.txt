Instructions on using the reordering tools.   (4/2022)

(I assume basic familiarity with BIGSTICK and access to 
the BIGSTICK manual.)

These tools will read interaction files for BIGSTICK, 
either in 'standard' isospin format or in XPN proton-neutron format, 
and convert from one single-particle space to another single-particle
space.  This can be simply changing the order of the single-particle 
orbitals, or it can be a truncation. For truncations, there is an 
option to sum over the excluded space to get effective single-particle
energies.

To reorder/truncate with isospin format files, use reorder.f90
To reorder/truncation with XPN proton-neutron format, use reorderpn.f90
(If you need to convert from isospin to XPN format, use convertiso2pn.v2.f90 
or later version.  You can also truncated in isospin format and then 
convert to XPN format.)

Compile with fortran, e.g.,

gfortran -o reorder.x reorder.f90

The code will ask for the intial and final .sps files. (They can be in 
'iso' format even when reordering XPN files.)

It will then ask for the name of the old .int file and the name for 
the new .int file you wish to create.

If you want to sum over a core, you must make sure there are different
W weights on the core orbitals.  For example, suppose you have a file 
in the s-p-sd-pf space and want to truncate to the sd-pf space, but 
with effective single-particle energies to account for the core. 

Then the INITIAL spsdpf.sps should look like this:

iso
         10 
  0.0  0.0  0.5   0
  0.0  1.0  0.5   0
  0.0  1.0  1.5   0
  1.0  0.0  0.5   1
  0.0  2.0  1.5   1
  0.0  2.0  2.5   1
  1.0  1.0  0.5   1
  1.0  1.0  1.5   1
  0.0  3.0  2.5   1
  0.0  3.0  3.5   1 

(It won't matter what the FINAL sdpf.sps file looks like, as long as one exists.)

IMPORTANT: You must turn off auto-scaling when reordering: if the first line looks like:

   -768  -3.164 1.647 -3.948 -0.497 -1.216 7.622 2.663  40.0000   42.0000 0.3000000
you need to change -768 (which signals autoscaling) to +768 (which turns off autoscaling).
Note that scaling generally won't work if you sum over core particles; however 
in most cases that won't be relevant anyway.

Here is a sample run:

  Program to reorder single-particle orbits in .int files 
  
  Enter name of INITIAL .sps file (leave off extension )
spsdpf
  Enter name of FINAL .sps file (leave off extension )
sdpf
  Do you want to sum over core for induced s.p.e.s (y/n)?
y
  Enter max W for inclusion in core 
  (You may have to check the initial .sps file)
  Suggested max W =            0
0
  These orbits are included in the core 
   N   L  2J 
   0   0   1
   0   1   1
   0   1   3
  Enter name of initial .int interaction file 
a40hcm
  successfully opened 
  
  Enter overall scaling for TBMEs (choose 1 if you dont want to scale) 
1
  
  Kept          198  matrix elements out of          448
  Energy of exclude core =   0.749999940    
  Enter name of final .int interaction file 
a40hcmsdpf

  
   
