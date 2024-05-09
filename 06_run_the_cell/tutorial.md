---
date: 2024-05-06
authors: Jan Stevens, Marieke Westendorp
goal: >
  Construct small model of a cell with a chromosome, an envelope, and cytosolic proteins and metabolites.
time: ?? minutes
software: GROMACS 2024.1, polyply, TS2CG, bentopy
optional software: VMD
version: beta
---

# Simulating a Martini Cell Model

Download tarball from:
https://drive.google.com/file/d/1BsWglcEDCO9TxwkZ5tdyVXMVaj-CoZXT/view?usp=sharing

## Compiling GROMACS

```
wget https://ftp.gromacs.org/gromacs/gromacs-2024.1.tar.gz
tar xfz gromacs-2024.1.tar.gz
cd gromacs-2024.1
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_DOUBLE=on -DCMAKE_INSTALL_PREFIX=...
make -j16
sudo make install
source .../bin/GMXRC
```

Replace the `...` by a path in your local directory.

## Setting up a conda environment

```
module load gcc anaconda3_cpu
conda init
```

disconnect from Delta

```
conda create -n martini_tutorial
conda activate martini_tutorial
```
Finally, install `MDAnalysis.`

```
pip3 install MDAnalysis
```


## Setting up the simulation

```
mkdir -p em-vac
gmx_d grompp -f mdp_files/em.mdp -c cell.gro -p topol.top -o em-vac/em.tpr -maxwarn 1
gmx_d mdrun -nt 16 -v -deffnm em-vac/em -pin on
```

```
mkdir -p solvate
gmx_d solvate -cp vacuum_minimization/em.gro -cs mdp_files/water.gro -radius 0.21, -o solvate/system_solvated.gro
```

```python
u = mda.Universe(system_solvated.gro')

radius = 8.0  # Specify your radius
# Find badly placed W beads
not_near_tail_solvent = u.select_atoms(f'not( (resname W) and around {radius} (name C1 C2A C2B C3A C3B C4A C4B D2A D2B) )')

#shuffle water beads in place
num_W = len(not_near_tail_solvent)

# Shuffle water beads to easily add ions later
temp = np.copy(not_near_tail_solvent.atoms.positions)
rng.shuffle(temp[-num_W:], axis=0)
shuffle_system.atoms.positions = temp
```

Add the solvent:
```python
# Determine the number of ions to add
conc_NaCl = 0.150
magic_factor = 55.5 / (2 * conc_NaCl) / 4 + 1
num_NA = num_CL = int(num_W / magic_factor // 2)
if system_charge > 0:
    num_CL += abs(system_charge)
else:
    num_NA += abs(system_charge)
num_W -= num_NA + num_CL

# Add antifreeze -> only for Martini2!!!!
num_WF = int(num_W /  10)
num_W -= num_WF

# Write ions and solvent to topology file
with open(topology_file, "a") as file_:
    file_.write(f"W {num_W}\n")
    file_.write(f"WF {num_WF}\n")
    file_.write(f"NA {num_NA}\n")
    file_.write(f"CL {num_CL}\n")
clean_directory(working_dir)
```

## Performing energy minimization

