---
date: 2024-05-06
authors: Jan Stevens, Marieke Westendorp
goal:
time: minutes
software: GROMACS 2024.1
optional_software: visualization software [VMD](https://www.ks.uiuc.edu/Research/vmd), Xmgrace
version: beta
---


# Simulating a Martini 3 protein model


In this tutorial, we will create a Martini 3 coarse-grain (cg) model based on a predicted atomistic structure of a peptide deformylase ([P47352][protein_uniprot]).
After creating the cg model, we will move on to energy minimization, solvation, equilibration, and finally a short production run.

### Programs

In this tutorial, we will use the following programs.

- `martinize2` (`pip3 install vermouth`).
- `gmx` (`source .../GMXRC`).
- `xmgrace` or some other means of viewing `xvg` files.
- Common command line utilities.

### Files

We will use the Martini 3 force field.
Download and unpack the [Martini3 parameter files][martini3_parameters].

```sh {execute}
wget http://md.chem.rug.nl/images/martini_v300.zip
mkdir martini_v3.0.0
unzip martini_v300.zip -d martini_v3.0.0
rm martini_v3.0.0  # Clean up.
```

## Preparing the input structure

First, we download the structure file from the AlphaFold Protein Structure Database.
For clarity and readability in this tutorial, we store the file as `protein.pdb`.

```sh {execute}
wget https://alphafold.ebi.ac.uk/files/AF-P47352-F1-model_v4.pdb -O protein.pdb
```

Since we are using a structure from the AlphaFold database, we are already in very good shape to immediately move on to mapping the structure to a coarse-grain model.
Other models may require some clean-up prior to that step, however.
Generally, cleaning the structure comes down to removing any non-protein atoms, such as waters.

## Create a coarse-grained structure using _martinize2_.

Now we are ready to map the atomistic structure to a coarse-grain model using _martinize2_.
We use the following options:

- `-f`: path to the (cleaned) structure.
- `-o`: output path for the topology file for the coarse-grained structure.
- `-x`: output path for the coarse-grained structure.
- `-p backbone`: instruct _martinize2_ to apply position restrains to the backbone beads.
- `-ff martini3001`: the target force field is Martini 3.
- `-scfix`: apply side chain corrections. **TODO:** What does that mean?
- `-cys auto`: create cysteine bonds. **TODO:** What does auto mean here?
- `-elastic`: write elastic bonds.
- `-ef 700.0`: set the elastic bond force constant.
- `-el 0.5`: elastic bond lower cutoff.
- `-eu 0.9`: elastic bond upper cutoff.
- `-ea 0`: elastic bond decay factor _a_.
- `-ep 0`: elastic bond decay factor _p_.

**TODO:** Actually explain what these values mean.

```sh {execute}
martinize2 -f protein.pdb -o topol.top -x protein_cg.pdb \
    -p backbone -ff martini3001 \
    -scfix -cys auto -elastic \
    -ef 700.0 -el 0.5 -eu 0.9 -ea 0 -ep 0
```

Three outputs files will be generated: the coarse-grained structure file (`protein_cg.pdb`), the associated main topology file (`topol.top`), and a protein topology file (`molecule_0.itp`).
List the contents of your working directory with `ls` to see whether they have been created successfully.

Since _martinize2_ creates a rather generically named `itp` file for us, we want to rename it to avoid confusion.
To match our naming strategy for this tutorial, we will rename the file to `protein.itp` and change the internal title as well.

```sh {execute}
mv molecule_0.itp protein.itp
```

With a text editor open `protein.itp`, and substitute the `[ moleculetype ]` from `molecule_0 1` to `protein 1`.

<details>
<summary>Is there a command for this?</summary>

We can also use something like a stream editor, such as `sed` to do this task.

```sh {execute}
sed -i 's/molecule_0/protein/' protein.itp
```
</details>

## Adjust the topology file

In the main topology file (`topol.top`), we need to reference all `itp` files that are required for setting up simulations with the structures that are defined in it.
Moreover, we changed the name of the `molecule_0.itp` file to `protein.itp`, so we also make that modification in the main topology file.

Earlier, we instructed _martinize2_ to create a topology and model for the Martini 3 force field.
Now, we add the appropriate `#include`s for this force field at the top of `topol.top`.
These force field `itp`s can be found in the `martini_v3.0.0` directory.
If we list its contents, we see the following items:

```sh
ls martini_v3.0.0
```

```
martini_v3.0.0_ions_v1.itp
martini_v3.0.0.itp
martini_v3.0.0_nucleobases_v1.itp
martini_v3.0.0_phospholipids_v1.itp
martini_v3.0.0_proteins
martini_v3.0.0_small_molecules_v1.itp
martini_v3.0.0_solvents_v1.itp
martini_v3.0.0_sugars_v1.itp
```

From these, we will need the base `martini_v3.0.0.itp` file as well as _ions_ and _solvents_ `itp`s, which will become relevant as we solvate our system in later steps.
You can add these with an editor you enjoy using, such that your `topol.top` is equivalent to our example below.
Notice that we set the title of the `[ system ]` to something appropriate, and that the `[ molecules ]` section features one instance of our renamed `protein`.

```conf
#include "martini_v3.0.0/martini_v3.0.0.itp"
#include "martini_v3.0.0/martini_v3.0.0_ions_v1.itp"
#include "martini_v3.0.0/martini_v3.0.0_solvents_v1.itp"
#include "protein.itp"

[ system ]
Peptide deformylase

[ molecules ]
protein 1
```

<details>
<summary>Do it with a command!</summary>

```sh {execute}
cat > topol.top <<EOF
#include "martini_v3.0.0/martini_v3.0.0.itp"
#include "martini_v3.0.0/martini_v3.0.0_ions_v1.itp"
#include "martini_v3.0.0/martini_v3.0.0_solvents_v1.itp"
#include "protein.itp"

[ system ]
Peptide deformylase

[ molecules ]
protein 1
EOF
```

</details>

## Prepare our box

Now that we have a coarse-grained structure, we will put it through a minimization run, solvate it, apply another minimization followed by equilibration steps, before finally doing a short production simulation.
To prepare for these steps, we will configure an appropriate box based on our cg structure.

We want to specify a box that is large enough to fit our protein model without risk of it interacting with itself over the periodic boundaries.
A dodecahedron box shape is appropriate here, since it will limit the volume around the protein while still providing adequate distance from the protein's periodic images.

- `-d 1.0`: set the minimum distance from the protein to any box side to 1.0 nm.
- `-bt dodecahedron`: select the box type.

```sh {execute}
gmx editconf -f protein_cg.pdb -d 1.0 -bt dodecahedron -o protein_cg.gro
```

The `gro` file that is produced in this step contains all beads that were defined in the cg `pdb` file that served as the input and also features the box vectors that were determined.
To inspect the box vectors, you can look at the last line in the `gro` file.

```sh
tail protein_cg.gro
```

## A short energy minimization

Before we move on to solvation, we will apply a short energy minimization in vacuum.

... **TODO** explain? ...

```sh {execute}
mkdir em_vac
gmx grompp -f mdp_files/em.mdp -p topol.top -c protein_cg.gro -r protein_cg.gro -o em_vac/em_vac.tpr
gmx mdrun -v -deffnm em_vac/em_vac
```

## Solvate the structure

### Introducing water

We are ready to solvate the structure.
First, we will introduce water into the system by repeatedly placing an equilibrated water box (`water.gro`).
It is important to set a reasonable Van der Waals radius (`-radius 0.21`), since the Martini water beads are larger than atomistic waters.

```sh {execute}
mkdir sol
gmx solvate -cp em_vac/em_vac.gro -cs water.gro -radius 0.21 -o sol/sol.gro
```

Now that we have embedded the structure in water, we need to update the topology to reflect that.
We add another entry to the `[ molecules ]` section to describe the water beads that are present in the solvated structure file.
The Martini water beads we inserted are denoted by the name `W`.
We can count the number of waters we added by searching for lines in the solvated structure that contain the bead name `W`.

```sh {execute}
echo W $(grep -c W sol/sol.gro) >> topol.top  # Append the W entry.
tail -n3 topol.top  # Check the last 3 lines to see if the W entry was appended correctly.
```

The `tail` command should show that the end of `topol.top` looks as follows (note that the exact number of `W` molecules may be different in your case).

```conf
[ molecules ]
protein 1
W 1933
```

**TODO**: In my run on md25, the W count was 1927, not 1933. This may not be deterministic (without a seed??). Would it be useful to use a seed for such a tutorial?

### Adding ions

After adding this water, we will need to insert ions to neutralize the charge of our protein and bring the final salt concentration to 0.15 M NaCl.

While running `gmx genion`, you will be prompted to select a molecule group.
As you may expect, we want to place the ions in the `W` group.

```sh
gmx grompp -f mdp_files/ions.mdp -c sol/sol.gro -p topol.top -o sol/ions.tpr -maxwarn 1
gmx genion -s sol/ions.tpr -p topol.top -pname NA -nname CL -neutral -conc 0.15 -o sol/sol_neutral.gro
# When prompted, select the solvent group by entering `W`.
```

<details>
<summary>Can the group selection be automated?</summary>

Sure! You could also run the `gmx genion` command below with `echo W` before it.
For example:

```sh {execute}
gmx grompp -f mdp_files/ions.mdp -c sol/sol.gro -p topol.top -o sol/ions.tpr -maxwarn 1
echo W | gmx genion -s sol/ions.tpr -p topol.top -pname NA -nname CL -neutral -conc 0.15 -o sol/sol_neutral.gro
```

This works in many such instances.
Embrace ~~laziness~~ automation!
</details>


<details>
<summary>What is that warning about?</summary>

The system has a +2 net charge.
Gromacs warns that the results may not be valid since Ewald electrostatics are used, which can produce artifacts in non-neutral systems.
But we want to carry on anyway, since it is this net charge that we intend to neutralize in this very step.
So, we can safely ignore the warning using the `-maxwarn 1` option.
</details>

## More energy minimization

Since we have introduced new beads into the system, we want to go through a short energy minimization process.

```sh {execute}
mkdir em
gmx grompp -f mdp_files/em.mdp -p topol.top -c sol/sol_neutral.gro -r sol/sol.gro -o em/em.tpr
gmx mdrun -v -deffnm em/em
```

## An equilibration run

**TODO** ... talk about position restraints? ...

```sh {execute}
mkdir eq
gmx grompp -f mdp_files/eq.mdp -p topol.top -c em/em.gro -r sol/sol.gro -o eq/eq.tpr -maxwarn 3
gmx mdrun -v -deffnm eq/eq
```

<details>
<summary>What are these warnings about?</summary>

You can expect to see three warnings regarding the use of a Berendsen thermostat and barostat, as well as a warning about pressure coupling with absolute position restraints.
We can ignore these, since the artifacts that we are warned about are not of concern to us in this phase.
In this step, we want to equilibrate the system rapidly _before_ we move on to the actual, proper molecular dynamics run.
In that final run, we will use an `mdp` recipe that addresses these concerns by selecting a more appropriate thermostat and barostat.
The warning about pressure coupling with absolute position restraints will not be relevant anymore, since we let go of those restraints by then.
</details>

## A short production simulation

Finally, we are ready to run a proper molecular dynamics simulation.
The setup will be straight-forward and the execution is quite familiar, by now.

But first, we will take a look at the `mdp` file.

```sh {execute}
less mdp_files/md.mdp
```

You can read that we will run the simulation for 2500000 steps, with a 20 fs time step.
From this, we can quickly calculate that the simulation will sample 2500000 &sdot; 20 fs = 5&times;10<sup>7</sup> fs = 50 ns.
The number of steps per compressed frame output (`nstxout-compressed`) is set to 50000, so we will store 50 frames.
The temperature coupling is set to a velocity-rescaling protocol, and a Parrinello-Rahman barostat is used.
Since the system can be considered to be symmetrical across all axes, an isotropic pressure coupling is the appropriate choice here.
The constraints are set to `none`, since we are interested in the unconstrained dynamics, here.

```sh {execute}
mkdir md
gmx grompp -f mdp_files/md.mdp -p topol.top -c eq/eq.gro -o md/md.tpr
gmx mdrun -v -deffnm md/md
```

## Visualization

- `-pbc whole`: make molecules that are broken up by the pbc whole.
- `-center`: place the molecules in the center of the box.

```sh {execute}
mkdir viz
gmx trjconv -s md/md.tpr -f md/md.xtc -o viz/whole.xtc -pbc whole -center
```

**TODO:** To come: vmd script to easily display it using vmd.

## Analysis

For clarity, we create an analysis directory to write our output files into.

```sh {execute}
mkdir analysis
```

### RMSD

```sh {execute}
gmx rms -s md/md.tpr -f viz/whole.xtc -o analysis/rmsd.xvg
```

You will be asked to provide two selections.
The first group you select represents the atoms over which the least squares fit is carried out.
For the second group, the root mean square deviation is determined.
Since we are interested in the RMSD of our protein beads, we select the _Protein_ group (1) in both cases.

### RMSF

```sh {execute}
gmx rmsf -s md/md.tpr -f viz/whole.xtc -o analysis/rmsf.xvg
```

**TODO** Consider adding instructions to project the rmsf as the backbone color through setting it as the B-factor, for example.

[protein_uniprot]: https://www.uniprot.org/uniprotkb/P47352/entry
[martini3_parameters]: http://md.chem.rug.nl/images/martini_v300.zip
