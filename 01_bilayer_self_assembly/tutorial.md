---
date: 2024-05-06
authors: Jan Stevens, Marieke Westendorp
goal:
time: minutes
software: GROMACS 2024.1
optional_software: visualization software [VMD](https://www.ks.uiuc.edu/Research/vmd), Xmgrace
version: beta
---

# Bilayer self-assembly

based on http://cgmartini.nl/index.php/2021-martini-online-workshop/tutorials/561-1-lipid-bilayers-i

The *Martini* coarse-grained (CG) model was initially developed for lipids[^marrink2004] [^marrink2007]. The underlying philosophy of *Martini* is to build an extendable CG model based on simple modular building blocks and to use only a few parameters and standard interaction potentials to maximize applicability and transferability. Martini 3 greatly expanded the number of possible interactions but retained this building-block approach[^souza2021]. Due to Martini's modularity, a large set of different lipid types has been parameterized. Their parameters are available with the [Martini 3 Release archive](../files/martini_v3.0.0). This tutorial will teach you how to set up lipid-water system simulations with lipids from this collection, focusing on bilayers.

The tutorial aims to create and study the properties of CG Martini models of lipid bilayer membranes. First, we will attempt to spontaneously self-assemble a POPC bilayer and check various properties that characterize lipid bilayers, such as the area per lipid, bilayer thickness, order parameters, and diffusion. Next, we will change the nature of the lipid head groups and tails and study how that influences the properties.

To start this tutorial, don't forget to navigate into the respective folder in the `martini-workshop` repository:

```bash
cd 01_bilayer_self_assembly
```

> [!TIP]
> You can download the worked examples of this tutorial [here](...). (GROMACS version 2024.1) 

### Programs

In this tutorial, we will use the following programs:

- `gmx` (`source .../GMXRC`).
- `xmgrace` or some other means of viewing `xvg` files.
- Common command line utilities.

## Preparing the starting structure

POPC (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine)

We will begin with self-assembling a dipalmitoyl-phosphatidylcholine (POPC) bilayer from a random configuration of lipids and water in the simulation box. The first step is to create a simulation box containing a random configuration of 128 POPC lipids. This can be done by starting from a file containing a single POPC molecule. You can find a martini model of POPC in the current directory.

<details>
<summary>feel free to open it in vmd and inspect</summary>

    ```sh
    Some code
    ```

</details>

The gromacs tool insert-molecules can take the `POPC.gro` this single-molecule conformation and attempt to place it in a simulation box multiple times at a random position and random orientation, each time checking that there are no overlaps between the consecutively placed molecules. For help on any gromacs tool, you can add the -h flag.

```sh {execute}
gmx insert-molecules -ci DOPC.gro -box 7.5 7.5 7.5 -nmol 128 -radius 0.21 -try 500 -o 128_DOPC.gro
```
The value of the flag -radius (default van der Waals radii) has to be increased from its default atomistic length (*0.105* nm) to a value reflecting the size of Martini CG beads.
Preparing the topology
To run a simulation we not only need a starting structure but also a topology. In the topology of our simulation, all the interactions (bonded and non-bonded) are defined. These interactions, combined with the positions of all the beads, specify the forces in our simulation and, as a consequence, what the system's dynamics will be. We will create the topology file ourselves by getting the topology for water and for the POPC lipid and organizing them as a topology (`.top`) file.

…

Note that the Martini 3 release is organized into several `.itp` files, each with the definitions for a class of molecules. For this tutorial, you won’t need all of the Martini 3 .itps, only the one where water is defined (hint: it’s a ‘solvent’) and the one where POPC is defined (hint: it’s a ‘phospholipid’). There is a third .itp you will need, which is the one with all the Martini 3 particle definitions (hint: it’s `martini_v3.0.0.itp`). The needed .itps should be placed in the tutorial directory.

Making the topology file

To create the .top file (we’ll call it `topol.top`) that describes the system topology to GROMACS, you can follow the template below. Semi-colons indicate comments, which are ignored, but hashtags aren’t: they’re preprocessing directives. Namely, it is the #include directive that allows us to bring into the .top the particle/molecule information in the .itps. Use your editor of choice (gedit/vi/other) to create a file `topol.top` and copy/paste the template below.

```sh
#include "martini_v3.0.0.itp" ; the particle definitions should be included first
#include "martini_v3.0.0_solvents_v1.itp" ; include here the relevant .itps defining the molecules to use
#include "martini_v3.0.0_phospholipids_v1.itp"

[ system ]
DOPC BILAYER SELF-ASSEMBLY

[ molecules ]
; Molecule types and their numbers are written in the order
; they appear in the structure file
DOPC 128
```
## Adding the water

Using the Gromacs tool solvate to add the water beads to the simulation box. `gmx solvate` needs to have the structure of an equilibrated water box to use as a template to fill the empty space in `128_POPC.gro`. An equilibrated water box gro file is provided in the `mdp_files` directory.

```sh {execute}
gmx solvate -cp 128_DOPC.gro -cs mdp_files/water.gro -radius 0.21 -p topol.top -o 128_DOPC_solvated.gro
```

The flag -radius value is used to reflect the size of Martini CG beads. A new file named `128_POPC_solvated.gro` is produced containing 128 lipids and added water beads.
Also if you open the file ‘topol.top’ you see that at the bottom of the file a line with is added of the form “W …”, indicating how many water beads have been added to the system.

<div align="center">
<img src="../figures/01_initial_structure.png" width="70%"/>
</div>
*__Figure 1: Bilayer__ Snapshot of the simulation after a short simulation.*

## Energy minimize the system

A short energy minimization
Now you will perform an energy minimization of the solvated system to get rid of high forces between beads that may have been placed quite close to each other. The settings file minimization.mdp is provided in the `mdp_files` directory.

```sh {execute}
gmx grompp -f mdp_files/minimization.mdp -c 128_DOPC_solvated.gro -p topol.top -o min.tpr
gmx mdrun -v -s min.tpr -c minimized.gro
```
## Running the production simulation

Now, you are ready to run the self-assembly MD simulation, by using the md.mdp run settings file and the just-minimized minimized.gro. *30* ns, or *1.5* million simulation steps at *20* fs per step, should suffice.

```sh {execute}
gmx grompp -f mdp_files/md.mdp -c minimized.gro -p topol.top -o md.tpr
gmx mdrun -v -s md.tpr -x md.xtc -c md.gro
```

This might take approximately *10* minutes on a single CPU but by default gmx mdrun will use all available CPUs on your machine. The -v option shows an estimate of the time to completion, and it is interesting to observe how the generations of desktop computers have sped up this *30* ns simulation over the years. See gmx mdrun’s help (-h) for instructions on how to tune the numbers of parallel threads gmx mdrun uses. You may want to check the progress of the simulation to see whether the bilayer has already formed before the end of the simulation.

<details>
<summary> TIP </summary>

To properly sample in an isothermal-isobaric ensemble, you should at this point switch to the Parrinello-Rahman barostat (12 ps is a typical tau-p value to use with it). The Parrinello-Rahman barostat is less robust than the Berendsen one, and may diverge (crash) if the system is far from equilibrium. As such, is usually used only on production runs, whereas Berendsen is used in preparation ones.

Because of potentially poor heat transfer across the membrane-water interface, it is recommended that the solvent and the membrane groups of particles each be coupled to their own thermostat, to prevent unequal heat accumulation. You can set that in your .mdp using the tc-grps option.

Buildup in numerical precision error may cause the system to gain overall momentum. This is undesirable because such translation will be interpreted as temperature by the thermostat, and result in an excessively cooled system. Such center-of-mass motion (COMM) is corrected using comm-mode = linear. When membranes are involved, it is also possible (even in the absence of precision errors, or when controlling for COMM) that the membrane phase gains momentum relative to the water phase. In this case, the COMM should be corrected for each phase separately, using the comm-grps option. In some applications, it may be needed to further correct for the COMM of each leaflet separately.

</details>

<div align="center">
<img src="../figures/01_bilayer.png" width="70%"/>
</div>
*__Figure 2: Bilayer__ Snapshot of the simulation after a short simulation.*

The initial and final snapshots should look similar to Fig. 1, at least if the self-assembly resulted in a bilayer in the allotted time, which is not guaranteed. You may have noticed, however, that there is relatively little water, and some part of the initial solvated box is actually devoid of water. This helps to more or less guarantee a bilayer, but does make the simulation a little less realistic. You can test by solvating the lipids with more solvent. You can do this by changing the  -maxsol flag on the  gmx solvate command.

## Visualization

Use VMD to inspect the progress of the simulation. VMD suffers from the fact that Martini bonds are usually not drawn because they are much longer than standard bonds, and the default visualization is not very informative.

The initial and final snapshots should look similar to Fig. 1, at least if the self-assembly resulted in a bilayer in the allotted time, which is not guaranteed.

In the meantime, have a closer look at the POPC section of the phospholipid .itp file. The
available bead types and their interactions are defined in `martini_v3.0.0.itp` and described in
the Martini 3 paper[^3]. Understanding how these files work will help you work with the Martini 3 model and define new molecules or refine existing models.




Figure 1 | POPC bilayer formation. A) 128 randomly localized POPC lipid molecules with water. B) After a x ns long simulation the POPC lipids have aggregated together and formed a bilayer.


When the simulation has finished, check whether you got a bilayer. If yes, check if the formed membrane is normal to the z-axis (i.e., the membrane is situated in the xy-plane). Have a look at the self-assembly process: can you recognize intermediate states, such as micelles, and do you see metastable structures, such as a water pore (water spanning the membrane interior) and/or a lipid bridge (lipid tail(s) spanning the water layer)?

## Analysis

## References
[^marrink2004]: Marrink, S. J., De Vries, A. H., and Mark, A. E. (2004) Coarse grained model for semiquantitative lipid simulations. J. Phys. Chem. B 108, 750–760. DOI:10.1021/jp036508g
[^marrink2007]: Marrink, S. J., Risselada, H. J., Yefimov, S., Tieleman, D. P., and De Vries, A. H. (2007) The MARTINI force field: coarse grained model for biomolecular simulations. J. Phys. Chem. B 111, 7812–7824. DOI:10.1021/jp071097f
[^souza2021]: Souza, P.C.T., Alessandri, R., Barnoud, J. et al. Martini 3: a general purpose force field for coarse-grained molecular dynamics. Nat Methods 18, 382–388 (2021). https://doi.org/10.1038/s41592-021-01098-3
