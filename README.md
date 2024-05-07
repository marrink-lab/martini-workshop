
# Workshop: Simulating a Minimal Bacterial Cell

<img align="right" width="365" height="365" src="./figures/cell.png">

In this workshop we aim to demonstrate how the Martini forcefield and associated software ecosystem can be used to study biological systems and ultimately simulate whole-cell models. The purpose is to introduce you to the basic and some more advanced aspects of using GROMACS to perform Martini coarse-grained simulations.<br>
A prior knowledge of basic concepts in molecular dynamics and familiarity with GROMACS is required. For an explanation of how to work with GROMACS and the specification of force field and run parameters, please refer to the GROMACS [user guid](https://manual.gromacs.org/current/user-guide/index.html) and [web pages](www.gromacs.org). In addition, excellent GROMACS tutorials are available [here](https://tutorials.gromacs.org/).

This tutorial series follows a logical path, where each tutorial introduces a new element of the Martini ecosystem, ultimately combining these different elements into a model of a bacterial cell.
<br>
<div align="center">

### Tutorials

[Tutorial I: Bilayer Self-assembly](01_bilayer_self_assembly/tutorial.md)

[Tutorial II: Protein Basics](02_protein_basics/tutorial.md)

[Tutorial III: Membranes and Vesicles](03_membranes_and_vesicles/tutorial.md)

[Tutorial IV: Polymers and DNA](04_polymers_and_DNA/tutorial.md)

[Tutorial V: Constructing Martini cell model](05_constructing_martini_cell/tutorial.md)

</div>

<br>

We start, historically accurate, by introducing the Martini model for lipids, followed by a simulation of a small protein. In Tutorial III, these two models will be combined to generate a lipid vesicle containing membrane proteins. In Tutorial IV, we will use Polyply to fill the vesicle we previously generated with a long fragment of single-stranded Martini2 DNA. To conclude, we will then move on to Tutorial V, where we will merge all these elements to create a toy model of a bacterial cell.

**To begin the tutorial**: 

Start by connecting to the Delta HPC using your account:

```
ssh $USERNAME@dt-login01.delta.ncsa.illinois.edu
```
> [!WARNING]
> Make sure to replace `$USERNAME` by your own.

Once connected to the Delta HPC, request a GPU node to complete the simulations for this tutorial:

```
salloc --account=bcuj-delta-gpu --partition=gpuA100x4 --time=05:00:00 --mem=64g --gpus-per-node=1 --ntasks=1 --cpus-per-task=16 --nodes=1
```

```
salloc: Pending job allocation 3583167
salloc: job 3583167 queued and waiting for resources
salloc: job 3583167 has been allocated resources
salloc: Granted job allocation 3583167
salloc: Waiting for resource configuration
salloc: Nodes gpua028 are ready for job
```

Above is an example of the terminal output from the Delta HPC. The assigned node is mentioned on the final line; you can now SSH to it:

```
ssh gpua028
```

> [!WARNING]
> Make sure to replace `gpua028` by the node assigned to you!

After connecting to the node, the required modules have to be loaded on Delta:
```
module load openmpi/4.1.5+cuda
module load gromacs/2024.1.cuda
module load python/3.11.6
module load rust
```

Finally, you need to clone the `martini-workshop` repository into your home directory:
```
git clone https://github.com/marrink-lab/martini-workshop.git
cd martini-workshop
```
We are now ready to start the tutorial! Each tutorial will be completed in their respective folders where all the files are prepared for use.

> [!TIP]
> If one of the simulations takes too long, you can download a worked example for each tutorial [here](https://github.com/marrink-lab/martini-workshop/...).

These tutorials were written for the *"Workshop: Simulations and Visualization of a Minimal Bacterial Cell"* in Champaign-Urbana by *Jan Stevens, Marieke Westendorp, Mert Bozoflu*. Parts of this workshop were based on the [*2021 Martini online workshop*](http://cgmartini.nl/index.php/2021-martini-online-workshop).
