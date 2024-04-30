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

# Constructing a Martini Cell Model

## Make Chromosome

`genome.ig` is given and contains a random sequence of 20.000 basepairs.

```sh {execute}
polyply gen_params -lib martini2 -o chromosome.itp -name chromosome -seqf genome.ig -dsdna
```

`coords.oxdna` is given and is generated from subsampeleing the mesoscale model.

```sh {execute}
polyply gen_coords -p topol.top -box 120 120 120 -o dsDNA.gro -lib martini2 -bm_fudge 1.0 -bm_mode by-frame -mc coords.oxdna
```

- Convert constraints to bonds.
- add elastic network to `.ITP` file.

## Make Envelope

Given is

- Directory with membrane proteins containing: ATP synthase, magnesium transporter, calcium transporter, potassium transporter.
- sphere.tsi -> membrane of radius 120nm with selected inclusions
- input.str input for TS2CG -> tweaked ALP
- Martini2.LIB (standard)


### Pointilise mesh

Subsample the mesh to have enough points for all the lipids

```sh {execute}
PLM -TSfile sphere.tsi -Mashno 3
```

### Backmack mesh to Martini2

Place the lipids and proteins on the correct vertices.

```sh {execute}
PCG -str input.str -Bondlength 0.2 -LLIB Martini2.LIB -bilayerThickness 2.0 -defout topol
```

## Make cytosol

Now that we have a structure for the chromosome and the cell's envelope, we can make the first steps in bringing the structures together.
Inside of the cell envelope, we want to model the cytosol.
We will place the chromosome into this compartment, and we want to fill the remaining space with protein and metabolites.
For these steps, we have developed a tool called **bentopy**.

**bentopy** is a tool for packing molecules in spaces.
Through a spectral space reduction scheme combined with a random-placement strategy, **bentopy** can quickly set up well-stirred systems of any number of input structures.
One of its goals is to enable packing of large spaces in a user-friendly and performance-conscious manner.

Its central subcommand is **bentopy pack**, which takes an input file in which a space and a list of structures are defined, and places these structures within the space.
The placements are stored in a placement list that associates the placed structures with their rotations and positions.
A space can be specified using voxel masks that define the regions where placement is allowed.
The structures are specified as a list of entries that include a name, path, and the desired number of placed instances.

The **bentopy render** subcommand creates a structure (and topology) file from the placement list that is written by **bentopy pack**.

To create the masks that define a space, **bentopy mask** is available, which can be used to identify and select different compartments in a provided structure.
The masks are represented as compressed boolean numpy arrays (`.npz`), which provide a very flexible interface for defining these spaces.
The `mask` subcommand merely serves as a convenient interface for setting up these compressed arrays.
But since any boolean numpy array can be saved as a valid mask, there are many other possibilities for setting these up for more specialized applications.

An additional convenient tool is provided by **bentopy grocat**.
As its name suggests, this subcommand concatenates `.gro` files.
Doing this by other means requires minimal effort, but this command provides an additional ability to replace the residue names of all particles in a provided file with some other identifier.
This new residue name is set by appending `:<name>` to the file path to apply this name to.
For mesoscale models, this can be a very useful or even necessary step in distinguishing between different large sets of particles and structures.

### Creating mask

First, we merge the chromosome model with the membrane model.
This allows us to define a space based on the inside of this merged model in the next step.
After the paths to the chromosome and membrane structures, we can set the residue names for each of the beads in those files with the `:<residue name>` notation.
For all atoms in the `chromosome.gro` structure, we set the residue name to `CHROM`, for `membrane.gro` we set them to `MEM`.
The output file is specified with the `-o`/`--output` flag.

```sh {execute}
bentopy grocat chromosome.gro:CHROM membrane.gro:MEM -o chromosome_membrane.gro
```

This merged model consists of the chromosome and the envelope around it.
With **bentopy mask**, we can select the _inside_ of this compartment.
When we use notions like inside and outside in conversation, we have a very strong sense of what that means.
Yet, it can be tough to define the distinction of such compartments in a larger space in a computational workflow.
The [mdvcontainment][mdvc] package provides a powerful and robust way of making this distinction, even in periodic systems.
The **mask** subcommand wraps mdvcontainment to quickly create masks that serve as input for the packing process.

We are interested in the space inside our vesicle, but want to exclude the chromosome.
To select this space from the merged `chromosome_membrane.gro` model, we process it with the **mask** subcommand.
This command operates in an interactive mode by default,[^bentopy-mask_interactive] which allows you to load and process the structure once and select compartments while it is held in memory.

To inspect the labels that are given to different regions of the space, we can give a path to the `-b`/`--inspect-labels-path` flag.
A `.gro` file in which beads are placed at the center of each voxel that is considered by mdvcontainment, will be written to that path.
Each of the beads is named according to its compartment label.
By selecting these different labels in a molecule viewer, you can find the exact space you are interested in.

In most cases, however, this is not required, since some inner compartment is desired from a simple system.
This also happens to be the case for us, in this example.
The `--autofill` flag instructs **mask** to automatically select all _innermost_ compartments.
To understand which spaces this flag selects, you can imagine a nested doll.
In a nested doll, the `--autofill` flag will select the space inside the innermost doll.

We write the resulting mask to `mask.npz`.
As mentioned above, this is simply a compressed boolean numpy array.
Any other means of creating such an array of the expected size would work as well---the **mask** subcommand merely serves as a convenient tool for setting these up for most cases.

TODO: Containment vs mask resolution? Seems ipmortant to actually mention the mask resolution, at least, since it does come up later in the _space_ section for the input file.

```sh {execute}
bentopy mask --inspect-labels-path labels.gro --autofill chromosome_membrane.gro mask.npz
```

Though we used the `--autofill` flag to pick the inner compartment automatically, it can still be helpful and interesting to inspect the `labels.gro` file.
When loaded into a molecule viewer, the different label groups can be selected according to their atom name.
For instance, in VMD, the selection `name "-1"` will show the inside of the vesicle in our case, the envelope can be selected with `name "1"`, and `name "-2"` selects the outside.
Note that the quotes are necessary for correctly selecting negative-numbered labels.

[^bentopy-mask_interactive]: For automated applications, the `--no-interactive`
flag can be set to require the specification of all parameters through command
line arguments and with no input at runtime.

### Pack the cell

With this mask ready, we can already move on to packing that space with protein and metabolites.
In order to specify what objects we want place in which space, we need to create an input file.

Let's learn how `cytosol_input.json` is set up.

#### The _space_ section

To define the space in which the placement procedure will take place, we must set the following fields:

- _size_: a three-integer list describing the dimensions of the space in nm, (_x_, _y_, _z_) order.
  Note that this _size_ must match the nm dimensions of the mask we just set up.
- _resolution_: the size of the voxels in nm that are used to represent the space internally.
  Note that the _resolution_ must match the _mask-resolution_ that was used in setting up the mask.
- _compartments_: a list of compartments. A compartment has the following structure:
    - _id_: a name for the compartment.
    - One of the following two options:
        - _voxels_: provide the _path_ to a voxel mask (any compressed boolean numpy array `.npz`).
        - _shape_: use an analytical function to set up the internal voxel mask.
          Takes the name of a shape from the following options: "spherical", "cuboid", "none".
          (Note that this option is very likely to change.)

```json
{
	"space": {
		"size": [120, 120, 120],
		"resolution": 0.5,
		"compartments": [
			{
				"id": "cytosol",
				"voxels": {
					"path": "mask.npz"
				}
			}
		]
	},
```

#### The _output_ section

In the _output_ section, we can specify what information should be associated with the placement, and where the output should be written to.
As is explained in detail later, the output of **bentopy pack** is not the final structure.
Instead, an intermediate instance-based list of placements is written out---a file we call the _placement list_.
It takes the following information:

- _title_: the name of this placement. This will become the title of the placement list.
- _dir_: the directory to place the output files into.
- _topol_includes_: a list of files that will be included into the topology (`.top`) file that is rendered from the placement list.

```json
	"output": {
		"title": "workshop_cell",
		"dir": "outputs",
		"topol_includes": [
			"martini_v2.1-dna.itp",
			"martini_v2.0_ions.itp",
			"chromosome.itp",
			"metabolites.itp",
			"cytosolic_proteins.itp",
			"membrane_proteins.itp",
			"lipids.itp"
		]
	},
```

#### The _segments_ section

The _segments_ section lists all structures to be placed.
Each segment definition has the following fields:

- _name_: a name for the structure that corresponds to its definition in the `.itp` files defined above in the _output_._topol_includes_ field. Making sure the _name_ is correct is important for generating correct topology files after the placement procedure.
- _path_: the path to the relevant structure file (`.pdb` or `.gro`).
- _number_: an integer value that represents the number of desired placements for this structure.
- _compartments_: a list of the compartments in which this structure may be placed. In our case, this is the `"cytosol"` compartment for all of these structures.

It is common to place this section last, since it tends to be rather long.

```json
	"segments": [
		{
			"name": "syn539",
			"path": "structures/proteins/syn539_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
        // More cg protein structures...
		{
			"name": "syn6",
			"path": "structures/proteins/syn6_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "GLU",
			"path": "structures/metabolites/AminoAcids/GLU.gro",
			"number": 944,
			"compartments": [ "cytosol" ]
		},
        // More amino acids...
		{
			"name": "ATP",
			"path": "structures/metabolites/NucleotideCofactors/ATP.gro",
			"number": 3850,
			"compartments": [ "cytosol" ]
		},
        // More nucleotide cofactors...
	]
}
```

<details>
<summary>The full `cytosol_input.json` file.</summary>

```json
{
	"space": {
		"size": [120, 120, 120],
		"resolution": 0.5,
		"compartments": [ { "id": "cytosol",
				"voxels": {
					"path": "mask.npz"
				}
			}
		]
	},
	"output": {
		"title": "workshop_cell",
		"dir": "outputs",
		"topol_includes": [
			"martini_v2.1-dna.itp",
			"martini_v2.0_ions.itp",
			"chromosome.itp",
			"metabolites.itp",
			"cytosolic_proteins.itp",
			"membrane_proteins.itp",
			"lipids.itp"
		]
	},
	"segments": [
		{
			"name": "syn539",
			"path": "structures/proteins/syn539_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn804",
			"path": "structures/proteins/syn804_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn637",
			"path": "structures/proteins/syn637_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn163",
			"path": "structures/proteins/syn163_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn451",
			"path": "structures/proteins/syn451_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn305",
			"path": "structures/proteins/syn305_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn353",
			"path": "structures/proteins/syn353_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn220",
			"path": "structures/proteins/syn220_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn407",
			"path": "structures/proteins/syn407_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn809",
			"path": "structures/proteins/syn809_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn447",
			"path": "structures/proteins/syn447_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn260",
			"path": "structures/proteins/syn260_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn308",
			"path": "structures/proteins/syn308_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn95",
			"path": "structures/proteins/syn95_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn661",
			"path": "structures/proteins/syn661_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn348",
			"path": "structures/proteins/syn348_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn79",
			"path": "structures/proteins/syn79_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn142",
			"path": "structures/proteins/syn142_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn297",
			"path": "structures/proteins/syn297_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "syn6",
			"path": "structures/proteins/syn6_cg.pdb",
			"number": 50,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "GLU",
			"path": "structures/metabolites/AminoAcids/GLU.gro",
			"number": 944,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "GLN",
			"path": "structures/metabolites/AminoAcids/GLN.gro",
			"number": 597,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "ASN",
			"path": "structures/metabolites/AminoAcids/ASN.gro",
			"number": 823,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "ARG",
			"path": "structures/metabolites/AminoAcids/ARG.gro",
			"number": 1155,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "HIS",
			"path": "structures/metabolites/AminoAcids/HIS.gro",
			"number": 2111,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "VAL",
			"path": "structures/metabolites/AminoAcids/VAL.gro",
			"number": 1364,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "SER",
			"path": "structures/metabolites/AminoAcids/SER.gro",
			"number": 350,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "LEU",
			"path": "structures/metabolites/AminoAcids/LEU.gro",
			"number": 1163,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "ASP",
			"path": "structures/metabolites/AminoAcids/ASP.gro",
			"number": 89,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "GLY",
			"path": "structures/metabolites/AminoAcids/GLY.gro",
			"number": 1115,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "ALA",
			"path": "structures/metabolites/AminoAcids/ALA.gro",
			"number": 884,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "LYS",
			"path": "structures/metabolites/AminoAcids/LYS.gro",
			"number": 966,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "PRO",
			"path": "structures/metabolites/AminoAcids/PRO.gro",
			"number": 1413,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "THR",
			"path": "structures/metabolites/AminoAcids/THR.gro",
			"number": 287,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "TRP",
			"path": "structures/metabolites/AminoAcids/TRP.gro",
			"number": 838,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "ILE",
			"path": "structures/metabolites/AminoAcids/ILE.gro",
			"number": 1103,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "CYS",
			"path": "structures/metabolites/AminoAcids/CYS.gro",
			"number": 2097,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "PHE",
			"path": "structures/metabolites/AminoAcids/PHE.gro",
			"number": 278,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "TYR",
			"path": "structures/metabolites/AminoAcids/TYR.gro",
			"number": 1110,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "MET",
			"path": "structures/metabolites/AminoAcids/MET.gro",
			"number": 2587,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "FAD",
			"path": "structures/metabolites/NucleotideCofactors/FAD.gro",
			"number": 13,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "RBFL",
			"path": "structures/metabolites/NucleotideCofactors/RBFL.gro",
			"number": 2,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "Alt_NAD",
			"path": "structures/metabolites/NucleotideCofactors/Alt_NAD.gro",
			"number": 560,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "NADH",
			"path": "structures/metabolites/NucleotideCofactors/NADH.gro",
			"number": 7,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "NADPH",
			"path": "structures/metabolites/NucleotideCofactors/NADPH.gro",
			"number": 10,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "ATP",
			"path": "structures/metabolites/NucleotideCofactors/ATP.gro",
			"number": 3850,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "AMP",
			"path": "structures/metabolites/NucleotideCofactors/AMP.gro",
			"number": 110,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "ADP",
			"path": "structures/metabolites/NucleotideCofactors/ADP.gro",
			"number": 459,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "FMN",
			"path": "structures/metabolites/NucleotideCofactors/FMN.gro",
			"number": 197,
			"compartments": [ "cytosol" ]
		},
		{
			"name": "Alt_NADP",
			"path": "structures/metabolites/NucleotideCofactors/Alt_NADP.gro",
			"number": 91,
			"compartments": [ "cytosol" ]
		}
	]
}
```

</details>

#### Packing

Now that we have our input file, we can move on to actually packing the cell.
We provide the following command line options:

- `--rearrange`: order the specified segments according to the voxel volume they occupy, placing the largest structures first and the smallest structures last.
  This rearrangement is a heuristic to improve the quality of the packing
  If a number of small structures are placed first, large structures may be hard to place even when enough voxel volume is available
  If this flag is not set, the user-provided order from the input file is respected.
- `--seed 5172`: set the random number generator seed. Since **bentopy pack** places according to a random sampling scheme, we must set a seed to make the placements deterministic.
- `--rotations 3`: the number of random rotations to sample for each of the listed segments.
  If the desired placement number for some segment is set to 30, and the number of rotations is set to 3, ten placements will be made for each of the three random rotations (unless the space is full and placement is obstructed).

```sh {execute}
bentopy pack --rearrange --seed 5172 --rotations 3 cytosol_input.json
```

As we specified in the _output_ section of the input configuration, the placements that are determined during the packing procedure will be written to the `outputs` directory as an instance-based list.
In this list, each structure is matched with each rotations and all its associated translations (structure center positions).

<details>
<summary>
Writing the placements to such an intermediate file format has a couple of advantages.
</summary>

- No time is spent waiting for costly formatting and disk writes during the packing procedure.
- Representing the placement of a larger structure as a simple translation and rotation minimizes memory overhead. Storing all particles for the placed structures can take up a considerable amount of memory at large scales.
- The light-weight placement list file is trivial to transfer while a rendered structure file may be very large and slow to send around.
- Inspecting a full structure file can be slow or even prohibitive. Rendering only a part of the beads involved or even only one bead per structure instance based on the placement list can be very helpful in such cases.
- The final structure can be rendered very quickly. **bentopy render**, which renders the placement list to a structure file is a stand-alone, optimized executable written in a very fast language.
</details>

### Render the cell to `.gro`

Convert the placements list to a `.gro` file

With the `-t`/`--topology` flag, we specify a path to write a topology file to, based on the placement list.
The _name_ field for each segment entry will be used as the identifier for each structure in the topology file.
We render the structure to `cytosol.gro`.

```sh {execute}
bentopy render -t topol.top outputs/workshop_cell_placements.json cytosol.gro
```

## Assemble Cell

We will use **bentopy grocat** a final time to combine our cytosol with the chromosome-membrane structure we created earlier.
Since the `chromosome_membrane.gro` structure was already labeled with the `CHROM` and `MEM` labels, we will not set a residue name for that structure.
To properly distinguish the atoms from `cytosol.gro` in the final structure, we label them with the `CYT` residue name.
The final structure is written to `cell.gro`.

```sh {execute}
bentopy grocat chromosome_membrane.gro cytosol.gro:CYT -o cell.gro
```

Setup topology file.

[mdvc]: https://github.com/BartBruininks/mdvcontainment
