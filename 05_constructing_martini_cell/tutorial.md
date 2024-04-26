---
date: 06/05/2024
authors: Jan Stevens, Marieke Westendorp
goal:
time: __ minutes
software: GROMACS 2024.1
optional software: VMD, Xmgrace
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

### Creating mask

Merge `.GRO`  files
```sh {execute}
bentopy grocat chromosome.gro:CHROM membrane.gro:MEM chromosome_membrane.gro
```

Convert `.GRO` files to voxel mask
```sh {execute}
bentopy mask --interactive --debug-labels debug-labels.gro chromosome_membrane.gro mask.npz
```

-> view labels in vmd
-> select ID: -1

### Pack the cell

Explain the `input.json`.

```sh {execute}
bentopy pack --rearrange --seed 1234 --rotations 3 cytosol_input.json
```

Explain that the result of bentopy is a placements list -> instanced based file format which
contains all the information about the packed system.

### Render the cell to .GRO

Convert the placements list to a `.GRO` file

```sh {execute}
bentopy render -t topol.top outputs/workshop_cell_placements.json cytosol.gro
```

## Assemble Cell

Merge `.GRO`  files
```sh {execute}
bentopy grocat chromosome_membrane.gro cytosol.gro:CYT -o cell.gro
```

Setup topology file.
