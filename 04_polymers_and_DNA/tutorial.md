---
date: 2024-05-06
authors: Jan Stevens, Marieke Westendorp
goal:
time: ?? minutes
software: GROMACS 2024.1, VMD, Polyply
version: beta
---

# polymers and DNA

In this tutorial we will discuss how to generate input files to perform simulations of complex
polymer systems. We will use **Polyply** to generate both the input parameter files aswell the starting
coordinates.[^polyply] Note that **Polyply** can be used to generate input files for all sorts of
polymers, ssDNA, and carbohydrates, this in both all-atom and coarse-grained resolution.

During this tutorial we will focus on coarse-grained single stranded DNA (ssDNA) simulations using
the Martini2 forcefield. [^martini2] The workflow showcased can be applied to other types of
polymer systems

We will show how **Polyply** can be used to model polymer melt, incorporaitng the stiffness of the polymer, by means of a persistence length to generate initial chain coordinates.
Aswell as how to generate closed  To conclude, we will generate a model that resembles the ssDNA enclosed inside the capsid of a DNA virus (not modelled
here).

This tutorial is based on the Polypy wiki[^wiki], more examples on how to use the suite look here.

## Generating .itp files

Generate .itp files with **Polyply**.

List all available forcefields in **polyply**

```{execute}
polyply -list-lib
```

- `-lib`:
- `-seq`:
- `-name`:
- `-seq`:
```{execute}
polyply gen_params -lib martini2 -o ssDNA.itp -name ssDNA -seq DT5:1 DT:23 DT3:1
```

For more information on the functionality of the `gen_params` subrouting, run `polyply gen_params -h`

```text
#include "martini_v2.1-dna.itp"
#include "martini_v2.0_ions.itp"
#include "ssDNA.itp"
[ system ]
ssDNA in capsid in water
[molecules]
ssDNA 100
```

## Generating starting coordinates

### Polymer melt

multiple molecules in a box given density

usual to define persistence length

[ molecule ]
ssDNA 0 100

[ persistence_length ]
; model lp first residue last residue
WCM 1.0 0 24


- `-p`:
- `-b`:
- `-name`:
- `-dens`:
- `-o`:
```{execute}
polyply gen_coords -p topol.top -b build_file.bld -name ssDNA -dens 250 -o output.gro
```

For more information on the functionality of the `gen_coords` subrouting, run `polyply gen_coords -h`

## Circular polymers

Can also do circular polymer, showcase with circular DNA

Use a `.ig` file to generate parameters for circular polymers
```
; Circular DNA
Random 25 bp sequence
TCCCGGCGAACTTAAAGTTGTAATG2
```

```{execute}
polyply gen_params -lib martini2 -o ssDNA.itp -name ssDNA -seqf sequence.ig
```

```{execute}
polyply gen_coords -p topol.top -b build_file.bld -name ssDNA -dens 250 -o output.gro -cycles ssDNA
```

## Confined polymers


For convenients generate parameters for a polythymine piece off ssDNA.

```{execute}
polyply gen_params -lib martini2 -o ssDNA.itp -name ssDNA -seq DT5:1 DT:1000 DT3:1
```

For this long piece of ssDNA reduce the
```text
#include "martini_v2.1-dna.itp"
#include "martini_v2.0_ions.itp"
#include "ssDNA.itp"
[ system ]
ssDNA in capsid in water
[molecules]
ssDNA 1
```


```text
[ volumes ]
DT 1

[ molecule ]
ssDNA 0 1
[ sphere ]
DT 0 1002 in 10.0 10.0 10.0 10.0

```

!confinement is defined per resname per molecule

```{execute}
polyply gen_coords -p topol.top -b build_file.bld -name ssDNA -box 20 20 20 -o output.gro
```

## Visualisation

## References
[^Polyply]: Grünewald, F., Alessandri, R., Kroon, P.C. et al. Polyply; a python suite for facilitating simulations of macromolecules and nanomaterials. Nat Commun 13, 68 (2022). https://doi.org/10.1038/s41467-021-27627-4
[^Martini2]: Uusitalo, J.J., Ingólfsson, H.I., Akhshi, P., Tieleman, D.P. and Marrink, S.J., 2015. Martini coarse-grained force field: extension to DNA. Journal of chemical theory and computation, 11(8), pp.3932-3945.
[^PolyplyWiki]: https://github.com/marrink-lab/polyply_1.0/wiki
