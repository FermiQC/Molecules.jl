<p align="center">
  <img src="images/logo.png" width="400" alt=""/>
</p>

<table align="center">
  <tr>
    <th>CI</th>
    <th>Coverage</th>
    <th>License</th>
  </tr>
  <tr>
    <td align="center">
      <a href=https://github.com/FermiQC/Molecules.jl/actions/workflows/CI.yml>
      <img src=https://github.com/FermiQC/Molecules.jl/actions/workflows/CI.yml/badge.svg>
      </a> 
    </td>
    <td align="center">
      <a href=https://codecov.io/gh/FermiQC/Molecules.jl>
      <img src=https://codecov.io/gh/FermiQC/Molecules.jl/branch/main/graph/badge.svg?token=NQDJ0QYLB0>
      </a> 
    </td>
    <td align="center">
      <a href=https://github.com/FermiQC/Molecules.jl/blob/main/LICENSE>
      <img src=https://img.shields.io/badge/License-MIT-blue.svg>
      </a>
    </td>
  </tr>
</table>

This package offers creating and manipulation of atomic ensamble for molecular simulations. 

- Parsing XYZ data
- Manipulation of atomic/geometric structure
- Properties such as center of mass, nuclear repulsion, nuclear dipoles and moments of inertia.

# Examples

Read a string and convert into a vector of `Atom` objects

```
julia> water = Molecules.parse_string("""
         O        1.2091536548      1.7664118189     -0.0171613972
         H        2.1984800075      1.7977100627      0.0121161719
         H        0.9197881882      2.4580185570      0.6297938830
""")

julia> water[1]
Oxygen Z = 8 M = 15.9990
 Center:    1.2091536548     1.7664118189    -0.0171613972
```

> Masses are taken from standard values, you can input the desired mass values as an extra column after the element symbol.

One can also create a `Molecule` object, where some information about the number of electrons, charge, and multiplicty are deduced.

```
julia> water = Molecule("""
         O        1.2091536548      1.7664118189     -0.0171613972
         H        2.1984800075      1.7977100627      0.0121161719
         H        0.9197881882      2.4580185570      0.6297938830
         """)
Molecule:

O    1.209153654800    1.766411818900   -0.017161397200
H    2.198480007500    1.797710062700    0.012116171900
H    0.919788188200    2.458018557000    0.629793883000


Charge: 0   Multiplicity: 1   
Nuclear repulsion:    8.8880641743

julia> Molecules.nuclear_dipole(water)
3-element Vector{Float64}:
 12.7914974341
 18.3870231709
  0.5046188773

julia> Molecules.center_of_mass(water)
3-element Vector{Float64}:
 1.2483188267782848
 1.8068607904101412
 0.020676111103880096
```

# Contribute

Checkout our issue [section](https://github.com/FermiQC/Molecules.jl/issues) for a list of desired features and milestone.