**Create a prim.json file in a folder containing the system of interest**
```
{
  "basis" : [
    {
      "coordinate" : [ 0.000000000000, 0.000000000000, 0.000000000000 ],
      "occupant_dof" : [ "Pb", "Na" ]
    }
  ],
  "coordinate_mode" : "Fractional",
  "description" : "Pb1",
  "lattice_vectors" : [
    [ 0.000000000000, 2.475000000000, 2.475000000000 ],
    [ 2.475000000000, 0.000000000000, 2.475000000000 ],
    [ 2.475000000000, 2.475000000000, 0.000000000000 ]
  ],
  "title" : "FCC"
}
```

Enumerate supercell data
```
casm enum --method ConfigEnumAllOccupations --min 16 --max 16 --filter 'eq(comp(a),0.125)'
```