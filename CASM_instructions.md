# 1. Create a prim.json file in a folder containing the system of interest
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
# 2. Create composition axes
```
casm composition --calc
casm composition --select 1
```

# 3. Enumerate supercell data
--min and --max sets minimum and maximum supercell size
--filter includes filtering according to some criteria
```
casm enum --method ConfigEnumAllOccupations --max 8
```
To enumerate for a specific composition:
```
casm enum --method ConfigEnumAllOccupations --min 16 --max 16 --filter 'eq(comp(a),0.125)'
```

# 4. Perform DFT calculations
Firstly create an arbritary calc.json file at training_data/settings/calctype.default/
```
{"software":"vasp",
"method":"relax"}
```
Select all configurations and set up structure files.
```
casm select --set-on -c ALL
casm-calc --setup
```
Second command will give an error, but this is okay as it will create all the structures in the training_data file
