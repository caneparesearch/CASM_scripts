# CASM v1.3.0 instructions
IMPORTANT: As of CASM v1.3.0, importing pre-calculated structures is not working. Do not use CASM v1.3.0 if you are importing pre-calculated structures.

## 1. Create a prim.json file in a folder containing the system of interest
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
## 2. Create composition axes
```
casm composition --calc
casm composition --select 1
```

## 3. Enumerate supercell data
--min and --max sets minimum and maximum supercell size
--filter includes filtering according to some criteria
```
casm enum --method ConfigEnumAllOccupations --max 8
```
To enumerate for a specific composition:
```
casm enum --method ConfigEnumAllOccupations --min 16 --max 16 --filter 'eq(comp(a),0.125)'
```

## 4. Perform DFT calculations
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
Now we can use pymatgen to generate the VASP input files in the training_data folder and transfer them to any cluster to submit the jobs.
Example of write_vasp_input.py script to be used in the training_data folder:
```
import glob, json, os
from pymatgen.io.vasp.sets import MITRelaxSet
from pymatgen.core.structure import Structure

incar_mod = {"NCORE":16, "ALGO":"Fast", "ISMEAR":1, "ISPIN":1, "ISYM":0, "ICHARG":2, "SIGMA":0.2, "EDIFF":1e-06, "SYMPREC":1e-09,
    "LREAL":False,"LSCALAPACK":False}

filenames = glob.glob("*/*/structure.json")

for i in filenames:
    with open(i, "r") as f:
        data = json.load(f)
    structure = Structure(lattice=data["lattice_vectors"],species=data["atom_type"],coords=data["atom_coords"],coords_are_cartesian=True)
    dir_name = os.path.dirname(i)
    relax = MITRelaxSet(structure,user_incar_settings=incar_mod, user_kpoints_settings={"grid_density":5000},force_gamma=True)
    relax.write_input(output_dir=dir_name, make_dir_if_not_present=False,include_cif=True)
```
Once calculations are done, transfer them back to the same folder. At this point, vasp_relax_report is broken so we can manually write a vasp_relax_report.py script to get the properties.calc.json file for each calculation:
```
from pymatgen.io.vasp.outputs import Outcar, Vasprun
import json, glob, os
import warnings
warnings.filterwarnings("ignore")

filenames = glob.glob("*/*/vasprun.xml.relaxed.gz")
for i in filenames:
    dir_name = os.path.dirname(i)
    print(f"Begin VASP relax report for {dir_name}")
    
    outcar = Outcar(f"{dir_name}/OUTCAR.relaxed.gz")
    vasprun = Vasprun(f"{dir_name}/vasprun.xml.relaxed.gz")
        
    data_new = {}
    data_new["global_properties"] = {} 
    data_new["global_properties"]["energy"] = {}
    data_new["global_properties"]["energy"]["value"] = outcar.final_energy_wo_entrp
    structure = vasprun.final_structure
    data_new["atom_type"] = [str(x) for x in structure.species]
    data_new["mol_type"] = [str(x) for x in structure.species]
    data_new["coordinate_mode"] = "direct"
    data_new["lattice_vectors"] = structure.lattice.matrix.tolist()
    data_new["atom_coords"] = structure.frac_coords.tolist()
    data_new["mol_coords"] = structure.frac_coords.tolist()
    
    with open(f"{dir_name}/calctype.default/properties.calc.json", "w") as f:
        json.dump(data_new, f, indent=4)
```
Use ```casm update``` to update the configurations to the master list.

## 5. Fit cluster expansion
Create the basis set in basis_sets/bset.default/bspecs.json. The recommended value for pairwise basis functions is the radius of the sphere that completely circumscribes the _largest_ supercell you have in your training data. In other words, the longest distance possible in the _largest_ supercell divided by 2.
```
{
  "basis_function_specs" : {
    "global_max_poly_order": 4,
    "dof_specs": {
      "occ": {
        "site_basis_functions" : "chebychev"
      }
    }
  },
  "cluster_specs": {
    "method": "periodic_max_length",
    "params": {
      "orbit_branch_specs": {
        "2" : {"max_length" : 12.0000},
	"3" : {"max_length" : 7.00000},
	"4" : {"max_length" : 6.00000}
      }
    }
  }
}
```
Create the basis set
```
casm bset -u
```
Create a new directory 'fit' and move into the directory. Select all configurations.
```
casm select --set is_calculated -o train
```
It is recommeneded to remove those structures with unstable relaxations for the cluster expansion model. To do so, go into the train file and deselect them. The list of unstable structures can be found in the reports folder. Then select these structures for fitting.
```
casm query -c train -k formation_energy corr -o casm_learn_input
```
Create a fit.json file with parameters for fitting
```
{
    "estimator": {
        "method": "Lasso",
        "kwargs": {
            "alpha": 0.00001,
            "max_iter": 1000000
        }
    },
    "feature_selection": {
        "method": "SelectFromModel",
        "kwargs": null
    },
    "problem_specs": {
        "data": {
            "y": "formation_energy",
            "X": "corr",
            "kwargs": null,
            "type": "selection",
            "filename": "casm_learn_input"
        },
        "weight": {
            "method": "wCustom"
        },
        "cv": {
            "penalty": 0,
            "method": "LeaveOneOut"
        }
    },
    "n_halloffame": 25
}
```
We are using custom weights, so put a weight field into casm_learn_input which contains all the weights. This function is implemented in tune_weight.py
* tune_weight.py enables you to change the weight of a specific structure
* energy.py {composition} enables you to query the error of the fit at a specific composition
* optimize_weight.py is a script to optimize a set of weights using Bayesian Optimization

## 6. Monte Carlo calculations
