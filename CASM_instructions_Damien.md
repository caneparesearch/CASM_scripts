# CASM v1.2.0 instructions
This file provides instructions for using CASM v1.2.0, which contains some unresolved issues.

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
*Skip this step if you want to import DFT data from somewhere else.*  
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

Then set chemical references for formation energy calculation using ```casm ref --set-auto```.

## 5. (Optional) Import VASP data
If you have pre-calculated VASP data, you can import the structures into CASM. Firstly run the vasp_relax_report.py code above to generate the properties.calc.json file for each calculation. Then write a list of the path of the properties.calc.json files using ```find $(readlink -f .) -name "properties.calc.json" > DFT_data```
Create a settings.json file:
```
{
  "data" : {
    "copy_additional_files" : true,
    "copy_structure_files" : true,
    "import_properties" : true,
    "overwrite" : true
  },
  "mapping" : {
    "cost_tol" : 0.000010000000,
    "fix_lattice" : false,
    "fix_volume" : false,
    "ideal" : false,
    "k_best" : 1,
    "lattice_weight" : 0.500000000000,
    "max_va_frac" : 0.500000000000,
    "max_vol_change" : 0.300000000000,
    "min_va_frac" : 0.000000000000,
    "primitive_only" : false,
    "robust" : false,
    "strict" : false
  }
}
```
Then run the following command:
```
casm import --batch DFT_data -s settings.json
```
CASM will attempt to map your files into a configuration. Ensure that the properties.calc.json file is written correctly.

## 6. Fit cluster expansion
Create the basis set in basis_sets/bset.default/bspecs.json. The recommended value for pairwise basis functions is the radius of the sphere that completely inscribes the _largest_ supercell you have in your training data. 
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

If casm gives an error saying eci.json does not exist, create the eci using:
```
casm settings --new-eci default
```

## 7. Monte Carlo calculations
Copy casm_learn_input to a new directory "monte" and move to the directory. Create a monte.json file:
```
{
{
    "comment" : "Built from example",
    "debug" : false,
    "ensemble" : "grand_canonical",
    "method" : "metropolis",
    "model" : {
      "formation_energy" : "formation_energy"
    },
    "supercell" : [
      [24, 0, 0],
      [0, 24, 0],
      [0, 0, 24]
    ],
    "data" : {
      "sample_by" : "pass",
      "sample_period" : 1,
      "min_pass" : 100,
      "max_pass" : 100,
      "confidence" : 0.95,
      "measurements" : [
        {
          "quantity" : "formation_energy",
          "precision" : 1e-3
        },
        {
          "quantity" : "potential_energy",
          "precision" : 1e-3
        },
        {
          "quantity" : "clex_hull_dist(casm_learn_input,comp)",
          "precision" : 1e-3
        },
        {
          "quantity" : "atom_frac"
        },
        {
          "quantity" : "site_frac"
        },
        {
          "quantity" : "comp",
          "precision" : 1e-3
        },
        {
          "quantity" : "comp_n"
        }
      ],
      "storage" : {
        "write_observations" : false,
        "write_trajectory" : false,
        "output_format" : ["csv", "json"]
      }
    },
    "driver" : {
      "dependent_runs": false,
      "mode" : "incremental",
      "motif" : {
        "configname" : "auto"
      },
      "initial_conditions" : {
        "param_chem_pot" : {
          "a" : -0.6,
          "b" : 0
        },
        "temperature" : 5,
        "tolerance" : 0.001
      },
      "final_conditions" : {
        "param_chem_pot" : {
          "a" : 0.7,
          "b" : 0
        },
        "temperature" : 5,
        "tolerance" : 0.001
      },
      "incremental_conditions" : {
        "param_chem_pot" : {
          "a" : 0.1,
          "b" : 0
        },
        "temperature" : 0,
        "tolerance" : 0.001
      }
    }
  }
```
The above file is for low T to verify that GCMC can reproduce the convex hull. Set $\mu$ by changing the initial conditions and final conditions of "a". Ensure the $\mu$ range can cover the whole composition range. Ensure that Monte Carlo does not produce any structure below the hull.

To obtain the $\mu$ range of each tangent line in the convex hull, change the setting ```"configname" : "auto"``` to ```"configname" : "restricted_auto"```. This ensures that the supercell will be constructed from the structures in your training set. After the calculation is complete, the results.json file contains the results at each value of $\mu$. From this, you can choose the range of $\mu$ which corresponds to each region in the 0K convex hull. For example, $\mu$=-0.6 to -0.4 corresponds to x=0 to 0.25, etc.

We follow the algorithm in this [paper](https://arxiv.org/pdf/2309.11761.pdf). To obtain the grand canonical potential energy at a specific T and $\mu$, an integration path must be constructed. This means starting from the ground state, raise the temperature at intervals of $\Delta T$ at constant $\mu$ (T_up). this is done by changing the conditions in the monte.json file accordingly:
```
"initial_conditions" : {
	"param_chem_pot" : {
	  "a" : -0.30,
	  "b" : 0
	},
	"temperature" : 5,
	"tolerance" : 0.001
	},
"final_conditions" : {
	"param_chem_pot" : {
	  "a" : -0.30,
	  "b" : 0
	},
	"temperature" : 705,
	"tolerance" : 0.001
	},
"incremental_conditions" : {
	"param_chem_pot" : {
	  "a" : 0,
	  "b" : 0
	},
	"temperature" : 5,
	"tolerance" : 0.001
```
To obtain the grand canonical free energies, we perform free energy integration:

$\beta\Phi(T,\mu)=\beta_0\Phi(T_0,\mu)+\int_{\beta_0}^{\beta} [E-\mu x] d\beta$

$\Phi$ is the grand canonical potential energy, $\Phi = E-TS-\mu x$

$E-\mu x$ is the <potential_energy> in results.json

Next, at each temperature interval, fix temperature and vary $\mu$ across the range determined before. Do both $\mu$ _up and $\mu$ _down scan. This means that if the $\mu$ range is from 0.2 to 0.4, do 2 separate scans from 0.2-0.4, and from 0.4-0.2. **IMPORTANT**: the starting structure must be from the final structure at the end of the T_up procedure. The starting structure for both $\mu$ _up and $\mu$ _down is different as the ground state at the lower $\mu$ and upper $\mu$ is different! This is done by changing this part of the monte.json file to the path which contains the final structure:
```
"motif" : {
        "configdof" : "../../T_up_0/chempot_-0.60/conditions.1/final_state.json"
      },
```
Set the initial and final conditions accordingly and run. Integrate both sets of data:

$\Phi(T,\mu)=\phi(T,\mu_0)-\int_{\mu_0}^{\mu} x d\mu$

Resolve hysteresis by taking the lower envelope of $\Phi$ at each $\mu$ value. Now you can plot the grand canonical energies on a single plot to determine the phase boundaries. The phase boundaries can be identified from the intersections of $\Phi$ for each phase in the $\mu$ space and from the discontinuities in $x$ vs. $\mu$ plots.
