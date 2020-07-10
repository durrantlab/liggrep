# Copyright 2020 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy
# of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations
# under the License.

import argparse
import os
import glob
import json
from .. import output
from .. import tests
import sys
from . import validate


def get_parameters(args):
    """
       Parse the name of the receptor and directory containing many docked
       ligands from the command line.

       :param args: A list of arguments from the command line.
       :type args: list
       :returns: The parameters, properly processed, with the specified mode.
       :rtype: dict
    """

    # First check if running in test mode.
    if "--test" in args or "-t" in args:
        tests.start()
        sys.exit(0)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
        description="""
LigGrep 1.0.0, a free, open-source program for identifying docked poses that
participate in user-specified receptor/ligand interactions. As input, LigGrep
accepts a protein receptor file (PDB, PDBQT), a directory containing many
docked-compound files (PDB, PDBQT, SDF), and a list of user-specified filters
(JSON). It evaluates each docked pose and outputs the names of the compounds
with poses that pass all filters.""",
        epilog="""
EXAMPLES OF USE:

1. Prepare a virtual library and save all 3D models to a single SDF file in
   the present directory:

python liggrep.py ./liggrep/examples/receptors/receptor.pdb \\
    "./liggrep/examples/ligands/sdf/*.sdf" ./liggrep/examples/filters.json

2. Assign bond orders to PDBQT-formatted ligands using SMI files in the same
   directory:

python liggrep.py ./liggrep/examples/receptors/receptor.pdbqt \\
    "./liggrep/examples/ligands/pdb/*.pdb" ./liggrep/examples/filters.json \\
    --mode SMILES

3. Or use Open Babel to assign bond orders:

python liggrep.py ./liggrep/examples/receptors/receptor.pdbqt \\
    "./liggrep/examples/ligands/pdb/*.pdb" ./liggrep/examples/filters.json \\
    --mode OPENBABEL --babel_exec /usr/local/bin/obabel

4. By default, LigGrep saves the names of the poses that pass all filters to
   "output.txt". You can optionally specify a different output file:

python liggrep.py ./liggrep/examples/receptors/receptor.pdb \\
    "./liggrep/examples/ligands/sdf/*.sdf" ./liggrep/examples/filters.json \\
    --file other_output.txt

5. You can tell LigGrep to output why each ligand/pose is accepted or
   rejected:

python liggrep.py ./liggrep/examples/receptors/receptor.pdb \\
    "./liggrep/examples/ligands/sdf/*.sdf" ./liggrep/examples/filters.json \\
    --verbose

6. By default, LigGrep runs in serial mode, meaning it runs on only one
   processor. You can also use multiple processors to speed the process.
   Requesting -1 processors means all processors will be used.

python liggrep.py ./liggrep/examples/receptors/receptor.pdb \\
    "./liggrep/examples/ligands/sdf/*.sdf" ./liggrep/examples/filters.json \\
    --job_manager multiprocessing --num_processors -1
""",
    )

    parser.add_argument(
        "receptor", type=str, help="PDBQT file containing receptor information."
    )
    parser.add_argument(
        "ligands",
        type=str,
        help="Directory of PDBQT files containing the docked-ligand files.",
    )
    parser.add_argument(
        "filters",
        type=str,
        help="JSON file containing filters, which are formatted as a list of dictionaries.",
    )
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        default="NONE",
        help="Optional user-specified bond-order mode. If OPENBABEL, LigGrep will assign bond orders using the open-babel "
        + "executable specified via the --babel_exec parameter. If SMILES, LigGrep will use "
        + "SMILES files to assign bond orders. These files must be in the same directory as "
        + "the docked-ligand files, and they must be similarly named (except with the .smi "
        + "extension). If NONE, LigGrep will not attempt to assign bond orders. The docked ligands "
        + "must either be in a format that itself specifies bond orders (e.g., SDF), or the "
        + "user-specified filters must not depend on bond orders. Default: NONE",
    )
    parser.add_argument(
        "-o",
        "--babel_exec",
        type=str,
        help="Optional path to the OpenBabel executable.",
    )
    parser.add_argument(
        "-f",
        "--file",
        type=str,
        default="output.txt",
        help='The name of the file were LigGrep analysis should be saved. Defaults to "output.txt".',
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_const",
        const=True,
        default=False,
        help="Indicate why molecules and poses are rejected (standard output).",
    )
    parser.add_argument(
        "--num_processors",
        "-p",
        type=int,
        metavar="N",
        default=1,
        help="Number of processors to use for parallel calculations. "
        + "Defaults to 1.",
    )
    parser.add_argument(
        "--job_manager",
        type=str,
        default="serial",
        choices=["serial", "multiprocessing", "mpi"],
        help="Determine what style of multiprocessing to use: serial, "
        + "mpi, or multiprocessing. Serial will override the "
        + "num_processors flag, forcing it to be one. MPI mode "
        + "requires mpi4py 2.1.0 or higher. Defaults to serial."
        # Old advice. Not sure still valid:
        # "and should be executed "
        # + "as: mpirun -n $NTASKS python -m mpi4py liggrep.py "
        # + "...-settings...",
    )
    parser.add_argument(
        "-t",
        "--test",
        action="store_const",
        const=True,
        default=False,
        help="Run optional tests to verify that code updates don't break functionality.",
    )
    parser.add_argument(
        "-i",
        "--internal_test",
        action="store_const",
        const=True,
        default=False,
        help="(Optional parameter used for internal testing.)",
    )
    params = vars(parser.parse_args(args))
    orig_params = {k: params[k] for k in params.keys()}

    if not params["internal_test"]:
        print("\nPARAMETERS:\n")
        print(json.dumps(params, indent=4))
        print("")

    params["debug_msg"] = ""
    params["debug_warn"] = ""
    params["debug_error"] = ""

    params = expand_and_validate_parameters(params)

    return params, orig_params


def expand_and_validate_parameters(params):
    """
       Process and load the receptor, directory containing many docked
       ligands, and filters. Ensure all necessary information is
       provided for the specified mode.

       :param params: A dictionary of the filters, with the specified mode.
       :type params: dict
       :returns: The parameters, properly processed and loaded.
       :rtype: dict
    """

    # Make sure receptor file exists
    if not validate.file_exists(params["receptor"], params):
        return params

    # Get the list of ligand files
    if os.path.isdir(params["ligands"]):
        output.msg(
            "The --ligands parameter is a directory, not a file mask. LigGrep will use "
            + 'all files in the specified directory that end in ".pdb", ".pdbqt", or ".sdf".',
            params,
        )
        ligand_filenames = (
            glob.glob(params["ligands"] + "/*.pdb")
            + glob.glob(params["ligands"] + "/*.pdbqt")
            + glob.glob(params["ligands"] + "/*.sdf")
        )
    else:
        # it must be a file mask like *.out
        ligand_filenames = glob.glob(params["ligands"])

    if not validate.found_ligands(ligand_filenames, params):
        return params

    params["ligands"] = ligand_filenames

    # Go through the ligands and collect all the extensions
    ligand_filename_exts = set([])
    for filename in ligand_filenames:
        _, file_extension = os.path.splitext(filename.lower())
        ligand_filename_exts.add(file_extension)
    params["ligand_exts"] = ligand_filename_exts

    # Get the filters
    if not validate.file_exists(params["filters"], params):
        return params

    with open(params["filters"]) as filter_data:
        params["filters"] = json.load(filter_data)

    # For each of the filters, if no exclude is specified, set it to false.
    for i in range(len(params["filters"])):
        if not "exclude" in params["filters"][i]:
            params["filters"][i]["exclude"] = False

    # Validate each filter.
    for filtr in params["filters"]:
        if not validate.proper_filter_keys(filtr, params):
            return params

    # Make sure the specified mode works.
    params["mode"] = params["mode"].upper()
    if not validate.obabel_specified_if_needed(params):
        return params

    if not validate.smi_files_exist_if_needed(params):
        return params

    if not validate.no_sdf_in_openbabel_smiles_mode(params):
        return params

    return params
