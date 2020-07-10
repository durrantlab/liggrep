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

import sys

if sys.version_info[0] < 3:
    raise Exception("LigGrep is designed to run on Python3.")

from liggrep.parameters import parameters
from liggrep.ligand import get_ligand_rdkit_mols
from liggrep.ligand import matches
from liggrep import receptor
from liggrep import output
from liggrep.parameters import validate
from liggrep.parallelizer import Parallelizer
from . import scoria
from rdkit import Chem
from rdkit.Chem import AllChem
import glob
import json
import os


def main(args):
    """The main definition. Runs LigGrep.

    :param args: The command-line arguments.
    :type args: list
    :return: A dictionary containing any log, warning, and error messages.
    :rtype: dict
    """

    params, orig_params = parameters.get_parameters(args)

    if params["debug_error"] != "":
        return params  # Useful for debugging

    if params["internal_test"] == False:
        f = open(params["file"], "w+")
        if params["verbose"]:
            f.write(
                "PARAMETERS\n==========\n"
                + "\n".join(
                    k.rjust(14) + "  " + str(orig_params[k]) for k in orig_params.keys()
                )
                + "\n\n"
            )

    receptor_atom_scoria_mols = receptor.get_query_points_as_scoria_dummy_atoms(params)

    if receptor_atom_scoria_mols is None:
        return params

    # Go through the ligands in the directory and process each one. First,
    # format the inputs for use in the parallelizer.
    smarts = matches.get_smarts_substructs(params)
    if smarts is None:
        # Useful for debugging. If it gets here, you're running in test mode.
        # Because otherwise matches.get_smarts_substructs() would have thrown
        # an error.
        return params

    # If any of the identified smarts are None, the smarts string was probably
    # not valid.
    if not validate.all_smarts_valid(smarts, params):
        # One of the smarts substructures must be invalid.
        return params

    inputs = tuple(
        [
            tuple([ligand_filename, receptor_atom_scoria_mols, smarts, params])
            for ligand_filename in params["ligands"]
        ]
    )

    # Remove "ligands" from the params because it could potentially take up a
    # lot of memory. Why pass it to all the processors?
    del params["ligands"]

    # Run on multiple processors
    results = Parallelizer(
        mode=params["job_manager"], num_procs=params["num_processors"]
    ).run(
        inputs,
        parallel_filter_pose_files,
        params["num_processors"],
        params["job_manager"],
    )

    # Collect all the messages from the various processors.
    msgs = {
        "debug_msg": " ".join(list(set([r["debug_msg"] for r in results]))),
        "debug_warn": " ".join(list(set([r["debug_warn"] for r in results]))),
        "debug_error": " ".join(list(set([r["debug_error"] for r in results]))),
    }

    if params["internal_test"] == False:
        # Reformat all the messages to write them to the file.
        new_msg = {}
        for key in ["debug_msg", "debug_warn", "debug_error"]:
            new_msg[key] = (
                "\n".join(
                    [l.split(":", 1)[1].strip() for l in msgs[key].strip().split("\n")]
                )
                if msgs[key] != ""
                else ""
            )
        outpt = (
            (
                "MESSAGES\n========\n" + new_msg["debug_msg"]
                if params["verbose"]
                else new_msg["debug_msg"].strip()
            )
            + (
                "\n\nWARNINGS\n========\n" + new_msg["debug_warn"]
                if new_msg["debug_warn"] != ""
                else ""
            )
            + (
                "\n\nERRORS\n======\n" + new_msg["debug_error"]
                if new_msg["debug_error"] != ""
                else ""
            )
        )
        f.write(outpt)
        f.close()
        print('\nOutput saved to "' + params["file"] + '"\n')

    return msgs  # Useful for debugging


def parallel_filter_pose_files(
    ligand_filename, receptor_atom_scoria_mols, smarts, params
):
    """Applies the filters to a given ligand. Designed to enable parallel
    processing.

    :param ligand_filename: The filename of the ligand file.
    :type ligand_filename: str
    :param receptor_atom_scoria_mols: A list of .Molecule objects.
    :type receptor_atom_scoria_mols: list
    :param smarts: The smarts strings, mapped to corresponding RDKit mol
        objects.
    :type smarts: dict
    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :return: A dictionary containing any log, warning, and error messages.
    :rtype: dict
    """

    # Make a copy of params so it's by value, not reference.
    params = {k: params[k] for k in params.keys()}

    # Go through the ligands in the directory and load each ligand.
    # Create the ligand molecule based using the user-specified mode.
    ligand_mols = get_ligand_rdkit_mols(params, ligand_filename)

    if not validate.mols_created(ligand_mols, ligand_filename, params):
        # A required smi file doesn't exist or there is otherwise some
        # error. Continue on to the next ligand. It's a bad ligand
        # somehow.
        if params["verbose"]:
            output.warn(
                "Could not create RDKit molecule object: "
                + os.path.basename(ligand_filename),
                params,
            )
        return {k: params[k] for k in ["debug_msg", "debug_warn", "debug_error"]}

    for ligand_mol in ligand_mols:
        # Start by assuming all filters are matched. Set to false if any one
        # of them fails.
        passes = True

        # Loop through each of the filters to verify that the ligand possess all
        # the substructures. If it doesn't, no reason to proceed. Note that if
        # exclude is True, doesn't mater whether it has the substruucture or not.
        bad_filtr_substruct = ""
        for filtr in params["filters"]:
            if filtr["exclude"] == False:
                # Check if ligand_mol contains the required substructures.
                smrts = filtr["ligandSubstructSMARTS"]
                match = ligand_mol.HasSubstructMatch(smarts[smrts])
                if match == False:
                    passes = False
                    bad_filtr_substruct = filtr["ligandSubstructSMARTS"]
                    break

        if passes == False:
            # One of the substrctures is not present in the ligand. Move on to
            # the next ligand.
            if params["verbose"]:
                output.warn(
                    'Molecule "{ligand_filename}" does not contain substructure "{substruct}".'.format(
                        ligand_filename=os.path.basename(ligand_filename),
                        substruct=bad_filtr_substruct,
                    ),
                    params,
                )
            # return {k: params[k] for k in ["debug_msg", "debug_warn", "debug_error"]}
            continue

        # If you get here, the ligand does have all the required substrctures.
        # But are those substructures close enough to the specified receptor
        # atom (for each filter)? Let's check.

        # Go through all filters. Are they all satisfied?
        for i, filtr in enumerate(params["filters"]):
            receptor_atom_mol = receptor_atom_scoria_mols[i]

            # Make a dummy molecule with the coordinates of the ligand atoms
            # in the matching substructure of this filter.
            ligand_atoms_mol = matches.get_ligand_selection_scoria_mol(
                ligand_mol, smarts[filtr["ligandSubstructSMARTS"]]
            )
            if ligand_atoms_mol.get_coordinates() is None:
                # This happens if the substructure doesn't exist on the
                # ligand. Only possible if filtr["exclude"] is True, because
                # otherwise you would have rejected this ligand up above.
                # Since exclude is True, we can continue on to the next filter
                # without further evaluation.
                continue

            # Get the minimum distance between the receptor and ligand
            # matches.
            dist = receptor_atom_mol.get_distance_to_another_molecule(
                ligand_atoms_mol, True
            )

            # If the distance is more than the cutoff, stop checking filters.
            # They are not collectively satisfied.
            if filtr["exclude"] == True:
                if dist <= filtr["distance"]:
                    passes = False
                    if params["verbose"]:
                        output.warn(
                            'Molecule "{ligand_filename}" (pose {pose}) has substructure "{substruct}" within {dist} Å of the filter query point.'.format(
                                ligand_filename=os.path.basename(ligand_filename),
                                substruct=filtr["ligandSubstructSMARTS"],
                                dist=str(filtr["distance"]),
                                pose=ligand_mol.GetProp("pose"),
                            ),
                            params,
                        )
                    break
            else:
                # exclude is False
                if dist > filtr["distance"]:
                    passes = False
                    if params["verbose"]:
                        output.warn(
                            'Molecule "{ligand_filename}" (pose {pose}) has substructure "{substruct}" more than {dist} Å from the filter query point.'.format(
                                ligand_filename=os.path.basename(ligand_filename),
                                substruct=filtr["ligandSubstructSMARTS"],
                                dist=str(filtr["distance"]),
                                pose=ligand_mol.GetProp("pose"),
                            ),
                            params,
                        )
                    break

        # If it matches all filters, print out the name of the ligand:
        # ligand_filename.
        if passes == True:
            output.msg(
                "Molecule {ligand_filename} (pose {pose}) passes all filters.".format(
                    ligand_filename=os.path.basename(ligand_filename),
                    pose=ligand_mol.GetProp("pose"),
                ),
                params,
            )

    return {k: params[k] for k in ["debug_msg", "debug_warn", "debug_error"]}


if __name__ == "__main__":
    main(sys.argv)
