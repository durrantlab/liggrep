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

import os
import json
from .. import output


def file_exists(path, params):
    """
       Check that specified file exists.

       :param path: A path to a specified file as a string object.
       :type path: str
       :param params: A dictionary of the user parameters and filters.
       :type params: dict
       :returns: 'True' if file exists and 'False' otherwise.
       :rtype: bool
    """

    if os.path.exists(path):
        return True
    output.error("File does not exist: " + path, params)
    return False


def found_ligands(ligand_filenames, params):
    """Whether the user-specified ligand-file mask matches any files. If so,
    adds those to the params object.

    :param ligand_filenames: A list of the filenames, if any.
    :type ligand_filenames: list
    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :return: 'True' if it validates, 'False' otherwise.
    :rtype: bool
    """

    if len(ligand_filenames) == 0:
        output.error('No ligand files match "' + params["ligands"] + '"', params)
        return False
    else:
        params["ligands"] = ligand_filenames
        return True


def proper_filter_keys(filtr, params):
    """Checks if a given filter contains the proper keys.

    :param filtr: The individual filter data.
    :type filtr: dict
    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :return: 'True' if it validates, 'False' otherwise.
    :rtype: bool
    """

    # If a filter has both "receptorAtom" and "coordinate" keys, throw an
    # error. Must be either or.
    if "receptorAtom" in filtr and "coordinate" in filtr:
        output.error(
            'One of your filters contains both a "receptorAtom" and '
            + '"coordinate" key. You must pick one or the other: '
            + json.dumps(filtr),
            params,
        )
        return False

    # But it must have one or the other.
    if not "receptorAtom" in filtr and not "coordinate" in filtr:
        output.error(
            'All filters must contain either a "receptorAtom" or a '
            + '"coordinate" key (but not both). One of your filters '
            + "has neither: "
            + json.dumps(filtr),
            params,
        )
        return False

    required_keys = [
        "ligandSubstructSMARTS",
        "distance",
    ]

    acceptable_keys = required_keys + ["exclude", "coordinate", "receptorAtom"]

    # Make sure all the keys of the filter are in the acceptable_keys list.
    for k in filtr:
        if not k in acceptable_keys:
            output.error(
                "One of your filters contains an invalid key: " + json.dumps(filtr),
                params,
            )
            return False

    # Make sure the filter contains all required filters.
    for k in required_keys:
        if not k in filtr.keys():
            output.error(
                'One of your filters does not contain the required key "'
                + k
                + '": '
                + json.dumps(filtr),
                params,
            )
            return False

    return True


def obabel_specified_if_needed(params):
    """Checks if the user has specified a valid path to the obabel executable.

    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :return: 'True' if it validates, 'False' otherwise.
    :rtype: bool
    """

    if params["mode"] == "OPENBABEL":
        if params["babel_exec"] == None:
            output.error(
                "You're running LigGrep in OPENBABEL mode, but you haven't specified "
                + "the open-babel executable via the --babel_exec parameter.",
                params,
            )
            return False
        elif not os.path.basename(params["babel_exec"].lower()).startswith("obabel"):
            output.error(
                "The open-babel executable specified via the --babel_exec paramter "
                + 'does not start with "obabel": "'
                + params["babel_exec"]
                + '"',
                params,
            )
            return False
        elif not os.path.exists(params["babel_exec"]):
            output.error(
                "You're running LigGrep in OPENBABEL mode, but the open-babel executable "
                + 'that you specified via the --babel_exec parameter does not exist: "'
                + params["babel_exec"]
                + '"',
                params,
            )
            return False
    return True


def smi_files_exist_if_needed(params):
    """Check if smi files exist for each ligand, if needed.

    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :return: 'True' if it validates, 'False' otherwise.
    :rtype: bool
    """

    if params["mode"] == "SMILES":
        for lig_filename in params["ligands"]:
            basename = os.path.splitext(lig_filename)[0]
            if not os.path.exists(lig_filename + ".smi") and not os.path.exists(
                basename + ".smi"
            ):
                output.error(
                    "You're running LigGrep in SMILES mode, but the ligand file \""
                    + lig_filename
                    + '"'
                    + ' does not have an associated SMILES file named "'
                    + lig_filename
                    + '.smi" '
                    + 'or "'
                    + basename
                    + '.smi".',
                    params,
                )
                return False
    return True


def receptor_filter_match_has_atoms(receptor_atom_scoria_mol, filtr, params):
    """Verifies that receptor-atom specifications match actual atoms.

    :param receptor_atom_scoria_mol: A .Molecule object of the receptor.
    :type receptor_atom_scoria_mol: .Molecule
    :param filtr: The filter data.
    :type filtr: dict
    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :return: 'True' if it validates, 'False' otherwise.
    :rtype: bool
    """

    if len(receptor_atom_scoria_mol.get_coordinates()) == 0:
        # Throw an error. Receptor selection matches no atoms.
        output.error(
            "No receptor atom matches filter: " + json.dumps(filtr), params,
        )
        return False
    return True


def appropriate_bond_orders(params, smrts_mol, smrts):
    """Checks if a SMARTS substring specification has appropriate bond orders
    given the user-specified mode.

    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :param smrts_mol: RDKit mol object of the SMARTS string.
    :type smrts_mol: RDKit mol object.
    :param smrts: The SMARTS string.
    :type smrts: str
    :return: 'True' if it validates, 'False' otherwise.
    :rtype: bool
    """

    # Test if double bonds are inappropriately specified.
    if params["mode"] == "NONE" and (
        ".pdb" in params["ligand_exts"] or ".pdbqt" in params["ligand_exts"]
    ):
        bond_orders = [b.GetBondTypeAsDouble() for b in smrts_mol.GetBonds()]
        bond_orders = [o for o in bond_orders if o != 1.0]
        if len(bond_orders) > 0:
            # So it has bonds with orders greater than 1
            output.error(
                "When processing PDB- and PDBQT-formatted ligands in NONE "
                + "mode, LigGrep ignores bond orders and simply "
                + "assumes that all appropriately juxtaposed atoms are "
                + "connected by single bonds. But one (or more) of your "
                + "filters describes a substructure with bonds of higher "
                + "orders: "
                + smrts,
                params,
            )
            return False
    return True


def no_sdf_in_openbabel_smiles_mode(params):
    """Make sure the user hasn't specified OPENBABEL mode when processing SDF
    files.

    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :return: 'True' if it validates, 'False' otherwise.
    :rtype: bool
    """

    # Do a little bit of validation here.
    if params["mode"] in ["OPENBABEL", "SMILES"] and ".sdf" in params["ligand_exts"]:
        # OPENBABEL/SMILES mode is not appropriate for SDF files, because
        # these files already include bond-order information. Throw error any
        # of the ligands are SDF files.
        output.error(
            "OPENBABEL/SMILES mode is not appropriate for SDF files, because "
            + "these files already include bond-order information. Use NONE "
            + "mode instead.",
            params,
        )
        return False
    return True


def mols_created(ligand_mols, ligand_filename, params):
    """Checks to make sure a mol has been successfully created.

    :param ligand_mols: A list of RDKit object of the potential ligand.
    :type ligand_mols: list
    :param ligand_filename: The file name of the ligand.
    :type ligand_filename: str
    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :return: 'True' if it validates, 'False' otherwise.
    :rtype: bool
    """

    if ligand_mols == None or len(ligand_mols) == 0:
        # A required smi file doesn't exist or there is otherwise some
        # error. Continue on to the next ligand. It's a bad ligand
        # somehow.
        output.warn(
            "Could not process ligand in " + ligand_filename + ". Skipping...", params,
        )
        return False
    return True


def all_smarts_valid(smarts, params):
    """[summary]

    :param smarts: The smarts strings, mapped to corresponding RDKit mol
        objects.
    :type smarts: dict
    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :return: 'True' if it validates, 'False' otherwise.
    :rtype: bool
    """

    for smrt in smarts:
        if smarts[smrt] is None:
            # User-specified SMARTS string is likely invalid.
            output.error(
                "User-specified SMARTS string is not valid: " + smrt, params,
            )
            return False
    return True
