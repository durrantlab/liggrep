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

from . import smiles
from . import openbabel
from . import none


def get_ligand_rdkit_mols(params, ligand_filename):
    """Loads a ligand into a rdkit mol object, per the mode specified in
    params["mode"] (SMILES, OPENBABEL, or NONE).

    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :param ligand_filename: The ligand filename.
    :type ligand_filename: str
    :return: A list of ligands (each conformation, as an RDKit mol object).
    :rtype: list
    """

    ligand_mols = None
    if params["mode"] == "SMILES":
        ligand_mols = smiles.get_ligand_mols(ligand_filename, params)
        # Note that ligand_mol has hydrogen atoms without coordinate.
        # Shouldn't be a problem.
    elif params["mode"] == "OPENBABEL":
        ligand_mols = openbabel.get_ligand_mols(
            ligand_filename, params["babel_exec"], params
        )
    elif params["mode"] == "NONE":
        ligand_mols = none.get_ligand_mols(ligand_filename, params)

    return ligand_mols
