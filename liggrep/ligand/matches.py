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

from .. import scoria
from rdkit import Chem
from ..parameters import validate


def get_ligand_selection_scoria_mol(ligand_mol, ligand_atom_mol_from_smarts):
    """
       Isolates only the instances of the existing substructure within the ligand.
       Removes all other atoms.

       :param ligand_mol: A RDKit molecule object, containing the appropriate atom
           information and hybridization states of the ligand.
       :type ligand_mol: RDKit molecule object
       :param ligand_atom_mol_from_smarts: SMARTS string of the substructure within
           the ligand.
       :type ligand_atom_mol_from_smarts: str
       :returns: A single .Molecule object, containing the coordinates of
           the atoms in the substructure of the ligand (for each conformer
           separately).
       :rtype: .Molecule
    """

    # Ligand does contain substructure. Make a dummy molecule with the
    # coordinates of the atoms in the substructure (for each conformer
    # separately).
    coors = ligand_mol.GetConformers()[0].GetPositions()
    substruc_coors = []

    matches = ligand_mol.GetSubstructMatches(ligand_atom_mol_from_smarts)

    for match in matches:
        match = list(match)
        substruc_coors.append(coors[match])

    ligand_atoms_mol = scoria.Molecule()

    for coors in substruc_coors:
        for coord in coors:
            ligand_atoms_mol.add_atom(coordinates=coord)

    return ligand_atoms_mol


def get_smarts_substructs(params):
    """Associate the SMARTS strings of each small-molecule substructure in
    params["filters"] with an RDKit mol object.

    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :return: A mapping of the strings to the objects.
    :rtype: dict
    """

    smarts = {}
    for filtr in params["filters"]:
        smrts = filtr["ligandSubstructSMARTS"]
        mol = Chem.MolFromSmarts(smrts)
        smarts[smrts] = mol

        # Test if double bonds are inappropriately specified.
        if mol is not None and not validate.appropriate_bond_orders(params, mol, smrts):
            return None

    return smarts
