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

# Functions for loading molecules.

import os
from rdkit import Chem
from rdkit.Chem import AllChem
from ..output import error
import re
import io
from .. import scoria
import copy

# pdb, pdbqt, sdf


def load_ligand_mols(ligand_filename, params):
    """
       Load ligand into an RDKit molecule.

       :param ligand_filename: File name of the ligand pdb or pdbqt file.
       :type ligand_filename: str
       :param params: A dictionary of the user parameters and filters.
       :type params: dict
       :returns: A list of RDKit molecule objects, containing the
           appropriate atom information and hybridization states of the ligand
           conformations.
       :rtype: list
    """

    filename_root, filename_ext = os.path.splitext(ligand_filename)
    ligand_mols = None

    if filename_ext.lower() == ".pdb":
        # rdkit can load PDB files directly.
        pdb_text = open(ligand_filename).read()
        pdb_text = pdb_text.replace(  # prevents openbabel problem.
            "ENDMDL\nEND", "ENDMDL"
        )
        ligand_mols = [AllChem.MolFromPDBBlock(t.strip()) for t in pdb_text.split("END") if "ATOM" in t or "HETATM" in t]
    elif filename_ext.lower() == ".sdf":
        # Each SDF file must contain the same ligand. Can be different
        # conformers, but the same ligand.
        suppl = Chem.SDMolSupplier(ligand_filename, sanitize=False)
        ligand_mols = []
        for i, m in enumerate(suppl):
            fix_formal_charges(m)

            # Need to sanitize here, because didn't above.
            Chem.SanitizeMol(m)

            if m is not None:
                # Add conformer to existing molecule list.
                ligand_mols.append(m)
    elif filename_ext.lower() == ".pdbqt":
        # Deal with pdbqt here.
        # Get the text of the vina ligand
        vina_text = open(ligand_filename).read()

        # Get multiple frames
        prts = re.sub(r"^END$", "ENDMDL", vina_text, flags=re.MULTILINE)
        prts = re.split(r"^ENDMDL$", prts, flags=re.MULTILINE)

        # Remove ones that are empty (usually last one).
        prts = [p for p in prts if p.strip() != ""]

        # Load it into a RDKit Mol object.
        pdb_texts = []
        for p in prts:
            # Convert PDBQT to a PDB file.
            ligand_file_obj = io.StringIO(p)
            ligand_mol_convert = scoria.Molecule()
            ligand_mol_convert.load_pdbqt_into_using_file_object(ligand_file_obj)
            ligand_mol_convert.assign_elements_from_atom_names()
            pdb_text = ligand_mol_convert.save_pdb(return_text=True)
            pdb_texts.append(pdb_text)

        ligand_mols = [AllChem.MolFromPDBBlock(t) for t in pdb_texts]
    else:
        error(
            "Unrecognized file format: "
            + filename_ext
            + ". Please use PDB, PDBQT, or SDF (+/- SMI).",
            params,
        )

    # Add hydrogen atoms always.
    if ligand_mols is not None:
        ligand_mols = [Chem.AddHs(m) for m in ligand_mols if m is not None]

        # Also add pose number
        for i, m in enumerate(ligand_mols):
            m.SetIntProp("pose", i + 1)

    return ligand_mols

def fix_formal_charges(mol):
    """Unfortunately, RDKit and Open Babel don't always play nice. The specific
    problem is OpenBabel SDF files do not assign partial charges for
    nitrogens, etc. These must be determined based on the connectivity. Code
    below improves compatibility. Adapted from Dimorphite-DL.

    :param mol: The RDKit mol object to modify (in place)
    :type mol: RDKit mol object
    """

    for atom in mol.GetAtoms():
        explicit_bond_order_total = sum(
            [b.GetBondTypeAsDouble() for b in atom.GetBonds()]
        )

        # Assign the protonation charge, with special care for nitrogens
        element = atom.GetAtomicNum()
        if element == 7:  # nitrogens
            if explicit_bond_order_total == 4:
                atom.SetFormalCharge(1)
            elif explicit_bond_order_total == 3:
                atom.SetFormalCharge(0)
            elif explicit_bond_order_total == 2:
                atom.SetFormalCharge(-1)
        elif element == 8 or element == 16: # O and S
            if explicit_bond_order_total == 2:
                atom.SetFormalCharge(0)
            elif explicit_bond_order_total == 1:
                atom.SetFormalCharge(-1)
            # Note that SO4 has explicit_bond_order_total, so not changed
            # here.

    mol.UpdatePropertyCache(strict=False)
