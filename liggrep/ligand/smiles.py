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

from .. import output
from . import load

from rdkit import Chem
from rdkit.Chem import AllChem
import os


def get_ligand_mols(ligand_filename, params):
    """
       Create a molecule of the ligand, using the smi file of the ligand to
       assign bond orders.

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

    template_smi = filename_root + ".smi"

    # If path to ligand.smi doesn't exist in the ligand directory, throw an
    # error.
    if (
        os.path.exists(template_smi) == False
        and os.path.exists(ligand_filename + ".smi") == False
    ):
        output.warn("Ligand SMI file does not exist: " + filename_root + ".smi", params)
        # all_filters_matched = False
        # ligand_mol = None
        return None

    mols = load.load_ligand_mols(ligand_filename, params)

    # Get the rdkit mol object.
    template_smi = [
        f for f in [template_smi, ligand_filename + ".smi"] if os.path.exists(f)
    ][0]
    smi_string = open(template_smi, mode="r").readlines()[
        0  # Because should only be one SMILES
    ]
    ref_mol = AllChem.MolFromSmiles(smi_string)
    if ref_mol is None:
        # output.msg("Trying to fix " + template_smi, params)
        smi_string2 = prevent_smiles_errors(smi_string)
        ref_mol = AllChem.MolFromSmiles(smi_string2)
        if ref_mol is not None:
            # output.msg("Successfully fixed " + template_smi, params)
            pass
        else:
            # output.msg("Trying to fix (again) " + template_smi, params)
            smi_string2 = prevent_smiles_errors2(smi_string)
            ref_mol = AllChem.MolFromSmiles(smi_string2)
            if ref_mol is not None:
                # output.msg("Successfully fixed " + template_smi, params)
                pass
            else:
                output.warn("Failed to load " + template_smi, params)

    try:
        ref_mol = Chem.RemoveHs(ref_mol)
        mols = [Chem.RemoveHs(m) for m in mols]

        # Note that the next line sometimes produces the following warning:
        # "WARNING: More than one matching pattern found - picking one". I
        # don't think this is a problem based on this link:
        # https://sourceforge.net/p/rdkit/mailman/message/31844451/ "If you
        # get this warning, it means that there is some symmetry in the
        # "all-single-bonds-stage" of your molecule. In your case, I guess
        # it's the carboxylic acids which can match two ways when there are
        # only single bonds."
        ligand_mols = [AllChem.AssignBondOrdersFromTemplate(ref_mol, m) for m in mols]
        ligand_mols = [Chem.AddHs(m) for m in ligand_mols]
    except:
        output.warn("ERROR processing " + ligand_filename + ". Skipping...", params)
        return None

    return ligand_mols


def prevent_smiles_errors(smiles):
    """
       Resolve valence problems from Open Babel generated smiles strings that
       are incompatible with RDKit.

       :param smiles: Open Babel generated smiles string.
       :type smiles: str
       :returns: Revised smiles string that is compatible with RDKit.
       :rtype: str
    """

    # I can't anticipate all smiles problems. Users are responsible for making
    # sure their smiles are properly formatted. But let's at least try to
    # avoid some of the more common problems...

    # These are to avoid valence problems that crop up in open-babel generated
    # smiles strings.
    smiles = smiles.replace("[NH]", "N")
    smiles = smiles.replace("[NH2]", "N")
    smiles = smiles.replace("[NH3]", "N")

    return smiles


def prevent_smiles_errors2(smiles):
    """
       Resolve amine protonation problems from Open Babel generated smiles strings
       that are incompatible with RDKit.

       :param smiles: Open Babel generated smiles string.
       :type smiles: str
       :returns: Revised smiles string that is compatible with RDKit.
       :rtype: str
    """

    # tertiary amines should be protonated.
    smiles = smiles.replace("[N]", "[N+]")
    smiles = smiles.replace("[N@]", "[N@+]")
    smiles = smiles.replace("[N@@]", "[N@@+]")
    smiles = smiles.replace("[nH]", "[nH+]")

    return smiles
