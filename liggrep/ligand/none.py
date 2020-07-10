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

from . import load


def get_ligand_mols(ligand_filename, params):
    """
       Create a molecule of the ligand without attempting to assign bond orders.
       Docked ligands must either be in a format that itself specifies bond
       orders or the user-specified filters must not depend on bond orders.

       :param ligand_filename: File name of the ligand pdb or pdbqt file.
       :type ligand_filename: str
       :param params: A dictionary of the user parameters and filters.
       :type params: dict
       :returns: A list of RDKit molecule objects, containing the
           appropriate atom information and hybridization states of the ligand
           conformations.
       :rtype: list
    """

    return load.load_ligand_mols(ligand_filename, params)
