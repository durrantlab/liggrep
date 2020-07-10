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
from . import scoria
from .parameters import validate

receptor_mol = None


def load_receptor(receptor):
    """
       Load receptor into a .Molecule object.

       :param receptor: File name of the receptor pdb or pdbqt file.
       :type receptor: str
       :returns: A single .Molecule object, containing the atomic
           information of the receptor protein.
       :rtype: .Molecule
    """

    receptor_mol = scoria.Molecule()
    root, ext = os.path.splitext(receptor)
    if ext == ".pdbqt":
        receptor_mol.load_pdbqt_into(receptor, False, False, False, False)
    elif ext == ".pdb":
        receptor_mol.load_pdb_into(receptor, False, False, False, False)

    return receptor_mol


def get_query_points_as_scoria_dummy_atoms(params):
    """Returns scoria molecule objects containing the query points, as
    .Molecule objects containing one atom.

    :param params: A dictionary of the user parameters and filters.
    :type params: dict
    :return: A list of .Molecule objects.
    :rtype: list
    """

    # First, get all the receptor atoms for each filter.
    receptor_atom_scoria_mols = []
    for filtr in params["filters"]:
        # Get a scoria molecule with the matching receptor atoms (or dummy
        # atoms if it's a coordinate filter).
        receptor_atom_scoria_mol = get_query_point_as_scoria_dummy_atom(params["receptor"], filtr)

        if not validate.receptor_filter_match_has_atoms(
            receptor_atom_scoria_mol, filtr, params
        ):
            return None

        receptor_atom_scoria_mols.append(receptor_atom_scoria_mol)

    return receptor_atom_scoria_mols


def get_query_point_as_scoria_dummy_atom(receptor, filtr):
    """
       Create a scoria Molecule with the specified receptor atoms or dummy atoms
       with the specified coordinates.

       :param receptor: File name of the receptor pdb or pdbqt file.
       :type receptor: str
       :param filtr: Dictionary containing receptor information, which can
           be a specified atom or coordinate, for each filter.
       :type filtr: dict
       :returns: A single .Molecule object, containing the specified atom
           or a dummy atom at the specified coordinate.
       :rtype: .Molecule
    """

    global receptor_mol

    receptor_atom_mol = None

    # Are any atom types? If so, then load the receptor (the actual protein).
    if "receptorAtom" in filtr:
        receptor_atom_inf = filtr["receptorAtom"]

        # Load receptor if it hasn't been already.
        if receptor_mol is None:
            receptor_mol = load_receptor(receptor)

        # Get the selection of the receptor atom.
        sel = {}

        if "chain" in receptor_atom_inf.keys():
            sel["chainid"] = receptor_atom_inf["chain"]
        if "resid" in receptor_atom_inf.keys():
            sel["resseq"] = receptor_atom_inf["resid"]
        if "atomname" in receptor_atom_inf.keys():
            sel["name"] = receptor_atom_inf["atomname"]

        receptor_atom_sel = receptor_mol.select_atoms(sel)

        # Make new molecule just containing that one atom.
        receptor_atom_mol = receptor_mol.get_molecule_from_selection(
            receptor_atom_sel, False, False
        )

    # Is it a coordinate type? If so, make a molecule with a dummy atom at the
    # coordinate.
    elif "coordinate" in filtr:
        coor = filtr["coordinate"]

        receptor_atom_mol = scoria.Molecule()
        receptor_atom_mol.add_atom(coordinates=[coor[0], coor[1], coor[2]])

    return receptor_atom_mol
